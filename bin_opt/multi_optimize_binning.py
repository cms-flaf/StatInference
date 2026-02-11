import bayes_opt
import json
import math
import numpy as np
import os
import queue
import random
import re
import ROOT
import shutil
import threading
import time
from sortedcontainers import SortedSet
import yaml
from optimize_binning import HistToNumpy, arrayToStr, strtobool, Yields, ExtractYields

input_binning_opt_config = os.path.join(os.environ["ANALYSIS_PATH"], "StatInference", "bin_opt", "bin_optimization.yaml")
with open(input_binning_opt_config, "r") as f:
    input_binning_opt_config_dict = yaml.safe_load(f)

if input_binning_opt_config_dict["input"]["num_original_bins"] == 1000:
    min_step = 0.001
    max_value_int = input_binning_opt_config_dict["input"]["num_original_bins"]
elif input_binning_opt_config_dict["input"]["num_original_bins"] == 5000:
    min_step = 0.0002
    max_value_int = input_binning_opt_config_dict["input"]["num_original_bins"]

class MultiBinning:
    def __init__(self, edges):
        self.edges = edges
        self.exp_limit = None
        self.input_datacard = None
    

    @staticmethod
    def fromEntry(entry):
        binning = MultiBinning(np.array(entry['bin_edges']))
        binning.exp_limit = float(entry['exp_limit'])
        binning.input_datacard = entry['input_datacard']
        return binning
    
    def isEquivalent(self, other):
        if len(self.edges) != len(other.edges):
            return False
        if self.input_datacard != other.input_datacard:
            return False
        return np.count_nonzero(self.edges == other.edges) == len(self.edges)

    def toPoint(self, n_thrs, n_cards):
        rel_thrs = self.getRelativeThresholds(n_thrs, n_cards)
        point = {}
        for idx in range(n_cards):
            for n in range(n_thrs):
                point[f'rel_thr_card{idx}_{n}'] = rel_thrs[idx][n]
        return point

    def isBetter(self, other):
        if other is None or other.exp_limit is None:
            return True
        if self.exp_limit != other.exp_limit:
            return self.exp_limit < other.exp_limit
        return len(self.edges) < len(other.edges)
        
    
    def getRelativeThresholds(self, n_thrs, n_cards):
        n_inner_edges = len(self.edges) - 2
        if n_thrs < n_inner_edges:
            raise RuntimeError(f"Invalid number of output thresholds: {n_thrs} should not be < {n_inner_edges}")
        rel_thrs = {}
        for i in range(n_cards):
            rel_thrs[i] = (np.random.rand(n_thrs) + min_step) / (1 - min_step)
            rel_thrs[i][0:n_inner_edges] = np.flip(self.edges[2:] - self.edges[1:-1])
            if n_inner_edges < n_thrs:
                rel_thrs[i][n_inner_edges] = max(self.edges[1], rel_thrs[i][n_inner_edges])
        return rel_thrs

    @staticmethod
    def fromRelativeThresholds(rel_thrs, bkg_yields):
        edges = [ max_value_int ]
        prev_yield = -1
        for rel_thr in rel_thrs:
            edge_up = edges[-1]
            edge_down = max(edge_up - int(round(rel_thr * step_int_scale)), 0)
            all_ok = False
            while edge_down >= 0:
                all_ok, new_yield = bkg_yields.test(edge_down, edge_up, prev_yield)
                if all_ok: break
                edge_down -= 1
            if all_ok:
                edges.append(edge_down)
                prev_yield = new_yield
            if edge_down == 0: break
        if len(edges) == 1:
            edges.append(0)
        edges.reverse()
        edges = np.array(edges, dtype=float)
        edges = edges / step_int_scale
        if edges[-1] != 1:
            edges[-1] = 1
        if edges[0] != 0:
            edges[0] = 0
        return Binning(edges)

class MultiBayesianOptimization:
    def __init__(self, max_n_bins, working_area, workers_dir, acq, kappa, xi,
                input_datacards_list, poi, bkg_yields_dict, input_queue_size, random_seed=12345):
        self.max_n_bins = max_n_bins
        self.binnings = []
        self.best_binning = None
        self.working_area = working_area
        self.workers_dir = workers_dir
        self.log_output = os.path.join(working_area, 'results_multibinning.json')
        self.poi = poi
        self.best_binning_split = False
        self.input_queue = queue.Queue(input_queue_size)
        self.open_requests = []
        self.binning_lock = threading.Lock()
        self.optimizer_lock = threading.Lock()
        self.print_lock = threading.Lock()
        self.bkg_yields_dict = bkg_yields_dict

        self.input_datacards = {}
        bounds = {}
        for idx,card in enumerate(input_datacards_list):
            self.input_datacards[idx] = card
            for n in range(max_n_bins):
                upper_bound = 1. if n > 0 else 1.5
                bounds[f'rel_thr_card{idx}_{n}'] = (min_step, upper_bound) 

        self.bounds_transformer = None
        self.optimizer = bayes_opt.BayesianOptimization(f=None, pbounds=bounds, random_state=random_seed, verbose=1)
        self.utilities = [
            bayes_opt.acquisition.UpperConfidenceBound(kappa=kappa),
            bayes_opt.acquisition.ExpectedImprovement(xi=xi),
            bayes_opt.acquisition.ProportionalIntegralDerivative()
        ]

        if not os.path.isdir(working_area):
            os.makedirs(working_area)
        if not os.path.isdir(workers_dir):
            os.makedirs(workers_dir)
        if os.path.isfile(self.log_output):
            with open(self.log_output, 'r') as f:
                prev_binnins = json.loads('[' + ', '.join(f.readlines()) + ']')
            for entry in prev_binnings:
                binning = MultiBinning.fromEntry(entry)
                if self.findEquivalent(binning) is None:
                    point = binning.toPoint(self.max_n_bins, len(self.input_datacards))
                    self.register(point, len(binning.edges), binning.exp_limit)
                    if binning.isBetter(self.best_binning):
                        self.best_binning = binning
                    self.binnings.append(binning)
        
        if self.best_binning is not None:
            self.print(f'The best binning from previous iterations: {arrayToStr(self.best_binning.edges)}, exp_limit = {self.best_binning.exp_limit}, datacard = {self.best_binning.input_datacard}')

        self.bounds_transformer = bayes_opt.SequentialDomainReductionTransformer()
        self.bounds_transformer.initialize(self.optimizer.space)

        suggestions_file = os.path.join(working_area, 'to_try.json')
        self.suggestions = []
        if os.path.isfile(suggestions_file):
            with open(suggestions_file, 'r') as f:
                suggestions = json.load(f)
            for edges in suggestions:
                self.addSuggestion(edges)

    def findEquivalent(self, binning, lock=True):
        equivalent_binning = None
        if lock:
            self.binning_lock.acquire()
        for prev_binning in self.binnings:
            if prev_binning.isEquivalent(binning):
                equivalent_binning = prev_binning
                break
        if lock:
            self.binning_lock.release()
        return equivalent_binning
    
    def register(self, point, n_edges, exp_limit):
        loss = -(exp_limit*1e6 + n_edges)
        self.optimizer_lock.acquire()
        self.optimizer.register(params=point, target=loss)
        self.updateBounds(n_edges - 2)
        self.optimizer_lock.release()
    
    def print(self, msg):
        self.print_lock.acquire()
        print(msg)
        self.print_lock.release()

    def updateBounds(self, n_points):
        if self.bound_transformer:
            new_bounds = self.bounds_transformer.transform(self.optimizer.space)
            prev_fixed = True
            for cat_i in range(len(self.input_datacards)):
                for n in range(self.max_n_bins):
                    key = f'rel_thr_card{cat_i}_{n}'
                    if prev_fixed:
                        prev_fixed = new_bounds[key][1] - new_bounds[key][0] < 0.1 and n < n_points
                    else:
                        new_bounds[key] = np.array([min_step, 1.])
            for key in new_bounds:
                if new_bounds[key][0] >= new_bounds[key][1]:
                    new_bounds[key] = np.array([new_bounds[key][1], new_bounds[key][0]])
                if new_bounds[key][0] < min_step:
                    new_bounds[key][0] = min_step
                
                upper_bound = 1. if key != f'rel_thr_card{cat_i}_0' else 1.5
                if new_bounds[key][1] > upper_bound:
                    new_bounds[key][1] = upper_bound
                if np.any(np.isnan(new_bounds[key])):
                    new_bounds[key] = np.array([min_step, upper_bound])
            self.optimizer.set_bounds(new_bounds)
    
    def suggest(self, utility_index):
        self.optimizer_lock.acquire()
        if self.suggestions is None or len(self.suggestions) == 0:
            self.optimizer._acquisition_function = self.utilities[utility_index]
            point = self.optimizer.suggest()
        else:
            point = self.suggestions.pop(0)
            self.suggestions.remove(point)
        self.optimize_lock.release()
        return point

    def addSuggestion(self, edges):
        suggested_binning = Binning(np.array(edges))
        rel_thrs = suggested_binning.getRelativeThresholds(self.max_n_bins, len(self.input_datacards))
        binning = Binning.fromRelativeThresholds(rel_thrs, self.bkg_yields)
        if self.findEquivalent(binning, False) is None:
            point = binning.toPoint(self.max_n_bins, len(self.input_datacards))
            self.suggestions.append(point)

###check below for compatibility with multi category optimization

    def addNewBinning(self, binning):
        self.binning_lock.acquire()
        if binning.isBetter(self.best_binning):
            self.best_binning = binning
            self.best_binning_split = False
            self.print('New best binning: {}, exp_limit = {}' \
                       .format(arrayToStr(binning.edges), binning.exp_limit))
        self.binnings.append(binning)
        self.binning_lock.release()

    def tryBestBinningSplit(self):
        self.binning_lock.acquire()
        if not self.best_binning_split:
            self.print('Splitting best binning...')
            n_best = len(self.best_binning.edges)
            self.optimizer_lock.acquire()
            if n_best - 2 < self.max_n_bins:
                for k in range(n_best - 1):
                    edges = np.zeros(n_best + 1)
                    edges[0:k+1] = self.best_binning.edges[0:k+1]
                    edges[k+2:] = self.best_binning.edges[k+1:]
                    for delta_edge in [ 0.1, 0.5, 0.9 ]:
                        edges[k+1] = self.best_binning.edges[k] \
                                     + (self.best_binning.edges[k+1] - self.best_binning.edges[k]) * delta_edge
                        self.addSuggestion(edges)
            if n_best > 2:
                for k in range(1, n_best):
                    edges = np.zeros(n_best - 1)
                    edges[0:k] = self.best_binning.edges[0:k]
                    edges[k:] = self.best_binning.edges[k+1:]
                    self.addSuggestion(edges)
            self.optimizer_lock.release()
            self.best_binning_split = True
        self.binning_lock.release()

    def findOpenRequest(self, binning):
        equivalent_binning = None
        self.binning_lock.acquire()
        for prev_binning in self.open_requests:
            if prev_binning.isEquivalent(binning):
                equivalent_binning = prev_binning
                break
        self.binning_lock.release()
        return equivalent_binning

    def addOpenRequest(self, binning):
        self.binning_lock.acquire()
        self.open_requests.append(binning)
        self.binning_lock.release()

    def clearOpenRequest(self, binning):
        self.binning_lock.acquire()
        to_remove = []
        for prev_binning in self.open_requests:
            if prev_binning.isEquivalent(binning):
                to_remove.append(prev_binning)
        for b in to_remove:
            self.open_requests.remove(b)
        self.binning_lock.release()

    def waitOpenRequestsToFinish(self):
        has_open_requests = True
        while has_open_requests:
            self.binning_lock.acquire()
            has_open_requests = len(self.open_requests) > 0
            self.binning_lock.release()
            time.sleep(1)

    def PointGenerator(self, n_eq_steps):
        n = 0
        utility_index = 0
        open_request_sleep = 1
        while n < n_eq_steps:
            point = self.suggest(utility_index)

            rel_thrs = np.zeros(len(point))
            for k in range(self.max_n_bins):
                rel_thrs[k] = point['rel_thr_{}'.format(k)]
            self.print('rel_thrs: {}'.format(arrayToStr(rel_thrs)))

            binning = Binning.fromRelativeThresholds(rel_thrs, self.bkg_yields)
            equivalent_binning = self.findEquivalent(binning)
            if equivalent_binning is None:
                n = 0
                open_request = self.findOpenRequest(binning)
                if open_request is None:
                    self.print('Next binning to probe: {}'.format(arrayToStr(binning.edges)))
                    self.input_queue.put(binning)
                    self.addOpenRequest(binning)
                    utility_index = 0
                    open_request_sleep = 1
                else:
                    self.print('Open request for binning found: {}'.format(arrayToStr(open_request.edges)))
                    time.sleep(open_request_sleep)
                    open_request_sleep = open_request_sleep * 2
                    utility_index = (utility_index + 1) % len(self.utilities)
            else:
                self.print('Equivalent binning found: {}'.format(arrayToStr(equivalent_binning.edges)))
                self.register(point, len(equivalent_binning.edges), equivalent_binning.exp_limit)
                n += 1
                if n >= 5:
                    self.tryBestBinningSplit()
                if n == n_eq_steps - 1:
                    self.print('Waiting for open requests to finish...')
                    self.waitOpenRequestsToFinish()

        self.input_queue.put(None)

    def JobDispatcher(self):
        while True:
            time.sleep(1)
            for worker_dir in os.listdir(self.workers_dir):
                worker_dir = os.path.join(self.workers_dir, worker_dir)
                if not os.path.isdir(worker_dir): continue
                task_file = os.path.join(worker_dir, 'task.txt')
                task_file_tmp = os.path.join(worker_dir, '.task.txt')
                if os.path.isfile(task_file):
                    result_file = os.path.join(worker_dir, 'result.txt')
                    if os.path.isfile(result_file):
                        with open(result_file, 'r') as f:
                            result = json.load(f)
                        if result['input_datacard'] == self.input_datacard:
                            binning = Binning.fromEntry(result)
                            self.clearOpenRequest(binning)
                            if self.findEquivalent(binning) is None:
                                point = binning.toPoint(self.max_n_bins)
                                self.addNewBinning(binning)
                                with open(self.log_output, 'a') as f:
                                    f.write(json.dumps(result) + '\n')
                                self.register(point, len(binning.edges), binning.exp_limit)

                        os.remove(task_file)
                        os.remove(result_file)
                else:
                    try:
                        binning = self.input_queue.get(True, 1)
                        if binning is None: return
                        task = {
                            'input_datacard': self.input_datacard,
                            'bin_edges': [ x for x in binning.edges ],
                            'poi': self.poi,
                            'other_datacards': self.other_datacards,
                        }
                        with open(task_file_tmp, 'w') as f:
                            json.dump(task, f)
                        shutil.move(task_file_tmp, task_file)
                    except queue.Empty:
                        pass

    def maximize(self, n_eq_steps):

        def generator():
            self.PointGenerator(n_eq_steps)

        def dispatcher():
            self.JobDispatcher()

        generator_thread = threading.Thread(target=generator)
        generator_thread.start()

        dispatcher_thread = threading.Thread(target=dispatcher)
        dispatcher_thread.start()

        generator_thread.join()
        dispatcher_thread.join()


###end of check for multi_bo

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Find optimal binning for multiple catgories that minimizes the expected limits.')
    parser.add_argument('--input', required=True, type=str, help="comma separated input datacards, e.g. datacard_cat1,datacard_cat2,..")
    parser.add_argument('--shape-file', required=True, type=str, help="input root file with shapes, comma separated")
    parser.add_argument('--output', required=True, type=str, help="output directory")
    parser.add_argument('--workers-dir', required=True, type=str, help="output directory for workers results")
    parser.add_argument('--max_n_bins', required=True, type=int, help="maximum number of bins per category")
    parser.add_argument('--poi', required=False, type=str, default='r', help="parameter of interest")
    parser.add_argument('--params', required=False, type=str, default=None,
                        help="algorithm parameters in format param1=value1,param2=value2,...")
    parser.add_argument('--verbosity', required=False, type=int, default=1, help="verbosity")

    args = parser.parse_args()

    param_dict = {}
    if args.params is not None:
        params = args.params.split(',')
        for param in params:
            p = param.split('=')
            if len(p) != 2:
                raise RuntimeError(f'invalid parameter definition {param}.')
            param_dict[p[0]] = p[1]
    
    ref_bkgs = {}
    for main_bkg in input_binning_opt_config_dict["background"]["main"]:
        ref_bkgs[f'{main_bkg}'] = re.compile(f'^{main_bkg}$')
    
    nonbkgg_regex = re.compile('(data_obs|^ggHH.*|^qqHH.*|^DY_[0-2]b.*)')
    ignore_unc_variations = re.compile('(CMS_bbtt_201[6-8]_DYSFunc[0-9]+|CMS_bbtt_.*_QCDshape)(Up|Down)')

    input_datacards = args.input.split(',')
    for idx,card in enumerate(input_datacards):
        input_datacards[idx] = os.path.abspath(card)

    bkg_yields = {}
    input_shape_files = args.shape_file.split(',')
    for idx,shape in enumerate(input_shape_files):
        input_shape_files[idx] = os.path.abspath(shape)
        print("Extracting yields for background processes", ref_bkgs, "and shape", shape)
        bkg_yields[idx] = ExtractYields(input_shape_files[idx], ref_bkgs, nonbkg_regex, ignore_unc_variations)
        bkg_yields[idx].printSummary()

        for name,value in param_dict.items():
            if not hasattr(bkg_yields[idx], name):
                raise RuntimeError(f'Unknown parameter {name}')
            def_value = getattr(bkg_yields[idx], name)
            def_value_type = type(def_value)
            if def_value_type == bool:
                new_value = bool(strtobool(value))
            else:
                new_value = def_value_type(value) #this might throw an error because def_value_type function is not defined.

            print(f"Setting {name} = {new_value}")
            setattr(bkg_yields[idx], name, new_value)

    multi_bo = MultiBayesianOptimization(max_n_bins=args.max_n_bins,
                                            working_area=args.output,
                                            workers_dir=args.workers_dir,
                                            acq='ucb', kappa=2.57, xi=0.1,
                                            input_datacards=input_datacards,
                                            poi=args.poi,
                                            bkg_yields=bkg_yields,
                                            input_queue_size=2, random_seed=None)

    multi_bo.maximize(20)

    print("Minimization finished.")
    print("Best binning per category:\n")
