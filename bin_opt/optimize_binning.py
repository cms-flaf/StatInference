import bayes_opt
import json
import math
import numpy as np
import os
import queue
from queue import Queue
import random
import re
import ROOT
import shutil
import threading
import time
from sortedcontainers import SortedSet
import distutils.util

min_step = 0.001
min_step_digits = -int(math.log10(min_step))
step_int_scale = 10 ** min_step_digits
max_value_int = step_int_scale

def HistToNumpy(hist):
    epsilon = 1e-7
    x = np.zeros(max_value_int)
    x_err2 = np.zeros(max_value_int)
    if hist.GetNbinsX() != max_value_int:
        raise RuntimeError("Inconsistent number of bins")
    for n in range(max_value_int):
        if abs(hist.GetBinLowEdge(n + 1) - n * min_step) > epsilon \
                or abs(hist.GetBinLowEdge(n + 2) - (n + 1) * min_step) > epsilon:
            raise RuntimeError("Inconsistent binning")
        x[n] = hist.GetBinContent(n + 1)
        x_err2[n] = hist.GetBinError(n + 1) ** 2
    return x, x_err2

def arrayToStr(a):
    return '[ ' + ', '.join([ str(x) for x in a ]) + ' ]'

class Yields:
    from sortedcontainers import SortedSet
    def __init__(self, ref_bkgs, n_bins):
        self.n_bins = n_bins
        self.ref_bkgs = ref_bkgs
        self.yields = {}
        self.processes = SortedSet([ 'total' ] + list(ref_bkgs.keys()))
        # print(f'debug: self.processes: {self.processes}')
        self.input_processes = SortedSet()
        # print(f'debug: self.input_processes: {self.input_processes}')
        self.unc_variations = SortedSet([''])
        self.ref_bkg_thr = 1e-7
        self.total_bkg_thr = 0.03 #0.18
        self.consider_non_central = True
        self.monotonicity_relax_thr = 0.25
        self.rel_err_thr = 1.
        self.rel_err_relax_thr = 1.
        self.rel_unc_variation_thr = 0.2

    def addProcess(self, process, yields, err2, unc_variation):
        names = [ 'total' ]
        for ref_bkg_name, ref_bkg_regex in self.ref_bkgs.items():
            # print(f'Debug: ref_bkg_name: {ref_bkg_name}, ref_bkg_regex: {ref_bkg_regex}, process: {process}')
            if ref_bkg_regex.match(process) is not None:
                names.append(ref_bkg_name)
        if process not in self.input_processes:
            self.input_processes.add(process)
            # print(f'debug: Added process {process} to input_processes: {self.input_processes}')
        for name in names:
            if unc_variation not in self.unc_variations:
                self.unc_variations.add(unc_variation)
            key = (name, unc_variation)
            # print(f'key: {key}')
            # print(f'self.yields {self.yields}')
            if key not in self.yields:
                # print(f'self.n_bins {self.n_bins} np.zeros(self.n_bins) {np.zeros(self.n_bins)}')
                self.yields[key] = [ np.zeros(self.n_bins), np.zeros(self.n_bins) ]
            self.yields[key][0] += yields
            self.yields[key][1] += err2

    def test(self, start, stop, yield_thr):
        total_yield = {}
        # print(f'Debug: self.processes {self.processes}')
        for process in self.processes:
            # print(f'Debug: process {process}')
            is_ref = process != 'total'
            # print(f'Debug: self.unc_variations {self.unc_variations}')
            for unc_variation in self.unc_variations:
                # print(f'Debug: unc_variation {unc_variation}')
                is_central = unc_variation == ''
                # print(f'Debug: is_central {is_central}, self.consider_non_central {self.consider_non_central}')
                if not is_central and not self.consider_non_central: continue
                key = (process, unc_variation)
                # print(f'Debug: key {key}')
                # print(f'Debug: self.yields.get(key, None) {self.yields.get(key, None)}')
                # print(f'Debug: self.yields[key][0] {self.yields[key][0]} ')
                # print(f'Debug: self.yields[key][1] {self.yields[key][1]}')
                if self.yields.get(key, None) is None: continue
                sum = np.sum(self.yields[key][0][start:stop])
                err2 = np.sum(self.yields[key][1][start:stop])
                err = math.sqrt(err2)

                if key not in self.yields:
                    return False, -1
                if is_ref:
                    thr = self.ref_bkg_thr
                elif is_central:
                    if err / sum <= self.monotonicity_relax_thr:
                        thr = self.total_bkg_thr
                    else:
                        thr = max(self.total_bkg_thr, yield_thr)
                else:
                    thr = self.total_bkg_thr
                if sum < thr:
                    return False, -1
                if not is_ref:
                    total_yield[unc_variation] = (sum, err)
        rel_err = total_yield[''][1] / total_yield[''][0]
        max_delta = 0
        if self.consider_non_central:
            for unc_variation in self.unc_variations:
                max_delta = max(max_delta, abs(total_yield[''][0] - total_yield[unc_variation][0]))
        if rel_err <= self.rel_err_thr or (self.consider_non_central and rel_err <= self.rel_err_relax_thr \
                                           and max_delta / total_yield[''][0] < self.rel_unc_variation_thr):
            return True, total_yield[''][0]
        return False, -1

    def printSummary(self):
        print('total = {}'.format(' + '.join(self.input_processes)))
        for unc_variation in self.unc_variations:
            unc_name = unc_variation if len(unc_variation) > 0 else 'central'
            print(unc_name)
            for process in self.processes:
                x = 0
                err = 0
                key = (process, unc_variation)
                if key in self.yields:
                    x = np.sum(self.yields[key][0])
                    err = math.sqrt(np.sum(self.yields[key][1]))
                print('\t{}: {:.3f} +/- {:.3f}'.format(process, x, err))

def ExtractYields(input_shapes, ref_bkgs, nonbkg_regex, ignore_variations_regex):
    yields = Yields(ref_bkgs, max_value_int)

    input_root = ROOT.TFile.Open(input_shapes)
    print(f'Extracting yields from {input_shapes}')
    hist_names = [ str(key.GetName()) for key in input_root.GetListOfKeys() ]
    nuis_name_regex = re.compile('(.*)_(CMS_.*(Up|Down))')
    for hist_name in sorted(hist_names):
        if nonbkg_regex.match(hist_name) is not None:
            continue
        nuis_match = nuis_name_regex.match(hist_name)
        if nuis_match is not None:
            process = nuis_match.groups()[0]
            unc_variation = nuis_match.groups()[1]
        else:
            process = hist_name
            unc_variation = ''
        if ignore_variations_regex.match(unc_variation):
            continue

        hist = input_root.Get(hist_name)
        hist_yield, hist_err2 = HistToNumpy(hist)
        yields.addProcess(process, hist_yield, hist_err2, unc_variation)
    input_root.Close()
    return yields

class Binning:
    def __init__(self, edges):
        self.edges = edges
        self.exp_limit = None

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

    @staticmethod
    def fromEntry(entry):
        binning = Binning(np.array(entry['bin_edges']))
        binning.exp_limit = float(entry['exp_limit'])
        return binning

    def getRelativeThresholds(self, n_thrs):
        n_inner_edges = len(self.edges) - 2
        if n_thrs < n_inner_edges:
            raise RuntimeError("Invalid number of output thresholds")
        rel_thrs = (np.random.rand(n_thrs) + min_step) / (1 - min_step)
        rel_thrs[0:n_inner_edges] = np.flip(self.edges[2:] - self.edges[1:-1])
        if n_inner_edges < n_thrs:
            rel_thrs[n_inner_edges] = max(self.edges[1], rel_thrs[n_inner_edges])
        return rel_thrs

    def toPoint(self, n_thrs):
        rel_thrs = self.getRelativeThresholds(n_thrs)
        point = {}
        for n in range(n_thrs):
            point['rel_thr_{}'.format(n)] = rel_thrs[n]
        return point

    def isEquivalent(self, other):
        if len(self.edges) != len(other.edges):
            return False
        return np.count_nonzero(self.edges == other.edges) == len(self.edges)
    def isBetter(self, other):
        if other is None or other.exp_limit is None:
            return True
        if self.exp_limit != other.exp_limit:
            return self.exp_limit < other.exp_limit
        return len(self.edges) < len(other.edges)


class BayesianOptimization:
    def __init__(self, max_n_bins, working_area, workers_dir, acq, kappa, xi, input_datacard, poi,
                 bkg_yields, input_queue_size, random_seed=12345, other_datacards=[]):
        self.max_n_bins = max_n_bins
        self.binnings = []
        self.best_binning = None
        self.working_area = working_area
        self.workers_dir = workers_dir
        self.log_output = os.path.join(working_area, 'results.json')
        self.input_datacard = input_datacard
        self.poi = poi
        self.bkg_yields = bkg_yields
        self.other_datacards = other_datacards
        self.best_binning_split = False
        self.input_queue = Queue(input_queue_size)
        self.open_requests = []
        self.binning_lock = threading.Lock()
        self.optimizer_lock = threading.Lock()
        self.print_lock = threading.Lock()

        bounds = {}
        for n in range(max_n_bins):
            upper_bound = 1. if n > 0 else 1.5
            bounds['rel_thr_{}'.format(n)] = (min_step, upper_bound)

        self.bounds_transformer = None
        self.optimizer = bayes_opt.BayesianOptimization(f=None, pbounds=bounds, random_state=random_seed, verbose=1)
        self.utilities = [
            # bayes_opt.util.UtilityFunction(kind='ucb', kappa=kappa, xi=xi),# kappa_decay=0.99), # supported in bayes_opt.__version__ <2.0.0
            # bayes_opt.util.UtilityFunction(kind='ei', kappa=kappa, xi=xi),
            # bayes_opt.util.UtilityFunction(kind='poi', kappa=kappa, xi=xi),
            bayes_opt.acquisition.UpperConfidenceBound(kappa=kappa),
            bayes_opt.acquisition.ExpectedImprovement(xi=xi),
            bayes_opt.acquisition.ProbabilityOfImprovement(xi=xi),
        ]

        if not os.path.isdir(working_area):
            os.mkdir(working_area)
        if not os.path.isdir(workers_dir):
            os.mkdir(workers_dir)
        if os.path.isfile(self.log_output):
            with open(self.log_output, 'r') as f:
                prev_binnings = json.loads('[' + ', '.join(f.readlines()) + ']')
            for entry in prev_binnings:
                binning = Binning.fromEntry(entry)
                if self.findEquivalent(binning) is None:
                    point = binning.toPoint(self.max_n_bins)
                    self.register(point, len(binning.edges), binning.exp_limit)
                    if binning.isBetter(self.best_binning):
                        self.best_binning = binning
                    self.binnings.append(binning)
        if self.best_binning is not None:
            self.print('The best binning from previous iterations: {}, exp_limit = {}' \
                       .format(arrayToStr(self.best_binning.edges), self.best_binning.exp_limit))

        self.bounds_transformer = bayes_opt.SequentialDomainReductionTransformer()
        self.bounds_transformer.initialize(self.optimizer.space)

        suggestions_file = os.path.join(working_area, 'to_try.json')
        self.suggestions = []
        if os.path.isfile(suggestions_file):
            with open(suggestions_file, 'r') as f:
                suggestions = json.load(f)
            for edges in suggestions:
                self.addSuggestion(edges)

    def print(self, msg):
        self.print_lock.acquire()
        print(msg)
        self.print_lock.release()

    def updateBounds(self, n_points):
        if self.bounds_transformer:
            new_bounds = self.bounds_transformer.transform(self.optimizer.space)
            prev_fixed = True
            for n in range(self.max_n_bins):
                key = 'rel_thr_{}'.format(n)
                if prev_fixed:
                    prev_fixed = new_bounds[key][1] - new_bounds[key][0] < 0.1 and n < n_points
                else:
                    new_bounds[key] = np.array([min_step, 1.])
            for key in new_bounds:
                if new_bounds[key][0] >= new_bounds[key][1]:
                    new_bounds[key] = np.array([new_bounds[key][1], new_bounds[key][0]])
                if new_bounds[key][0] < min_step:
                    new_bounds[key][0] = min_step
                upper_bound = 1. if key != 'rel_thr_0' else 1.5
                if new_bounds[key][1] > upper_bound:
                    new_bounds[key][1] = upper_bound
                if np.any(np.isnan(new_bounds[key])):
                    new_bounds[key] = np.array([min_step, upper_bound])
            self.optimizer.set_bounds(new_bounds)

    def register(self, point, n_edges, exp_limit):
        loss = -(exp_limit*1e6 + n_edges)
        self.optimizer_lock.acquire()
        self.optimizer.register(params=point, target=loss)
        self.updateBounds(n_edges - 2)
        self.optimizer_lock.release()

    def suggest(self, utility_index):
        self.optimizer_lock.acquire()
        if self.suggestions is None or len(self.suggestions) == 0:
            # point = self.optimizer.suggest(self.utilities[utility_index]) # bayes_opt.__version__ >=2.0.0 does not support passing utility functions
            point = self.optimizer.suggest()
        else:
            point = self.suggestions[0]
            self.suggestions.remove(point)
        self.optimizer_lock.release()
        return point

    def addSuggestion(self, edges):
        suggested_binning = Binning(np.array(edges))
        rel_thrs = suggested_binning.getRelativeThresholds(self.max_n_bins)
        binning = Binning.fromRelativeThresholds(rel_thrs, self.bkg_yields)
        if self.findEquivalent(binning, False) is None:
            point = binning.toPoint(self.max_n_bins)
            self.suggestions.append(point)

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
            print(f'PointGenerator: n {n}, n_eq_steps {n_eq_steps}')
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
                self.print('Equivalent binnig found: {}'.format(arrayToStr(equivalent_binning.edges)))
                self.register(point, len(equivalent_binning.edges), equivalent_binning.exp_limit)
                n += 1
                if n >= 5:
                    self.tryBestBinningSplit()
                if n == n_eq_steps - 1:
                    self.print('Waiting for open requests to finish...')
                    self.waitOpenRequestsToFinish()
        print('PointGenerator')
        self.input_queue.put(None)

    def JobDispatcher(self):
        while True:
            print('JobDispatcher')
            time.sleep(1)
            print(f'self.workers_dir {self.workers_dir}')
            for worker_dir in os.listdir(self.workers_dir):
                worker_dir = os.path.join(self.workers_dir, worker_dir)
                # print(f'os.path.isdir(worker_dir) {os.path.isdir(worker_dir)} {worker_dir}')
                if not os.path.isdir(worker_dir): continue
                task_file = os.path.join(worker_dir, 'task.txt')
                task_file_tmp = os.path.join(worker_dir, '.task.txt')
                print(f'os.path.isfile(task_file) {os.path.isfile(task_file)}')
                if os.path.isfile(task_file):
                    result_file = os.path.join(worker_dir, 'result.txt')
                    print(f'os.path.isfile(result_file) {os.path.isfile(result_file)}')
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
                        print('task.txt doesnt exist')
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
                            print(f'Task written to {task_file_tmp}')
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

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Find optimal binning that minimises the expected limits.')
    parser.add_argument('--input', required=True, type=str, help="input datacard")
    parser.add_argument('--output', required=True, type=str, help="output directory")
    parser.add_argument('--workers-dir', required=True, type=str, help="output directory for workers results")
    parser.add_argument('--max_n_bins', required=True, type=int, help="maximum number of bins")
    parser.add_argument('--poi', required=False, type=str, default='r', help="parameter of interest")
    parser.add_argument('--params', required=False, type=str, default=None,
                        help="algorithm parameters in format param1=value1,param2=value2 ...")
    parser.add_argument('--verbosity', required=False, type=int, default=1, help="verbosity")
    parser.add_argument('other_datacards', type=str, nargs='*',
                        help="list of other datacards to be combined together with the current target")

    args = parser.parse_args()

    input_datacard = os.path.abspath(args.input)
    input_name = os.path.splitext(os.path.basename(input_datacard))[0]
    input_shapes = os.path.splitext(input_datacard)[0] + '.input.root'

    other_datacards = [ os.path.abspath(p) for p in args.other_datacards ]

    param_dict = {}
    if args.params is not None:
        params = args.params.split(',')
        for param in params:
            p = param.split('=')
            if len(p) != 2:
                raise RuntimeError('invalid parameter definition "{}"'.format(param))
            param_dict[p[0]] = p[1]

    ref_bkgs = {
        'DY': re.compile('^DY$'),
        'TT': re.compile('^TT$'),
    }

    nonbkg_regex = re.compile('(data_obs|^ggHH.*|^qqHH.*|^DY_[0-2]b.*)')
    ignore_unc_variations = re.compile('(CMS_bbtt_201[6-8]_DYSFunc[0-9]+|CMS_bbtt_.*_QCDshape)(Up|Down)')
    print("Extracting yields for background processes...")
    bkg_yields = ExtractYields(input_shapes, ref_bkgs, nonbkg_regex, ignore_unc_variations)
    bkg_yields.printSummary()

    for name,value in param_dict.items():
        if not hasattr(bkg_yields, name):
            raise RuntimeError('Unknown parameter "{}"'.format(name))
        def_value = getattr(bkg_yields, name)
        def_value_type = type(def_value)
        if def_value_type == bool:
            new_value = bool(distutils.util.strtobool(value))
        else:
            new_value = def_value_type(value)
        print("Setting {} = {}".format(name, new_value))
        setattr(bkg_yields, name, new_value)

    bo = BayesianOptimization(max_n_bins=args.max_n_bins,
                              working_area=args.output,
                              workers_dir=args.workers_dir,
                              acq='ucb', kappa=2.57, xi=0.1,
                              input_datacard=input_datacard,
                              poi=args.poi,
                              bkg_yields=bkg_yields,
                              input_queue_size=2, random_seed=None,
                              other_datacards=other_datacards)
    print("here")
    bo.maximize(20)
    print("here here")

    print('Minimization finished.')
    print('Best binning: {}'.format(arrayToStr(bo.best_binning.edges)))
    print('95% CL expected limit: {}'.format(bo.best_binning.exp_limit))
