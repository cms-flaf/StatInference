#!/usr/bin/env python3
"""rateParam uncertainties (dc_make.uncertainty.RateParamUncertainty).

* fromConfig: `value` defaults to 1 and `range` to [0, 10], both read as numbers; a range that
  is not [low, high] with low < high around the value, a name with a CombineHarvester
  placeholder, or a key a rateParam does not take raises.
* configureParameters sets the range on the card's parameter of that name, and does nothing
  when the uncertainty matched no bin.
* With CombineHarvester available: one rateParam added to two bins is one parameter, written
  once per bin with its range.

CombineHarvester is stubbed where it is not installed, so the parsing tests run anywhere.
"""

import importlib
import os
import sys
import tempfile
import types
import unittest
from unittest import mock

STATINFERENCE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(STATINFERENCE))


def _has_combine_harvester():
    try:
        importlib.import_module("CombineHarvester.CombineTools.ch")
        return True
    except Exception:
        return False


HAS_CH = _has_combine_harvester()


class _SystMap:
    def __call__(self, *values):
        return values


def _stub_modules():
    ch = types.ModuleType("CombineHarvester.CombineTools.ch")
    ch.SystMap = _SystMap
    tools = types.ModuleType("CombineHarvester.CombineTools")
    tools.ch = ch
    root = types.ModuleType("CombineHarvester")
    root.CombineTools = tools
    return {
        "CombineHarvester": root,
        "CombineHarvester.CombineTools": tools,
        "CombineHarvester.CombineTools.ch": ch,
    }


def _import_uncertainty():
    name = os.path.basename(STATINFERENCE) + ".dc_make.uncertainty"
    if HAS_CH:
        return importlib.import_module(name)
    with mock.patch.dict(sys.modules, _stub_modules()):
        sys.modules.pop(name, None)
        return importlib.import_module(name)


unc_mod = _import_uncertainty()
Uncertainty = unc_mod.Uncertainty
RateParamUncertainty = unc_mod.RateParamUncertainty


def entry(**kwargs):
    e = {"name": "CMS_rate_TT", "type": "rateParam", "processes": ["TT"]}
    e.update(kwargs)
    return e


class FromConfigTest(unittest.TestCase):
    def test_defaults(self):
        unc = Uncertainty.fromConfig(entry())
        self.assertIsInstance(unc, RateParamUncertainty)
        self.assertEqual(unc.value, 1.0)
        self.assertEqual(unc.value_range, [0.0, 10.0])
        self.assertFalse(unc.needShapes)
        self.assertFalse(unc.canIgnore(1e9))

    def test_numbers(self):
        # PyYAML reads 1e1 (no dot) as a string; it is still a number here
        unc = Uncertainty.fromConfig(entry(value=2, range=[0, "1e1"]))
        self.assertEqual((unc.value, unc.value_range), (2.0, [0.0, 10.0]))
        self.assertIsInstance(unc.value, float)

    def test_bad_ranges(self):
        for bad in ([], [5.0], [0.0, 1.0, 2.0], 5, [0.0, "x"], [2.0, 0.0], [2.0, 5.0]):
            with self.subTest(range=bad), self.assertRaises(RuntimeError):
                Uncertainty.fromConfig(entry(range=bad))

    def test_bad_value(self):
        for bad in ([0.9, 1.1], "one", 20.0):
            with self.subTest(value=bad), self.assertRaises(RuntimeError):
                Uncertainty.fromConfig(entry(value=bad))

    def test_placeholder_name(self):
        with self.assertRaises(RuntimeError):
            Uncertainty.fromConfig(entry(name="CMS_rate_TT_$BIN"))

    def test_unsupported_key(self):
        with self.assertRaises(RuntimeError):
            Uncertainty.fromConfig(entry(hasMultipleValues=True))


class ConfigureParametersTest(unittest.TestCase):
    def test_sets_range(self):
        param = mock.Mock()
        cb = mock.Mock()
        cb.GetParameter.return_value = param
        Uncertainty.fromConfig(entry(range=[0, 5])).configureParameters(
            cb, "CMS_rate_TT"
        )
        cb.GetParameter.assert_called_once_with("CMS_rate_TT")
        param.set_range.assert_called_once_with(0.0, 5.0)

    def test_absent_parameter(self):
        cb = mock.Mock()
        cb.GetParameter.return_value = None
        Uncertainty.fromConfig(entry()).configureParameters(cb, "CMS_rate_TT")


@unittest.skipUnless(HAS_CH, "CombineHarvester not available")
class CardTest(unittest.TestCase):
    def test_written_once_per_bin_with_range(self):
        import CombineHarvester.CombineTools.ch as ch

        cb = ch.CombineHarvester()
        cats = [(0, "SR"), (1, "CR"), (2, "other")]
        cb.AddObservations(["*"], ["ana"], ["era"], ["chn"], cats)
        cb.AddProcesses(["*"], ["ana"], ["era"], ["chn"], ["TT"], cats, False)
        cb.AddProcesses(["*"], ["ana"], ["era"], ["chn"], ["sig"], cats, True)
        unc = Uncertainty.fromConfig(entry(range=[0, 5]))
        # a scope short of every bin, as a control region and its signal region (all bins
        # would be written as one "*" line)
        cb.cp().process(["TT"]).bin(["SR", "CR"]).AddSyst(
            cb, unc.name, "rateParam", unc.valueToMap()
        )
        unc.configureParameters(cb, unc.name)
        cb.ForEachProc(lambda p: p.set_rate(1.0))
        cb.ForEachObs(lambda o: o.set_rate(1.0))
        with tempfile.TemporaryDirectory() as tmp:
            card = os.path.join(tmp, "card.txt")
            cb.WriteDatacard(card)
            with open(card) as f:
                lines = [l.split() for l in f if "rateParam" in l]
        self.assertEqual(len(lines), 2)
        for parts in lines:
            self.assertEqual(parts[:2], ["CMS_rate_TT", "rateParam"])
            self.assertEqual(parts[3], "TT")
            self.assertEqual(parts[-1], "[0,5]")
        self.assertEqual({parts[2] for parts in lines}, {"SR", "CR"})


if __name__ == "__main__":
    unittest.main()
