from __future__ import division
import unittest

import collections
import json
import pathlib
import tempfile
import random
import string
import numbers

import numpy as np
import scipy.stats
import fit2x


x = np.linspace(0, 20, 32)
irf_position = 2.0
irf_width = 0.15
irf_y = scipy.stats.norm.pdf(x, loc=irf_position, scale=irf_width)


def fill_random(d: fit2x.Decay):
    letters = string.ascii_letters
    for k, v in d.parameter.items():
        print(k, v)
        if isinstance(v, str):
            d.set(k, ''.join(random.choice(letters) for _ in range(4)))
        elif isinstance(v, np.ndarray):
            d.set(k, np.random.random(4))
        elif isinstance(v, bool):
            d.set(k, bool(random.getrandbits(1)))
        elif isinstance(v, numbers.Number):
            d.set(k, random.random())


def compare_parameter(parameter, other):
    r = list()
    for k, v in parameter.items():
        if isinstance(v, collections.Iterable) and not isinstance(v, str):
            r.append(np.allclose(v, other[k]))
        else:
            r.append(v == other[k])
    return np.alltrue(r)


class Tests(unittest.TestCase):

    def test_parameter_json(self):
        data_settings = fit2x.DataSettings(
            time_axis=x,
            data=irf_y,
            irf_histogram=irf_y
        )
        lifetime_settings = fit2x.LifetimeSettings(
            lifetime_spectrum=[1, 4]
        )
        convolution_settings = fit2x.ConvolutionSettings()
        scatter_settings = fit2x.PatternSettings()
        background_settings = fit2x.PatternSettings()
        pileup_settings = fit2x.PileupSettings()
        linearization_settings = fit2x.LinearizationSettings()
        scaling_settings = fit2x.ScalingSettings()
        score_settings = fit2x.ScoreSettings()

        dc = fit2x.Decay(
            data_settings=data_settings,
            lifetime_settings=lifetime_settings,
            scatter_settings=scatter_settings,
            convolution_settings=convolution_settings,
            background_settings=background_settings,
            pileup_settings=pileup_settings,
            linearization_settings=linearization_settings,
            scaling_settings=scaling_settings,
            score_settings=score_settings
        )
        dc.update()

        # Test filling, reading and writing. Fill with random data
        fill_random(dc)
        # JSON parameter comparison
        json_str = dc.to_json()
        p_json = json.loads(json_str)
        p_obj = dc.parameter
        self.assertEqual(compare_parameter(p_json, p_obj), True)

        # Write to JSON file and restore Decay object from JSON
        fh, fn = tempfile.mkstemp(".json")
        fp = open(fh, "w")
        fp.write(json_str)
        fp.close()

        decay_from_json_str = fit2x.Decay.from_json(json_str)
        decay_from_json_fn = fit2x.Decay.from_json(pathlib.Path(fn))

        self.assertEqual(compare_parameter(decay_from_json_str.parameter, p_obj), True)
        self.assertEqual(compare_parameter(decay_from_json_fn.parameter, p_obj), True)

    def test_file_io(self):
        dc = fit2x.Decay()
        dc.update()
        fill_random(dc)
        # JSON parameter comparison
        json_str = dc.to_json()
        # Write to JSON file and restore Decay object from JSON
        fh, fn = tempfile.mkstemp(".json")
        fp = open(fh, "w")
        fp.write(json_str)
        fp.close()

        decay_from_json_str = fit2x.Decay.from_json(json_str)
        decay_from_json_fn = fit2x.Decay.from_json(pathlib.Path(fn))

        self.assertEqual(compare_parameter(decay_from_json_str.parameter, dc.parameter), True)
        self.assertEqual(compare_parameter(decay_from_json_fn.parameter, dc.parameter), True)

    def test_decay_convolution_lifetime(self):
        np.random.seed(42)

        x = np.linspace(0, 20, 32)
        irf_position = 2.0
        irf_width = 0.25
        irf_y = scipy.stats.norm.pdf(x, loc=irf_position, scale=irf_width)

        data_settings = fit2x.DataSettings(
            time_axis=x,
            data=irf_y * 100000,
            irf_histogram=irf_y
        )
        lifetime_settings = fit2x.LifetimeSettings(
            lifetime_spectrum=[1, 2]
        )
        convolution_settings = fit2x.ConvolutionSettings()
        background_settings = fit2x.PatternSettings()
        pileup_settings = fit2x.PileupSettings()
        linearization_settings = fit2x.LinearizationSettings()
        scaling_settings = fit2x.ScalingSettings()
        score_settings = fit2x.ScoreSettings()

        dc = fit2x.Decay(
            data_settings=data_settings,
            lifetime_settings=lifetime_settings,
            convolution_settings=convolution_settings,
            background_settings=background_settings,
            pileup_settings=pileup_settings,
            linearization_settings=linearization_settings,
            scaling_settings=scaling_settings,
            score_settings=score_settings
        )
        dc.update()
        dc.data = np.random.poisson(1000 / max(dc.y) * dc.y)
        scores = list()
        lts = np.linspace(0.0, 4, 11)
        for lt in lts:
            dc.lifetime_spectrum = [1, lt]
            scores.append(dc.score)

        scores_ref = [3641.4913153320936,
                      3023.436739258023,
                      1648.7875155218765,
                      731.3338936915842,
                      190.31400668610308,
                      23.229701466698423,
                      377.664429861769,
                      1529.1599392280614,
                      3840.3130365516517,
                      7698.247948184009,
                      13453.340913686912]
        self.assertEqual(np.allclose(scores, scores_ref), True)
    #
    # def test_decay_convolution_score(self):
    #     for score_type in ['chi2', 'poisson', 'gauss']:
    #         scores = list()
    #         dc.score_type = score_type
    #         lts = np.linspace(0.0, 4, 11)
    #         for lt in lts:
    #             dc.lifetime_spectrum = [1, lt]
    #             scores.append(dc.score)
    #         plt.plot(lts, scores, label=score_type)
    #     plt.legend()
    #     plt.show()
    #
    #     # scatter
    #     dc.lifetime_spectrum = [1, 4]
    #     dc.scatter_fraction = 0.8
    #     dc.update()
    #     plt.semilogy(dc.x, dc.y, label="sc 0.8")
    #     dc.scatter_fraction = 0.2
    #     dc.update()
    #     plt.semilogy(dc.x, dc.y, label="sc 0.2")
    #     plt.legend()
    #     plt.show()
    #
    #     # data noise
    #     dc.lifetime_spectrum = [1, 4]
    #     dc.scatter_fraction = 0.0
    #     dc.update()
    #     plt.semilogy(dc.x, dc.data, label="data")
    #     plt.semilogy(dc.x, dc.data_noise, label="noise")
    #     plt.legend()
    #     plt.show()
    #
    #     scores = list()
    #     scatter = np.linspace(0, 1, 10)
    #     for sc in scatter:
    #         dc.scatter_fraction = sc
    #         scores.append(dc.score)
    #     plt.plot(scatter, scores)
    #     plt.legend()
    #     plt.show()
    #
