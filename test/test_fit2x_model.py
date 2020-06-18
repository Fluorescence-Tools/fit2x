from __future__ import division

import unittest
import numpy as np
import scipy.stats

import fit2x


def model_irf(
        n_channels: int = 256,
        period: float = 32,
        irf_position_p: float = 2.0,
        irf_position_s: float = 18.0,
        irf_width: float = 0.25
):
    time_axis = np.linspace(0, period, n_channels * 2)
    irf_np = scipy.stats.norm.pdf(time_axis, loc=irf_position_p,
                                  scale=irf_width) + \
             scipy.stats.norm.pdf(time_axis, loc=irf_position_s,
                                  scale=irf_width)
    return irf_np, time_axis


class Tests(unittest.TestCase):

    def test_LVDoubleArray(self):
        lv_array = fit2x.CreateLVDoubleArray(10)
        self.assertEqual(
            len(lv_array), 10
        )
        lv_array[0] = 10
        vec = np.linspace(1, 10., 10, dtype=np.float)
        for i, v in enumerate(vec):
            lv_array[i] = v
        l = [v for v in lv_array]
        self.assertListEqual(
            list(lv_array), l
        )

    @unittest.expectedFailure
    def test_LVDoubleArray_slice(self):
        lv_array = fit2x.CreateLVDoubleArray(10)
        vec = np.linspace(1, 10., 10, dtype=np.float)
        # will fail slicing not supported
        lv_array[:] = vec

    def test_LVI32Array(self):
        lv_array = fit2x.CreateLVI32Array(10)
        self.assertEqual(
            len(lv_array), 10
        )
        lv_array[0] = 10
        vec = np.linspace(1, 10, 10, dtype=np.int)
        # will fail slicing not supported
        # lv_array[:] = vec
        for i, v in enumerate(vec):
            lv_array[i] = int(v)
        l = [v for v in lv_array]
        self.assertListEqual(
            list(lv_array), l
        )

    @unittest.expectedFailure
    def test_LVI32Array_slice(self):
        lv_array = fit2x.CreateLVI32Array(10)
        vec = np.linspace(1, 10, 10, dtype=np.int)
        lv_array[:] = vec

    def test_MParam_2(self):
        # create filled MParam structure
        n_corrections = 5
        corrections_np = np.zeros(n_corrections, dtype=np.float)
        irf_np, timeaxis = model_irf()
        dt = timeaxis[1] - timeaxis[0]
        bg_np = np.zeros_like(irf_np)
        parameter_group = fit2x.CreateMParam(
            irf=irf_np,
            background=bg_np,
            corrections=corrections_np,
            dt=dt,
        )
        irf = parameter_group.get_irf()
        bg = parameter_group.get_background()
        corrections = parameter_group.get_corrections()
        self.assertEqual(len(irf_np), len(irf))
        self.assertEqual(len(bg_np), len(bg))
        irf[0] = 123
        self.assertEqual(irf[0], 123)
        bg[0] = 123
        self.assertEqual(bg[0], 123)
        corrections[0] = 11
        self.assertEqual(corrections[0], 11)

    def test_lv_param(self):
        for i in range(200):
            lv_array = fit2x.CreateLVDoubleArray(10)
            a = np.ones(111, dtype=np.float)
            lv_array.set_data(a)
            x = [x for x in lv_array]
            self.assertListEqual(list(a), x)

        for i in range(200):
            lv_array = fit2x.CreateLVI32Array(10)
            a = np.ones(111, dtype=np.int32)
            lv_array.set_data(a)
            x = [x for x in lv_array]
            self.assertListEqual(list(a), x)

        for i in range(1000):
            dt = 0.032
            a = np.ones(111, dtype=np.float)
            parameter_group = fit2x.CreateMParam(
                irf=a,
                background=a,
                corrections=a,
                dt=dt
            )
            irf = parameter_group.get_irf()
            irf.set_data(a)
            x = [x for x in irf]
            self.assertListEqual(list(a), x)

            bg = parameter_group.get_background()
            bg.set_data(a)
            x = [x for x in bg]
            self.assertListEqual(list(a), x)

            corrections = parameter_group.get_corrections()
            corrections.set_data(a)
            x = [x for x in corrections]
            self.assertListEqual(list(a), x)
