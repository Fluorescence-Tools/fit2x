import numpy as np
import fit2x


class Fit23(object):

    def __init__(
            self,
            dt: float,
            irf: np.ndarray,
            background: np.ndarray,
            period: float,
            g_factor: float = 1.0,
            l1: float = 0.0,
            l2: float = 0.0,
            convolution_stop: int = -1,
            soft_bifl_scatter_flag: bool = True,
            p2s_twoIstar_flag: bool = False
    ):
        """
        :param dt: time difference between microtime bins
        :param irf: counting histogram of the instrument response function in
        Jordi format
        :param background: background counting histogram in Jordi format
        :param period: excitation period of the light source
        :param g_factor: g-factor
        :param l1: factor correcting mixing between parallel and perpendicular
        detection channel
        :param l2: factor correcting mixing between parallel and perpendicular
        detection channel
        :param convolution_stop: maximum micro time channel for comvolution. If
        no value is provided the length of the IRF is used
        :param soft_bifl_scatter_flag: if set to True the returned Istar value
        is reduced by the background photon contribution (background photons do
        no inform on the fluorescence lifetime)
        :param p2s_twoIstar_flag: If this is set to True the sum decay composed
        by P + 2S (P - Parallel, S - Perpendicular) is optimized. Otherwise
        (default) the decays P and S are optimized individually in a global fit.
        """
        if len(irf) != len(background):
            raise ValueError("The IRF and the background differ in size")
        if len(irf) % 2 != 0:
            raise ValueError("The length of the input arrays is not divisible "
                             "by two. Inputs need to be in Jordi format.")
        if convolution_stop < 0:
            convolution_stop = len(irf) // 2 - 1
        self._bifl_scatter = -1 if soft_bifl_scatter_flag else 0
        self._p_2s_flag = p2s_twoIstar_flag
        self._corrections = np.array(
            [period, g_factor, l1, l2, convolution_stop]
        )
        self._irf = irf
        self._background = background
        self._dt = dt
        self._m_param = fit2x.CreateMParam(
            irf=self._irf,
            background=self._background,
            corrections=self._corrections,
            dt=self._dt
        )

    @property
    def model(self):
        return np.array([x for x in self._m_param.get_model()])

    @property
    def data(self):
        return np.array([x for x in self._m_param.get_data()])

    def __call__(
            self,
            data: np.ndarray,
            initial_values: np.ndarray,
            fixed: np.ndarray = None,
            include_model: bool = False
    ) -> dict:
        """

        :param data: counting histogram containing experimental data
        :param initial_values: initial values of the model parameters that can
        be optimized. [tau, gamma, r0, rho]
        :param fixed: optional array of short (16bit) integers that specifies if
        a parameter is fixed. Parameters that are fixed are not optimized.
        :param include_model: if set to True (default is False) the realization
        of the model that corresponds to the optimized parameters is included in
        the returned dictionary.
        :return: dictionary containing a quality parameter (key: "Istar"), the
        corresponding optimized model parameter values (key: "x"), and an array
        which parameters were fixed (key: "fixed").
        """
        if len(initial_values) < 4:
            raise ValueError(
                "Provide initial values for all for all 4 fitting "
                "parameters."
            )
        if fixed is None:
            # lifetime free
            fixed = np.array([0, 1, 1, 1], dtype=np.int16)
        elif isinstance(fixed, np.ndarray):
            if len(fixed) < 4:
                raise ValueError(
                    "The fixed array is too short. Specify for all 4 fitting "
                    "parameters if they are fixed."
                )
        else:
            raise ValueError(
                "The fixed array is of the wrong type. Use an numpy array of "
                "length 4 to specify the fixed state for all 4 model "
                "parameters."
            )
        r = dict()
        d = self._m_param.get_data()
        d.set_data(data.astype(dtype=np.int32))
        x = np.zeros(8, dtype=np.float64)
        x[:4] = initial_values
        x[4] = self._bifl_scatter
        x[5] = self._p_2s_flag
        # the other x values are used as outputs
        fixed = fixed.astype(dtype=np.int16)
        Istar = fit2x.fit23(x, fixed, self._m_param)
        r['x'] = x
        r['fixed'] = fixed
        r['Istar'] = Istar
        if include_model:
            r['model'] = self.model
        return r

