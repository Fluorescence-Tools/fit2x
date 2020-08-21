
@property
def model(self):
    return self.get_model()


@property
def data(self):
    return self.get_data()


@data.setter
def data(self, v):
    self.set_data(v)


@property
def irf(self):
    return self.get_irf()


@irf.setter
def irf(self, v):
    self.set_irf(v)


@property
def corrected_irf(self):
    return self.get_corrected_irf()


@property
def weighted_residuals(self):
    return self.get_weighted_residuals()


@property
def time_axis(self):
    return self.get_time_axis()


@time_axis.setter
def time_axis(self, v):
    return self.set_time_axis(v)


@property
def lifetime_spectrum(self):
    return self.get_lifetime_spectrum()


@lifetime_spectrum.setter
def lifetime_spectrum(self, v):
    self.set_lifetime_spectrum(v)


@property
def chi2(self):
    return self.get_score()


@property
def score_range(self):
    return self.get_score_range()


@score_range.setter
def score_range(self, v):
    self.set_score_range(v[0], v[1])


@property
def linearization(self):
    return self.get_linearization()


@linearization.setter
def linearization(self, v):
    self.set_linearization(v)


@property
def mean_lifetime(self):
    return self.get_mean_lifetime()


@property
def n_scoring_channels(self):
    """Number of channels used for scoring
    """
    score_range = self.score_range
    return score_range[1] - score_range[0]


@property
def convolution_range(self):
    return self.get_convolution_range()


@convolution_range.setter
def convolution_range(self, v):
    self.set_convolution_range(v)


@property
def parameter(self) -> dict:
    re = {
        'irf_background_counts': self.irf_background_counts,
        'irf_shift_channels': self.irf_shift_channels,
        'scatter_fraction': self.scatter_fraction,
        'constant_offset': self.constant_offset,
        'lifetime_spectrum': np.copy(self.lifetime_spectrum), # make a copy as this is otherwise a view
        'number_of_photons': self.number_of_photons,
        'acquisition_time': self.acquisition_time,
        'instrument_dead_time': self.instrument_dead_time,
        'convolution_range': self.convolution_range,
        'excitation_period': self.excitation_period,
        'scale_model_to_data': self.scale_model_to_data,
        'score_range': self.score_range,
        'use_corrected_irf_as_scatter': self.use_corrected_irf_as_scatter,
        'amplitude_threshold': self.amplitude_threshold,
        'use_amplitude_threshold': self.use_amplitude_threshold,
        'use_pile_up_correction': self.add_pile_up,
        'convolution_method': self.convolution_method,
        'use_linearization': self.use_linearization
    }
    return re


def __repr__(self):
    imin, imax = self.score_range[0], self.score_range[1]
    n_lifetime_parameter = len(self.lifetime_spectrum)
    used_channels = abs(imax - imin)
    nu = used_channels - n_lifetime_parameter
    s = 'DECAY \n'
    s += '-- Lifetime spectrum: %s \n' % self.lifetime_spectrum
    s += '-- n_lifetime_parameter: %s \n' % n_lifetime_parameter
    s += '-- Constant background / counts: %s \n' % self.constant_offset
    s += '-- Areal scatter fraction: %s \n' % self.scatter_fraction
    s += '-- Irf shift / channels: %s \n' % self.irf_shift_channels
    s += '-- Irf background / counts: %s \n' % self.irf_background_counts
    s += '-- Score range: (%d, %d) \n' % (imin, imax)
    s += '-- Score: %s \n' % self.chi2
    s += '-- Score / n_channels: %s \n' % (self.get_score() / used_channels)
    s += '-- Chi2 / (n_channels - n_lifetime_parameter): %s ' % (
            self.get_score(score_type='normal') / nu
    )
    return s
