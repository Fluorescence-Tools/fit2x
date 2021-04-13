# Scoring
######################
@property
def weighted_residuals(self):
    return self.get_weighted_residuals()

@property
def score_range(self):
    return self.get_score_range()

@score_range.setter
def score_range(self, v):
    self.set_score_range(v)

@property
def chi2(self):
    return self.evaluate()

@property
def n_scoring_channels(self):
    """Number of channels used for scoring
    """
    score_range = self.score_range
    return score_range[1] - score_range[0]

# @property
# def parameter(self) -> dict:
#     re = {
#         'irf_background_counts': self.irf_background_counts,
#         'irf_shift_channels': self.irf_shift_channels,
#         'scatter_fraction': self.scatter_fraction,
#         'constant_offset': self.constant_offset,
#         'lifetime_spectrum': np.copy(self.lifetime_spectrum), # make a copy as this is otherwise a view
#         'number_of_photons': self.number_of_photons,
#         'acquisition_time': self.acquisition_time,
#         'instrument_dead_time': self.instrument_dead_time,
#         'convolution_range': self.convolution_range,
#         'excitation_period': self.excitation_period,
#         'scale_model_to_data': self.scale_model_to_data,
#         'score_range': self.score_range,
#         'use_corrected_irf_as_scatter': self.use_corrected_irf_as_scatter,
#         'amplitude_threshold': self.amplitude_threshold,
#         'use_amplitude_threshold': self.use_amplitude_threshold,
#         'use_pile_up_correction': self.use_pile_up_correction,
#         'convolution_method': self.convolution_method,
#         'use_linearization': self.use_linearization
#     }
#     return re

def __repr__(self):
    imin, imax = self.score_range[0], self.score_range[1]
    used_channels = abs(imax - imin)
    s = 'DECAY SCORE \n'
    s += '-- Score range: (%d, %d) \n' % (imin, imax)
    s += '-- Score: %s \n' % self.chi2
    s += '-- Score / n_channels: %s \n' % (self.get_score() / used_channels)
    return s
