
    # A list of all accessible parameters
    _parameter_keys = [
            # Data
            'data',
            'acquisition_time',
            'data_weights',
            'time_axis',
            'irf_histogram',
            # Lifetime spectrum
            'lifetime_spectrum',
            'use_amplitude_threshold',
            'abs_lifetime_spectrum',
            'amplitude_threshold',
            # Convolution
            'convolution_range',
            'use_corrected_irf_as_scatter',
            'scatter_fraction',
            'convolution_method',
            'excitation_period',
            'irf_shift_channels',
            'irf_background_counts',
            # Background
            'constant_offset',
            # Pile up
            'pile_up_model',
            'instrument_dead_time',
            'use_pile_up_correction',
            # Scaling
            'scale_model_to_data',
            'number_of_photons',
            # Linearization
            'linearization_table',
            'use_linearization',
            # Scoring
            'score_range',
            'score_type'
    ]

    # Score
    #######################
    @property
    def score_range(self):
        return self.get_score_range()

    @score_range.setter
    def score_range(self, v):
        self.set_score_range(v[0], v[1])

    # Lifetime spectrum
    #######################
    @property
    def lifetime_spectrum(self):
        return self.get_lifetime_spectrum()

    @lifetime_spectrum.setter
    def lifetime_spectrum(self, v):
        self.set_lifetime_spectrum(v)

    # CONVOLUTION
    #######################
    @property
    def convolution_range(self):
        return self.get_convolution_range()

    @convolution_range.setter
    def convolution_range(self, v):
        self.set_convolution_range(v)

    @property
    def irf(self):
        return self.get_irf()

    @irf.setter
    def irf(self, v):
        self.set_irf(v)

    @property
    def irf_histogram(self):
        return self.get_irf()

    # For backwards compatability
    @irf_histogram.setter
    def irf_histogram(self, v):
        self.set_irf(v)

    @property
    def corrected_irf(self):
        return self.get_corrected_irf()

    # LINEARIZATION
    #######################
    @property
    def linearization_table(self):
        return self.get_linearization_table()

    @linearization_table.setter
    def linearization_table(self, v):
        self.set_linearization_table(v)

    # DATA
    #######################
    @property
    def data(self):
        return self.get_data()

    @data.setter
    def data(self, v):
        self.set_data(v)

    @property
    def data_weights(self):
        return self.get_data_weights()

    @data_weights.setter
    def data_weights(self, v):
        return self.set_data_weights(v)

    @property
    def time_axis(self):
        return self.get_time_axis()

    @time_axis.setter
    def time_axis(self, v):
        return self.set_time_axis(v)

    # Model
    ######################
    @property
    def model(self):
        return self.get_model()

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
        return self.get_score()

    @property
    def n_scoring_channels(self):
        """Number of channels used for scoring
        """
        score_range = self.score_range
        return score_range[1] - score_range[0]

    def set(self, **kwargs):
        for k in kwargs:
            if k in Decay._parameter_keys:
                v = kwargs[k]
                if v is not None:
                    self.__setattr__(k, kwargs[k])
            else:
                raise AttributeError("The parameter '%s' cannot be set." % k)

    @property
    def parameter(self):
        return dict(
            zip(
                self._parameter_keys,
                [self.__getattribute__(k) for k in Decay._parameter_keys]
            )
        )

    @parameter.setter
    def parameter(self, v):
        self.set(v)

    @property
    def mean_lifetime(self):
        return self.get_mean_lifetime()


