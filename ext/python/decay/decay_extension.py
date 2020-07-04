
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

def __repr__(self):
    s = 'DECAY \n'
    s += '-- Lifetime spectrum: %s \n' % self.lifetime_spectrum
    s += '-- Constant background / counts: %s \n' % self.background
    s += '-- Areal scatter fraction: %s \n' % self.scatter_fraction
    s += '-- Irf shift / channels: %s \n' % self.irf_shift
    s += '-- Irf background / counts: %s \n' % self.irf_background
    s += '-- Chi2: %s \n' % self.chi2
    return s
