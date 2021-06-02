
import numpy as np
import scipy.stats
import fit2x

x = np.linspace(0, 20, 32)
irf_position = 2.0
irf_width = 0.15
irf_y = scipy.stats.norm.pdf(x, loc=irf_position, scale=irf_width)

data_settings = fit2x.DataSettings(
    time_axis=x,
    data=irf_y * 100000,
    irf_histogram=irf_y
)
lifetime_settings = fit2x.LifetimeSettings(
    lifetime_spectrum=[1, 4]
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

# lifetime
import pylab as plt
dc.lifetime_spectrum = [1, 4]
plt.semilogy(dc.x, dc.y, label="lt 4")
dc.lifetime_spectrum = [1, 2]
dc.update()
plt.semilogy(dc.x, dc.y, label="lt 2")
plt.legend()
plt.show()

# scatter
dc.lifetime_spectrum = [1, 4]
dc.scatter_fraction = 0.8
dc.update()
plt.semilogy(dc.x, dc.y, label="sc 0.8")
dc.scatter_fraction = 0.2
dc.update()
plt.semilogy(dc.x, dc.y, label="sc 0.2")
plt.legend()
plt.show()

# data noise
dc.lifetime_spectrum = [1, 4]
dc.scatter_fraction = 0.0
dc.update()
plt.semilogy(dc.x, dc.data, label="data")
plt.semilogy(dc.x, dc.data_noise, label="noise")
plt.legend()
plt.show()

scores = list()
scatter = np.linspace(0, 1, 10)
for sc in scatter:
    dc.scatter_fraction = sc
    scores.append(dc.score)
plt.plot(scatter, scores)
plt.legend()
plt.show()
print(scores)


dc.scatter_fraction = 0.01
dc.lifetime_spectrum = [1, 2]
dc.update()
y = dc.y
dc.data = np.random.poisson(y / max(y) * 10000)

scores = list()
lts = np.linspace(0.01, 4, 10)
for lt in lts:
    dc.lifetime_spectrum = [1, lt]
    scores.append(dc.score)
    plt.semilogy(dc.y)
plt.show()

plt.plot(lts, scores)
plt.legend()
plt.show()
print(scores)

plt.semilogy(dc.x, dc.y)
plt.semilogy(dc.x, dc.data)
plt.show()

dc.lifetime_spectrum = [1, 2]
dc.update()
m = dc.y
d = dc.data
e = dc.data_noise

ssw = ((d-m)/e)
plt.plot(ssw)
plt.show()

