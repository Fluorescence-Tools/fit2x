"""
===============================
Pile-up
===============================

"""
import pylab as p
import numpy as np
import scipy.stats
import fit2x


def model_irf(
        n_channels: int = 256,
        period: float = 32,
        irf_position_p: float = 2.0,
        irf_position_s: float = 2.0,
        irf_width: float = 0.25
):
    time_axis = np.linspace(0, period, n_channels * 2)
    irf_np = scipy.stats.norm.pdf(time_axis, loc=irf_position_p, scale=irf_width) + \
             scipy.stats.norm.pdf(time_axis, loc=irf_position_s, scale=irf_width)
    return irf_np, time_axis


# setup some parameters
n_channels = 128
n_corrections = 5
n_photons = 120
irf_position_p = 2.0
irf_position_s = 18.0
irf_width = 0.25
period, g, l1, l2, conv_stop = 32, 1.0, 0.1, 0.1, n_channels // 2 - 1
tau, gamma, r0, rho = 2.0, 0.01, 0.38, 1.2
np.random.seed(0)

# compute a irf
irf_np, time_axis = model_irf(
    n_channels=n_channels,
    period=period,
    irf_position_p=irf_position_p,
    irf_position_s=irf_position_s,
    irf_width=irf_width
)
dt = time_axis[1] - time_axis[0]
conv_stop = min(len(time_axis), conv_stop)
param = np.array([tau, gamma, r0, rho])
corrections = np.array([period, g, l1, l2, conv_stop])

# compute a model function
model = np.zeros_like(time_axis)
bg = np.zeros_like(time_axis)
fit2x.modelf23(param, irf_np, bg, dt, corrections, model)
n_photons = 5e5
model *= n_photons

pile_up_model = np.copy(model)
fit2x.add_pile_up_to_model(
    model=pile_up_model,
    data=model, # the model is modified in-place. Thus make a copy,
    repetition_rate=1./period * 1000,
    dead_time=120.0,
    measurement_time=0.1,
    pile_up_model="coates"
)
p.semilogy(model, label='no pileup')
p.semilogy(pile_up_model, label='with pileup')
#p.plot(irf_np)
p.legend()
p.show()

