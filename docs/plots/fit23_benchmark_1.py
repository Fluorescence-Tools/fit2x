"""
Generate n_samples random decays in range (tau_min, tau_max) with
n_photons photons and fits lifetime. Compares recovered lifetime with
fitted lifetime
"""
import fit2x
import numpy as np
import pylab as p
np.random.seed(0)

# setup some parameters
tau_min = 0.1
tau_max = 5.0

n_photons_min = 20
n_photons_max = 120
n_photon_step = 2

n_samples = 2000

n_channels = 32
irf_position_p = 2.0
irf_position_s = 18.0
irf_width = 0.25
period, g_factor, l1, l2, conv_stop = 32, 1.0, 0.1, 0.1, 31
dt = 0.5079365079365079
np.random.seed(0)

irf_np = np.array([0, 0, 0, 260, 1582, 155, 0, 0, 0, 0,
                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 0, 22, 1074, 830, 10, 0, 0,
                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                   0, 0, 0, 0], dtype=np.float64)

bg = np.zeros_like(irf_np)

p.show()
fit23 = fit2x.Fit23(
    dt=dt,
    irf=irf_np,
    background=bg,
    period=period,
    g_factor=g_factor,
    l1=l1, l2=l2
)

tau, gamma, r0, rho = 2.0, 0.01, 0.38, 1.22
x0 = np.array([tau, gamma, r0, rho])
fixed = np.array([0, 1, 1, 0])


n_photon_dict = dict()
for n_photons in range(n_photons_min, n_photons_max, n_photon_step):
    tau_sim = list()
    tau_recov = list()
    n_photons = int(n_photons)
    for i in range(n_samples):
        tau = np.random.uniform(tau_min, tau_max)
        # n_photons = int(5./tau * n_photon_max)
        param = np.array([tau, gamma, r0, rho])
        corrections = np.array([period, g_factor, l1, l2, conv_stop])
        model = np.zeros_like(irf_np)
        bg = np.zeros_like(irf_np)
        fit2x.modelf23(param, irf_np, bg, dt, corrections, model)
        model *= n_photons / np.sum(model)
        data = np.random.poisson(model)
        r = fit23(
            data=data,
            initial_values=x0,
            fixed=fixed
        )
        # print("tau_sim: %.2f, tau_recov: %s" % (tau, r['x'][0]))
        tau_sim.append(tau)
        tau_recov.append(r['x'][0])
        n_photon_dict[n_photons] = {
                'tau_simulated': np.array(tau_sim),
                'tau_recovered': np.array(tau_recov)
        }

devs = list()
for k in n_photon_dict:
    tau_sim = n_photon_dict[k]['tau_simulated']
    tau_recov = n_photon_dict[k]['tau_recovered']
    dev = (tau_recov - tau_sim) / tau_sim
    devs.append(dev)


fig, ax = p.subplots(nrows=1, ncols=3)
ax[0].semilogy([x for x in fit23._m_param.get_data()], label='Data')
ax[0].semilogy([x for x in fit23._m_param.get_irf()], label='IRF')
ax[0].semilogy([x for x in fit23._m_param.get_model()], label='Model')
ax[0].set_ylim((0.1, 10000))
ax[0].legend()
k = list(n_photon_dict.keys())[0]
tau_sim = n_photon_dict[k]['tau_simulated']
tau_recov = n_photon_dict[k]['tau_recovered']
dev = (tau_recov - tau_sim) / tau_sim
ax[1].plot(tau_sim, dev, 'o', label='#Photons: %s' % k)
ax[1].set_ylim((-1.5, 1.5))
ax[0].legend()
sq_dev = np.array(devs)**2
ax[2].plot(list(n_photon_dict.keys()), np.sqrt(sq_dev.mean(axis=1)))
p.show()
