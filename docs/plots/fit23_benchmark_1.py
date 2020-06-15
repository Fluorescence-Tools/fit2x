import fit2x
import numpy as np
import pylab as p

# setup some parameters
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
data = np.array(
    [0, 0, 0, 1, 9, 7, 5, 5, 5, 2, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 2, 2, 2, 3, 0, 1, 0,
     1, 1, 1, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
)

p.show()
fit23 = fit2x.Fit23(
    dt=dt,
    irf=irf_np,
    background=bg,
    period=period,
    g_factor=g_factor,
    l1=l1, l2=l2
)

tau, gamma, r0, rho = 2.2, 0.01, 0.38, 1.22
x0 = np.array([tau, gamma, r0, rho])
fixed = np.array([0, 1, 1, 0])

np.random.seed(0)
tau_sim = list()
tau_recov = list()
for i in range(500):
    tau = np.random.uniform(0.1, 5.0)
    n_photon_max = 120
    n_photons = n_photon_max
    # n_photons = int(5./tau * n_photon_max)
    param = np.array([tau, gamma, r0, rho])
    corrections = np.array([period, g_factor, l1, l2, conv_stop])
    model = np.zeros_like(irf_np)
    bg = np.zeros_like(irf_np)
    fit2x.modelf23(param, irf_np, bg, dt, corrections, model)
    model *= n_photons / sum(model)
    data = np.random.poisson(model)
    r = fit23(
        data=data,
        initial_values=x0,
        fixed=fixed
    )
    # print("tau_sim: %.2f, tau_recov: %s" % (tau, r['x'][0]))
    tau_sim.append(tau)
    tau_recov.append(r['x'][0])

p.semilogy([x for x in fit23._m_param.get_data()])
p.semilogy([x for x in fit23._m_param.get_irf()])
p.semilogy([x for x in fit23._m_param.get_model()])
p.ylim((0.1, 10000))
p.show()
tau_sim = np.array(tau_sim)
tau_recov = np.array(tau_recov)
p.plot(tau_sim, (tau_recov - tau_sim) / tau_sim, 'o')
p.ylim((-0.6, 0.6))
p.show()
