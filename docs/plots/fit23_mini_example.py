import fit2x
import numpy as np
import pylab as p

irf = np.array(
    [0, 0, 0, 260, 1582, 155, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 22, 1074, 830, 10, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0], dtype=np.float64
)

data = np.array(
    [
        0, 0, 0, 1, 9, 7, 5, 5, 5, 2, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 2, 2, 2, 3, 0, 1, 0,
        1, 1, 1, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    ]
)

settings = {
    'dt': 0.5079365079365079,
    'g_factor': 1.0,
    'l1': 0.1,
    'l2': 0.2,
    'conv_stop': 31,
    'irf': irf,
    'background': np.zeros_like(irf)
}

fit23 = fit2x.Fit23(**settings)

tau, gamma, r0, rho = 2.2, 0.01, 0.38, 1.22
x0 = np.array([tau, gamma, r0, rho])
fixed = np.array([0, 1, 1, 0])
r = fit23(
    data=data,
    initial_values=x0,
    fixed=fixed
)
data = fit23.data
model = fit23.model
