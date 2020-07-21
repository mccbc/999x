from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
rc('text', usetex=True)
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-p', type=float, nargs='?', default=4.02e16)
parser.add_argument('-t', type=float, nargs='?', default=2.88e7)

args = parser.parse_args()

Pc, Tc = args.p, args.t
print(Pc, Tc)

### CONSTANTS ###
G = 6.67e-8                    # Gravitational constant
c = 2.99e10                    # Speed of light
a = 7.5646e-15                 # Radiation density constant
M_sun = 1.99e33                # Solar mass
L_sun = 3.839e33               # Solar luminosity
R_sun = 6.96e10                # Solar radius
m_H = 1.67e-24                 # Hydrogen mass
k = 1.38e-16                   # Boltzmann constant

### PARAMETERS ###
M_star = 8. * M_sun                       # Stellar mass
r_bounds = (1e-2, 8. * R_sun)             # Bounds of integration in r
ivs = [Pc, Tc, 0., 0.]                    # Initial conditions
X, Y, Z = (0.7, 0.27, 0.03)               # Mass fraction of H, He, Heavy
X_CN = 0.01                               # Mass fraction of CN nuclei

### DERIVED QUANTITIES ###
m = m_H/(2.*X + 0.75*Y + 0.5*Z)

ad_array, rad_array, r_array = np.array([]), np.array([]), np.array([])

# Equations of stellar structure
def struct(r, (P, T, M, L)):
    global ad_array, rad_array, r_array

    # Derived variable-dependent quantities
    rho = m * (P - (1./3.) * a * T**4.) / (k * T)
    kappa = 3e25 * Z * (1. + X) * rho * T**-3.5 + 0.2 * (1. + X)

    theta = a * T**3. * m / rho / k
    grad_ad = (18. + 30.*theta + 8.*theta**2.)/(45. + 120.*theta + 32.*theta**2.)
    grad_rad = 3. * P * kappa * L / (16. * np.pi * a * c * T**4. * G * M)
    grad = min(grad_ad, grad_rad)

    ad_array, rad_array, r_array = (np.append(ad_array, grad_ad), 
                                    np.append(rad_array, grad_rad),
                                    np.append(r_array, r))

    tau_pp = 33.8 * (T * 1e-6)**(-1./3.)
    tau_CN = 152.3 * (T * 1e-6)**(-1./3.)
    epsilon = X * rho * (2.12e3 * X * tau_pp**2. * np.exp(-tau_pp)
                         + 3.7e23 * X_CN * tau_CN**2. * np.exp(-tau_CN))

    # Differential equations
    return [-G * M * rho / (r**2.), 
            -grad * G * M * rho * T / (r**2.*P),
            4. * np.pi * r**2. * rho,
            4. * np.pi * r**2. * rho * epsilon]

# Solve initial value problems
atol, rtol = 3e-14, 3e-14
sol = solve_ivp(struct, r_bounds, ivs, atol=atol, rtol=rtol)
while len(sol.t) == 1:
    atol = atol * 10.
    rtol = rtol * 10.
    sol = solve_ivp(struct, r_bounds, ivs, atol=atol, rtol=rtol)

plt.figure()
plt.plot(r_array/R_sun, ad_array, label=r'$\nabla_{ad}$')
plt.plot(r_array/R_sun, rad_array, label=r'$\nabla_{rad}$')
plt.xlabel(r'R / R$_\odot$')
plt.legend()
plt.show()

# Plot aesthetics
axis_scales = [1., 1e6, M_sun, L_sun]
y_labels = [r'P (dyn/cm$^3$)', r'T / $10^6$ K', r'M/M$_{\odot}$ ', r'L/L$_{\odot}$']
fig, axs = plt.subplots(2, 2, sharex=True, figsize=(16, 9))
axs = np.ndarray.flatten(axs)

for i in range(len(axs)):
  axs[i].plot(sol.t/R_sun, sol.y[i]/axis_scales[i])
  axs[i].set_ylabel(y_labels[i])
  axs[i].set_xlabel(r'R/R_{\odot}')

print(sol.y[1][-1]/1e6, sol.y[2][-1]/M_sun)

plt.tight_layout()
plt.show()

