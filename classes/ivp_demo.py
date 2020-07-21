from scipy.integrate import solve_ivp
import numpy as np

### CONSTANTS ###
# These can be defined or imported --- they don't have to be global to be
# used by the functions below
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
ivs = [4.02e16, 2.88e7, 0., 0.]           # Initial conditions
X, Y, Z = (0.7, 0.27, 0.03)               # Mass fraction of H, He, Heavy
X_CN = 0.01                               # Mass fraction of CN nuclei

### DERIVED QUANTITIES ###
m = m_H / (2. * X + 0.75 * Y + 0.5 * Z)

# Equations of stellar structure
def struct(r, y, *args):

    # Unpack the dependent variables
    (P, T, M, L) = y

    # Do something with args
    arg1, arg2, arg3 = args
    print(arg1 / arg2 + arg3)

    # Derived variable-dependent quantities. Keep in mind that more lines
    # added here will greatly slow down the solver --- use as few operations
    # as possible.
    rho = m * (P - (1. / 3.) * a * T**4.) / (k * T)
    kappa = 3e25 * Z * (1. + X) * rho * T**-3.5 + 0.2 * (1. + X)

    theta = a * T**3. * m / rho / k
    grad_ad = (18. + 30. * theta + 8. * theta**2.) / \
        (45. + 120. * theta + 32. * theta**2.)
    grad_rad = 3. * P * kappa * L / (16. * np.pi * a * c * T**4. * G * M)
    grad = min(grad_ad, grad_rad)

    tau_pp = 33.8 * (T * 1e-6)**(-1. / 3.)
    tau_CN = 152.3 * (T * 1e-6)**(-1. / 3.)
    epsilon = X * rho * (2.12e3 * X * tau_pp**2. * np.exp(-tau_pp) +
                         3.7e23 * X_CN * tau_CN**2. * np.exp(-tau_CN))

    # Differential equations
    return [-G * M * rho / (r**2.),
            -grad * G * M * rho * T / (r**2. * P),
            4. * np.pi * r**2. * rho,
            4. * np.pi * r**2. * rho * epsilon]

# Relative and absolute tolerances --- local error estimates are kept less
# than atol+rtol*abs(r)
atol, rtol = 3e-14, 3e-14

# Arguments can be passed to the "struct" function. These are just dummy
# variables for demonstration
arg1, arg2, arg3 = (1., 2., 3.)

# Solve the initial value problem, passing args to the diffeq function
sol = solve_ivp(struct, r_bounds, ivs, atol=atol, rtol=rtol,
                args=(arg1, arg2, arg3))

# If the solver fails, increase the tolerance and try again
while len(sol.t) == 1:
    atol = atol * 10.
    rtol = rtol * 10.
    sol = solve_ivp(struct, r_bounds, ivs, atol=atol, rtol=rtol, 
                    args=(arg1, arg2, arg3))

print(sol)

# Here are your solutions as a function of r!
r = sol.t
(P, T, M, L) = sol.y
