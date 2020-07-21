from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.rc('text', usetex=True)
import pdb
import pickle

### CONSTANTS ###
G = 6.67e-8                    # Gravitational constant
c = 2.99e10                    # Speed of light
a = 7.5646e-15                 # Radiation density constant
M_sun = 1.99e33                # Solar mass
L_sun = 3.839e33               # Solar luminosity
R_sun = 6.96e10                # Solar radius
m_H = 1.67e-24                 # Hydrogen mass
k = 1.38e-16                   # Boltzmann constant


def shoot(Pc, Tc, plot=False):

    global aout
    ### PARAMETERS ###
    M_star = 8. * M_sun                       # Stellar mass
    r_bounds = (1e-2, 8. * R_sun)             # Bounds of integration in r
    ivs = [Pc, Tc, 0., 0.]                    # Initial conditions
    X, Y, Z = (0.7, 0.27, 0.03)               # Mass fraction of H, He, Heavy
    m = m_H/(2.*X + 0.75*Y + 0.5*Z)           # Mean atomic mass of gas
    X_CN = 0.01                               # Mass fraction of CN nuclei


    # Equations of stellar structure
    def struct(r, (P, T, M, L)):
        global aout

        # Derived variable-dependent quantities
        rho = m * (P - (1./3.) * a * T**4.) / (k * T)
        kappa_kramer = 3e25 * Z * (1. + X) * rho * T**-3.5
        kappa_es = 0.2 * (1. + X)
        kappa = kappa_kramer + kappa_es

        theta = a * T**3. * m / rho / k
        grad_ad = (18. + 30.*theta + 8.*theta**2.)/(45. + 120.*theta + 32.*theta**2.)
        grad_rad = 3. * P * kappa * L / (16. * np.pi * a * c * T**4. * G * M)
        grad = min(grad_ad, grad_rad)

        tau_pp = 33.8 * (T * 1e-6)**(-1./3.)
        tau_CN = 152.3 * (T * 1e-6)**(-1./3.)
        epsilon_pp = X * rho * (2.12e3 * X * tau_pp**2. * np.exp(-tau_pp))
        epsilon_CN = X * rho * (3.7e23 * X_CN * tau_CN**2. * np.exp(-tau_CN))
        epsilon = epsilon_pp + epsilon_CN

        gas_pressure = rho * k * T / m 
        rad_pressure = a * T**4. / 3.

        # Output some variables to an array for plotting later
        output_vars = [r, grad_ad, grad_rad, epsilon_pp, epsilon_CN, 
                       kappa_kramer, kappa_es, gas_pressure, rad_pressure]
        try:
            aout = np.vstack((aout, np.array(output_vars)))
        except:
            aout = np.array(output_vars)

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

    pickle.dump(sol, open('sol.p', 'wb'))
    pickle.dump(aout.T, open('aout.p', 'wb'))

    if plot==True:
        axis_scales = [1., 1e6, M_sun, L_sun]
        y_labels = [r'P (dyn/cm$^2$)', r'T / $10^6$ K', 
                    r'M/M_{\odot} ', r'L/L_{\odot}']

        fig, axs = plt.subplots(2, 2, sharex=True, figsize=(16, 9))
        axs = np.ndarray.flatten(axs)

        for i in range(len(axs)):
            axs[i].plot(sol.t/R_sun, sol.y[i]/axis_scales[i], 'c')
            axs[i].set_ylabel(y_labels[i])
            axs[i].set_xlabel(r'R/R_{\odot}')

        plt.suptitle(r'P$_c$ = {}, T$_c$ = {}'.format(Pc, Tc))
        plt.tight_layout()
        plt.subplots_adjust(top=0.95)
        plt.savefig('ivp_vars.pdf')

        aout = aout.T
        fig, axs = plt.subplots(2, 2, sharex=True, figsize=(11, 8.5))
        axs = np.ndarray.flatten(axs)

        linelabels = [r'$\nabla_{ad}$', r'$\nabla_{rad}$', r'$\epsilon_{pp}$', 
                      r'$\epsilon_{CN}$', r'$\kappa_{Kramer}$', 
                      r'$\kappa_{es}$', r'P$_{gas}$', r'P$_{rad}$']

        ylabels = [r'$\nabla$', r'$\epsilon$ [erg s$^{-1}$]', 
                   r'$\kappa$ [cm$^2$ g$^{-1}$]', r'P [dyn cm$^{-2}$]']

        for j in range(len(axs)):
            axs[j].plot(aout[0]/R_sun, aout[(2*j)+1], 'c', label=linelabels[2*j])
            axs[j].plot(aout[0]/R_sun, aout[(2*j)+2], 'm', label=linelabels[2*j+1])
            axs[j].set_xlabel(r'R/R_{\odot}')
            axs[j].set_ylabel(ylabels[j])
            axs[j].legend()

        axs[1].set_yscale('log')
        axs[3].set_yscale('log')

        plt.tight_layout()
        plt.savefig('aout_vars.pdf')

    return sol

def extract_TM(solution):
    return solution.y[1][-1]/1e6, solution.y[2][-1]/M_sun

if __name__ == "__main__":

    refine = False
    Pguess, Tguess = 4.020018591101419e+16, 28822175.80485528
    #Pguess, Tguess = 4.02e16, 2.88e7
    if refine:
        Pstep, Tstep = 1e2, 1e-3
        T, M = extract_TM(shoot(Pguess, Tguess))

        while (T > 1. or T < 1e-2) or M > 8.1:

            if T > 1.:
                Tguess += Tstep
                T, M = extract_TM(shoot(Pguess, Tguess))
                print(Pguess, Tguess, T, M)
            elif T < 1e-2:
                Tguess -= Tstep
                T, M = extract_TM(shoot(Pguess, Tguess))
                print(Pguess, Tguess, T, M)
            else:
                Tstep = Tstep * 0.1

            if M > 8.:
                Pguess += Pstep
                T, M = extract_TM(shoot(Pguess, Tguess))
                print(Pguess, Tguess, T, M)
            else:
                Pstep = Pstep * 0.1

        print(Pguess, Tguess, T, M)
    sol = shoot(Pguess, Tguess, plot=True)


    # Pick a fitting radius and calculate stellar radius, effective temperature
    R = sol.t
    P, T, M, L = sol.y

    R_ = np.linspace(1., 8., 200.) * R_sun
    afit = np.array([])

    for j in range(len(R_)):
        i = (np.abs(R - R_[j])).argmin()

        R_fit = R[i]
        P_fit, T_fit, M_fit, L_fit = P[i], T[i], M[i], L[i]

        X, Y, Z = (0.7, 0.27, 0.03)               # Mass fraction of H, He, Heavy
        m = m_H/(2.*X + 0.75*Y + 0.5*Z)           # Mean atomic mass of gas

        R_star = (1./R_fit - 4.*k*T_fit/G/M[-1]/m)**-1.
        T_eff = (L[-1] / np.pi / R_star**2. / c / a)**0.25

        try:
            afit = np.vstack((afit, [R_fit, R_star, T_eff]))
        except:
            afit = np.array([R_fit, R_star, T_eff])

    afit = afit.T
    
    fig, ax1 = plt.subplots(1, 1, figsize=(8, 6))
    ax1.plot(afit[0]/R_sun, afit[1]/R_sun, 'c', label=r'R$_\star$')
    ax1.plot(afit[0]/R_sun, afit[2]/(10.**4.), 'm', label=r'T$_{eff}$ / 10$^4$ K')
    ax1.set_xlabel(r'R$_{fit}$/R$_{\odot}$')
    ax1.set_title('Stellar Radius and Effective Temperature v. Fitting Radius')
    ax1.legend()

    plt.tight_layout()
    plt.savefig('stellar_radius.pdf')

    R_fit = afit[0, 52]
    R_star = afit[1, 52]
    T_eff = afit[2, 52]
    print("R_fit={}, R_star={}, T_eff={}".format(R_fit/R_sun, R_star/R_sun, T_eff))

    # Fitting mass integral
    rho_0 = np.pi * c * G * M[-1] * m * a / 48. / aout[6, -1] / L[-1] / k * (G * M[-1] * m / k / R_star)**3.
    x = R_fit / R_star
    Mstar_minus_Mfit = 4 * np.pi * rho_0 * R_star**3. * (-np.log(x) + 3.*x - 3.*x**2./2. + x**3./3. - 11./6.)
    print(Mstar_minus_Mfit/M_sun)

