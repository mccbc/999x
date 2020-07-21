import numpy as np
import matplotlib.pyplot as plt
from read_profile import read_profile
import pickle
import matplotlib
matplotlib.rc('text', usetex=True)
import pdb
from glob import glob
import argparse
from labellines import *
import matplotlib.cm as cm

R_sun = 6.96e10                # Solar radius
M_sun = 1.99e33                # Solar mass
L_sun = 3.839e33               # Solar luminosity
Pc, Tc = 4.020018591101419e+16, 28822175.80485528

parser = argparse.ArgumentParser(description='Plot a comparison of MESA and model gradients')
parser.add_argument('profiles', type=int, nargs='*', help='Profile number(s) to use with the MESA model')
parser.add_argument('-b', '--beginning', type=int)
parser.add_argument('-e', '--end', type=int)
parser.add_argument('-s', '--step', type=int)
args = parser.parse_args()

if args.profiles == []:
    profiles = list(np.arange(args.beginning, args.end, args.step))
else:
    profiles = args.profiles

profilepath = '/home/bcm2vn/Programs/mesa/star/test_suite/7M_prems_to_AGB/LOGS'

sol = pickle.load(open('sol.p', 'rb'))

# Plot the simulation variables
axis_scales = [1., 1e6, M_sun, L_sun]
y_labels = [r'P (dyn/cm$^2$)', r'T / $10^6$ K', 
            r'M/M_{\odot} ', r'L/L_{\odot}']
ylims = [8e16, 36, 9, 4250]
fig, axs = plt.subplots(2, 2, sharex=True, figsize=(16, 9))
axs = np.ndarray.flatten(axs)
colors = cm.rainbow(np.linspace(0, 1, len(profiles)))

for i in range(len(axs)):
    axs[i].plot(sol.t/R_sun, sol.y[i]/axis_scales[i], 'k', lw=2., label="Model")
    for profile in profiles:
        mesaout = read_profile(profilepath+'/profile'+str(profile)+'.data')
        axs[i].plot(mesaout[2][::-1]/R_sun, mesaout[i+3][::-1]/axis_scales[i], label='MESA profile 110', alpha=0.75, color=colors[profiles.index(profile)])
    axs[i].set_ylabel(y_labels[i])
    axs[i].set_xlabel(r'R/R_{\odot}')
    axs[i].set_xlim((0, 3.25))
    axs[i].set_ylim(top=ylims[i])
    if i in [0, 1]:
        axs[i].legend(loc='upper right')
    else:
        axs[i].legend(loc='lower right')
   # labelLines(axs[i].get_lines(), align=False, fontsize=12)

plt.suptitle(r'P$_c$ = {}, T$_c$ = {}'.format(Pc, Tc))
plt.tight_layout()
plt.subplots_adjust(top=0.95)
plt.show()
#plt.savefig('zams_compare.pdf')
