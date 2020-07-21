import numpy as np
import matplotlib.pyplot as plt
from read_profile import read_profile
import pickle
import matplotlib
matplotlib.rc('text', usetex=True)
import argparse
from labellines import *
import pdb
import matplotlib.cm as cm

R_sun = 6.96e10                # Solar radius
M_sun = 1.99e33                # Solar mass
L_sun = 3.839e33               # Solar luminosity

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
aout = pickle.load(open('aout.p', 'rb'))

colors = cm.rainbow(np.linspace(0, 1, len(profiles)))
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
ax.plot(aout[0]/R_sun, aout[1], 'k--', lw=2, label=r'Model $\nabla_{ad}$', alpha=0.75)
ax.plot(aout[0]/R_sun, aout[2], 'k', lw=2, label=r'Model $\nabla_{rad}$', alpha=0.75)
ax.axvspan(0, 0.7444, alpha=0.25, color='yellow', label='Convective in model')
ax.axvspan(0, 0.572842, alpha=0.25, color='blue', label='Convective in MESA')
ax.legend()

for profile in profiles:
  mesaout = read_profile(profilepath+'/profile{}.data'.format(profile))
  ax.plot(mesaout[2][::-1]/R_sun, mesaout[8][::-1], linestyle='dashed', lw=1, label=profile, alpha=0.75, color=colors[profiles.index(profile)])
  ax.plot(mesaout[2][::-1]/R_sun, mesaout[7][::-1], lw=1, label=profile, alpha=0.75, color=colors[profiles.index(profile)])
ax.set_xlabel(r'R/R_{\odot}')
ax.set_xlim((0, 3))
labelLines(plt.gca().get_lines()[2:], align=False, fontsize=12, xvals=np.linspace(0.1, 2, 2*len(profiles)))

plt.tight_layout()
plt.show()
#plt.savefig('aout_vars.pdf')
