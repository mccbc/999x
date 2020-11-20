import sys, os
sys.path.append('/home/bcm2vn/research/wind/vis/python/')
cwd = os.getcwd()
from plot_1d import *
from glob import glob
sys.path.append(cwd)

plot_var_anim(sorted(glob('eddingtonwind.block0.out2.*')), [1, [8, 11], 3, 1], ylim=[None, (1e-6, 1), None, None], xlim=[None, None, (1, 10), None], xscale='log', yscale='log', dt=1e-4, sparsity=10, trace=5)
plot_var_anim(sorted(glob('eddingtonwind.block0.out2.*')), [1, [8, 11], 3, 1], xscale='log', yscale='log', dt=1e-4, sparsity=10, trace=5, outfile="highspatialres_plot.mp4")
