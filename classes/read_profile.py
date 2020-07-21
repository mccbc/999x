import numpy as np

def read_profile(filename):

    # read the data in as lines, each line a string
    with open(filename, 'r') as f:
        data = f.readlines()
        f.close()

    # split the header lines containing a number, string and value
    hname = data[1].split()
    hval = data[2].split()

    # make the header info into a dictionary
    hdata = {}
    for name, val in zip(hname, hval):
        try:
            hdata[name] = float(val)
        except:
            pass

    # read the data into an ndarray using genfromtxt
    data = np.genfromtxt(filename, skip_header=5, names=True)

    r = data['radius']*6.96e10
    T = data['temperature']
    P = data['pressure']
    M = data['q'] * hdata['star_mass']*1.99e33
    L = data['luminosity']*3.839e33
    del_rad = data['gradr']
    del_ad = data['grada']

    return hdata, data, r, P, T, M, L, del_rad, del_ad
