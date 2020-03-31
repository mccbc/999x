import athena_read
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import scipy.optimize as opt
import scipy.integrate as integrate
import pdb

class container:
    """Simple class to store input"""
    pass

def plot_files(ivars,range,base='eddingtonwind.block0.out2',yscale='linear',step=1):
    
    athfiles = []
    for i in xrange(range[0],range[1]+1,step):
        file = base+'.'+str(i).zfill(5)+'.tab'
        athfiles.append(file)

    plot_var(athfiles,ivars,yscale=yscale)


def plot_var_anim(athfiles, ivars, outfile=None, xscale='log',
                  yscale='log', xlim=None, 
                  ylim=None, norm=1, sparsity=1, trace=1, 
                  dt=None, title=None):

    n_images = len(ivars)
    xplots = int(np.ceil(np.sqrt(n_images)))
    yplots = int(np.around(np.sqrt(n_images)))
    fig, axes = plt.subplots(ncols=yplots, nrows=xplots, figsize=(12, 12))

    # Make single ivars elements iterable, for compatibility with multi-line plot functions
    ivars = [[x] if type(x) is not list else x for x in ivars]

    # From sparsity and number of files, figure out how many lines to make
    num_plots = int(np.floor(len(athfiles)/sparsity))
    
    # Load in the athena data for the files we're using
    athdata = [athena_read.tab(athfiles[k*sparsity]) for k in range(num_plots)]
    num_points = len(athdata[0][0, 0, :, 0])
    print("Number of files: {}\nNumber of plots: {}".format(len(athfiles), num_plots))

    line_attrs = []
    plots = []
    texts = []

    for j, ax in enumerate(np.ndarray.flatten(axes)):

        nested_line_attrs = [] # store lines in nested lists so they can be plotted appropriately

        for l in range(len(ivars[j])): # if more than one var is being plotted on the same plot
            # Set up an empty data structure that contains line data and color
            line_attr = np.zeros(num_plots, dtype=[('x1v', float, num_points), ('var', float, num_points), ('color', float, 4)])

            # Initialize the lines with the athena data
            for i, line in enumerate(line_attr):
                line['x1v'] = athdata[i][0, 0, :, 0]
                line['var'] = athdata[i][0, 0, :, ivars[j][l]]
            nested_line_attrs.append(line_attr)

        line_attrs.append(nested_line_attrs)

        # Plot aesthetics, labels, etc
        if ylim is not None:
            if ylim[j] is not None:
                ax.set_ylim(ylim[j])
        if xlim is not None: 
            if xlim[j] is not None:
                ax.set_xlim(xlim[j])

        ax.set_xscale(xscale)
        ax.set_yscale(yscale)

        first_plot_id = int(athfiles[0].split(".tab")[0][-5:])
        last_plot_id = int(athfiles[-1].split(".tab")[0][-5:])
        plot_id_range = "{}-{}".format(first_plot_id, last_plot_id)

        labels = ['x1v', 'dens', 'energy/pressure', 'v1', 'v2', 'v3', 'Er', 'Er0', 
                  'Fr1', 'Fr2', 'Fr3', 'Fr01', 'Fr02', 'Fr03', 'Pr11', 'Pr22', 
                  'Pr33', 'Pr', 'Pr', 'Pr', 'pr', 'Pr', 'Pr', 'sigma_s', 'sigma_a',
                  'sigma_p']

        ax.set_title(', '.join(labels[z] for z in ivars[j])[:-1]+' vs. x1v')
        ax.set_ylabel(', '.join(labels[z] for z in ivars[j])[:-1])
        ax.set_xlabel('x1v')

        # Draw all the lines, which will be updated during animation
        nested_plots = []
        for l in range(len(ivars[j])):
            line_attr = nested_line_attrs[l]
            plot = ax.plot(line_attr['x1v'].T, line_attr['var'].T, color=(0, 0, 0, 0))
            negplot = ax.plot(line_attr['x1v'].T, -line_attr['var'].T, '--', color=(0, 0, 0, 0))
            nested_plots.append([plot, negplot])
        plots.append(nested_plots)
        text = ax.text(.05, .05, "INIT", transform=ax.transAxes)
        texts.append(text)

    # Make each element of these lists iterable, if not already
    line_attrs = [[x] if type(x) is not list else x for x in line_attrs]
    plots = [[x] if type(x) is not list else x for x in plots]
    texts = [[x] if type(x) is not list else x for x in texts]

    def update(frame_number):

        for u in range(len(ivars)):
            for v in range(len(ivars[u])):
                for w in [0, 1]:
                    line_attr = line_attrs[u][v]
                    plot = plots[u][v][w]
                    text = texts[u][0]

                    # Loop over all the lines up to the current frame
                    for k, line in enumerate(line_attr[:frame_number]):

                        # Only have n lines visible, fade toward edge (n = trace)
                        stepfunc = 1 - np.heaviside(frame_number - trace - k - 1, 1)
                        alpha = abs(stepfunc * (frame_number - trace - k - 1) / trace)

                        # Set the color in the attribute array

                        colors = [
                            (alpha, 0, 1-alpha, alpha), # red to blue
                            (1-alpha, 0, alpha, alpha), # blue to red
                            (1-alpha, alpha, 0, alpha), # green to red
                            (alpha, 1-alpha, alpha, alpha), # magenta to green
                            (1-alpha, alpha, alpha, alpha), # cyan to red
                            (alpha, alpha, 1-alpha, alpha), # yellow to blue
                        ]

                        line['color'] = colors[v]
                    
                        # Set the line's color according to the attribute array, set the text to the current frame number
                        plot[k].set_color(line['color'])

                        # Display either a step ticker or a time ticker, depending on if dt is specified
                        if dt is None:
                            text.set_text("Step {}".format(frame_number*sparsity))
                        else:
                            text.set_text("t = {:.3f}s".format(frame_number*sparsity*dt))

                    # If the animation is complete, reset all lines to alpha of 0
                    if frame_number + 1 == num_plots:
                        line_attr['color'] = np.zeros((num_plots, 4))

    anim = animation.FuncAnimation(fig, update, num_plots, interval=40, blit=False)

    if title is not None:
        plt.suptitle(title)

    plt.tight_layout()

    if outfile is None:
        plt.show()
    else:
        if outfile.split('.')[1] == 'gif':
            writer = 'imagemagick'
        else:
            writer = None
        anim.save(outfile, dpi=300, fps=30./sparsity, writer=writer)


    # 0: x1v
    # 1: dens
    # 2: energy/pressure
    # 3: v1
    # 4: v2
    # 5: v3
    # 6: Er
    # 7: Er0
    # 8: Fr1
    # 9: Fr2
    #10: Fr3
    #11: Fr01
    #12: Fr02
    #13: Fr03
    #14: Pr11
    #15: Pr22
    #16: Pr33
    #17: Pr
    #18: Pr
    #19: Pr
    #20: pr
    #21: Pr
    #22: Pr
    #23: sigma_s
    #24: sigma_a
    #25: sigma_p


def plot_intens(athfiles,ivars,outfile='intens.pdf',xscale='log',yscale='linear',xlim=None,
ylim=None):



    #mu=[1.160841e-01,2.681522e-01,4.137792e-01,5.494671e-01,6.719567e-01,7.783057e-01,8.659595e-01,9.328128e-01,9.772599e-01]
    #r = np.empty(9)
    #yx = np.empty(9)
    #for i in xrange(9):
    #    r[i] = 1./(1.-mu[i])*1.e-4
    #    yx[i] = 0.5
    #print  r
    plt.xscale(xscale)
    plt.yscale(yscale)
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    for file in athfiles:
        data = athena_read.tab(file)

        #nx3=data.shape[0]
        #nx2=data.shape[1]
        nx1=data.shape[2]

        x1v = np.empty(nx1)
        var = np.empty(nx1)
        for j in ivars:
            for i in xrange(nx1):
                x1v[i] = data[0,0,i,0]
                if (data[0,0,i,8] > 1.e-15):
                    #var[i] = data[0,0,i,j]
                    var[i] = data[0,0,i,j]/data[0,0,i,8]
                else:
                    var[i] = 0.
            plt.plot(x1v,var)     
            #print var

    #plt.plot(r,yx,'o')
    plt.savefig(outfile)
    plt.close()


def intens_sphere(mu,r):

    B= 10.
    kapa = 10.
    R = 1.
    if (r < R):
        s = r*mu+R*(1-(r/R)**2*(1.-mu**2))**0.5
    else:
        if (mu > (1-(R/r)**2)**0.5):
            s=2*R*(1-(r/R)**2*(1.-mu**2))**0.5
        else:
            s=0.
    #print r,s,mu
    return (1.-np.exp(-s*kapa))*B/kapa

def flux_sphere(mu,r):

    return mu*intens_sphere(mu,r)

def get_er(r):

    nr = len(r)
    er = np.empty(nr)
    i=0
    nmu=20
    for rad in r:
        dmu=2./float(nmu)
        sum=0.
        mu=-1.
        for k in range(nmu):
            sum=sum+integrate.quad(intens_sphere, mu, mu+dmu, args=(rad))[0]
            mu=mu+dmu        
        er[i]=sum
        i=i+1
    return er/2

def get_fr(r):

    nr = len(r)
    fr = np.empty(nr)
    i=0
    nmu=20
    for rad in r:
        dmu=2./float(nmu)
        sum=0.
        mu=-1.
        for k in range(nmu):
            sum=sum+integrate.quad(flux_sphere, mu, mu+dmu, args=(rad))[0]
            mu=mu+dmu        
        fr[i]=sum
        i=i+1
    return fr/2

def get_intens(r,ml,mh):

    nr = len(r)
    intens = np.empty(nr)
    i=0
    for rad in r:
        intens[i]=integrate.quad(intens_sphere, ml, mh, args=(rad))[0]
        i=i+1
    return intens/(mh-ml)

def read_angles(infile="Rad_angles.txt"):


    muc = []
    wmu = []
    nmu = 19
    delimiter = " "
    with open(infile) as f:
        k = 0
        for line in f:
            k = k+1
            if (k > 10+nmu):
                muc.append(float(line.split()[1]))
                wmu.append(float(line.split()[4]))
    
    muf=np.empty(nmu+1)
    muf[0] = 0.
    for i in xrange(1,nmu):
        muf[i] = muf[i-1]+wmu[i-1]*2
    muf[nmu]= 1.
    return muf

def comp_er(athfiles,outfile='er_comp.pdf',xscale='linear',yscale='log',xlim=None,
ylim=None):

    
    plt.xscale(xscale)
    plt.yscale(yscale)
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    
    for file in athfiles:
        data = athena_read.tab(file)

        nx1=data.shape[2]

        x1v = np.empty(nx1)
        er = np.empty(nx1)
        mer = np.empty(nx1)
        fr = np.empty(nx1)
        mfr = np.empty(nx1) 

        for i in xrange(nx1):
            x1v[i] = data[0,0,i,0]
            er[i] = data[0,0,i,6]
            fr[i] = data[0,0,i,8]
        plt.plot(x1v,er)
        plt.plot(x1v,fr)
        mer = get_er(x1v)
        mfr = get_fr(x1v)
        plt.plot(x1v,mer,':')
        plt.plot(x1v,mfr,':')
        #print var
        #print mod
    plt.xlabel("r")
    plt.ylabel("Er, Fr")
    plt.savefig(outfile)
    plt.close()

def comp_intens(athfiles,outfile='intens_comp.pdf',xscale='linear',yscale='log',xlim=None,
ylim=None):

    #ivars=[26,45]
    #ivars=[42,43,44,45]
    #ivars=[26,27,28,29]
    #ivars=[34,35,36,37]
    #ivars=[38,39,40,41]
    ivars=[44,39,34,29]
    #mu=[1.160841e-01,2.681522e-01,4.137792e-01,5.494671e-01,6.719567e-01,7.783057e-01,8.659595e-01,9.328128e-01,9.772599e-01,9.982377e-01]
    
    muf=read_angles()
    plt.xscale(xscale)
    plt.yscale(yscale)
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)

    
    for file in athfiles:
        data = athena_read.tab(file)

        #nx3=data.shape[0]
        #nx2=data.shape[1]
        nx1=data.shape[2]

        x1v = np.empty(nx1)
        var = np.empty(nx1)
        mod = np.empty(nx1)
 
        for k in xrange(0,len(ivars)):
            j = ivars[k]
            jm = j-26+20
            print muf[jm],muf[jm+1]
            print 'lim:',1./np.sqrt(1-muf[jm]**2),1./np.sqrt(1-muf[jm+1]**2)
            for i in xrange(nx1):
                x1v[i] = data[0,0,i,0]
                var[i] = data[0,0,i,j]

            plt.plot(x1v,var)
            mod = get_intens(x1v,muf[jm],muf[jm+1])
            #print mod
            #print var
            plt.plot(x1v,mod,':')

    plt.xlabel("r")
    plt.ylabel("I")
    #plt.plot(r,yx,'o')
    plt.savefig(outfile)
    plt.close()
