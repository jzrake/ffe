#!/usr/bin/env python

import sys
import h5py
import matplotlib.pyplot as plt
import numpy as np



def pspec_plot(ax, k, P, comp=0, window=0,
               match_kstar=False,
               match_Pstar=False,
               plot_kwargs={ }):

    # Smooth the spectrum data
    # ----------------------------------------------------------
    if window != 0:
        pass
        #W = window
        #P = smooth(P, 2*W+1)[W:-W]

    # Offset to peak wavenumber and normalize to peak power
    # ----------------------------------------------------------
    if match_kstar: k /= k[np.argmax(P)]
    if match_Pstar: P /= P.max()

    # Set the x-axis label
    # ----------------------------------------------------------
    if match_kstar:
        ax.set_xlabel(r"$k / k^*$")
    else:
        ax.set_xlabel(r"$k / k_0$")

    # Set the y-axis label and compensation index
    # ----------------------------------------------------------
    if type(comp) == int:
        c = float(comp)
        if comp == 1:
            ax.set_ylabel(r"$k \ P(k)$")
        elif comp == 0:
            ax.set_ylabel(r"$P(k)$")
        else:
            ax.set_ylabel(r"$k^{%d} P(k)$" % comp)
            
    elif comp[0] == 0:
        c = 0.0
        ax.set_ylabel(r"$P(k)$" % comp)
    else:
        c = float(comp[0]) / comp[1]
        ax.set_ylabel(r"$k^{%d/%d} P(k)$" % comp)

    # Do the plot
    # ----------------------------------------------------------
    ax.semilogx(k, P * k**c, **plot_kwargs)



def imshow_field(ax, which, component, directory, checkpoint):
    h5f = h5py.File("%s/chkpt.%04d.h5" % (directory, checkpoint))
    dset = '%c%d' % ({'m':'B','e':'E'}[which[0]], component)
    data = h5f[which][dset]
    if len(data.shape) == 1: B = h5f[which][dset][:]
    if len(data.shape) == 2: B = h5f[which][dset][:,:]
    if len(data.shape) == 3: B = h5f[which][dset][:,:,0]
    h5f.close()

    m = ax.imshow(B, interpolation='nearest')
    plt.colorbar(m, ax=ax)
    ax.set_title(which)



def show_log(ax, directory, **kwargs):

    ffe_dat = np.loadtxt("%s/ffe.dat" % directory)
        
    if len(ffe_dat) == 0: # it was empty
        return

    t                 = ffe_dat[:,0]
    electric_energy   = ffe_dat[:,1]
    magnetic_energy   = ffe_dat[:,2]
    magnetic_helicity = ffe_dat[:,3]
    magnetic_monopole = ffe_dat[:,4]

    ax.semilogy(t, electric_energy,   c='g', ls='-', label='E energy',   **kwargs)
    ax.semilogy(t, magnetic_energy,   c='b', ls='-', label='B energy',   **kwargs)
    ax.semilogy(t, magnetic_monopole, c='r', ls='-', label='B monopole', **kwargs)



def summarize_run():

    directory = sys.argv[1]
    checkpoint_numbers = [int(n) for n in sys.argv[2:]]
    nrows = len(checkpoint_numbers)

    fig = plt.figure(figsize=[12,10])

    ax0 = fig.add_subplot(nrows + 1, 1, nrows + 1)
    show_log(ax0, directory)


    for cpn, arg in enumerate(checkpoint_numbers):

        ax1 = fig.add_subplot(nrows + 1, 2, 2 * cpn + 1)
        ax2 = fig.add_subplot(nrows + 1, 2, 2 * cpn + 2)

        imshow_field(ax1, 'magnetic', 3, directory, arg)
        imshow_field(ax2, 'electric', 3, directory, arg)

    plt.show()



def plot_pspec():

    fig = plt.figure(figsize=[12,10])

    ax1 = fig.add_subplot(3, 1, 1)
    ax2 = fig.add_subplot(3, 1, 2)
    ax3 = fig.add_subplot(3, 1, 3)

    directory = sys.argv[1]
    h5f = h5py.File("%s/analysis.h5" % directory)


    for ax, which in zip([ax1, ax2, ax3], ['magnetic', 'electric', 'helicity-real']):
        for g in ['spectra-000400', 'spectra-000600', 'spectra-000800']:
            try:
                k = h5f[g][which]['binlocX'][:]
                P = h5f[g][which]['binval'][:]
            except KeyError as e:
                print "something happened:", e, g, which
                continue

            try:
                pspec_plot(ax, k, P)
            except ValueError as e:
                print "something happened:", e, g, which
                print P
                continue

        ax.set_title(which)

    ax1.set_xlabel('')
    ax2.set_xlabel('')
    plt.show()



def compare_logs():

    fig = plt.figure(figsize=[12,10])
    ax0 = fig.add_subplot(1, 1, 1)
    num = len(sys.argv[1:])

    for n, arg in enumerate(sys.argv[1:]):
        show_log(ax0, arg)

    ax0.legend(loc='best')
    plt.show()



if __name__ == "__main__":
    plot_pspec()
    #summarize_run()
    #compare_logs()
