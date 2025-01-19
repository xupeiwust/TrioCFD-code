#!/bin/python3
# -*- coding: utf-8 -*-

import glob
import os
import sys

import numpy as np

proj = os.getenv("project_directory")
dir_module = os.path.join(proj, "share", "PyTools3")
sys.path.append(dir_module)
dir_module = os.path.join(dir_module, "commons")
sys.path.append(dir_module)

import commons.DNSTools as dtool
from matplotlib.pyplot import *
from follow_calcs.thplot import get_list_of_files

import readline
readline.parse_and_bind('tab:complete')

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.font_manager as font_manager
print(mpl.rcParams['font.family'])
# mpl.rc('font', family='sans-serif') 
# mpl.rc('font', serif='Helvetica Neue') 
mpl.rc('text', usetex='false') 

# regle un pb de mpl pour enregistrer des figures trop lourdes
matplotlib.rcParams['agg.path.chunksize'] = 100000


def parse_args(jdd):
    if len(sys.argv) >= 2:
        figname_prefix = sys.argv[1]
    else:
        raise Exception("usage: Missing header string to name the plots")
    if len(sys.argv) >= 3:
        tdeb = float(sys.argv[2])
    else:
        tdeb = 0.
    if len(sys.argv) == 4:
        repr_file = sys.argv[3]
        if not os.path.isfile(repr_file):
            raise Exception("usage: Provided file %s does not exist" % repr_file)
        else:
            print("repr_file= %s [given]" % repr_file)
    else:
        repr_file = dtool.getParam(jdd, 'nom_reprise', string=True)
        print("repr_file= %s [read]" % repr_file)
    return figname_prefix, tdeb, repr_file


def read_acceleration_out(head=None, tdeb=0., tmax=10**10, path=None, for_each=1):
    cwd = os.getcwd()
    if path is not None:
        os.chdir(os.path.join(cwd, path))

    fic = "%s_acceleration.out"
    fic_list = get_list_of_files(fic, head=head)
    print(fic_list)
    print("tstep\ttime\tVx\trhoVx\ttauw\tda/dt\tNewT\tacceleration")
    if len(fic_list) == 0:
        os.chdir(cwd)
        return 'none', np.zeros(())

    fict = fic_list[0]
    mat = np.genfromtxt(fict, usecols=(0, 1, 2, 3, 4, 5, 6, 7))
    tdeb = max(tdeb, mat[0, 1])
    tmaxx = min(tmax, mat[-1, 1])
    itdeb = np.argmin(abs(mat[:, 1] - tdeb))
    itmax = np.argmin(abs(mat[:, 1] - tmaxx))
    print("iterations selected: ", itdeb, " : ", itmax)
    mat = mat[itdeb:itmax, :]
    mat = mat[::for_each, :]

    os.chdir(cwd)
    return fict, mat


if __name__ == "__main__":
    print("Provide a name for headers... as in defo180_Eo2.34_Db0213 or defo250_Eo2.34_Db03")
    tmax = 5.e6
    jdd = glob.glob('*.data')[0]

    fig, tdeb, repr_file = parse_args(jdd)

    print("provided header for name_fig is %s." % fig)
    head = jdd.rstrip('.data')

    if os.path.isdir('OUT'):
        out_dir = 'OUT/'
    else:
        out_dir = './'

    # Ici on récupère tous les chemins des ..._acceleration.out dans des dossiers parents
    # fics=glob.glob(fic)
    # p=subprocess.Popen("find . -maxdepth 15 -name "+fic, stdout=subprocess.PIPE, shell=True)
    # lis1=p.stdout.read().split()
    # HACK: ne retient que le premier fichier en fait

    rhol, rhov, mul, muv, alv, beta, sigma, Lz, rhom, Eo, g, Svx, vb, nb, tauw, Db, Lx, Lyi, nx, ny, nz = \
        dtool.get_prop_all(jdd=jdd, repr_file=repr_file, Ret=0.)

    try:
        facv = 1. / (alv * (rhol - rhov))
        facl = 1. / ((1. - alv) * (rhol - rhov))
    except:
        facv = 0.
        facl = 1.

    print("rhom = ", rhom)
    # print("tstep\ttime\tVx\trhoVx\ttauw\tda/dt\tNewT\tacceleration")

    fic, mat = read_acceleration_out(head, tdeb, tmax, path=out_dir)
    t = mat[:, 1]
    #
    f0 = figure(0)
    plot(t, rhol * np.sqrt(abs(mat[:, 4]) / rhol) * Lz / 2. / mul, label=fic.replace('_', '-'))
    ##########
    f1 = figure(1)
    plot(t, mat[:, 3], label=fic.replace('_', '-'))
    ##########
    f11 = figure(11)
    plot(t, mat[:, 4], label=fic.replace('_', '-'))
    ##########
    f2 = figure(2)
    # plot(x, [0.104070213285, 0.104070213285], '--', label="Ret=250")
    # plot(x, [1.00828349999999997e-01,1.00828349999999997e-01], '--', label="Ret=180")
    plot(t, mat[:, 7], label=fic.replace('_', '-'))
    legend()
    ##########
    f3 = figure(3)
    ul = facl * (mat[:, 3] - rhov * mat[:, 2])
    uv = facv * (rhol * mat[:, 2] - mat[:, 3])
    plot(mat[:, 1], ul[:], '-', label='ul')
    plot(mat[:, 1], uv[:], '-', label='uv')
    nlast = min(len(t),5000)
    uvm = uv[-nlast:].mean()
    ulm = ul[-nlast:].mean()
    urm = uvm - ulm
    plot([t[-nlast], t[-1]], [uvm, uvm], '--', label='uvm[-%d:] = %f' % (nlast, uvm))
    print("uv-last %d: %f" % (nlast, uvm))
    print("ur-last %d: %f" % (nlast, urm))
    try:
        print("Reb = %g" % (rhol * urm * Db / mul))
    except:
        print('Reb = db not found, monophasic simulation ?')
    print("rhom*uv[-%d:] = %f" % (nlast, rhom * uvm))
    legend(loc=0)
    ##########
    f4 = figure(4)
    debl = (1. - alv) * rhol * ul
    debv = alv * rhov * uv
    debt = debl + debv
    plot(t, debl[:], '-', label=r'$\alpha_l\,\rho_l\,u_l$')
    plot(t, debv[:], '-', label=r'$\alpha_v\,\rho_v\,u_v$')
    print("debl [min:max]=", debl[:].min(), debl[:].max())
    print("debv [min:max]=", debv[:].min(), debv[:].max())
    print("debt [min:max]=", debt[:].min(), debt[:].max())
    plot(t, debt[:], '-', label=r'$\rho\,u$')
    xlabel(r'$t$ [s]')
    ylabel(r'$\rho\,u\:\mathrm{[kg.m}^{-2}\mathrm{.s}^{-1}]$')
    legend()
    ##########
    f5 = figure(5)
    plot(t, mat[:, 4], '-', label='tauw')
    legend()
    # import pdb; pdb.set_trace()
    ##########i
    f6 = figure(6)
    plot(mat[:, 1], uv[:] - ul[:], '-', label='ur')
    legend(loc=0)

    if os.path.isdir('FIGURES'):
        fig_dir = 'FIGURES'
    else:
        fig_dir = ""

    disable_savefig = False
    if not disable_savefig:
        f0.gca().grid(True)
        f0.savefig(os.path.join(fig_dir, fig + '_Retau.png'))

        f1.gca().grid(True)
        f1.savefig(os.path.join(fig_dir, fig + '_rhoub.png'))

        f11.gca().grid(True)
        f11.savefig(os.path.join(fig_dir, fig + '_ub.png'))

        f2.gca().grid(True)
        f2.savefig(os.path.join(fig_dir, fig + '_S.png'))

        f3.gca().grid(True)
        f3.savefig(os.path.join(fig_dir, fig + '_ul_uv.png'))

        f4.gca().grid(True)
        f4.savefig(os.path.join(fig_dir, fig + '_debits.png'))

        f5.gca().grid(True)
        f5.savefig(os.path.join(fig_dir, fig + '_tauw.png'))

        f6.gca().grid(True)
        f6.savefig(os.path.join(fig_dir, fig + '_ur.png'))

    fic = os.path.join(out_dir, f"{head}_bulles_external_force_every_0.out")
    if os.path.isfile(fic):
        mat = np.genfromtxt(fic)
        tdeb = float(sys.argv[2])
        tdeb = max(tdeb, mat[0, 1])
        tmaxx = min(tmax, mat[-1, 0])
        itdeb = np.argmin(abs(mat[:, 0] - tdeb))
        t = mat[itdeb:, 0]
        Fmean = mat[:, 0:].mean(axis=1)
        f7 = figure(7)
        plot(t, Fmean[itdeb:], '-', label='Fmean')
        legend(loc=0)
        if not disable_savefig:
            f7.gca().grid(True)
            f7.savefig(os.path.join(fig_dir, fig + '_force.png'))
        close(f7)
        n = min(len(t), 1000)
        print("Mean force last %d step : %g" % (n, Fmean[:-n].mean()))

    fic = out_dir + "%s_bulles_centre_x.out" % head

    if os.path.isfile(fic):
        for w in ["x", "y", "z"]:
            mat = np.genfromtxt(fic.replace("_x.out", "_%s.out" % w))
            # mat=np.genfromtxt(fic)
            tdeb = float(sys.argv[2])
            tdeb = max(tdeb, mat[0, 0])
            tmaxx = min(tmax, mat[-1, 0])
            itdeb = np.argmin(abs(mat[:, 0] - tdeb))
            t = mat[itdeb:, 0]
            xmean = mat[:, 1:].mean(axis=1)
            f7 = figure(7)
            plot(t, xmean[itdeb:], '-', label='%smean' % w)
            legend(loc=0)
            if not disable_savefig:
                f7.gca().grid(True)
                f7.savefig(os.path.join(fig_dir, fig + '_%sb.png' % w))
            close(f7)
            n = min(len(t)-10, 1000)
            print("Mean %s position bubble last %d step : %g" % (w, n, xmean[:-n].mean()))
