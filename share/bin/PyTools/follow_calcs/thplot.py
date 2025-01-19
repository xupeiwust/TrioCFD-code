#!/bin/python
import argparse
import glob
import os
import sys

#import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# regle un pb de mpl pour enregistrer des figures trop lourdes
#matplotlib.rcParams['agg.path.chunksize'] = 100000


# ce script va definir des fonctions a utiliser dans la console et peut s'executer avec son main

def read_temperature_out(head=None, tdeb=0., tfin=10**12, path=None):
    """Read the temperature_*.out files and stock them separately into pandas dataframe.

    Args:
        head: name of the data file
        tdeb: starting time for statistics
        tfin: end time for statistics

    Returns:
        list of pandas DataFrame containing ['tstep', 'time', 'theta_adim', 'nu', 'rho_cp_u'] for each temperature_*.out
         file.

    """
    cwd = os.getcwd()
    if path is not None:
        os.chdir(os.path.join(cwd, path))
    fict = "%s*_temperature_*.out"
    listt = get_list_of_files(fict, head=head)
    print("columns : ['tstep', 'time', 'theta_adim', 'nu', 'rho_cp_u']")

    matt_list = []
    for i, fict in enumerate(listt):
        print(fict)
        matt = pd.read_csv(fict, sep=' ', names=['tstep',  'time', 'theta_adim', 'nu', 'rho_cp_u'])
        matt = matt.query('@tfin >= time >= @tdeb')
        matt.drop_duplicates(subset='time', inplace=True)
        matt.sort_values('time', inplace=True)
        matt_list.append(matt)

    if path is not None:
        os.chdir(cwd)
    return matt_list


def read_temperature_swarm_out(head=None, tdeb=0., tfin=10**12, path=None, for_each=1):
    """Read the temperature_*bulles.out files and stock them separately into pandas dataframe.

    Args:
        head: name of the data file
        tdeb: starting time for statistics
        tfin: end time for statistics

    Returns:
        list of pandas DataFrame containing ["tstep", "time", "Tl", "Tv", "Etot_a", "Etot_h", "Elpu", "Elta", "Elth",
        "Evpu", "Evta", "Evth", "Ema", "Emh"] for each temperature_*.out
         file.

    """
    cwd = os.getcwd()
    if path is not None:
        os.chdir(os.path.join(cwd, path))
    fict = "%s_source_temperature*bulles.out"
    listt = get_list_of_files(fict, head=head)
    print("columns : [tstep, time, Tl, Tv, Etot_a, Etot_h, Elpu, Elta, Elth, Evpu, Evta, Evth, Ema, Emh]")

    matt_list = []
    for i, fict in enumerate(listt):
        print(fict)
        matt = pd.read_csv(fict, sep=' ', names=["tstep", "time", "Tl", "Tv", "Etot_a", "Etot_h", "Elpu",
                                                 "Elta", "Elth", "Evpu", "Evta", "Evth", "Ema", "Emh"])
        matt = matt.query('@tfin >= time >= @tdeb')
        matt.drop_duplicates(subset='time', inplace=True)
        matt = matt.loc[::for_each]
        matt = get_nussselt_for_swarm(matt)

        matt_list.append(matt)

    if path is not None:
        os.chdir(cwd)
    return matt_list


def get_list_of_files(regex, head=None, strip=None):
    if head is None:
        jdd = glob.glob('*.data')
        head = jdd[0].split('.')[0]
    if strip is None:
        strip = regex % head
    fict = regex % head
    print(fict)
    print(os.getcwd())
    listt = glob.glob(fict)
    if len(listt) > 1:
        listt.sort(key=lambda item: float(item.strip(strip)))
    print('file list to read : ', listt)
    return listt


def load_matt_and_sort_time(fict, tdeb, tfin, for_each):
    matt = pd.read_csv(fict, sep=' ', header=None)
    matt.rename({0: 'time'}, axis='columns', inplace=True)
    matt = matt.query('@tfin >= time >= @tdeb')
    matt.drop_duplicates(subset='time', inplace=True)
    matt = matt.loc[::for_each]
    matt.set_index('time', inplace=True)
    matt.dropna(axis='columns', inplace=True)
    return matt


def read_Ti_qi_Si_out(head=None, tdeb=0., tfin=10**12, path=None, for_each=1):
    """Read the *_bulles_Ti_*.out and *_bulles_qi_*.out files and stock them into pandas dataframe.

    Args:
        head: name of the data file
        tdeb: starting time for statistics
        tfin: end time for statistics

    Returns:
        list of pandas DataFrame containing ["time", "Ti_0", "Ti_1", ...]
         and list of pandas DataFrame containing ["time", "qi bulle0", "qi bulle1", ...] for each *bulles_*i_*.out.out
         file.

    """
    cwd = os.getcwd()
    if path is not None:
        os.chdir(os.path.join(cwd, path))
    fict_Ti = "%s_bulles_Ti_*.out"
    listt = get_list_of_files(fict_Ti, head=head)
    print("columns : [time, Ti_0, ..., Ti_bullen]")

    matt_list_T = []
    for i, fict in enumerate(listt):
        matt = load_matt_and_sort_time(fict, tdeb, tfin, for_each)
        matt_list_T.append(matt)

    fict_qi = "%s_bulles_phin_*.out"
    listt = get_list_of_files(fict_qi, head=head)
    print("columns : [time, qi_0, ..., qi_bullen]")

    matt_list_q = []
    for i, fict in enumerate(listt):
        matt = load_matt_and_sort_time(fict, tdeb, tfin, for_each)
        matt_list_q.append(matt)

    fict_si = "%s_bulles_surface.out"
    listt = get_list_of_files(fict_si, head=head)
    print("columns : [time, si_0, ..., si_bullen]")

    matt_list_s = []
    for i, fict in enumerate(listt):
        matt = load_matt_and_sort_time(fict, tdeb, tfin, for_each)
        matt_list_s.append(matt)

    if path is not None:
        os.chdir(cwd)
    return matt_list_T, matt_list_q, matt_list_s


def plot_temperature_stat(matt_list, fig_name='DNS', savefig=True, show=True):
    """Creates figures from the temperature stat from ``matt_list``.

    Args:
        matt_list: list of pandas DataFrame containing ['tstep', 'time', 'theta_adim', 'nu', 'rho_cp_u']
        fig_name: the prefix name for the figures saved
        savefig (bool): flag to save the figures
        show (bool): flag to show figures

    .. jupyter-execute::

       from follow_calcs.thplot import *
       matt_list = read_temperature_out(path='./stat_file/Dabiri/', head='DNS')
       plot_temperature_stat(matt_list)

    """

    ftemp = []
    for i, matt in enumerate(matt_list):
        ftemp.append(plt.figure())
        plt.plot(matt['time'], matt['theta_adim'], '-', label='Theta_adim_moy')
        plt.title('temperature %d' % i)
        plt.legend()

        ftemp.append(plt.figure())
        plt.plot(matt['time'], matt['nu'], '-', label='Nu')
        plt.title('temperature %d' % i)
        plt.legend()

        ftemp.append(plt.figure())
        plt.plot(matt['time'], matt['rho_cp_u'], '-', label='rho_cp_u')
        plt.title('temperature %d' % i)
        plt.legend()

    if savefig:
        for i, f in enumerate(ftemp):
            f.gca().grid(True)
            f.savefig(fig_name+'_T%d.png' % i)

    if show:
        plt.show()


def plot_temperature_swarm(matt_list, fig_name='DNS', savefig=True, show=True):
    """Creates figures from the temperature stat from ``matt_list``.

    Args:
        matt_list: list of pandas DataFrame containing ["tstep", "time", "Tl", "Tv", "Etot_a", "Etot_h", "Elpu", "Elta",
        "Elth", "Evpu", "Evta", "Evth", "Ema", "Emh"]
        fig_name: the prefix name for the figures saved
        savefig (bool): flag to save the figures
        show (bool): flag to show figures

    """

    ftemp = []
    for i, matt in enumerate(matt_list):
        ftemp.append(plt.figure())
        plt.plot(matt['time'], matt['Tl'], '-', label=r'$T_l$')
        plt.plot(matt['time'], matt['Tv'], '-', label=r'$T_v$')
        plt.title('temperature %d' % i)
        plt.legend()

        ftemp.append(plt.figure())
        plt.plot(matt['time'], matt['Nu'], '-', label='Nu')
        plt.title('temperature %d' % i)
        plt.legend()

    if savefig:
        for i, f in enumerate(ftemp):
            f.gca().grid(True)
            f.gca().minorticks_on()
            f.gca().grid(which='minor', alpha=0.2)
            f.savefig(fig_name+'_T%d.png' % i)

    if show:
        plt.show()

# def read_temperature_bulles_out(head=None, path=None, tdeb=0., tfin=10**12):
#     """Cette méthode vise à lire les fichiers *temperature_adim_bulles.out
#
#     Args:
#         path (str): nom du jdd
#
#     Return:
#
#     """
#     if path is not None:
#         cwd = os.getcwd()
#         os.chdir(os.path.join(cwd, path))
#     if head is None:
#         jdd = glob.glob('*.data')
#         head = jdd[0].split('.')[0]
#     fict = "%s_source_temperature_*bulles*.out" % head
#     print(fict)
#     print(os.getcwd())
#     listt = glob.glob(fict)
#     listt.sort()
#     print('file list to read : ', listt)
#
#     matt_list = []
#     for i, fict in enumerate(listt):
#         print(fict)
#         hd = get_header_out(fict)
#         matt = pd.read_csv(fict, sep=' ', names=hd, skiprows=1)
#         # matt = pd.read_csv(fict, sep=' ', names=['tstep', 'T_l', 'T_v'], skiprows=1)
#         matt = matt.query('@tfin >= tstep >= @tdeb')
#         matt_list.append(matt)
#
#     if path is not None:
#         os.chdir(cwd)
#     return matt_list


def get_header_out(fict):
    """
    Cette fonction récupère le header du fichier out demandé

    Args:
        fict (str): file name

    Returns:
        (str): the header list
    """
    with open(fict, 'r') as f:
        hd = f.readline()
        hd = hd.split('\t')
        hd[0] = hd[0].split(' ')[-1]
        hd[-1] = hd[-1].split('\n')[0]
    return hd


def get_nussselt_for_swarm(df):
    df['Nu'] = np.gradient(df.Tl, df.time) / (df.Tl - df.Tv)
    return df


if __name__ == '__main__':
    # ici on lance les methodes qui lisent les fichiers temperature.out puis on trace les statistiques selon les
    # arguments qui sont passes en entree.
    print(os.getcwd())
    jdd = glob.glob('*.data')[0]
    jdd = jdd.split('.')[0]
    parser = argparse.ArgumentParser(description="Ce programme peut enregistrer des images de l'évolution de la "
                                                 "thermique du cacul")
    parser.add_argument("--data", "-d", help="name for the header of _temperature_*.out", default=jdd)
    parser.add_argument("--fig_name", help="the prefix name of the figures", default="fig")
    parser.add_argument("--no_plot", action="store_true", help="plot the temperature statistics")
    parser.add_argument("--tdeb", "-s", type=float, help="stating time used for the plot", default=0.)
    parser.add_argument("--tfin", "-e", type=float, help="end time used for the plot", default=5.e6)
    parser.add_argument("--for_each", type=int, help="only consider one for n rows", default=1)
    args = parser.parse_args()

    print("provided header for name_fig is %s." % args.fig_name)

    a = glob.glob('*_temperature_*.out')
    if len(a) > 0:
        if not a[0].endswith('bulles.out'):
            matt_list = read_temperature_out(args.data, args.tdeb, args.tfin)
            plot_temperature_stat(matt_list, fig_name=args.fig_name, show=not args.no_plot)
        else:
            matt_list = read_temperature_swarm_out(args.data, args.tdeb, args.tfin, for_each=args.for_each)
            plot_temperature_swarm(matt_list, args.fig_name, show=not args.no_plot)
    else:
        print('There is no file to analyse')
