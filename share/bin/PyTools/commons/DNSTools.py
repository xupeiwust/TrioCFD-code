"""
Test docstrings echelle module
"""

#
# Boite a outils pour les post-traitement du IJK...
#
#
# from numpy import *
import numpy as np

# from math import *
import math
import glob, os, re

from IPython.display import display_latex
from commons.Tools import Field
from commons.TrustFiles import DTEVFile  # , BuildFromPath
import json

# import commons.BubbleTools as ob

# MODIF GAB
# import FormatTools as ofrm
# import BubbleTools as ob
# FIN MODIF GAB

alwidth = 0.7
mewidth = 0.75
errlw = 0.25
lwidth = 1.5
pltmewidth = 0
msize = 7.5


#### Utilisation de LaTeX

# try:
#     rc('text', usetex = True)
#     #rcParams['text.latex.preamble']=[r"\usepackage{gensymb}",r"\usepackage{txfonts}",r"\usepackage{nicefrac}",r"\usepackage{amssymb}",r"\usepackage{sistyle}"]
#     rcParams['text.latex.preamble']=[r"\usepackage{txfonts}",r"\usepackage{amssymb}"]
#     rc('legend', fontsize='medium',numpoints=2)
# except:


class SimuFieldEncoder(json.JSONEncoder):
    """
    Cette classe permet d'encoder des Fields au format json.
    """

    def default(self, obj):
        # if isinstance(obj, Simu):
        #     return obj.__dict__
        if isinstance(obj, Field):
            return obj.to_dict()
        return json.JSONEncoder.default(self, obj)


def simu_decoder(obj):
    if "value" in obj.keys():
        return Field.from_dict(obj)
    elif "je_suis_une_simu" in obj.keys():
        return Simu.get_from_dict(obj)
    else:
        return obj


def cd(fold):
    """

    Args:
        fold: a special explanation jsut for test

    Returns:
        nothing at all

    """
    os.chdir(fold)
    return


def pwd():
    return os.getcwd()


def getJddName():
    l = glob.glob("*data")
    if len(l) > 1:
        l = [s for s in l if "repr" not in s]
    jdd = l[0][:-5]
    return jdd


# Recuperation des cdg des bulles :
# t : Le temps
# nb : Le nombre de bulles
# coords[:,it,ib] : Les coords du cdg de la bulle ib a l'instant it
def getBarys(head=None):
    """
    Récupère les CDG des bulles.

    Args:
        Rien, ou bien le nom du jdd

    Returns:
        coords les coordonnées des cdg
        t les instants ou sont lues ces coordonnees
        nb le nombre de bulles

    """
    # En Python3 , on ne peut plus d�clarer directement des nouvelles
    # variables dans un exec. On a vu sur stackoverflow que passer par
    # un dico �a fait le taff.

    drc = {"x": 0, "y": 0, "z": 0}

    if head == None:
        head = getJddName()
        pass
    fics = glob.glob(head.replace(".data", "") + "_bulles_centre_[xyz].out")
    print(fics)
    if len(fics) != 3:
        raise Exception("Error in getBarys: Missing files :", fics)
        pass
    for fic in fics:
        ax = fic[-5]
        print(ax)
        drc[ax] = np.loadtxt(fic)
        drc[ax] = drc[ax][:, 1:]
        pass
    # bon a ce stade fic correspond a un des 3 fic [x,y ou z], on veut le temps
    t = np.loadtxt(fic)[:, 0]
    x, y, z = drc["x"], drc["y"], drc["z"]

    # ~ print ('in getBary : lenx,np.shapex0',len(x), np.shape(x)[0])
    lmin = min(min(np.shape(x)[0], np.shape(y)[0]), np.shape(z)[0])  # -1
    # ~ print ('in getBary : lmin', lmin)
    coords = np.array([x[:lmin, :], y[:lmin, :], z[:lmin, :]])
    _, nb = np.shape(x)

    ### ON NETTOIE LA MEMOIRE ##########################################
    del (drc, fics, ax, x, y, z)
    return (coords, t[:lmin], nb)


# GAB A JOUE ICI #######################################################
def get_R_Vol_Bulle(jdd=None):
    """
    Renvoie le rayon et le volume de la bulle

    Args:
        rien ou le nom du jdd

    Returns:
        rb, vb respectivement le rayon et le volume de la bulle

    """
    if jdd == None:
        jdd = getJddName()
        pass
    vb = getValue("vol_bulle_monodisperse", jdd)
    rb = ((3.0 * vb) / (4.0 * math.pi)) ** (1.0 / 3.0)
    return (rb, vb)


def get_dt(jdd=None):
    """
    Récupère différents pas de temps relatifs à la simulation

    Args:
        rien ou le nom du jdd

    Returns:
        pas de temps de post-traitement; de -//- pour les statistiques de bulle; -//- par plan;
        de sauvegarde. Renvoie aussi le nombre de pas de temps effectués

    """

    if jdd == None:
        jdd = getJddName() + ".data"
        pass
    dt_p = getValue("dt_post ", jdd)
    dt_psb = getValue("dt_post_stats_bulles ", jdd)
    dt_psp = getValue("dt_post_stats_plans ", jdd)
    dt_s = getValue("dt_sauvegarde ", jdd)
    nT = getValue("nb_pas_dt_max ", jdd)
    return (dt_p, dt_psb, dt_psp, dt_s, nT)


# TODO : Gab : Ecrire la fonction get_RelVel
# def get_RelVel(jdd_None):
#    """
#    ####################################################################
#    # On voudrai obtenir la vitesse relative de chaque bulle, pour #####
#    # orienter correctement le frame autour des bulles #################
#    # --> chooper expression_v[x,y,z]_init dans .data ##################
#    # --> se servir de getVelocity (vitesse d'une bulle) ###############
#    # --> ATTN : il faudra faire en sorte que l'orientation du frame ###
#    #            change pas trop vite / trop fort : prendre des vites- #
#    #            ses relatives MOYENNES, enfin les lisses quoi ( sur ###
#    #            dt_post par exemple) ##################################
#    ####################################################################
#    """
#    pass


def getNbElem(jdd=None):
    """
    Récupère le nombre de cellules dans chaque dimensions
    Args :
        Rien ou le nom du jdd
    Returns :
       Nombre de cellules dans chaque dimensions
    """
    if jdd == None:
        jdd = getJddName() + ".data"
        pass
    Nx = getValue("nbelem_i", jdd)
    Ny = getValue("nbelem_j", jdd)
    Nz = getValue("nbelem_k", jdd)
    return np.array([Nx, Ny, Nz])


# FIN JEU GAB ##########################################################


def getDimensions(jdd=None):
    """
    Récupère la longueur du domaine dans chaque direction

    Args :
        Rien ou le nom du jdd

    Returns :
       longueur du domaine dans chaque direction

    """
    if jdd == None:
        jdd = getJddName() + ".data"
    Lx = getValue("uniform_domain_size_i", jdd)
    Ly = getValue("uniform_domain_size_j", jdd)
    Lz = getValue("uniform_domain_size_k", jdd)
    return np.array([Lx, Ly, Lz])


def getVelocity(coords, time, nb, jdd=None):
    if jdd == None:
        jdd = getJddName() + ".data"
    Lx = getValue("uniform_domain_size_i", jdd)
    Ly = getValue("uniform_domain_size_j", jdd)
    Lz = getValue("uniform_domain_size_k", jdd)
    L = np.array([Lx, Ly, Lz])
    dt = time[1:] - time[:-1]
    tvel = (time[1:] + time[:-1]) / 2.0
    ndim, nt, nb = np.shape(coords)
    nt -= 1
    velocity = np.zeros((ndim, nt, nb))
    for direction in [0, 1, 2]:
        LL = L[direction]
        dx = coords[direction, 1:, :] - coords[direction, :-1, :]
        # Correction du dx si traversee de frontiere perio :
        nt, nb = np.shape(dx)
        for i in range(nt):
            for j in range(nb):
                x = dx[i, j]
                if abs(x) > LL * 0.5:
                    print(
                        "Crossing periodic boundary ",
                        direction,
                        " at iter=",
                        i,
                        " at t=",
                        time[i],
                    )
                    print("Old ", x)
                    # Un deplacement superieur a la moitie de la taille du domaine,
                    # C'est forcement un passage par une frontiere perio :
                    while x - LL * 0.5 > 0.0:
                        x -= LL
                    while x + LL * 0.5 < 0.0:
                        x += LL
                    print(" New ", x)
                dx[i, j] = x
        # dx = [ sign(x)*(x % L[direction]) for x in dx ]
        velocity[direction, :, :] = dx / dt[:, np.newaxis]
    print("len(tvel), velocity.shape")
    print(len(tvel), velocity.shape)
    return tvel, velocity


#
# Retourne : temps, valeurs, nombre de bulles
# val[it, ib] : valeur de var au temps it pour la bulle ib
# time, val, nbulles = getOnBulles("centre_x")
def getOnBulles(var, head=None):
    if head == None:
        head = getJddName()
    mat = np.loadtxt(head + "_bulles_" + var + ".out")
    time = mat[:, 0]
    val = mat[:, 1:]
    _, nbulles = np.shape(val)
    return time, val, nbulles


def getTempsIntegrationFile(fic):
    f = open(fic, "r")
    lines = f.readlines()
    tintegration = lines[0].split()[2]
    f.close()
    return float(tintegration)


def getTempsIntegrationList(fics):
    lt = []
    for fic in fics:
        lt.append(getTempsIntegrationFile(fic))
    return np.array(lt)


def buildDicoColonnes(fic):
    d = {}
    f = open(fic, "r")
    lines = f.readlines()
    for val, key in [
        (int(s.split()[2]) - 1, s.split()[4]) for s in lines if "# colonne" in s
    ]:
        d[key] = val
    f.close()
    return d


def getStatsInstant():
    return getStats("ins")


def getStatsMoy():
    return getStats()


# caract = moy or ins : instantan� ou moyen.
# Retourne une matrice avec les stats :
# resu[i,j,k] : i --> temps
#                    j --> Position en z
#                    k --> var.
#
# ltimes : liste des temps.
# lz      : liste des coords en z.
# lvar    : liste des variables stockees.
# val[it, iz, ivar] : La matrice 3D contenant tout ca.
#
# tintegration : le temps d'integration des resulats si caract == moy.
def getStats(caract="moy"):
    if caract == "moy":
        fic_stats = glob.glob("statistiques_*.txt")
        if fic_stats == []:
            fic_stats = glob.glob("*phasique_statistiques_*.txt")
        fic_stats.sort(
            key=lambda item: float(item.strip("monodiphasique_statistiques_.txt"))
        )
        ltimes = np.array(
            [float(f.strip("monodiphasique_statistiques_.txt")) for f in fic_stats]
        )
    else:
        fic_stats = glob.glob("moyenne_spatiale_*.txt")
        if fic_stats == []:
            fic_stats = glob.glob("*phasique_moyenne_spatiale_*.txt")
        fic_stats.sort(
            key=lambda item: float(item.strip("monodiphasique_moyenne_spatiale_.txt"))
        )
        ltimes = np.array(
            [float(f.strip("monodiphasique_moyenne_spatiale_.txt")) for f in fic_stats]
        )
    if fic_stats == []:
        raise Exception("La liste fic_stats est restee vide... ")
    lvar = buildDicoColonnes(fic_stats[0])
    nt = len(ltimes)
    if caract == "moy":
        tintegration = getTempsIntegrationList(fic_stats)
    else:
        tintegration = [0]
    for i, fic in enumerate(fic_stats):
        mat = np.loadtxt(fic)
        t = ltimes[i]
        lz = mat[:, 0]
        if i == 0:
            # Initialise la taille du resu :
            nz, nvar = np.shape(mat)
            resu = np.zeros((nt, nz, nvar))
        resu[i, :, :] = mat[:, :]
    return ltimes, lz, lvar, resu, tintegration


def getValue(key, fic, pre="", default=None):
    f = open(fic, "r")
    lines = f.readlines()
    f.close()
    nb = len(lines)
    #
    if "liquid" in key:
        rliq = re.compile(".*fluide1.*")
        for i, st in enumerate(lines):
            m = rliq.match(st)
            if m:
                for j in range(6):
                    lkey=key.split("_")[0]
                    li = lines[i+2+j]
                    if lkey in li:
                        return float(li.split()[3])
                    if "}" in li:
                        raise Exception("too bad")
        raise Exception(f"should have found {key} in {fic}")
        
    #
    #
    if key.endswith("_k"):
        key = key[:-2]
        rank = 2
        for i, st in enumerate(lines):
            if key in st:
                return float(st.split()[rank+1])
        raise Exception("too bad, not found")
    #
    rc = re.compile(
        pre + ".*" + key + ".*(?P<value>[\-]?[\d]*.?[\d]*[eE]?[+\-]?[\d]*)"
    )
    for i, st in enumerate(lines):
        m = rc.match(st)
        if m:
            val = float(m.group("value"))
            return val
            break
    if default != None:
        return default
    print("On a rien trouve dans ", fic)
    raise Exception("Etonnant, non?")


def getValues(key, fic, pre=""):
    f = open(fic, "r")
    lines = f.readlines()
    f.close()
    nb = len(lines)
    rc = re.compile(pre + key + "\s*(?P<value>[\-]?[\d]*.?[\d]*[eE]?[+\-]?[\d]*)")
    for i, st in enumerate(lines):
        it = re.finditer(rc, st)
        val = []
        for match in it:
            val.append(float(match.group("value")))
        if len(val):
            return np.array(val)
    print("On a rien trouve dans ", fic)
    raise Exception("Etonnant, non?")


def getSondesCoords(fic):
    x = getValues("x=", fic)
    y = getValues("y=", fic)
    z = getValues("z=", fic)
    return np.array([x, y, z]).T, len(x)


# En z+ et En z- :
def evaluateRetau(U, z, jdd=None):
    if jdd == None:
        jdd = getJddName() + ".data"
    rhol = getValue("rho_liquide", jdd)
    mul = getValue("mu_liquide", jdd)
    Lz = getValue("size_dom_k", jdd)
    h = Lz / 2.0
    tauwp = mul * U[:, 0] / z[0]
    tauwm = mul * U[:, -1] / (Lz - z[-1])
    #
    tauw = (tauwp + tauwm) / 2.0
    utau = np.array([np.sqrt(abs(x) / rhol) for x in tauw])
    Retau = rhol * utau * h / mul
    return Retau, tauwp, tauwm


def getFloatingAverage(tintegration, Um, n):
    if n >= len(tintegration):
        print("Too long averaging ", n, " required!")
        return tintegration.mean(), Um.mean()  # Whatever, in the good range
    tf = tintegration[n:]  # temps de fin
    ti = tintegration[:-n]  # temps de debut
    # print("getFloatingAverage:: Attention : il faut prendre ltintegration et pas tm")
    Uf = Um[n:]
    Ui = Um[:-n]
    dt = tf - ti
    t = (tf + ti) / 2.0
    var = (Uf * tf - Ui * ti) / dt
    return t, var


def getFloatBtwIters(tab, ltintegration, it, it2):
    dt = ltintegration[it2] - ltintegration[it]
    res = (ltintegration[it2] * tab[it2, :] - ltintegration[it] * tab[it, :]) / dt
    return res


# v : NumPy np.array
# dz : pas du maillage. On le suppose constant et on suppose que les bords sont a dz/2 du premier et dernier noeuds
#                              (ce qui correspond a un cell_center regulier).
# bc_type = 0 : Pour la vitesse, suppose une valeur nulle sur les bords. neuman_homogene
#              2 : Suppose un gradient normal nul a la paroi, par exemple pour la pression
#              1 : Ne suppose rien, decentre interieur... (ordre 1)
def cell_to_cell_gradient(v, dz, bc_type):
    # if ( (z[1:]-z[:-1]).max() -  (z[1:]-z[:-1]).min() ) > 1e-11): raise Exception("Fonction pour pas uniforme!")
    # Formule centree (ordre 2) pour pas cste dans le domaine :
    grad = np.zeros(len(v))
    grad[1:-1] = (v[2:] - v[:-2]) / (2.0 * dz)
    # Aux bords :
    # Au bord en z- :
    Uc = v[0]
    Up = v[1]
    if bc_type == 0:
        Um = 0.0
        grad[0] = (-4 * Um + 3 * Uc + Up) / (3 * dz)
    elif bc_type == 1:
        grad[0] = (Up - Uc) / dz
    else:
        Um = Uc + (Uc - Up) / 8.0
        grad[0] = (-4 * Um + 3 * Uc + Up) / (3 * dz)
    #
    # Au bord en z+ :
    Uc = v[-1]
    Um = v[-2]
    if bc_type == 0:
        Up = 0.0
        grad[-1] = (-Um - 3 * Uc + 4 * Up) / (3 * dz)
    elif bc_type == 1:
        grad[-1] = (Uc - Um) / dz
    else:
        Up = Uc + (Uc - Um) / 8.0
        grad[-1] = (-Um - 3 * Uc + 4 * Up) / (3 * dz)
    #
    return grad


def deriv(v, z, Lz):
    raise Exception("Tentative pour maillage a pas variable incorrecte sur les bords!")
    xf = (z[1:] + z[:-1]) / 2.0  # Centre des faces
    xf = np.append(np.append(0.0, xf), Lz)
    dx = xf[1:] - xf[:-1]  # dx_i = x_i+1 - x_i
    # L'evaluation du dx entre 2 cellules :
    h = (dx[:-1] + dx[1:]) / 2.0  # contient (Nsommets - 2) valeurs = Ncellules - 1
    h1 = h[:-1]
    h2 = h[1:]
    u_pl = v[2:]
    u_m = v[:-2]
    u_c = v[1:-1]
    # Formule centree (ordre 2) pour pas variable dans le domaine :
    grad = np.zeros(len(v))
    grad[1:-1] = (
        h1 / h2 * u_pl - h2 / h1 * u_m + (h2**2 - h1**2) / (h1 * h2) * u_c
    ) / (h1 + h2)
    # Formule decentree (ordre 2) pour pas variable sur le bord gauche :
    h1 = h[0]
    h2 = h[1]
    u_g = v[0]
    u_pl = v[1]
    u_pl2 = v[2]
    grad[0] = (1 / (h1 * h2 * (h1 + h2))) * (
        -(2 * h1 * h2 + h2**2) * u_g + ((h1 + h2) ** 2) * u_pl - (h1**2) * u_pl2
    )
    # Formule decentree (ordre 2) pour pas variable sur le bord droit :
    h1 = h[-2]
    h2 = h[-1]
    u_d = v[-1]
    u_m = v[-2]
    u_m2 = v[-3]
    grad[-1] = (
        1
        / (h1 * h2 * (h1 + h2))
        * (
            (h2**2) * u_m2
            - ((h1 + h2) ** 2) * u_m
            - ((h2**2) - (h1 + h2) ** 2) * u_d
        )
    )
    return grad


# v : NumPy np.array
# dz : pas du maillage. On le suppose constant et on suppose que les bords sont a dz/2 du premier et dernier noeuds
#                              (ce qui correspond a un cell_center regulier).
# bc_type = 0 : Pour la vitesse, suppose une valeur nulle sur les bords. neuman_homogene
#              2 : Suppose un gradient normal nul a la paroi, par exemple pour la pression
#              1 : Ne suppose rien, decentre interieur... (ordre 1)
def cell_to_cell_second_gradient(v, dz, bc_type):
    # Formule centree (ordre 2) pour pas cste dans le domaine :
    dz2 = dz * dz
    grad2 = np.zeros(len(v))
    # L'erreur est d'ordre 2
    grad2[1:-1] = (v[:-2] - 2 * v[1:-1] + v[2:]) / dz2
    # Aux bords :
    # Au bord en z- :
    Uc = v[0]
    Up = v[1]
    if bc_type == 0:
        Um = 0.0
        # L'erreur est d'ordre 1
        grad2[0] = 4.0 / 3.0 * (2 * Um - 3 * Uc + Up) / dz2
    elif bc_type == 1:
        Upp = v[2]
        # L'erreur est d'ordre 1
        grad2[0] = (-2 * Up + Uc + Upp) / (2 * dz2)
    else:
        Um = Uc + (Uc - Up) / 8.0
        # L'erreur est d'ordre 1, mais plus grande qu'avec bc_type=0
        grad2[0] = 4.0 / 3.0 * (2 * Um - 3 * Uc + Up) / dz2
    #
    # Au bord en z+ :
    Uc = v[-1]
    Um = v[-2]
    if bc_type == 0:
        Up = 0.0
        grad2[-1] = 4.0 / 3.0 * (Um - 3 * Uc + 2 * Up) / dz2
    elif bc_type == 1:
        Umm = v[-3]
        grad2[-1] = (-2 * Um + Uc + Umm) / (2 * dz2)
    else:
        Up = Uc + (Uc - Um) / 8.0
        grad2[-1] = 4.0 / 3.0 * (Um - 3 * Uc + 2 * Up) / dz2
    #
    return grad2


def test_cell_to_cell_gradient(z, Lz):
    # Test d'une fonction nulle aux bords :
    u = np.array([np.sin(2.0 * np.pi * zz / Lz) for zz in z])
    nz = len(u)
    dz = Lz / nz
    dudz = cell_to_cell_gradient(u, dz, 0)  # bc_type = champ nul au bord
    dudz_ana = np.array([2.0 * np.pi / Lz * np.cos(2.0 * np.pi * zz / Lz) for zz in z])
    err_max = np.fabs(dudz - dudz_ana).max()
    tol = 2 * np.pi * np.pi / (nz * nz * dz) * 2 * np.pi / nz
    if err_max > tol:
        raise Exception(
            "Failing function u test_deriv error = %g > %g = tol" % (err_max, tol)
        )
    print("Error max u : ", err_max, "(tol=", tol, ")")
    #
    # Test d'une fonction de derivee normale nulle aux bords :
    w = np.array([np.cos(2.0 * np.pi * zz / Lz) + 10 for zz in z])
    nz = len(w)
    dz = Lz / nz
    dwdz = cell_to_cell_gradient(
        w, dz, 2
    )  # bc_type = derivee normale du champ nulle au bord
    dwdz_ana = np.array([-2.0 * np.pi / Lz * np.sin(2.0 * np.pi * zz / Lz) for zz in z])
    err_max = np.fabs(dwdz - dwdz_ana).max()
    tol = 2 * np.pi * np.pi / (nz * nz * dz) * 2 * np.pi / nz
    if err_max > tol:
        raise Exception(
            "Failing function w test_deriv error = %g > %g = tol" % (err_max, tol)
        )
    print("Error max w : ", err_max, "(tol=", tol, ")")
    #
    # Test d'une fonction sans propriete aux bords :
    v = np.array([np.sin(2.0 * np.pi * zz / Lz) + 10 for zz in z])
    nz = len(v)
    dz = Lz / nz
    dvdz = cell_to_cell_gradient(v, dz, 1)  # bc_type = quelconque
    dvdz_ana = np.array([2.0 * np.pi / Lz * np.cos(2.0 * np.pi * zz / Lz) for zz in z])
    err_max = np.fabs(dvdz - dvdz_ana).max()
    tol = 2 * np.pi * np.pi / (nz * nz * dz)
    if err_max > tol:
        raise Exception(
            "Failing function v test_deriv error = %g > %g = tol" % (err_max, tol)
        )
    print("Error max v : ", err_max, "(tol=", tol, ")")
    #
    return True


def test_cell_to_cell_second_gradient(z, Lz):
    # Test d'une fonction nulle aux bords :
    u = np.array([np.sin(2.0 * np.pi * zz / Lz) + zz * (Lz - zz) * 1000.0 for zz in z])
    nz = len(u)
    dz = Lz / nz
    dduddz = cell_to_cell_second_gradient(u, dz, 0)  # bc_type = champ nul au bord
    dduddz_ana = np.array(
        [
            -((2.0 * np.pi / Lz) ** 2) * np.sin(2.0 * np.pi * zz / Lz) - 2000.0
            for zz in z
        ]
    )
    err_max = np.fabs(dduddz - dduddz_ana).max()
    tol = 1.0 / 6.0 * dz * (2 * np.pi / Lz) ** 3
    if err_max > tol:
        raise Exception(
            "Failing function u test_second_gradient error = %g > %g = tol"
            % (err_max, tol)
        )
    print("Error max u : ", err_max, "(tol=", tol, ")")
    #
    # Test d'une fonction de derivee normale nulle aux bords :
    w = np.array([np.cos(2.0 * np.pi * zz / Lz) + 1000 for zz in z])
    nz = len(w)
    dz = Lz / nz
    ddwddz = cell_to_cell_second_gradient(
        w, dz, 2
    )  # bc_type = derivee normale du champ nulle au bord
    ddwddz_ana = np.array(
        [-((2.0 * np.pi / Lz) ** 2) * np.cos(2.0 * np.pi * zz / Lz) for zz in z]
    )
    err_max = np.fabs(ddwddz - ddwddz_ana).max()
    tol = (1.0 / 6.0 - 1.0 / 8.0) * dz * (2 * np.pi / Lz) ** 3
    if err_max > tol:
        raise Exception(
            "Failing function w test_second_gradient error = %g > %g = tol"
            % (err_max, tol)
        )
    print("Error max w : ", err_max, "(tol=", tol, ")")
    #
    # Test d'une fonction sans propriete aux bords :
    v = np.array([np.sin(2.0 * np.pi * zz / Lz) + 10 + 1000.0 * zz for zz in z])
    nz = len(v)
    dz = Lz / nz
    ddvddz = cell_to_cell_second_gradient(v, dz, 1)  # bc_type = quelconque
    ddvddz_ana = np.array(
        [-((2.0 * np.pi / Lz) ** 2) * np.sin(2.0 * np.pi * zz / Lz) for zz in z]
    )
    err_max = np.fabs(ddvddz - ddvddz_ana).max()
    tol = 1.0 / 2.0 * dz * (2 * np.pi / Lz) ** 3
    if err_max > tol:
        raise Exception(
            "Failing function v test_second_gradient error = %g > %g = tol"
            % (err_max, tol)
        )
    print("Error max v : ", err_max, "(tol=", tol, ")")
    #
    return True


# v : NumPy np.array
# Derivee seconde (Laplacien)
# On suppose un pas uniforme
def deriv2(v, z, Lz):
    raise Exception(
        "Tentative de derivee seconde pour maillage a pas variable non validee!"
    )
    xf = (z[1:] + z[:-1]) / 2.0  # Centre des faces
    xf = np.append(np.append(0.0, xf), Lz)
    dx = xf[1:] - xf[:-1]  # dx_i = x_i+1 - x_i
    # L'evaluation du dx entre 2 cellules :
    h = (dx[:-1] + dx[1:]) / 2.0  # contient (Nsommets - 2) valeurs = Ncellules - 1
    if h.max() - h.min() > 0.01 * h.min():
        raise Exception("Bad use of deriv2 : varying step size!")
    #
    # On prend le h moyen
    h1 = h.mean()
    u_pl = v[2:]
    u_m = v[:-2]
    u_c = v[1:-1]
    # Formule centree (ordre 2) pour pas cste dans le domaine :
    laplacien = np.zeros(len(v))
    laplacien[1:-1] = (u_pl - 2 * u_c + u_m) / (h1 * h1)
    # Formule decentree (ordre 1) pour pas cste sur le bord gauche :
    u_g = v[0]
    u_pl = v[1]
    u_pl2 = v[2]
    laplacien[0] = (u_g - 2 * u_pl + u_pl2) / (h1 * h1)
    # Formule decentree (ordre 1) pour pas cste sur le bord droit :
    u_d = v[-1]
    u_m = v[-2]
    u_m2 = v[-3]
    laplacien[-1] = (u_d - 2 * u_m + u_m2) / (h1 * h1)
    return laplacien


def integrate(vec, dz):
    nz = len(vec)
    integral = np.zeros(nz)
    old = 0
    for i, val in enumerate(vec):
        old += val * dz
        integral[i] = old
    return integral


def getParam(fic, name, vec=False, compo=0, string=False):
    f = open(fic)
    lines = f.readlines()
    f.close()
    for i, line in enumerate(lines):
        if line.find("#") >= 0:
            lp = line.find("#")
            rp = line.rfind("#")
            if rp > lp:
                s = line[: lp - 1] + line[rp + 1 :]
                # print(lines[i], "replaced by ", s)
                lines[i] = s
    goodline = [line for line in lines if name in line.split()]
    # print(goodline)
    # Suppression des lignes contenant un commentaire :
    # goodline=[line for line in goodline if "#" not in line]
    # print(goodline)
    if len(goodline) == 1:
        if vec:
            val = goodline[0].split()[2 + compo]
        else:
            val = goodline[0].split()[1]
        try:
            return float(val)
        except:
            if not string:
                raise Exception("Trying to return a string? %s" % val)
            return val
    else:
        raise Exception(
            "ohoh, looking for name=%s commented line: %s" % (name, goodline)
        )


def get_thermal_param(jdd, name, ind=0):
    """
    Gets the thermal param in file ``jdd`` corresponding to ``name`` and ``ind``

    Args:
        jdd (str): file location
        name (str): param name
        ind (int): the thermal field number

    Returns:
        str: the param corresponding to ``ind`` and ``name``

    """
    f = open(jdd)
    text = f.read()
    f.close()
    text = re.sub(r"#.*#", r"", text)
    deb = text.find("thermique")
    if deb < 0:
        raise(Exception("Il n'a a pas de bloc thermique dans %s" % jdd))
    text = text[deb + 9 :]
    th_field = get_inside_brackets(text)[ind]
    param = th_field[th_field.find(name) + len(name) + 1 :]
    param = param[: param.find("\n")]
    return param


def get_inside_brackets(st):
    """
    Separates the different fields inside a set a brackets

    Args:
        st (str):

    Returns:
        list:

    Examples:
        >>>get_inside_brackets('{ {a}, bla {b {}} }')
        ['{a}', '{b {}}']

    """
    st = st[st.find("{") + 1 :]
    list_champs = []
    field = ""
    opened = 1
    ind = 0
    while opened > 0:
        ch = st[ind]
        ind += 1
        if ch == "{":
            opened += 1
        if opened > 1:
            field += ch
        if ch == "}":
            opened -= 1
            if opened == 1:
                list_champs.append(field)
                field = ""

    return list_champs


def get_thermal_prop(jdd="DNS.sauv", ind=0):
    ll = float(get_thermal_param(jdd, "lambda_liquid", ind))
    lv = float(get_thermal_param(jdd, "lambda_vapor", ind))
    cpl = float(get_thermal_param(jdd, "cp_liquid", ind))
    cpv = float(get_thermal_param(jdd, "cp_vapor", ind))
    try:
        qw = (
            float(get_thermal_param(jdd, "flux_impose_kmin", ind))
            + float(get_thermal_param(jdd, "flux_impose_kmax", ind))
        ) / 2
    except:
        qw = "Not fixed"
    return ll, lv, cpl, cpv, qw


def get_prop_all(jdd="DNS.data", repr_file=None, Ret=0.0, nb=57):
    Lx = getParam(jdd, "uniform_domain_size_i")
    Ly = getParam(jdd, "uniform_domain_size_j")
    Lz = getParam(jdd, "uniform_domain_size_k")
    mul = getParam(jdd, "mu_liquide")
    rhol = getParam(jdd, "rho_liquide")
    nx = getParam(jdd, "nbelem_i")
    ny = getParam(jdd, "nbelem_j")
    nz = getParam(jdd, "nbelem_k")
    try:
        muv = getParam(jdd, "mu_vapeur")
        rhov = getParam(jdd, "rho_vapeur")
    except:
        muv = mul
        rhov = rhol
    try:
        vb = getParam(jdd, "vol_bulle_monodisperse")
        # TODO : chopper le nombre de bulles (depuis le data ou depuis les out...)
        # nb = getParam(repr_file, 'bubble_groups')
        sigma = getParam(jdd, "sigma")
    except:
        print("Bubble assumed volume 0")
        vb = 0.0
        nb = 0.0
        sigma = 1.0
    try:
        g = getParam(jdd, "gravite", vec=True, compo=0)
    except:
        g = 0

    db = pow(vb * 6.0 / np.pi, 1.0 / 3.0)
    alv = nb * vb / (Lx * Ly * Lz)
    h = Lz / 2.0
    rhom = alv * rhov + (1 - alv) * rhol
    utau = Ret * mul / (rhol * h)
    tauw = rhol * utau * utau
    Eo = rhol * abs(g) * db**2 / sigma
    beta = tauw / h - rhom * g
    try:
        str_Svx = getParam(jdd, "expression_variable_source_x")
        print("Variable source read as: ", str_Svx)
        Svx = float(str_Svx)
    except:
        # The default value is 0.
        Svx = 0.0
    #
    verb = 1
    if verb:
        print("rhol= ", rhol)
        print("rhov= ", rhov)
        print("mul= ", mul)
        print("muv= ", muv)
        print("vb= ", vb)
        print("nb= ", nb)
        print("sigma= ", sigma)
        print("g= ", g)
        print("alv= ", alv)
        print("Ret= ", Ret)
        print("rhom ", rhom)
        print("beta ", beta)
        print("Svx ", Svx)
        print("tauw ", tauw)
        print("Db : ", db)
        print("Eo : ", Eo)
        print("nx : ", nx)
        print("ny : ", ny)
        print("nz : ", nz)
    return (
        rhol,
        rhov,
        mul,
        muv,
        alv,
        beta,
        sigma,
        Lz,
        rhom,
        Eo,
        g,
        Svx,
        vb,
        nb,
        tauw,
        db,
        Lx,
        Ly,
        nx,
        ny,
        nz,
    )


class Simu:
    def __init__(self, jdd=None, repr_file=None, Ret=0.0, thermal=None, dynamique=None):
        """
        | One class to rule them all,
        | one class to find them,
        | One class to bring them all
        | and in the darkness bind them

        This class possess some usefull arguments :

        * :attr:`dim_prop` : the list of dimensional non thermal properties
        * :attr:`thermal_prop` : the list of dimensional thermal properties
        * :attr:`adim_prop` : the list of addimentional properties

        Args:
            jdd (str): name of the data file
            repr_file (str): name of the sauv file
            Ret: the value of :math:`Re_\\tau`
            thermal (bool): if loading thermal fields
            indicatrice (bool): if loading indicatrice field
            velocity (bool): if loading velocity field
            pressure (bool): if loading pressure field

        """
        self.ret = Ret
        self.thermal = thermal
        self.jdd = jdd
        self.repr_file = repr_file
        self.je_suis_une_simu = True

        self.Ret = Field(self.ret, name=r"Re_{\tau}")
        # TODO: intégrer get_prop_all comme une méthode de classe, et réucpérer aussi les dimensions du maillage (nx, ny, nz)
        if jdd is not None:
            (
                rhol,
                rhov,
                mul,
                muv,
                alv,
                beta,
                sigma,
                Lz,
                rhom,
                Eo,
                g,
                Svx,
                vb,
                nb,
                tauw,
                Db,
                Lx,
                Ly,
                nx,
                ny,
                nz,
            ) = get_prop_all(jdd, repr_file, Ret)
        else:
            (
                rhol,
                rhov,
                mul,
                muv,
                alv,
                beta,
                sigma,
                Lz,
                rhom,
                Eo,
                g,
                Svx,
                vb,
                nb,
                tauw,
                Db,
                Lx,
                Ly,
                nx,
                ny,
                nz,
            ) = [np.nan] * 21
        (
            self.rhol,
            self.rhov,
            self.mul,
            self.muv,
            self.alv,
            self.beta,
            self.sigma,
            self.Lz,
            self.rhom,
            self.Eo,
            self.g,
            self.Svx,
            self.vb,
            self.nb,
            self.tauw,
            self.Db,
            self.Lx,
            self.Ly,
            self.nx,
            self.ny,
            self.nz,
        ) = (
            Field(rhol, name=r"\rho_l"),
            Field(rhov, name=r"\rho_v"),
            Field(mul, name=r"\mu_l"),
            Field(muv, name=r"\mu_v"),
            Field(alv, name=r"\alpha_v"),
            Field(beta, name=r"\beta"),
            Field(sigma, name=r"\sigma"),
            Field(Lz, name=r"2\delta"),
            Field(rhom, name=r"\rho_m"),
            Field(Eo, name=r"Eo"),
            Field(g, name=r"g"),
            Field(Svx, name=r"Sv_x"),
            Field(vb, name=r"V_b"),
            Field(nb, name=r"N_b"),
            Field(tauw, name=r"\tau_w"),
            Field(Db, name=r"D_b"),
            Field(Lx, name=r"L_x"),
            Field(Ly, name=r"L_y"),
            Field(nx, name=r"n_x"),
            Field(ny, name=r"n_y"),
            Field(nz, name=r"n_z"),
        )

        # Attention, si tout n'est pas défini au début ces valeurs ne seront pas bonnes !
        self.thermique = False

        if thermal is not None:
            self.thermique = True
            lambdal, lambdav, cpl, cpv, qw = get_thermal_prop(repr_file, thermal)
            self.lambdal, self.lambdav, self.cpl, self.cpv, self.qw = (
                Field(lambdal, name=r"\lambda_l"),
                Field(lambdav, name=r"\lambda_v"),
                Field(cpl, name=r"Cp_l"),
                Field(cpv, name=r"Cp_v"),
                Field(qw, name=r"q_w"),
            )
            self.thermal_prop = [
                self.lambdav,
                self.lambdal,
                self.cpv,
                self.cpl,
                self.qw,
            ]

            # Adim thermal prop
            self.Prl = self.mul * self.cpl / self.lambdal
            self.Prl.name = r"Pr_l"
            self.Prv = self.muv * self.cpv / self.lambdav
            self.Prv.name = r"Pr_v"
            self.rhocpfrac = (self.rhol * self.cpl) / (self.rhov * self.cpv)
            self.Effl = np.sqrt(self.rhol * self.cpl * self.lambdal)
            self.Effv = np.sqrt(self.rhov * self.cpv * self.lambdav)
            self.Diffl = self.lambdal / (self.rhol * self.cpl)
            self.Diffv = self.lambdav / (self.rhov * self.cpv)
            self.Effl.name = r"\textrm{Eff}_l"
            self.Effv.name = r"\textrm{Eff}_v"
            self.Diffl.name = r"\textrm{Diff}_l"
            self.Diffv.name = r"\textrm{Diff}_v"

        # MODIF GAB
        self.PostTrtVersion = "0"  # cf beau_final
        self.it0 = 0
        self.itf = 1
        self.listeTemps = []
        self.jdd = jdd
        if not (dynamique == None):
            self.champ = Field(dynamique, name=r"champ")

        # FIN MODIF GAB

        # la liste des fichiers utilisés (lus) dans la simulation en plus des fichiers .data et .sauv comme les
        # fichiers dt_ev
        self.__TRUST_FILES_CACHE = {}

    @property
    def all(self):
        return Field(1 - self.alv, name=r"\alpha_l")

    @property
    def dim_prop(self):
        return [
            self.rhol,
            self.rhov,
            self.mul,
            self.muv,
            self.all,
            self.alv,
            self.beta,
            self.sigma,
            self.Lz,
            self.rhom,
            self.g,
            self.Svx,
            self.vb,
            self.nb,
            self.tauw,
            self.Db,
        ]

    @property
    def rholrhov(self):
        return self.rhol / self.rhov

    @property
    def mulmuv(self):
        return self.mul / self.muv

    @property
    def adim_prop(self):
        if self.thermique:
            return [
                self.Ret,
                self.Eo,
                self.rholrhov,
                self.mulmuv,
                self.Prl,
                self.Prv,
                self.rhocpfrac,
                self.Diffl,
                self.Diffv,
                self.Effl,
                self.Effv,
            ]
        else:
            return [self.Ret, self.Eo, self.rholrhov, self.mulmuv]

    def getObject(self, filename):
        fileObj = self.__TRUST_FILES_CACHE.get(filename, None)
        if fileObj is None:
            fileObj = DTEVFile(filename, None)
            self.__TRUST_FILES_CACHE[filename] = fileObj
            pass
        return fileObj

    def getEntries(self, filename):
        fObj = self.getObject(filename)
        dvar = fObj.getEntries()
        # print("Available variables are : ", dvar)
        return dvar

    # TODO: à supprimer car non utilisée
    def printEntries(self, filename):
        dvar = self.getEntries(filename)
        print("Available variables are : ", dvar)
        return dvar

    @classmethod
    def get_from_json(cls, file_path):
        with open(file_path, "r") as fc:
            dic = json.load(fc, object_hook=simu_decoder)
        simu = Simu.get_from_dict(dic)
        return simu

    @classmethod
    def get_from_dict(cls, dic: dict):
        simu = Simu()
        for key in dic.keys():
            simu.__setattr__(key, dic[key])
        return simu

    def get_fields_from_stat(self, fstat):
        fstat = fstat.replace("_temperature_{}", "")
        verb = 0
        if verb:
            print("Load file %s waiting ....................." % fstat)
        fileObj = self.getObject(fstat)
        entries = fileObj.getEntries()
        if verb:
            print("end.")
        if not os.path.isfile(fstat):
            raise Exception("File %s missing to load stats" % fstat)
        self.I = Field.LoadFromFile(
            fileObj, ["coordonnee_K"], [r"I"], r"\chi", "z", None
        )
        Iinv = 1.0 / self.I
        Iv = 1 - self.I
        Ivinv = 1.0 / Iv
        # self.I.print()
        self.G = Field.initgravity([self.g, 0, 0], "g", self.I)
        self.UI = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            ["UI", "VI", "WI"],
            r"\overline{\mathbf{U}_l \chi}",
            r"z",
            None,
        )

        self.UI_v = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            ["UIv", "VIv", "WIv"],
            r"\mathbf{U}_v \chi_v",
            "z",
            None,
        )
        self.pression_l = (
            Field.LoadFromFile(
                fileObj, ["coordonnee_K"], ["PI"], r"\overline{P \chi}", "z", None
            )
            * Iinv
        )
        self.pression_v = (
            Field.LoadFromFile(
                fileObj,
                ["coordonnee_K"],
                ["PIv"],
                r"\overline{P (1 - \chi)}",
                "z",
                None,
            )
            * Ivinv
        )
        self.pression_l.settexname(r"\overline{P}^l")
        self.pression_v.settexname(r"\overline{P}^v")
        self.ai = Field.LoadFromFile(
            fileObj, ["coordonnee_K"], ["AI"], r"\overline{a_i}", "z", None
        )
        if ("PaiNx" in entries) and ("PaiNy" in entries) and ("PaiNz" in entries):
            self.PaiN = Field.LoadFromFile(
                fileObj,
                ["coordonnee_K"],
                ["PaiNx", "PaiNy", "PaiNz"],
                r"f_interf",
                r"$z$",
                r"$\sigma n\kappa\delta^i$",
            )  # GB 2022 : je doute du nom mis par Antoine...
        if "P_LIQ_I" in entries:
            self.pression_l_ext = Field.LoadFromFile(
                fileObj, ["coordonnee_K"], ["P_LIQ_I"], "P_LIQ_I", "z", None
            )
        if "P_VAP_Iv" in entries:
            self.pression_v_ext = Field.LoadFromFile(
                fileObj, ["coordonnee_K"], ["P_VAP_Iv"], "P_VAP_Iv", "z", None
            )
        if (
            ("P_LIQ_aiNx" in entries)
            and ("P_LIQ_aiNy" in entries)
            and ("P_LIQ_aiNz" in entries)
        ):
            self.p_l_ext = Field.LoadFromFile(
                fileObj,
                ["coordonnee_K"],
                ["P_LIQ_aiNx", "P_LIQ_aiNy", "P_LIQ_aiNz"],
                "p_l",
                r"$z$",
                r"$p_l\nabla\chi_l$",
            )
        if (
            ("P_VAP_aiNx" in entries)
            and ("P_VAP_aiNy" in entries)
            and ("P_VAP_aiNz" in entries)
        ):
            self.p_v_ext = Field.LoadFromFile(
                fileObj,
                ["coordonnee_K"],
                ["P_VAP_aiNx", "P_VAP_aiNy", "P_VAP_aiNz"],
                "p_v",
                r"$z$",
                r"$p_v\nabla\chi_v$",
            )
        self.uu_l = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            ["UUI", "UVI", "UWI", "UVI", "VVI", "VWI", "UWI", "VWI", "WWI"],
            r"\overline{\mathbf{U} \mathbf{U} \chi}",
            "z",
            None,
        )
        self.uu_v = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            ["UUIv", "UVIv", "UWIv", "UVIv", "VVIv", "VWIv", "UWIv", "VWIv", "WWIv"],
            "UU_v",
            "z",
            r"\langle(u\times u)\rangle)_v",
        )
        self.aii = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            ["kaiNx", "kaiNy", "kaiNz"],
            "f_interf",
            "z",
            r"$\kappa n\delta^i$",
        )
        self.kaiN = self.aii
        self.aiN = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            ["aiNx", "aiNy", "aiNz"],
            "f_interf",
            "z",
            r"$n\delta^i$",
        )
        self.vdvdxaiN = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            [
                "UdUdxaiNx",
                "UdVdxaiNx",
                "UdWdxaiNx",
                "UdVdxaiNx",
                "VdVdxaiNx",
                "VdWdxaiNx",
                "UdWdxaiNx",
                "VdWdxaiNx",
                "WdWdxaiNx",
            ],
            "vdvdxaiN",
            "z",
            r"$\langle(u\times u)\rangle)_l$",
        )
        self.vdvdyaiN = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            [
                "UdUdyaiNy",
                "UdVdyaiNy",
                "UdWdyaiNy",
                "UdVdyaiNy",
                "VdVdyaiNy",
                "VdWdyaiNy",
                "UdWdyaiNy",
                "VdWdyaiNy",
                "WdWdyaiNy",
            ],
            "UU_l",
            "z",
            r"$\langle(u\times u)\rangle)_l$",
        )
        self.vdvdzaiN = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            [
                "UdUdzaiNz",
                "UdVdzaiNz",
                "UdWdzaiNz",
                "UdVdzaiNz",
                "VdVdzaiNz",
                "VdWdzaiNz",
                "UdWdzaiNz",
                "VdWdzaiNz",
                "WdWdzaiNz",
            ],
            "UU_l",
            "z",
            r"$\langle(u\times u)\rangle)_l$",
        )
        self.IdU = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            [
                "IdUdx",
                "IdVdx",
                "IdWdx",
                "IdUdy",
                "IdVdy",
                "IdWdy",
                "IdUdz",
                "IdVdz",
                "IdWdz",
            ],
            "UU_l",
            "z",
            r"$ \nabla\cdot u_l$",
        )
        if (
            ("dUdxaiNx" in entries)
            and ("dUdyaiNx" in entries)
            and ("dUdzaiNx" in entries)
        ):
            self.duaiNx = Field.LoadFromFile(
                fileObj,
                ["coordonnee_K"],
                ["dUdxaiNx", "dUdyaiNx", "dUdzaiNx"],
                "UU_l",
                "z",
                r"$\nabla\cdot\u_l$",
            )
        if (
            ("dUdxaiNy" in entries)
            and ("dUdyaiNy" in entries)
            and ("dUdzaiNy" in entries)
        ):
            self.dvaiNy = Field.LoadFromFile(
                fileObj,
                ["coordonnee_K"],
                ["dVdxaiNy", "dVdyaiNy", "dVdzaiNy"],
                "UU_l",
                "z",
                r"$\nabla\cdot\v_l$",
            )
        if (
            ("dUdxaiNz" in entries)
            and ("dUdyaiNz" in entries)
            and ("dUdzaiNz" in entries)
        ):
            self.dwaiNz = Field.LoadFromFile(
                fileObj,
                ["coordonnee_K"],
                ["dWdxaiNz", "dWdyaiNz", "dWdzaiNz"],
                "UU_l",
                "z",
                r"$\nabla\cdot\w_l$",
            )
        self.dvdxaiN = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            ["dUdxaiNx", "dVdxaiNx", "dWdxaiNx"],
            "UU_l",
            "z",
            r"$\langle(u\times u)\rangle)_l$",
        )
        self.dvdyaiN = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            ["dUdyaiNy", "dVdyaiNy", "dWdyaiNy"],
            "UU_l",
            "z",
            r"$\langle(u\times u)\rangle)_l$",
        )
        self.dvdzaiN = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            ["dUdzaiNz", "dVdzaiNz", "dWdzaiNz"],
            "UU_l",
            "z",
            r"$\langle(u\times u)\rangle)_l$",
        )
        self.UPI = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            ["UPI", "VPI", "WPI"],
            r"UpI",
            r"$z$",
            r"$\alpha_l up_l$",
        )
        self.dudxdudx_l = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            [
                "IdUdxdUdx",
                "IdUdxdVdx",
                "IdUdxdWdx",
                "IdUdxdVdx",
                "IdVdxdVdx",
                "IdVdxdWdx",
                "IdUdxdWdx",
                "IdVdxdWdx",
                "IdWdxdWdx",
            ],
            "UU_l",
            "z",
            r"$\langle(u\times u)\rangle)_l$",
        )
        self.dudydudy_l = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            [
                "IdUdydUdy",
                "IdUdydVdy",
                "IdUdydWdy",
                "IdUdydVdy",
                "IdVdydVdy",
                "IdVdydWdy",
                "IdUdydWdy",
                "IdVdydWdy",
                "IdWdydWdy",
            ],
            "UU_l",
            "z",
            r"$\langle(u\times u)\rangle)_l$",
        )
        self.dudzdudz_l = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            [
                "IdUdzdUdz",
                "IdUdzdVdz",
                "IdUdzdWdz",
                "IdUdzdVdz",
                "IdVdzdVdz",
                "IdVdzdWdz",
                "IdUdzdWdz",
                "IdVdzdWdz",
                "IdWdzdWdz",
            ],
            "UU_l",
            "z",
            r"$\langle(u\times u)\rangle)_l$",
        )
        self.uuu_l = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            [
                "UUUI",
                "UUVI",
                "UUWI",
                "UUVI",
                "UVVI",
                "UVWI",
                "UUWI",
                "UVWI",
                "UWWI",
                "UUVI",
                "UVVI",
                "UVWI",
                "UVVI",
                "VVVI",
                "VVWI",
                "UVWI",
                "VVWI",
                "VWWI",
                "UUWI",
                "UVWI",
                "UWWI",
                "UVWI",
                "VVWI",
                "VWWI",
                "UWWI",
                "VWWI",
                "WWWI",
            ],
            "UUU",
            "z",
            r"\langle(u\times u\times u)\rangle)_v",
        )
        self.pdu_l = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            [
                "IPdUdx",
                "IPdUdy",
                "IPdUdz",
                "IPdVdx",
                "IPdVdy",
                "IPdVdz",
                "IPdWdx",
                "IPdWdy",
                "IPdWdz",
            ],
            "UU_l",
            "z",
            r"$\langle(u\times u)\rangle)_l$",
        )
        self.vaiN = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            [
                "UaiNx",
                "UaiNy",
                "UaiNz",
                "VaiNx",
                "VaiNy",
                "VaiNz",
                "WaiNx",
                "WaiNy",
                "WaiNz",
            ],
            "UU_l",
            "z",
            r"$\langle(u\times u)\rangle)_l$",
        )
        self.vPaiN = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            [
                "UPaiNx",
                "UPaiNy",
                "UPaiNz",
                "VPaiNx",
                "VPaiNy",
                "VPaiNz",
                "WPaiNx",
                "WPaiNy",
                "WPaiNz",
            ],
            "UU_l",
            "z",
            r"$\langle(u\times u)\rangle)_l$",
        )
        self.vvaiN = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            [
                "UUaiNx",
                "UUaiNy",
                "UUaiNz",
                "UVaiNx",
                "UVaiNy",
                "UVaiNz",
                "UWaiNx",
                "UWaiNy",
                "UWaiNz",
                "UVaiNx",
                "UVaiNy",
                "UVaiNz",
                "VVaiNx",
                "VVaiNy",
                "VVaiNz",
                "VWaiNx",
                "VWaiNy",
                "VWaiNz",
                "UWaiNx",
                "UWaiNy",
                "UWaiNz",
                "VWaiNx",
                "VWaiNy",
                "VWaiNz",
                "WWaiNx",
                "WWaiNy",
                "WWaiNz",
            ],
            "UU_l",
            "z",
            r"$\langle(u\times u)\rangle)_l$",
        )

    # MODIF GAB

    def get_field_from_med(self, fmed, dual_Mesh=None, give_Mesh=False):

        """
        Loads the field, and its attached time from a given .med file

        Args:
            fmed i   : nom du fichier med à post traiter (avec chemin en dur)
            NomRef   : nom du champ tel qu'il est ecrit dans le .lata
            NomCompo : nom des composantes du champ tels qu'elles apparaîtrons dans les plots et
                       dans Paraview
            NomBeau  : nom du champ tels qu'elles apparaîtrons dans les plots et
                         dans Paraview

        Returns:
            En vrai on return rien, mais on modifie :
            self.Champ      : tantôt la vitesse, l'indicatrice, la pression, ...
            self.listeTemps : Liste des instants sauvegardés par le code

        """

        self.champ.listeTemps = self.champ.LoadFromMED(
            fmed, self.it0, self.itf, self.nb, dual_Mesh=None, give_Mesh=False
        )

    # A changer ca ...
    def get_BxMesh(self):
        """
        Récupères le maillage MED, autour de la bulle.
        (DIBox=[0,0,0], NomCompo=[""], NomChamp="".)

        Args:
           self.Champ,[self.listeTemps[0],0,0],Mesh,
                                             FieldName, jddName, DIBox, NomCompo, NomChamp, False,
                                             'perioX, perioY, perioZ'

        Returns:
            Le maillage MED. ATTN, ce n'est pas au format numpy.

        """
        # _,_,BxMesh = ob.BubbleFieldWithNumpyI(self.champ,[self.listeTemps[0],0,0],Mesh,
        #                                      FieldName, jddName, DIBox, NomCompo, NomChamp, False,
        #                                      'perioX, perioY, perioZ')

        # return(BxMesh.getMesh())
        pass

    # FIN MODIF GAB

    def get_fields_from_therm_stat(self, fstat, file_number):
        fstati = fstat.format(str(file_number))
        verb = 0
        if verb:
            print("Load file %s waiting ....................." % fstati)
        fileObj = self.getObject(fstati)
        if verb:
            print("end.")
        self.ITU = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            ["ITU", "ITV", "ITW"],
            r"\overline{T\mathbf{U}\chi}",
            "z",
            None,
        )
        try:
            self.ITbulles = Field.LoadFromFile(
                fileObj,
                ["coordonnee_K"],
                ["ITbulles"],
                r"\overline{T_{\textrm{adim " r"bulles}} \chi}",
                "z",
                None,
            )
        except:
            print("ITbulles not found")
        try:
            self.TTI = Field.LoadFromFile(
                fileObj, ["coordonnee_K"], ["TTI"], r"\overline{T T \chi}", r"z", None
            )
        except:
            print("TTI was not registered")
        self.TI = Field.LoadFromFile(
            fileObj, ["coordonnee_K"], ["TI"], r"\overline{T\chi}", "z", None
        )
        self.ITUU = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            ["ITUU", "ITUV", "ITUW", "ITUV", "ITVV", "ITVW", "ITUW", "ITVW", "ITWW"],
            r"\overline{T\mathbf{UU}\chi}",
            "z",
            None,
        )
        self.ITdP = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            ["ITdPdx", "ITdPdy", "ITdPdz"],
            "ITdP",
            "z",
            r"ITdP",
        )
        self.ITddU = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            [
                "ITddUdxdx",
                "ITddUdydy",
                "ITddUdzdz",
                "ITddVdxdx",
                "ITddVdydy",
                "ITddVdzdz",
                "ITddWdxdx",
                "ITddWdydy",
                "ITddWdzdz",
            ],
            "ITddU",
            "z",
            r"ITddU",
        )
        self.IUddT = Field.LoadFromFile(
            fileObj,
            ["coordonnee_K"],
            [
                "IUddTdxdx",
                "IUddTdydy",
                "IUddTdzdz",
                "IVddTdxdx",
                "IVddTdydy",
                "IVddTdzdz",
                "IWddTdxdx",
                "IWddTdydy",
                "IWddTdzdz",
            ],
            "IUddT",
            "z",
            r"IUddT",
        )

    def print_prop(self):
        """
        Todo: Remplacer par get_tex

        Returns:

        """
        display_latex(r"$\rho_l = %s $" % self.rhol, raw=True)
        display_latex(r"$\rho_v = %s $" % self.rhov, raw=True)
        display_latex(r"$\mu_l = %s $" % self.mul, raw=True)
        display_latex(r"$\mu_v = %s $" % self.muv, raw=True)
        display_latex(r"$v_b = %s $" % self.vb, raw=True)
        display_latex(r"$n_b = %s $" % self.nb, raw=True)
        display_latex(r"$\sigma = %s $" % self.sigma, raw=True)
        display_latex(r"$g = %s $" % self.g, raw=True)
        display_latex(r"$\alpha_v = %s $" % self.alv, raw=True)
        display_latex(r"$\alpha_l = %s $" % self.all, raw=True)
        display_latex(r"$\rho_m = %s $" % self.rhom, raw=True)
        display_latex(r"$\beta = %s $" % self.beta, raw=True)
        display_latex(r"$Sv_x = %s $" % self.Svx, raw=True)
        display_latex(r"$\tau_w = %s $" % self.tauw, raw=True)
        display_latex(r"$Db = %s $" % self.Db, raw=True)

    def print_prop_termique(self):
        display_latex(r"$\lambda_v = %s $" % self.lambdav, raw=True)
        display_latex(r"$\lambda_l = %s $" % self.lambdal, raw=True)
        display_latex(r"$Cp_v = %s $" % self.cpv, raw=True)
        display_latex(r"$Cp_l = %s $" % self.cpl, raw=True)
        display_latex(r"$q_w = %s $" % self.qw, raw=True)

    def print_prop_adim(self):
        display_latex(r"$Re_{\tau} = %s $" % self.Ret, raw=True)
        display_latex(r"$Eo = %s $" % self.Eo, raw=True)
        rholrhov = self.rhol / self.rhov
        display_latex(r"$\rho_l / \rho_v = %s $" % rholrhov, raw=True)
        mulmuv = self.mul / self.muv
        display_latex(r"$\mu_l / \mu_v = %s $" % mulmuv, raw=True)
        if self.thermique:
            display_latex(r"$Pr_l = %s $" % self.Prl, raw=True)
            rhocpfrac = (self.rhol * self.cpl) / (self.rhov * self.cpv)
            display_latex(
                r"${\rho_l Cp_l} / {\rho_v Cp_v} = %s $" % rhocpfrac, raw=True
            )


if __name__ == "__main__":
    f = open("../stat_file/monop.sauv")
    s = f.read()
    st = re.sub(r"#.*#", r"", s)
    st = st[st.find("thermique") + 9 :]
    test = get_inside_brackets(st)
    print(test)
