#Script de recuperation des proprietes physiques et de la geometrie

import os, sys, math

def readPropertiesData(nomFic):
    #initialisation
    properties = {}
    properties['mu'] = -1
    properties['rho_he'] = -1
    properties['lambda_he'] = -1
    properties['Cp_he'] = -1
    properties['beta_th'] = -1
    properties['rho_c'] = -1
    properties['lambda_c'] = -1
    properties['Cp_c'] = -1
    properties['Lx'] = -1
    properties['Ly'] = -1
    properties['Lz'] = -1
    properties['V'] = -1
    properties['nx'] = -1
    properties['nz'] = -1
    properties['DHY'] = -1
    # ouverture des fichiers
    fic = open(nomFic,'r')
#
#
#
    for ligne in fic:
        ligne = ligne.strip().lower()
        tLigne = ligne.split()
        if ligne.find('champ_uniforme')>-1:
            if ligne.startswith('mu'):
                properties['mu'] = float(tLigne[-1])
            elif ligne.startswith('rho champ_uniforme 1 4'):
                properties['rho_he'] = float(tLigne[-1])
            elif ligne.startswith('lambda champ_uniforme 1 0'):
                properties['lambda_he'] = float(tLigne[-1])
            elif ligne.startswith('cp champ_uniforme 1 5'):
                properties['Cp_he'] = float(tLigne[-1])
            elif ligne.startswith('beta_th'):
                properties['beta_th'] = float(tLigne[-1])
            elif ligne.startswith('rho champ_uniforme 1 2'):
                properties['rho_c'] = float(tLigne[-1])
            elif ligne.startswith('lambda champ_uniforme 1 6'):
                properties['lambda_c'] = float(tLigne[-1])
            elif ligne.startswith('cp champ_uniforme 1 6'):
                properties['Cp_c'] = float(tLigne[-1])
            elif ligne.startswith('modele_turbulence prandtl'):
                properties['DHY'] = float(tLigne[-3])
        elif ligne.find('longueurs')>-1:
            properties['Lx'] = float(tLigne[-3])
            properties['Ly'] = float(tLigne[-2])
            properties['Lz'] = float(tLigne[-1])
        elif ligne.find('frontiere_ouverte_vitesse_imposee champ_front_uniforme')>-1:
            properties['V'] = abs(float(tLigne[-1]))
        elif ligne.find('nombre_de_noeuds')>-1:
            properties['nx'] = float(tLigne[-3])
            properties['nz'] = float(tLigne[-1])

    fic.close()
    return properties


def ecritureFichier(properties):
    #ecriture du fichier
    nomFic = 'propertiesGeometry.dat'
    fichier = open(nomFic, 'w')

    fichier.write('%18.7f %18.3f %18.3f %18.3f %18.7f %18.3f %18.3f %18.3f %18.5f %18.5f %18.5f %18.3f %18.5f %18.5f %18.5f\n' % ( properties['mu'], properties['rho_he'], properties['lambda_he'], properties['Cp_he'], properties['beta_th'],  properties['rho_c'], properties['lambda_c'], properties['Cp_c'], properties['Lx'], properties['Ly'], properties['Lz'], properties['V'], properties['nx'], properties['nz'], properties['DHY'] ))
    fichier.close()

def getPropertiesFromdat():
    properties = {}
    nomFichier = 'propertiesGeometry.dat'
    if os.path.isfile(nomFichier):
        #recupere les donnees du fichier
        f = open(nomFichier, 'r')
        lignes = f.readlines()
        f.close()
    else:
        print('Erreur getPropertiesFromdat : fichier %s non trouve !' % (nomFichier))
        sys.exit()
    ligne = (lignes[0]).strip()
    tabLigne = ligne.split()
    ind = 0
    try:
        properties['mu'] = float(tabLigne[ind])
        ind += 1
        properties['rho_he'] = float(tabLigne[ind])
        ind += 1
        properties['lambda_he'] = float(tabLigne[ind])
        ind += 1
        properties['Cp_he'] = float(tabLigne[ind])
        ind += 1
        properties['beta_th'] = float(tabLigne[ind])
        ind += 1
        properties['rho_c'] = float(tabLigne[ind])
        ind += 1
        properties['lambda_c'] = float(tabLigne[ind])
        ind += 1
        properties['Cp_c'] = float(tabLigne[ind])
        ind += 1
        properties['Lx'] = float(tabLigne[ind])
        ind += 1
        properties['Ly'] = float(tabLigne[ind])
        ind += 1
        properties['Lz'] = float(tabLigne[ind])
        ind += 1
        properties['V'] = float(tabLigne[ind])
        ind += 1
        properties['nx'] = float(tabLigne[ind])
        ind += 1
        properties['nz'] = float(tabLigne[ind])
        ind += 1
        properties['DHY'] = float(tabLigne[ind])
    except IndexError:
        print('Erreur getPropertiesFromdat : lecture element %d pour 0-%d elements...' % (ind, len(tabLigne)-1))
        sys.exit()
    except ValueError:
        print('Erreur getPropertiesFromdat : lecture element %d n\'est pas un float (%s)...' % (ind, tabLigne[ind]))
        sys.exit()
    return properties



if __name__ == '__main__':

    #recuperation du fichier data
    import glob
    #derniere ligne du ls
    #ficLS = os.popen('ls *.data')
    #lignes = ficLS.readlines()
    #Ligne = lignes[0]
    #suppression du \n en fin de nom
    #nomFic = Ligne[:len(Ligne)-1]
    listFics = glob.glob('*.data')
    if len(listFics)>0:
        nomFic = listFics[0]
        properties = readPropertiesData(nomFic)

        #ecriture du fichier
        ecritureFichier(properties)
    else:
        print('Erreur propertiesGeometry : pas de fichier data trouve !')