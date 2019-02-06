#!/usr/bin/env python3

__author__ = "Veselin Kolev"
__license__ = "GPLv2"
__version__ = "2019020600"
__maintainer__ = __author__
__email__ = "vesso.kolev@gmail.com"
__status__ = "Production"


def read_residue_topology_dicts_as_pkl_files(argument):
    '''
    Reads the residue definitions dictionaries, already stored as Pickle \n
    storage containers.
    '''

    argument.map_residues = \
        argument.pickle.load(open(argument.pickle_map_residues_file, 'rb'))

    argument.map_sequence = \
        argument.pickle.load(open(argument.pickle_map_sequence_file, 'rb'))

    argument.sequence = \
        argument.pickle.load(open(argument.pickle_sequence_file, 'rb'))

    argument.pro_residues = \
        argument.pickle.load(open(argument.pickle_pro_residues_file, 'rb'))

    argument.excluded_pro_residues = \
        argument.pickle.load(open(argument.pickle_excluded_pro_residues_file, 'rb'))

    argument.axis_dict = \
        argument.pickle.load(open(argument.pickle_axis_dict_file, 'rb'))

    argument.branches_to_rotate_dict = \
        argument.pickle.load(open(argument.pickle_branches_to_rotate_dict_file, 'rb'))

    argument.nb_neigh_dict = \
        argument.pickle.load(open(argument.pickle_nb_neigh_dict_file, 'rb'))

    argument.dihedrals_dict = \
        argument.pickle.load(open(argument.pickle_dihedrals_dict_file, 'rb'))

    argument.bonds_dict = \
        argument.pickle.load(open(argument.pickle_bonds_dict_file, 'rb'))

    argument.bond_angles_dict = \
        argument.pickle.load(open(argument.pickle_bond_angles_dict_file, 'rb'))


    return None


if __name__ == "__main__":

    import sys

    class Dummy():
        pass

    import pickle

    import sqlite3

    import os

    import time

    S = time.time()

    argument = Dummy()

    argument.sys = sys

    argument.os = os

    argument.pickle = pickle

    argument.sqlite3 = sqlite3

    argument.pickle_storage_dir = '../pickle/'

    argument.pickle_map_residues_file = argument.pickle_storage_dir+\
        '/map_residues.pkl'
        
    argument.pickle_map_sequence_file = argument.pickle_storage_dir+\
        '/map_sequence.pkl'

    argument.pickle_sequence_file = argument.pickle_storage_dir+\
        '/sequence.pkl'

    argument.pickle_pro_residues_file = argument.pickle_storage_dir+\
        '/pro_residues.pkl'

    argument.pickle_excluded_pro_residues_file = argument.pickle_storage_dir+\
        '/excluded_pro_residues.pkl'

    argument.pickle_axis_dict_file = argument.pickle_storage_dir+\
        '/axes.pkl'

    argument.pickle_branches_to_rotate_dict_file = argument.pickle_storage_dir+\
        '/branches.pkl'

    argument.pickle_nb_neigh_dict_file = argument.pickle_storage_dir+\
        '/nonbonded.pkl'

    argument.pickle_dihedrals_dict_file = argument.pickle_storage_dir+\
        '/dihedrals.pkl'

    argument.pickle_bonds_dict_file = argument.pickle_storage_dir+\
        '/bonds.pkl'

    argument.pickle_bond_angles_dict_file = argument.pickle_storage_dir+\
        '/bond_angles.pkl'

    # Load the definitions:

    read_residue_topology_dicts_as_pkl_files(argument)

    for i in argument.pro_residues:

        if not i in argument.excluded_pro_residues:

            print('\n\n***********************************')
            print('\nResidue    : ',i)
            print('\n***********************************\n')

            aname = [j for j in argument.map_residues[i].keys()]

            arr = [argument.map_residues[i][j][0] for j in argument.map_residues[i]]

            str_ = ['   atom', '   atom', '      vdW', '         vdW',\
                    '        role', '   name', '   charge', '    sigma',\
                    '       epsilon']

            print('\n      Atoms in the resudue and their input parameters:')

            print('\n      ------------------------------------------------')
            print('   ' + ' '.join(str_[0:5]))
            print('   ' + ' '.join(str_[5:]))
            print('      ------------------------------------------------\n')

            counter = 0

            for j in set(arr):

               counter += 1

               str_ = ['  ' + '{:3d}'.format(counter)]

               str_ += [aname[arr.index(j)]]

               str_ += ['{:12.5e}'.format(argument.map_residues[i][aname[arr.index(j)]][2])]

               str_ += ['{:12.5e}'.format(argument.map_residues[i][aname[arr.index(j)]][4])]

               str_ += ['{:12.5e}'.format(argument.map_residues[i][aname[arr.index(j)]][5])]

               str_ += ['{:2d}'.format(argument.map_residues[i][aname[arr.index(j)]][-1])]

               print(' '.join(str_))

            print('\n   Bonds      : ',argument.bonds_dict[i])
            print('\n   Bond angles: ',argument.bond_angles_dict[i])
            print('\n   Axes       : ',argument.axis_dict[i])
            print('\n   Branches   : ',argument.branches_to_rotate_dict[i])
            print('\n   Dihedrals  : ',argument.dihedrals_dict[i])
            print('\n   NB pairs   : ',argument.nb_neigh_dict[i])


    print('\nTime for execution:', '{:8.5e}'.format(time.time()-S), 'seconds\n')

    sys.exit(0)
