#!/usr/bin/env python3

__author__ = "Veselin Kolev"
__license__ = "GPLv2"
__version__ = "2019020600"
__maintainer__ = __author__
__email__ = "vesso.kolev@gmail.com"
__status__ = "Production"

__about__ = '''
This program code creates those of the dictionaries and lists,
related to the process of mutation of the the side chain of
protein molecules described in PDB files, whose role is to
bring a ligand into the topology.
'''

def error_0001(argument, exit):

    '''
    Prints an error message (code 0001).
    '''

    str_ = '\nFATAL ERROR: The folder\n\n'
    str_ += argument.ligand_defs_folder
    str_ += '\n\nwhich is supposed to contain files with ligand '
    str_ += 'definitions does not exist, it is not accessible or '
    str_ += 'the specified folder path points to a file.\n'

    print(str_)

    if exit:

        argument.sys.exit(1)

    return None


def error_0002(argument, i, exit):

    '''
    Prints an error message (code 0002).
    '''

    str_ = '\nFATAL ERROR: The ligand part "'
    str_ += i
    str_ += '" of the file name\n\n'
    str_ += argument.ligand_defs_folder+i+argument.ligand_file_prefix
    str_ += '\n\ncontains more than 3 symbols. '
    str_ += 'That length is against the PDB specifications. Correct '
    str_ += 'the file name and run this program again!\n'

    print(str_)

    if exit:

        argument.sys.exit(1)

    return None


def error_0003(argument, i, exit):

    '''
    Prints an error message (code 0003).
    '''

    str_ = '\nFATAL ERROR: The ligand part "'
    str_ += i
    str_ += '" of the file name\n\n'
    str_ += argument.ligand_defs_folder+i+argument.ligand_file_prefix
    str_ += '\n\ncontains symbols which cannot be part of a '
    str_ += 'residue name. The set of allowed symbols in the residue '
    str_ += 'name is\n\n'
    str_ += argument.allowed_ligand_name_letters
    str_ += '\n\nCorrect the file name and run this program again!\n'

    print(str_)

    if exit:

        argument.sys.exit(1)

    return None


def error_0004(argument, i, line, exit):

    '''
    Prints an error message (code 0004).
    '''

    str_ = '\nFATAL ERROR: An error has been found during processing '
    str_ += 'the file\n\n'
    str_ += argument.ligand_defs_folder+i+argument.ligand_file_prefix
    str_ += '\n\nThe line:\n\n'
    str_ += line.rstrip()
    str_ += '\n\ndoes not match the expected format! '
    str_ += 'Every line in that file should contain an atom '
    str_ += 'description of four columns, separated by a space.\n'

    print(str_)

    if exit:

        argument.sys.exit(1)

    return None

def error_0005(argument, file_name, line, tmp, exit):

    '''
    Prints an error message (code 0005).
    '''

    str_ = '\nFATAL ERROR: A problem with one of atom name '
    str_ += 'definitions has been found when reading the ligand '
    str_ += ' atom definitions from the file\n\n'
    str_ += file_name
    str_ += '\n\nThe line:\n\n'
    str_ += line.rstrip()
    str_ += '\n\ndoes not provide properly defined atom name. '
    str_ += 'because the atom name detected there, "'
    str_ += tmp[0]
    str_ += '", contains more than 4 symbols!\n'

    print(str_)

    if exit:

        argument.sys.exit(1)

    return None


def error_0006(argument, file_name, line, tmp, alnl, exit):

    '''
    Prints an error message (code 0006).
    '''

    str_ = '\nFATAL ERROR: A problem with one of the atom name definitions'
    str_ += 'has been detected when reading the ligand atom definitions '
    str_ += 'from the file\n\n'
    str_ += file_name
    str_ += '\n\nThe line:\n\n'
    str_ += line.rstrip()
    str_ += '\n\ndoes not provide properly defined atom name because '
    str_ += 'the atom name string detected there, "'
    str_ += tmp[0]
    str_ += '", contains symbols which are not included in the list '
    str_ += 'of allowed symbols:\n\n'
    str_ += alnl
    str_ += '\n\n'

    print(str_)

    if exit:

        argument.sys.exit(1)

    return None


def error_0007(argument, file_name, line, tmp, exit):

    '''
    Prints an error message (code 0007).
    '''

    str_ = '\nFATAL ERROR: A problem with one of the atom name definitions'
    str_ += 'has been detected when reading the ligand atom definitions '
    str_ += 'from the file\n\n'
    str_ += file_name
    str_ += '\n\nThe line:\n\n'
    str_ += line.rstrip()
    str_ += '\n\ndoes not provide properly defined atomic charge for the '
    str_ += 'atom ' + tmp[0] + '. The atomic charge record is expected '
    str_ += 'to appear as floating point number, but the supplied one\n\n'
    str_ += str(tmp[1]) + '\n\ndoes not match the floating point number ' 
    str_ += 'format.\n'

    print(str_)

    if exit:

        argument.sys.exit(1)

    return None


def error_0008(argument, file_name, line, tmp, exit):

    '''
    Prints an error message (code 0008).
    '''

    str_ = '\nFATAL ERROR: A problem with one of the atom name definitions'
    str_ += 'has been detected when reading the ligand atom definitions '
    str_ += 'from the file\n\n'
    str_ += file_name
    str_ += '\n\nThe line:\n\n'
    str_ += line.rstrip()
    str_ += '\n\ndoes not provide properly defined 12-6 LJ sigma value '
    str_ += 'for the atom ' + tmp[0] + '. The 12-6 LJ sigma value ' 
    str_ += 'should appear as floating point number, but the supplied one\n\n'
    str_ += str(tmp[2]) + '\n\ndoes not match the floating point number ' 
    str_ += 'format.\n'

    print(str_)

    if exit:

        argument.sys.exit(1)

    return None


def error_0009(argument, file_name, line, tmp, exit):

    '''
    Prints an error message (code 0009).
    '''

    str_ = '\nFATAL ERROR: A problem with one of the atom name definitions'
    str_ += 'has been detected when reading the ligand atom definitions '
    str_ += 'from the file\n\n'
    str_ += file_name
    str_ += '\n\nThe line:\n\n'
    str_ += line.rstrip()
    str_ += '\n\ndoes not provide properly defined 12-6 LJ epsilon value '
    str_ += 'for the atom ' + tmp[0] + '. The 12-6 LJ sigma value ' 
    str_ += 'should appear as floating point number, but the supplied one\n\n'
    str_ += str(tmp[3]) + '\n\ndoes not match the floating point number ' 
    str_ += 'format.\n'

    print(str_)

    if exit:

        argument.sys.exit(1)

    return None


def error_0010(argument, file_name, exit):

    '''
    Prints an error message (code 0010).
    '''

    str_ = '\nERROR: It seesm that the file:\n\n'
    str_ += _a_.ligand_defs_folder+i+_a_.ligand_file_prefix
    str_ += '\n\ndoes not contain '
    str_ += 'any ligand declarations. Check the file content.\n'

    print(str_)

    if exit:

        argument.sys.exit(1)

    return None

def get_ligand_defs_file(_a_):

    '''
    This function reads the content of the folder with ligands defintions\n
    and selects the files containing the records. The result is a list\n
    of files.
    '''

    import os

    #
    # Check if the folder containing the ligands exists.
    #
    if not os.path.isdir(_a_.ligand_defs_folder):

        error_0001(_a_, True)

        _a_.sys.exit(1)

    x = os.listdir(_a_.ligand_defs_folder)

    f1 = lambda x: [f for f in x if f.endswith(_a_.ligand_file_prefix)]

    y = lambda x : x.split(_a_.ligand_file_prefix)[:-1]

    _a_.ligand_names = [''.join(y(f)) for f in f1(x)]

    #
    # Extract the ligand names from the file names and check if the
    # collected names are really ligand names.
    #
    for i in _a_.ligand_names:

        if len(i) > 3:

            error_0002(_a_, i, True)

            _a_.sys.exit(1)

        else:

            x= [t.upper() in _a_.allowed_ligand_name_letters for t in i]

            if sum(x) != len(i):

                error_0003(_a_, i, exit)

                argument.sys.exit(1)

    return None


def read_ligand_definitions(_a_):

    '''
    This function reads the ligand atoms defined inside a text
    file. The list of text files is taken by calling the
    function "get_ligand_defs_file".
    '''

    for i in _a_.ligand_names:

        _a_.ligand[i.upper()] = {}

        counter = 0

        file_name = _a_.ligand_defs_folder+i+_a_.ligand_file_prefix

        with open(file_name, 'r') as f_obj:

            for line in f_obj:

                tmp = line.rstrip().split()

                #
                # If the line starts with "#" that is a comment line
                # and its content must not be examined. Any other
                # line will pass the examination listed bellow.
                #
                if tmp[0][0] != '#':

                    if len(tmp) < 4:

                        error_0004(_a_, i, line, exit)

                        _a_.sys.exit(1)

                    #
                    # Now examine the line columns in details. Rise an
                    # error and exit the execution of some column
                    # content does not match the expected format.
                    #

                    #
                    # Examine the atom name (first column in every line):
                    #
                    if len(tmp[0]) > 4:

                        #
                        # Examine if the atom name is more than 4 symbols
                        # long. The atom name have to follow the PDB
                        # specifications. Which means it should be
                        # non-empty string containing no more than 4
                        # symbols. 
                        #
                        error_0005(_a_, file_name, line, tmp, True)

                    if sum([j in _a_.allowed_ligand_name_letters\
                             for j in tmp[0]]) != len(tmp[0]):

                        #
                        # The next check scans the symbols included in the
                        # ligand file names. If they are not included in the
                        # list _a_.allowed_ligand_name_letters the examined
                        # string cannot be accepted as an atom name.
                        #
                        error_0006(_a_, file_name, line, tmp, \
                                   _a_.allowed_ligand_name_letters, True)

                    aname = '{:<4}'.format(tmp[0]).upper()

                    counter += 1

                    _a_.ligand[i.upper()][aname] = [counter]

                    #
                    # Examine the atomic charge:
                    #
                    try:

                        charge = float(tmp[1])

                    except ValueError:

                        error_0007(argument, file_name, line, tmp, True)

                    _a_.ligand[i.upper()][aname] += [charge]

                    #
                    # Examine the 12-6 LJ sigma:
                    #
                    try:

                        lj_sigma = float(tmp[2])

                    except ValueError:

                        error_0008(argument, file_name, line, tmp, True)

                    _a_.ligand[i.upper()][aname] += [lj_sigma]

                    #
                    # Examine the 12-6 LJ epsilon:
                    #
                    try:

                        lj_epsilon = float(tmp[3])

                    except ValueError:

                        error_0009(argument, file_name, line, tmp, True)

                    _a_.ligand[i.upper()][aname] += [lj_epsilon]

        if counter == 0:

            error_0010(argument, file_name, True)

    return None


def save_the_ligand_defs_as_pickles(_a_):

    '''
    This function saves the dictionary with the ligand definitions
    as a pickle binary container.
    '''

    import pickle

    with open(_a_.pickle_folder+_a_.pickle_file_name, 'wb') as f_obj:
        #
        pickle.dump(_a_.ligand, f_obj, pickle.HIGHEST_PROTOCOL)

    return None


if __name__ == '__main__':

    import pickle

    import sys

    class Dummy():

        pass

    argument = Dummy()

    argument.sys = sys

    argument.ligand_defs_folder = '../ligands/'

    argument.pickle_folder = '../pickle/'

    argument.pickle_file_name = argument.pickle_folder + '/ligands.pkl'

    argument.ligand_file_prefix = '.lig'

    argument.allowed_ligand_name_letters = 'ABCDEFGHIJKLMNOPQR'\
                                           'STUVWXYZ0123456789'

    if argument.ligand_defs_folder[-1] != '/':

        argument.ligand_defs_folder += '/'

    if argument.pickle_folder[-1] != '/':

        argument.pickle_folder += '/'

    argument.ligand_names = []

    argument.ligand = {}

    import time

    S = time.time()

    get_ligand_defs_file(argument)

    read_ligand_definitions(argument)

    save_the_ligand_defs_as_pickles(argument)

    total_time = '{:8.3e}'.format(time.time()-S)

    print('\n[OK] Time for performing the executon:',total_time,'s\n')
    
    sys.exit(0)
    
