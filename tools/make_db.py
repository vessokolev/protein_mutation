#!/usr/bin/env python3

__author__ = "Veselin Kolev"
__license__ = "GPLv2"
__version__ = "2019020600"
__maintainer__ = __author__
__email__ = "vesso.kolev@gmail.com"
__status__ = "Production"


__about__ = '''
This program code creates the dictionaries and lists required to
mutate the side chain of protein molecules described in PDB files.\n
You might need to supply manually the protein residue names as
elements of the list "pro_residues" bellow, as well as to mention
which residues to be excluded from the mutation process, by enlist
them as elements of the list "exluded_pro_residues".\n
Each protein residue must have be described as a separate file as
it shown in the folder "residues". Refer to the README file there
for explanations regarding the used syntax.
'''


############ LOOP FINDING BEGIN ############

def get_atom_connectivity(bonds, atom_serial_number):

    '''
    This function parses the list of bonds and finds\n
    all atoms connected to the atom with a serial number\n
    equal to "atom_serial_number".
    '''

    return [[j for j in i if j != atom_serial_number][0] for i in bonds\
             if i[0] == atom_serial_number or i[1] == atom_serial_number]


def get_collision_pos(_list_):

    '''
    Finds if two or more elements of the list "_list_"
    are not unique. Returns the positions of the collision
    backwardly (storing the problematic positions starting
    from the end of the array).
    '''

    if len(set(_list_)) != len(_list_):

        return [[-i - 1][0] for i in range(len(_list_)) for j in range(i + 1, len(_list_))\
                 if _list_[i] == _list_[j]][0]

    else:

        return 0


def update_atom_connectivity_tree(_a_, bonds, old_tree):

    '''
    This is a helper function to get the atom connectivity tree.
    It updates the branches of the connectivity that have been
    found so far.
    '''

    updated_tree = []

    for i in old_tree:

        if i[-1] > 0:

            dummy = get_collision_pos(i)

            if dummy < 0:

                tmp = i[:]

                tmp.append(dummy)

                updated_tree.append(tmp)

            else:
                connectivity = get_atom_connectivity(bonds, i[-1])

                if len(connectivity) == 1:

                    tmp = i[:]

                    if connectivity[0] == i[-2]:

                        tmp.append(0)

                    else:

                        tmp.append(connectivity[0])

                    updated_tree.append(tmp)

                else:

                    for j in connectivity:

                        tmp = i[:]

                        if j == i[-2]:
    
                            tmp.append(0)

                        else:

                            tmp.append(j)

                        updated_tree.append(tmp)
        else:

            tmp = i[:]

            tmp.append(0)

            updated_tree.append(tmp)


    # Mark the branches of the topology tree which are dead ends.
    # The dead ends are atoms connected to one atom only. Any of
    # those atoms is considered an end of a branch.

    for i in [j for j in range(len(updated_tree))]:

        if updated_tree[i][-1] == 0:

            indx_ = updated_tree[i].index(0)

            if updated_tree[i][indx_-1] in _a_.dead_ends:

                updated_tree[i][indx_] = -1

    return updated_tree


def get_atom_connectivity_tree(_a_, atom_serial_number,\
                               bond_id=0, revert_bond=False):

    '''
    This function finds all branches of the connectivity
    tree which pass through the bond "bond_id" or the atom
    with serial number "atom_serial_number".
    '''

    tree = []

    if bond_id == 0:

        connectivity = get_atom_connectivity(_a_.bond_list_numeric, atom_serial_number)

        tree = [[atom_serial_number, i] for i in connectivity]
    else:

        if revert_bond:

            tree = [[_a_.bond_list_numeric[bond_id][1], _a_.bond_list_numeric[bond_id][0]]]

        else:

            tree = [_a_.bond_list_numeric[bond_id]]

    flag = True

    while flag:

        tree = update_atom_connectivity_tree(_a_, _a_.bond_list_numeric, tree)

        #
        # Reduce the tree dimension by removing those of the
        # branches containing zero as last element of the
        # branch list, if the last non-zero number in the
        # list is greater than zero.
        #
        tree = [tree[i] for i in [i for i in range(len(tree)) if tree[i][-1] != 0 or\
                                 tree[i][tree[i].index(0)-1] < 0]]


        control = [i[-1] for i in tree]

        flag = len([x for x in control if x <= 0]) != len(control)

    return tree


def discover_dead_ends(_a_):

    '''
    Dead ends are the atoms from the side chain which are involved
    in only one chemical bond. That means they are dead ends of the
    side chain atom tree. Do call the function "collect_the_side_chain_atoms"
    in advance to get the list with relative atom serial numbers of the
    side chain atoms.
    '''

    #
    # Search if any of the side chain atoms appear twice in the bond list.
    #
    _a_.dead_ends = []

    #
    # Flatten the bond list to achieve fast speed of search.
    #
    bond_flat_list = [i for j in _a_.bond_list_numeric for i in j]

    for i in _a_.sidechain_atoms:

        flag = False

        try:

            flag = bond_flat_list.count(i) == 1

        except:

            str_ = '\nFATAL ERROR: It seems there is a serious problem'

            str_ += 'with either the topology or bond definitions for '

            str_ += 'the residue "'+_a_.tmp_resname+'"\n'

            print(str_)

            _a_.sys.exit(1)

        if flag:

            _a_.dead_ends.append(i)

    return None


def get_atom_connectivity_tree_loops(_a_, atom_serial_number,\
                                     bond_id=0, revert_bond=False):

    '''
    This function discovers loops in the molecular topology,
    based on investigating the bonds between the atoms. Some of
    the residues contains aromatic rings, which are considered
    loops.
    '''

    # Use a private functon here to examine each and every branch of the discovered
    # connectivity tree of the side chain.

    def check_tree_branch(tree):

        # Every branch needs to be terminated with a negative number - either -3 or -1.
        #
        # Every branch description have to contain only one negative number. If it contains
        # no negative number that means the branch is not terminated. Such a branch must
        # not be passed to this scope of the code! If there is more than one negative
        # number in the branch list, that is a general error and possible bug.

        check = sum([i < 0 for i in tree])

        if check > 1:

            # That means two or more numbers in the branch are negative. Terminate the execution
            # by rising error message!

            print('\nFATAL ERROR: When processing the branch', tree, 'detected for the ' +\
                  'side chain of the residue "' + _a_.tmp_resname + '", two or more termination ' +\
                  'numbers (numbers less than zero) have been detected. That is a general error and ' +\
                  'have to be reported to the developer of this code.\n')

            _a_.sys.exit(1)

        elif check == 0:

            # This is the case when the branch has no termination. That is also a general error
            # and have to be reported to the developer.

            print('\nFATAL ERROR: When processing the branch', tree, 'detected for the ' +\
                  'side chain of the residue "' + _a_.tmp_resname + '", no termination ' +\
                  'number (number less than zero) has been detected. That is a general error and ' +\
                  'have to be reported to the developer of this code.\n')

            _a_.sys.exit(1)

        # Find if there are any zeros positioned left of the termination number. That is a general
        # error and have to be reported.
        #
        # Get the index of the terminating number. If that number is not -3 and -1 that is a general
        # error of bug. Terminate the program.

        term_index = [i < 0 for i in tree].index(True)

        # The first two records in the branch list are the atom serial numbers (with respect to the
        # numbers defined in the res file for the processed residue) belong to C-alpha and C-beta
        # atoms. Those atoms are not part of the side chain. Therfore, if the term_index is less
        # than 2 that is a general error. It means that the bond selected for building the tree of
        # bonded atoms is not the one of CA-CB!

        if term_index < 2:

            print('\nFATAL ERROR: The first two atoms in the branch', tree, 'of the residue "', +\
                  _a_.tmp_resname + '" do not correspond to C-alpha C-betta bond! That ' +\
                  'means the bond selected for the root of the connectivity tree is not ' +\
                  'properly choosen!\n')

            _a_.sys.exit(1)

        if tree[term_index] != -3 and tree[term_index] != -1:

            print('\nFATAL ERROR: The number terminating the branch', tree, 'of the residue ' +\
                  '"' + _a_.tmp_resname+'" is not in the list of the expected numbers [-3,-1]. ' +\
                  'That is a general error and should be reported to the developer of this code.\n')

            _a_.sys.exit(1)

        elif sum([i == 0 for i in tree[:term_index]]) > 0:

            print('\nFATAL ERROR: At least one zero has been detected in the branch', tree, 'left ' +\
                  'the terminating number', tree[term_index], 'when processing the records for ' +\
                  'the residue "' + _a_.tmp_resname + '". That is a general error and should ' +\
                  'be reported to the developer of this code.\n')

            _a_.sys.exit(1)

        # Check if the atom numbers in tree[:term_index]] correspond to side chain atoms.

        check = [i in _a_.sidechain_atoms for i in tree[2:term_index]]

        if sum(check) != len(tree[2:term_index]):

            print('\nFATAL ERROR: At least on atom in the branch', tree[:term_index], 'of ' +\
                  'the residue "' + _a_.tmp_resname + '" is not a side chain atom. That ' +\
                  'is a general error and should be reported to the developer of this code.\n')

            _a_.sys.exit(1)

        return term_index

    # Discover dead ends in the connectivity in advance:

    discover_dead_ends(_a_)

    # Get the connectivity tree of the side chain:

    tree = get_atom_connectivity_tree(_a_, atom_serial_number, bond_id, revert_bond)

    # TODO: This is bug fix. In the next release this bug should be fixed completely.
    # Now only some workaround is implemented!

    for i in range(len(tree)):

        if tree[i][-1] < -1:

            tree[i][-1] = -3

        elif tree[i][-1] == 0:

            # Find the last non-zero element using backward search:

            indx = [j for j in range(len(tree[i]) - 1, -1, -1) if tree[i][j] < -1]

            if len(indx) > 0:

                tree[i][indx[0]] = -3

    # Examine each branch of the connectivity tree (each member of the list "tree"), and
    # get the index of the terminating number (that number should be either -3 or -1).

    term_index = [check_tree_branch(i) for i in tree]

    # Mark the loops in the list members of "tree". If any of them contains a loop, that loop
    # is the shortest possible for that branch (primitve loop). No loops including other loops
    # are stored there.

    loop = [i for i in range(len(tree)) if tree[i][term_index[i]] == -3]
    open = [i for i in range(len(tree)) if tree[i][term_index[i]] == -1]

    # For each branch (each element in the list "tree") get both starting and ending indexes
    # to isolate the loop as the slice "tree[i][ll_loop_indx:ul_loop_indx]".

    ll_loop_indx = [tree[i][:term_index[i] - 1].index(tree[i][term_index[i] - 1]) for i in loop]
    ul_loop_indx = [term_index[i] - 1 for i in loop]

    # By using the loops create list of axes, which are not suitable for being used as
    # axes of rotation.

    axis_to_avoid = []

    # The lambda "f1" is required to check if the examined axis is not already
    # listed in "axis_to_avoid". Note that "f1" should be defined and used only
    # of loops are detected.

    if len(loop) > 0:

        f1 = lambda k: (not [tree[i][k], tree[i][k + 1]] in axis_to_avoid and\
                        not [tree[i][k + 1], tree[i][k]] in axis_to_avoid)

    axis_to_avoid += [[tree[loop[i]][k], tree[loop[i]][k + 1]] for i in range(len(loop))\
                     for k in range(ll_loop_indx[i], ul_loop_indx[i]) if f1(k)]

    _a_.axis = [_a_.bond_list_numeric[bond_id]]

    # "f2" helps to check if a certain axis of rotation is already in the list
    # "_a_.axis".

    f2 = lambda k: (not [tree[i][k], tree[i][k + 1]] in _a_.axis and\
                    not [tree[i][k + 1], tree[i][k]] in _a_.axis)

    # "f3" helps to discover if at least one atom of the examined axis is a
    # dead end in the connectivity topology.

    f3 = lambda k: not tree[i][k] in _a_.dead_ends and\
                   not tree[i][k + 1] in _a_.dead_ends

    # Define the lambda "f" depending on the conditions. If no loops are detected
    # in the connectivity tree, there is no need to call the lambda "f1".

    if len(loop) > 0:

        f = lambda k: f1(k) and f2(k) and f3(k)

    else:

        f = lambda k: f2(k) and f3(k)

    # Proceed with the branches which does not contain loops (that does not mean
    # they do not include part of some loop - they just do not contain the whole
    # loop). It is not easy to get rid of the for loop in expanded for given
    # bellow, because the lambda "f" checks the presence of the axis in the list
    # _a_.axis which is a subject of management inside the loop.

    for i in open:

        _a_.axis += [[tree[i][k], tree[i][k + 1]] for k in range(1, term_index[i] - 1)\
                         if f(k)]

    # It is possible now to collect the branches of atoms which should be rotate around
    # each of the discovered axes. The branches to parse (those in the list "tree") are
    # composed in a way providing optimization - if one select to rotate the part of the
    # branche appearing right of the axis, the minimum number of atoms will be rotated.
    # In other words, if a rotation in the side chain of one protein is required it will
    # rotate only part of the side chain with respect to the protein and not the rest of
    # the protein atoms with respect to that part of the chosen side chain.

    _a_.branches_to_rotate = []

    for i in _a_.axis:

        tmp = []

        for j in range(len(tree)):

            indx = [k + 2 for k in range(term_index[j]-1)\

                    if tree[j][k] == i[0] and tree[j][k + 1] == i[1]]

            if len(indx) == 1:

                if j in loop:

                    m=term_index[j] - 1

                else:

                    m=term_index[j]

                tmp += [k for k in tree[j][indx[0]:m] if not k in tmp and not k in i]

        tmp.sort()

        _a_.branches_to_rotate.append(tmp)

    return None


############ LOOP FINDING END ############


def save_residue_topology_dicts_as_pkl_files(_a_):

    '''
    Saves the collected lists and dictionaries as Pickle storage containers.
    By specifying bellow pickle.HIGHEST_PROTOCOL, a binary files are created.
    Note the "wb" argument specified in "open". It tells Python to write a
    binary content to the file. Loading lists and dictionaries from within
    a binary file is always faster then reading them from text files and
    parsing them afterwards.
    '''

    proto = _a_.pickle.HIGHEST_PROTOCOL

    with open(_a_.pickle_map_residues_file, 'wb') as f_obj:
        _a_.pickle.dump(_a_.map_residues, f_obj, proto)

    with open(_a_.pickle_map_sequence_file, 'wb') as f_obj:
        _a_.pickle.dump(_a_.map_sequence, f_obj, proto)

    with open(_a_.pickle_sequence_file, 'wb') as f_obj:
        _a_.pickle.dump(_a_.sequence, f_obj, proto)

    with open(_a_.pickle_pro_residues_file, 'wb') as f_obj:
        _a_.pickle.dump(_a_.pro_residues, f_obj, proto)

    with open(_a_.pickle_excluded_pro_residues_file, 'wb') as f_obj:
        _a_.pickle.dump(_a_.excluded_pro_residues, f_obj, proto)

    with open(_a_.pickle_axis_dict_file, 'wb') as f_obj:
        _a_.pickle.dump(_a_.axis_dict, f_obj, proto)

    with open(_a_.pickle_branches_to_rotate_dict_file, 'wb') as f_obj:
        _a_.pickle.dump(_a_.branches_to_rotate_dict, f_obj, proto)

    with open(_a_.pickle_nb_neigh_file, 'wb') as f_obj:
        _a_.pickle.dump(_a_.nb_neigh_dict, f_obj, proto)

    with open(_a_.pickle_dihedrals_file, 'wb') as f_obj:
        _a_.pickle.dump(_a_.dihedrals_dict, f_obj, proto)

    with open(_a_.pickle_bonds_file, 'wb') as f_obj:
        _a_.pickle.dump(_a_.bonds, f_obj, proto)

    with open(_a_.pickle_bond_angles_file, 'wb') as f_obj:
        _a_.pickle.dump(_a_.bond_angles, f_obj, proto)

    return None


def collect_residue_topology_dicts(_a_):

    '''
    This function assembles all dictionaries and lists that will
    be stored on the disk as pickle containers in binary format.
    '''

    seq_counter = 0

    for i in _a_.pro_residues:

        file_name = _a_.res_raw_definitions_dir+'/'+i.lower()+'.res'

        _a_.map_residues[i] = {}

        _a_.map_sequence[i] = {}

        seq = []

        with open(file_name, 'r') as f_obj:

            counter = 1

            for line in f_obj:

                tmp = line.rstrip().split()

                if len(tmp) >= 6:

                    _a_.map_residues[i]['{:<4}'.format(tmp[0])] = \
                        [\
                            counter, '{:<2}'.format(tmp[2]), \
                            float(tmp[3]), float(tmp[4]), float(tmp[5]), \
                            float(tmp[6]), float(tmp[7]), float(tmp[8]), \
                            int(tmp[9])\
                        ]

                    seq.append('{:<4}'.format(tmp[0]))

                    counter += 1

        _a_.resname_seq.append(i)

        _a_.sequence.append(seq)

        _a_.map_sequence[i] = seq_counter

        seq_counter += 1

    return None


def create_dihedral_sqlite3_db(_a_):
    '''
    Create a table in the in-memory SQLite3 database for storing there
    the force field parameters related to the computation of the energy
    contribution due to the proper dihedral rotation. Also create
    indexes for the columns of interest to bust the search.\n
    TODO: Do not supply to the table the atom types as strings. Use
    integer numbers instead to point to every atom type. Combined with
    indexes that speeds up the SELECT based SQLite3 backward searches
    enormously.
    '''

    #
    # Create table named "dihedrals":
    #
    query = 'CREATE TABLE dihedrals(type_1 VARCHAR(2), type_2 VARCHAR(2), '+\
                                 'type_3 VARCHAR(2), type_4 VARCHAR(2), '+\
                                 'phi_0 REAL, kd REAL, n INT)'

    _a_.cursor.execute(query)

    #
    # Read the force field file line by line, strip the parameters and
    # insert them into the table "dihedrals":
    #
    with open(_a_.ff_file, 'r') as f_obj:

        for line in f_obj:

            tmp = line.rstrip().split()

            if len(tmp) >= 8 and not tmp[0] in ['#', ';', '*', '/']:

                #
                # Note that the check not tmp[0] in ['#',';','*','/'] assures
                # that the line is not a comment one.
                #
                query = 'INSERT INTO dihedrals VALUES(?,?,?,?,?,?,?)'

                _a_.cursor.execute(query, ('{:<2}'.format(tmp[0]),\
                                            '{:<2}'.format(tmp[1]),\
                                            '{:<2}'.format(tmp[2]),\
                                            '{:<2}'.format(tmp[3]),\
                                            float(tmp[5])*_a_.deg_2_rad,\
                                            float(tmp[6]),\
                                            int(tmp[7]),))
    #
    # Create the indexes to support fast SELECT searches:
    #
    query = 'CREATE INDEX type_1_indx ON dihedrals(type_1)'
    _a_.cursor.execute(query)

    query = 'CREATE INDEX type_2_indx ON dihedrals(type_2)'
    _a_.cursor.execute(query)

    query = 'CREATE INDEX type_3_indx ON dihedrals(type_3)'
    _a_.cursor.execute(query)

    query = 'CREATE INDEX type_4_indx ON dihedrals(type_4)'
    _a_.cursor.execute(query)

    return None


def check_dihedral_against_sqlite3_db(_a_):

    '''
    This function checks if a given dihedral, supplied as a Python
    list assigned to the variable _a_.dihedral_to_check, and containing
    the atom types of the i-th, j-th, k-th, and l-th atoms, exits in the
    database. If exists, the function returns the type 9 dihedral parameters
    as t_pointer. First, the function searches for a perfect match, which
    means that all for members of the _a_.dihedral_to_check are listed
    as a record(s). If no march is found as a result of that search, the
    search is repeating only for the j-th and k-th atom types pair.
    '''

    #
    # First check for a full (perfect) match:
    #
    query = 'SELECT * FROM dihedrals WHERE type_1=? AND type_2=? AND type_3=? AND type_4=?'

    _a_.cursor.execute(query, (_a_.dihedral_to_check[0],\
                               _a_.dihedral_to_check[1],\
                               _a_.dihedral_to_check[2],\
                               _a_.dihedral_to_check[3],))

    t_pointer = _a_.cursor.fetchall()

    if len(t_pointer) == 0:

        _a_.cursor.execute(query, (_a_.dihedral_to_check[3],\
                                   _a_.dihedral_to_check[2],\
                                   _a_.dihedral_to_check[1],\
                                   _a_.dihedral_to_check[0],))

        t_pointer = _a_.cursor.fetchall()

    if len(t_pointer) == 0:

        _a_.cursor.execute(query, ('X ',\
                                   _a_.dihedral_to_check[1],\
                                   _a_.dihedral_to_check[2],\
                                   'X ',))

        t_pointer = _a_.cursor.fetchall()

    if len(t_pointer) == 0:

        _a_.cursor.execute(query, ('X ',\
                                   _a_.dihedral_to_check[2],\
                                   _a_.dihedral_to_check[1],
                                   'X ',))

        t_pointer = _a_.cursor.fetchall()

    if len(t_pointer) == 0:

        return None

    else:

        return t_pointer


def check_if_record_is_bond(_a_):

    '''
    This function checks the format of the BND (bond) records.\n
    They must consist of two atom name, that belong to the\n
    corresponding residue of interest.
    '''

    # Check if the input is an bond record:

    if len(_a_.tmp_split_line) >= 3:

        # That might be a bond record. Check the fields:

        bond = ['{:<4}'.format(_a_.tmp_split_line[i].upper()) for i in range(1, 3)]

    else:

        # Rise an error flag if the bond record is not properly
        # formatted and exit the execution process.

        str_ = '\nFATAL ERROR: The bond definition "' + _a_.tmp_line.rstrip()
        str_ += '" made in the topology file "' + _a_.current_topology_file
        str_ += '", does not match the expected format\n'

        print(str_)

        _a_.sys.exit(1)

    return bond


def get_bond_angles(bond):

    '''
    This function discovers the bond angles based on a bond list,
    which is supplied as a Python list as the only input of the
    function. The function returns two lists: (i) the  list of the
    bond angles, and (ii) the list of the indexes of the bonds
    involved into each bond angle (the indexes are those of the
    bonds supplied by the input list).
    '''

    bond_angle = []

    bond_angle_composition = []

    for i in range(len(bond)-1):

        for j in range(i+1, len(bond)):

            factors = [[i, 0, j, 0, i, 1, i, 0, j, 1], \
                       [i, 0, j, 1, i, 1, i, 0, j, 0], \
                       [i, 1, j, 0, i, 0, i, 1, j, 1], \
                       [i, 1, j, 1, i, 0, i, 1, j, 0]]

            for factor in factors:

                if bond[factor[0]][factor[1]] == bond[factor[2]][factor[3]]:

                    addition = [bond[factor[4]][factor[5]], bond[factor[6]][factor[7]], \
                                bond[factor[8]][factor[9]]]

                    if not (addition in bond_angle or addition.reverse() in bond_angle):

                        # Note that the reverse() function reverses the bond candidate
                        # during the check above. But we do not really care if the bond
                        # angle is supplied as A-B-C or C-B-A. Otherwise the check must
                        # be done like this:
                        #
                        #  if not (addition in bond_angle\
                        #     or [k for k in reversed(addition)] in bond_angle)
                        #
                        # since reversed function returns an iterator object and not list.

                        bond_angle.append(addition)

                        bond_angle_composition.append([i, j])

    return bond_angle, bond_angle_composition


def get_proper_dihedral_angles(bond, bond_angle, bond_angle_composition):

    '''
    The function discovers the proper dihedrals based on the
    lists with defined bonds, bond angles, and the indexes of
    the bonds involved in the bond angles (see the function
    "get_bond_angles" to check how these indexes are taken).
    '''

    dihedral = []

    for i in range(len(bond)):

        for j in range(len(bond_angle)):

            if not i in bond_angle_composition[j]:

                if (bond[i][0] in bond_angle[j]) or (bond[i][1] in bond_angle[j]):

                    # The factors list makes is a list of permutations that speed up
                    # the process of checking the bond presence in given bond angle
                    # when constructing the proper dihedrals.

                    factors = [[i, 0, j, 0, i, 1, i, 0, j, 1, j, 2],\
                               [i, 0, j, 2, i, 1, i, 0, j, 1, j, 0],\
                               [i, 1, j, 0, i, 0, i, 1, j, 1, j, 2],\
                               [i, 1, j, 2, i, 0, i, 1, j, 1, j, 0]]

                    for factor in factors:

                        if bond[factor[0]][factor[1]] == bond_angle[factor[2]][factor[3]]:

                            addition = [bond[factor[4]][factor[5]],\
                                        bond[factor[6]][factor[7]],\
                                        bond_angle[factor[8]][factor[9]],\
                                        bond_angle[factor[10]][factor[11]]]

                            if not (addition in dihedral or addition.reverse() in dihedral):

                                # Note that the reverse() function reverses the dihedral candidate
                                # during the check above. But we do not really care if the dihedral
                                # is supplied as A-B-C-D or D-C-B-A. Otherwise the check must
                                # be done like this:
                                #
                                #  if not (addition in dihedral\
                                #     or [k for k in reversed(addition)] in dihedral)
                                #
                                # since reversed function returns an iterator object and not list.

                                dihedral.append(addition)

    return dihedral


def read_topology_file(_a_):
    '''
    This function reads the residue structure definitions
    from the residue topology file and collects the bond records.
    '''

    # Check the file size. If it is bigger than _a_.topo_file_max_size
    # then terminate the program!

    if _a_.os.stat(_a_.current_topology_file).st_size > _a_.topo_file_max_size:

        print('\nFATAL ERROR: The size of the topology file "' + _a_.current_topology_file +\
              '" is bigger than the maximum allowed size (' + _a_.topo_file_max_size + ' bytes)!\n')

        _a_.sys.exit(1)

    counter = 0

    # Assign empty lists to the storages for collecting the bonds, branches
    # and axes:

    _a_.bond_list = []

    # Do not load the topology file in the memory before reading it.
    # Read it directly line by line instead.

    with open(_a_.current_topology_file, 'r') as f_obj:

        for _a_.tmp_line in f_obj:

            counter += 1

            if counter <= _a_.max_lines_to_read:

                _a_.tmp_split_line = _a_.tmp_line.rstrip().split()

                # Make the parsing efficient. First, check _a_.tmp_split_line[0]
                # and depending on its value, send _a_.tmp_split_line to the
                # corresponding subroutine for further processing:

                if len(_a_.tmp_split_line) >= 2:

                    if _a_.tmp_split_line[0].upper() == 'BND':

                        # This is a bond. Process it accordingly. If the bond is already
                        # in the list, do not append it.

                        bond = check_if_record_is_bond(_a_)

                        # Check if the bond consists of two identical atoms names. That is a mistake:

                        if bond[0] == bond[1]:

                            print('\nFATAL ERROR: The bond record "' + _a_.tmp_line.rstrip() +\
                                  '" in the topology file "' + _a_.current_topology_file +\
                                  '" consists of two identical atom names!\n')

                            _a_.sys.exit(1)

                        else:

                            # Check if the atom names are part of the topology description of the
                            # residue:

                            flag = [i in _a_.map_residues[_a_.tmp_resname].keys() for i in bond]

                            for i in range(2):

                                if not flag[i]:

                                    print('\nFATAL ERROR: The bond record "' + _a_.tmp_line.rstrip() +\
                                          '" in the topology file "' + _a_.current_topology_file +\
                                          '" contains the atom name "' + bond[i] + '", which is not part of ' +\
                                          'the topology description of the residue "' + _a_.tmp_resname + '"!\n')

                                    _a_.sys.exit(1)

                            if not bond in _a_.bond_list and not bond.reverse() in _a_.bond_list:

                                _a_.bond_list.append(check_if_record_is_bond(_a_))

            else:

                print('\nFATAL ERROR: The topology file "' + _a_.current_topology_file +\
                      '" contains more lines than the maximum allowed number of lines (' +\
                      _a_.max_lines_to_read + ' lines)!\n')

                _a_.sys.exit(1)

    return None


def get_non_bonded_neighs_sidechain(_a_):

    # TODO!!! Too much nested loops. Use lambdas.

    '''
    This function selects the pairs of atoms which are
    recognized non-bonded neighbors. Each pair have to
    contains at least one atom that belongs to the side
    chain or the examined residue. The basic principle
    to follow when selecting the pair of non-bonded
    neighbors is to find any pair of atoms which do not
    belong to any bond, bond angle, or proper dihedral
    angle.
    '''
    #
    # If none side chain atoms are declared just skip the checks above!
    #
    if len(_a_.sidechain_atoms) > 0: # add bond and dihedral lengths. if no dihedrals then skip
        pass
    #
    # The relative atoms serial numbers (valid only when numbering the atoms in one
    # residue) are defined as a sequence. That in turn means that they start from 1
    # and go up by increment of 1. If they do not obey that schema of numbering,
    # rise an error flag and message!
    #
    # First, generate list with the sequence from 1 to the total number of atoms in
    # the residue.
    #
    len_ = len(_a_.sequence[_a_.map_sequence[_a_.tmp_resname]])
    #
    seq = [i + 1 for i in range(len_)]
    #
    if sum([i in seq for i in _a_.sidechain_atoms]) != len(_a_.sidechain_atoms):
        #
        # All atoms from sidechain_atoms list have to be members of seq list.
        # Otherwise the enumeration of the atom serial numbers is not done
        # properly in the dictionaries. That must be very unlike to happen
        # error. Nevertheless, keep the check.
        #
        print('\nFATAL ERROR (in get_non_bonded_neighs_sidechain): Some members '+\
              'of the list "sidchain_atoms" are not found in the atom serial numbers '+\
              'list "seq"!\n')
        print('\nThe content of the list "sidechain_atoms" is', _a_.sidechain_atoms)
        print('The content of the list "seq" is', seq, '\n')
        #
        _a_.sys.exit(1)
        #
    #
    # Start checking for non-bonded neighbors. To check them against the bond angle list
    # the bond list need to be re-arranged to allow the search to take place only for he
    # ending atoms of the angle, since the bonds are also checked and the bond angles are
    # composed by combining two bonds.
    #
    bond_angle_ends = [[i[0], i[2]] for i in _a_.bond_angle]
    #
    _a_.nb_neigh = []
    #
    for i in _a_.sidechain_atoms:
        #
        for j in range(len_):
            #
            if i != j + 1 and not [j + 1, i] in _a_.nb_neigh:
                #
                # Define [i,j+1] and check if matches any of the bonds
                #
                if not [i, j + 1] in _a_.bond_list and not [j + 1, i] in _a_.bond_list:
                    #
                    # i and j+1 are potential non-bonded neighbors. Let's check if they
                    # don't belong to any bond angle.
                    #
                    if not [i, j + 1] in bond_angle_ends and not [j + 1, i] in bond_angle_ends:
                        #
                        # If i and j+1 are not involved into the same bond or bond angle,
                        # check if they are not part of some of the proper dihedral angles.
                        # Use an advancing method. Iterate until file a match.
                        #
                        flag = True
                        #
                        for k in _a_.dihedrals:
                            #
                            if i in k and j + 1 in k:
                                #
                                flag = False
                                #
                                break
                                #
                        if flag:
                            #
                            _a_.nb_neigh.append([i, j + 1])
                            #
    return None


def collect_the_side_chain_atoms(_a_):

    '''
    This function collects the relative atom serial numbers
    of the side chain atoms based on the records provided by
    the residue topology definitions.
    '''

    _a_.sidechain_atoms += [_a_.map_residues[_a_.tmp_resname][i][0]\

         for i in _a_.sequence[_a_.map_sequence[_a_.tmp_resname]]\

            if _a_.map_residues[_a_.tmp_resname][i][8] == 1]

    return None


def get_the_branches_root(_a_):

    '''
    This function selects the bond between C-alpha and C-beta
    atoms, which is the root element in the connectivity topology.
    All branches of atoms in the side chain that might be rotated
    around bond (that bond is also part of the side chain) are
    discovered by starting from the C-alpha - C-beta bond.
    '''

    #
    # We are looking for the bond whose members are with roles -1 (C-alpha) and
    # 3 (C-beta).
    #
    t_pointer = _a_.map_residues[_a_.tmp_resname]

    counter = 0

    flag = True

    found = False

    for i in _a_.bond_list:

        if t_pointer[i[0]][8] == -1 and t_pointer[i[1]][8] == 3 or\
           t_pointer[i[1]][8] == -1 and t_pointer[i[0]][8] == 3:

            found = True

            if t_pointer[i[0]][8] == -1:

                flag = False

            break

        counter += 1

    if not found:

        print('\nFATAL ERROR: Cannot find the C-alpha - C-beta bond from the bond '+\
              'definitions made for the residue "' + _a_.tmp_resname + '"\n')

        _a_.sys.exit(1)

    return counter, flag


def analyze_bond_records(_a_):

    '''
    This function examines bond records, helps deriving the
    definitions of the proper dihedral angles, and read the
    force field parameters assigned to any of those angles.
    '''

    #
    # To process the bond list effectively the atom names need to be converted into
    # atom numbers, according to the records made in _a_.map_residues.
    #
    t_pointer = _a_.map_residues[_a_.tmp_resname]

    _a_.bond_list_numeric = [[t_pointer[i[0]][0], t_pointer[i[1]][0]]\
                                  for i in _a_.bond_list]

    #
    # Construct the proper dihedrals based on bond information. Bond angles are required
    # as intermediate data.
    #

    bond_angles, bond_angle_composition = get_bond_angles(_a_.bond_list_numeric)

    _a_.dihedrals = get_proper_dihedral_angles(_a_.bond_list_numeric,\
                                                   bond_angles,\
                                                   bond_angle_composition)

    if not _a_.tmp_resname in _a_.bond_angles.keys():

        _a_.bond_angles[_a_.tmp_resname] = bond_angles

    if not _a_.tmp_resname in _a_.bonds.keys():

        _a_.bonds[_a_.tmp_resname] = _a_.bond_list_numeric

        print(_a_.bonds)

        print()

    #
    # To get the force field parameters related to the proper dihedral angles,
    # the members of each dihedral need to be turned into atom type names.
    # Therefore, the members of the list _a_.dihedrals should be first
    # resolved to atom names and then to atom types.
    #
    _a_.dihedrals_atypes = []

    _a_.dihedrals_ff_params = []

    t_pointer = _a_.sequence[_a_.map_sequence[_a_.tmp_resname]]

    for i in _a_.dihedrals:

        _a_.dihedrals_atypes += [[None, None, None, None]]

        _a_.dihedrals_atypes[-1][:] =\
            [_a_.map_residues[_a_.tmp_resname][t_pointer[i[j]-1]][1]\
                for j in [0, 1, 2, 3]]

        dummy = []

        _a_.dihedral_to_check = _a_.dihedrals_atypes[-1]

        _a_.dihedrals_ff_params += [[list(j[4:])\
            for j in check_dihedral_against_sqlite3_db(_a_)]]

    return None


def create_binary_topology(_a_):

    '''
    This function organizes the discovery of the topology.
    '''

    _a_.bond_list = []

    _a_.bond_list_axes = []

    _a_.bond_angle = []

    _a_.sidechain_atoms = []

    _a_.nb_neigh = []

    _a_.dead_ends = []

    _a_.dihedral_defs = {}

    _a_.axis = []

    _a_.dihedrals = []

    _a_.branches_to_rotate = []

    _a_.sidechain_atoms = []
    
    _a_.dihedrals_to_save = []

    read_topology_file(_a_)

    analyze_bond_records(_a_)

    read_topology_file(_a_)

    analyze_bond_records(_a_)

    collect_the_side_chain_atoms(_a_)

    get_non_bonded_neighs_sidechain(_a_)

    resolved_topology[_a_.tmp_resname] = {}

    bond_id, reverse_bond = get_the_branches_root(_a_)
    
    get_atom_connectivity_tree_loops(_a_, 0, bond_id)
    
    return None


if __name__ == "__main__":


    import sys
    #
    # Check the Python version. Be sure that at least 3.4 is available.
    #
    version = sys.version_info

    if version.major < 3 or (version.major == 3 and version.minor < 4):

        print('\nRunning this script requires Python version 3.4 or higher.\n')

        sys.exit(1)

    #
    # Describe here all the protein residue names as a list of strings with
    # length 3. The residue names not listed in "pro_residues" will not be
    # recognized as a part of the protein molecule.
    #
    pro_residues = ['ALA', 'GLY', 'SER', 'THR', 'LEU', 'ILE', 'VAL', 'ASN',\
                    'GLN', 'ARG', 'HIS', 'HIE', 'HIP', 'TRP', 'PHE', 'TYR',\
                    'GLU', 'ASP', 'LYS', 'PRO', 'CYS', 'CYI', 'MET']

    #
    # Point to the residue name(s) which need to be excluded from the process
    # of mutation. For instance, the 'PRO' residue is not a subject of mutation.
    # One may exclude any of the residue names listed in "pro_residues" by just
    # specifying the index of the residue string in the list "pro_residues".
    #
    excluded_pro_residues = [pro_residues[1], pro_residues[19]]

    #
    # The 3 letter names of the protein residues must consist of capital letters
    # only. Therefore, convert the strings supplied by the user as elements of
    # the array "pro_residues" into ones strings containing capital letters
    # only. Do the same to the elements of "excluded_pro_residues".
    #
    pro_residues = [i.upper() for i in pro_residues]

    excluded_pro_residues = [i.upper() for i in excluded_pro_residues]

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

    argument.topo_file_max_size = 100000 # The size of the topology file in B

    argument.max_lines_to_read = 1000

    argument.res_raw_definitions_dir = '../residues/'

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

    argument.pickle_nb_neigh_file = argument.pickle_storage_dir+\
                                    '/nonbonded.pkl'

    argument.pickle_dihedrals_file = argument.pickle_storage_dir+\
                                     '/dihedrals.pkl'

    argument.pickle_bonds_file = argument.pickle_storage_dir+\
                                     '/bonds.pkl'

    argument.pickle_bond_angles_file = argument.pickle_storage_dir+\
                                     '/bond_angles.pkl'

    argument.pro_residues = pro_residues

    argument.excluded_pro_residues = excluded_pro_residues

    argument.map_residues = {}

    argument.map_sequence = {}

    argument.sequence = []

    argument.resname_seq = []

    argument.bonds = {}

    argument.bond_angles = {}

    argument.deg_2_rad = 0.017453292519943295

    argument.ff_file = '../forcefield/amber99.ff'

    collect_residue_topology_dicts(argument)

    #
    # Create in-memory SQLite3 database:
    #
    connection = argument.sqlite3.connect(':memory:')

    #
    # and assign the pointer to the cursor object:
    #
    argument.cursor = connection.cursor()

    #
    # Create the SQLite 3 table database containing the force field
    # parameters for the type 9 proper dihedral angles:
    #
    create_dihedral_sqlite3_db(argument)

    #
    # Define a variable (dictionary) to add there all
    #
    resolved_topology = {}

    argument.axis_dict = {}

    argument.branches_to_rotate_dict = {}

    argument.nb_neigh_dict = {}

    argument.dihedrals_dict = {}

    #
    # Collect the garbage here:
    #
    import gc

    gc.collect()

    for argument.tmp_resname in pro_residues:

        if not argument.tmp_resname in excluded_pro_residues:

            argument.current_topology_file = '../residues/' + argument.tmp_resname.lower() + '.top'

            create_binary_topology(argument)

            argument.axis_dict[argument.tmp_resname] = argument.axis

            argument.branches_to_rotate_dict[argument.tmp_resname] = argument.branches_to_rotate

            argument.nb_neigh_dict[argument.tmp_resname] = argument.nb_neigh

            f1 = lambda m: argument.dihedrals_ff_params[m][0][1] != 0.0

            f2 = lambda m: sum([l[1] for l in argument.dihedrals_ff_params[m]]) != 0.0

            argument.dihedrals_dict[argument.tmp_resname] = []

            for i in range(len(argument.dihedrals)):

                tmp_1=[]

                if [argument.dihedrals[i][1], argument.dihedrals[i][2]] in argument.axis\
                    or [argument.dihedrals[i][2], argument.dihedrals[i][1]] in argument.axis:

                    for j in range(len(argument.dihedrals_ff_params[i])):

                        if argument.dihedrals_ff_params[i][j][1] != 0.0:

                            tmp_1 += [argument.dihedrals_ff_params[i][j]]

                    if len(tmp_1) >= 1:

                        argument.dihedrals_dict[argument.tmp_resname] += [[argument.dihedrals[i], tmp_1]]

    save_residue_topology_dicts_as_pkl_files(argument)

    total_time = '{:8.5e}'.format(time.time()-S)

    print('\n[OK] Time for performing the executon:',total_time,'s\n')

    sys.exit(0)
