"""
.. module:: manipulate_matrices
    :synopsis: Merge or reduce matrices

.. moduleauthor:: Benjamin Audren
"""
import numpy as np
import sys


def merge_matrices(path1, path2):

    # Recover the names of the first matrix, and load it
    with open(path1, 'r') as file1:
        names1 = [elem.strip() for elem in
                  file1.readline().strip()[1:].split(',')]
    matrix1 = np.loadtxt(path1)

    # Recover the names of the second matrix, and load it
    with open(path2, 'r') as file2:
        names2 = [elem.strip() for elem in
                  file2.readline().strip()[1:].split(',')]
    matrix2 = np.loadtxt(path2)

    print
    print 'Here are the two matrices to merge: '
    print '-->', path1
    print '-->', path2
    print
    print 'List of parameters from ', path1
    for index, elem in enumerate(names1):
        print index+1, elem
    print
    print 'Give the list of parameters you want to keep from these'
    print 'You can specify either all single parameters, separated by commas,'
    print 'or specify ranges, like omega_b:tau_reio, or 1:6 (6 included)'

    answer = None
    while not answer:
        answer = raw_input('Specify parameters to keep: ')

    indices1 = extract_indices(answer, names1)

    final_list = []
    print
    print r' /!\ Keeping the following parameters from the first matrix:'
    for index, name in enumerate(names1):
        if index in indices1:
            print ' --> %s' % name
            final_list.append(name)

    print '--------------------------------------------------------'
    print 'Choose now the parameters from the second matrix'
    print
    print 'List of parameters from ', path2
    for index, elem in enumerate(names2):
        print index+1, elem
    print

    answer = raw_input('Specify parameters to keep (can be empty): ')

    try:
        indices2 = extract_indices(answer, names2)
        print
        print ' /!\ Keeping the following parameters from the second matrix:'
        for index, name in enumerate(names2):
            if index in indices2:
                print ' --> %s' % name
                final_list.append(name)
    except ValueError:
        indices2 = []
        print
        print ' /!\ No parameters selected for the second matrix'

    # Fill in the first submatrix
    N = len(indices1)
    # Define a mask
    mask1 = np.zeros(np.shape(matrix1))
    for i in range(np.shape(matrix1)[0]):
        for j in range(np.shape(matrix1)[0]):
            if i in indices1 and j in indices1:
                mask1[i, j] = 1
    submat1 = matrix1[mask1 == 1].reshape((N, N))

    # Fill in the second submatrix
    M = len(indices2)
    mask2 = np.zeros(np.shape(matrix2))
    for i in range(np.shape(matrix2)[0]):
        for j in range(np.shape(matrix2)[0]):
            if i in indices2 and j in indices2:
                mask2[i, j] = 1
    submat2 = matrix2[mask2 == 1].reshape((M, M))

    # Join the two submatrices
    final_mat = np.zeros((N+M, N+M), 'float64')
    final_mat[:N, :N] = submat1
    final_mat[N:, N:] = submat2

    print 'Final list of parameters:'
    for index, name in enumerate(final_list):
        print index+1, name
    print

    result = None
    while not result:
        result = raw_input('writing the resulting covariance matrix to file:\n')

    with open(result, 'w') as output:
        output.write('# '+', '.join(['%16s' % e for e in final_list]))
        output.write('\n')
        for i in range(len(final_list)):
            for j in range(len(final_list)):
                if final_mat[i][j] >= 0:
                    output.write('% .5e\t' % final_mat[i][j])
                else:
                    output.write('% .5e\t' % final_mat[i][j])
            output.write('\n')

    print '\ndone !'


def extract_indices(answer, names):

    # At this stage, elements could contain names or indices (+1)
    elements = [e.strip() for e in answer.split(',')]
    # the following list will contain all the indices to keep
    indices = []
    for elem in elements:
        # Try first if it contains semi colon
        if elem.find(':') != -1:
            lhs, rhs = elem.split(':')
            try:
                first_index = int(lhs)-1
                last_index = int(rhs)-1
            except ValueError:
                first_index = names.index(lhs)
                last_index = names.index(rhs)
            indices.extend([i for i in range(first_index, last_index+1)])
        else:
            try:
                index = int(elem)-1
            except ValueError:
                index = names.index(elem)
            indices.append(index)
    return indices


if __name__ == '__main__':
    if len(sys.argv) == 3:
        path1 = sys.argv[-2]
        path2 = sys.argv[-1]

        merge_matrices(path1, path2)
    elif len(sys.argv) == 2:
        raise NotImplemented(
            "Please input two matrices and select no parameters"
            "in the second one")
    else:
        print 'You should specify two matrices to manipulate'
