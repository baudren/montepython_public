"""
.. module:: mcmc
   :synopsis: Monte Carlo procedure
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>

This module defines one key function, :func:`chain`, that handles the Markov
chain. So far, the code uses only one chain, as no parallelization is done.

This function in turn calls several other routines. These are called just once:

* :func:`get_covariance_matrix`
* :func:`read_args_from_chain`
* :func:`read_args_from_bestfit`

Their usage is pretty straightforward, and detailled below anyway. On the
contrary, these routines are called at every step:

* :func:`compute_lkl` is called at every step in the Markov chain, returning
  the likelihood at the current point in the parameter space.
* :func:`get_new_position` returns a new point in the parameter space,
  depending on the proposal density.

The arguments of these functions will often contain **data** and/or **cosmo**.
They are both initialized instances of respectively :class:`data` and the
cosmological class. They will thus not be described for every function.
"""

import os
import sys
import math
import random as rd
import numpy as np
import io_mp
import data
import scipy.linalg as la
import unittest


def compute_lkl(cosmo, data):
    """
    Compute the likelihood, given the current point in parameter space.

    This function now performs a test before calling the cosmological model
    (**new in version 1.2**). If any cosmological parameter changed, the flag
    :code:`data.need_cosmo_update` will be set to :code:`True`, from the
    routine :func:`check_for_slow_step <data.data.check_for_slow_step>`.

    :Returns:
        - **loglike** (`float`) - the log of the likelihood
          (:math:`\\frac{-\chi^2}2`) computed from the sum of the likelihoods
          of the experiments specified in the input parameter file.

          This function returns :attr:`data.boundary_loglkie
          <data.data.boundary_loglike>`, defined in the module :mod:`data` if
          *i)* the current point in the parameter space has hit a prior edge,
          or *ii)* the cosmological module failed to compute the model. This
          value is chosen to be extremly small (large negative value), so that
          the step will always be rejected.


    """

    # If the cosmological module has already been called once, and if the
    # cosmological parameters have changed, then clean up, and compute.
    if (cosmo.state is True and data.need_cosmo_update is True):
        cosmo._struct_cleanup(set(["lensing", "nonlinear", "spectra",
                                    "primordial", "transfer", "perturb",
                                    "thermodynamics", "background", "bessel"]))

    # If the data needs to change, then do a normal call to the cosmological
    # compute function. Note that, even if need_cosmo update is True, this
    # function must be called if the jumping factor is set to zero. Indeed,
    # this means the code is called for only one point, to set the fiducial
    # model.
    if ((data.need_cosmo_update) or
            (not cosmo.state) or
            (data.jumping_factor == 0)):

        # Prepare the cosmological module with the new set of parameters
        cosmo.set(data.cosmo_arguments)

        # Compute the model, keeping track of the errors

        # In classy.pyx, we made use of two type of python errors, to handle
        # two different situations.
        # - AttributeError is returned if a parameter was not properly set
        # during the initialisation (for instance, you entered Ommega_cdm
        # instead of Omega_cdm).  Then, the code exits, to prevent running with
        # imaginary parameters. This behaviour is also used in case you want to
        # kill the process.
        # - NameError is returned if Class fails to compute the output given
        # the parameter values. This will be considered as a valid point, but
        # with minimum likelihood, so will be rejected, resulting in the choice
        # of a new point.
        try:
            cosmo.compute(["lensing"])
        except NameError:
            return data.boundary_loglike
        except (AttributeError, KeyboardInterrupt):
            io_mp.message("Something went terribly wrong with CLASS", "error")

    # For each desired likelihood, compute its value against the theoretical
    # model
    loglike = 0
    flag_wrote_fiducial = 0

    for likelihood in data.lkl.itervalues():
        if likelihood.need_update is True:
            value = likelihood.loglkl(cosmo, data)
            # Storing the result
            likelihood.backup_value = value
        # Otherwise, take the existing value
        else:
            value = likelihood.backup_value
        loglike += value
        # In case the fiducial file was written, store this information
        if value == 1:
            flag_wrote_fiducial += 1

    # Compute the derived parameters if relevant
    if data.get_mcmc_parameters(['derived']) != []:
        try:
            cosmo.get_current_derived_parameters(data)
        except NameError:
            print('Terminating now')
            exit()
    for elem in data.get_mcmc_parameters(['derived']):
        data.mcmc_parameters[elem]['current'] /= \
            data.mcmc_parameters[elem]['scale']

    # If fiducial files were created, inform the user, and exit
    if flag_wrote_fiducial > 0:
        if flag_wrote_fiducial == len(data.lkl):
            print '--> Fiducial file(s) was(were) created,',
            print 'please start a new chain'
            exit()
        else:
            print '--> Some previously non-existing fiducial files ',
            print 'were created, but potentially not all of them'
            print '--> Please check now manually on the headers ',
            print 'of the corresponding that all'
            print '--> parameters are coherent for your tested models'
            exit()

    return loglike


def read_args_from_chain(data, chain):
    """
    Pick up the last accepted values from an input chain as a starting point

    Function used only when the restart flag is set. It will simply read the
    last line of an input chain, using the tail command from the extended
    :class:`io_mp.File` class.

    :Parameters:
        * **chain** (`str`) - name of the input chain provided with the command
          line.

    .. warning::
        That method was not tested since the adding of derived parameters. The
        method :func:`read_args_from_bestfit` is the prefered one.

    .. warning::
        This method works because of the particular presentation of the chain,
        and the use of tabbings (not spaces). Please keep this in mind if you
        are having difficulties

    """
    Chain = io_mp.File(chain, 'r')
    parameter_names = data.get_mcmc_parameters(['varying'])

    i = 1
    for elem in parameter_names:
        data.mcmc_parameters[elem]['last_accepted'] = float(
            Chain.tail(1)[0].split('\t')[i])
        i += 1


def read_args_from_bestfit(data, bestfit):
    """
    Deduce the starting point either from the input file, or from a best fit
    file.

    :Parameters:
        * **bestfit** (`str`) - name of the bestfit file from the command line.
    """

    parameter_names = data.get_mcmc_parameters(['varying'])
    bestfit_file = open(bestfit, 'r')
    for line in bestfit_file:
        if line.find('#') != -1:
            bestfit_names = line.strip('#').replace(' ', '').\
                replace('\n', '').split(',')
            bestfit_values = np.zeros(len(bestfit_names), 'float64')
        else:
            line = line.split()
            for i in range(len(line)):
                bestfit_values[i] = line[i]

    print
    print('\nStarting point for rescaled parameters:')
    for elem in parameter_names:
        if elem in bestfit_names:
            data.mcmc_parameters[elem]['last_accepted'] = \
                bestfit_values[bestfit_names.index(elem)] / \
                data.mcmc_parameters[elem]['scale']
            print 'from best-fit file : ', elem, ' = ',
            print bestfit_values[bestfit_names.index(elem)] / \
                data.mcmc_parameters[elem]['scale']
        else:
            data.mcmc_parameters[elem]['last_accepted'] = \
                data.mcmc_parameters[elem]['initial'][0]
            print 'from input file    : ', elem, ' = ',
            print data.mcmc_parameters[elem]['initial'][0]


def get_covariance_matrix(data, command_line):
    """
    Compute the covariance matrix, from an input file or from an existing
    matrix.

    Reordering of the names and scaling take place here, in a serie of
    potentially hard to read methods. For the sake of clarity, and to avoid
    confusions, the code will, by default, print out a succession of 4
    covariance matrices at the beginning of the run, if starting from an
    existing one. This way, you can control that the paramters are set
    properly.

    .. note::

        The set of parameters from the run need not to be the exact same
        set of parameters from the existing covariance matrix (not even the
        ordering). Missing parameter from the existing covariance matrix will
        use the sigma given as an input.

    """

    # Setting numpy options in terms of precision (useful when writing to files
    # or displaying a result, but does not affect the precision of the
    # computation).
    np.set_printoptions(precision=2, linewidth=150)
    parameter_names = data.get_mcmc_parameters(['varying'])
    i = 0

    # if the user provides a .covmat file
    if command_line.cov is not None:
        cov = open('{0}'.format(command_line.cov), 'r')
        for line in cov:
            if line.find('#') != -1:
                # Extract the names from the first line
                covnames = line.strip('#').replace(' ', '').\
                    replace('\n', '').split(',')
                # Initialize the matrices
                M = np.zeros((len(covnames), len(covnames)), 'float64')
                rot = np.zeros((len(covnames), len(covnames)))
            else:
                line = line.split()
                for j in range(len(line)):
                    M[i][j] = np.array(line[j], 'float64')
                i += 1

        # First print out
        print('\nInput covariance matrix:')
        print(covnames)
        print(M)
        # Deal with the all problematic cases.
        # First, adjust the scales between stored parameters and the ones used
        # in mcmc
        scales = []
        for elem in covnames:
            if elem in parameter_names:
                scales.append(data.mcmc_parameters[elem]['scale'])
            else:
                scales.append(1)
        scales = np.diag(scales)
        # Compute the inverse matrix, and assert that the computation was
        # precise enough, by comparing the product to the identity matrix.
        invscales = np.linalg.inv(scales)
        np.testing.assert_array_almost_equal(
            np.dot(scales, invscales), np.eye(np.shape(scales)[0]),
            decimal=5)

        # Apply the newly computed scales to the input matrix
        M = np.dot(invscales.T, np.dot(M, invscales))

        # Second print out, after having applied the scale factors
        print('\nFirst treatment (scaling)')
        print(covnames)
        print(M)

        # Rotate M for the parameters to be well ordered, even if some
        # names are missing or some are in extra.
        # First, store the parameter names in temp_names that also appear in
        # the covariance matrix, in the right ordering for the code (might be
        # different from the input matri)
        temp_names = [elem for elem in parameter_names if elem in covnames]

        # If parameter_names contains less things than covnames, we will do a
        # small trick. Create a second temporary array, temp_names_2, that will
        # have the same dimension as covnames, and containing:
        # - the elements of temp_names, in the order of parameter_names (h index)
        # - an empty string '' for the remaining unused parameters
        temp_names_2 = []
        h = 0
        not_in = [elem for elem in covnames if elem not in temp_names]
        for k in range(len(covnames)):
            if covnames[k] not in not_in:
                temp_names_2.append(temp_names[h])
                h += 1
            else:
                temp_names_2.append('')

        # Create the rotation matrix, that will put the covariance matrix in
        # the right order, and also assign zeros to the unused parameters from
        # the input. These empty columns will be removed in the next step.
        for k in range(len(covnames)):
            for h in range(len(covnames)):
                try:
                    if covnames[k] == temp_names_2[h]:
                        rot[h][k] = 1.
                    else:
                        rot[h][k] = 0.
                except IndexError:
                    # The IndexError exception means that we are dealing with
                    # an unused parameter. By enforcing the corresponding
                    # rotation matrix element to 0, the resulting matrix will
                    # still have the same size as the original, but with zeros
                    # on the unused lines.
                    rot[h][k] = 0.
        M = np.dot(rot, np.dot(M, np.transpose(rot)))

        # Third print out
        print('\nSecond treatment (partial reordering and cleaning)')
        print(temp_names_2)
        print(M)

        # Final step, creating a temporary matrix, filled with 1, that will
        # eventually contain the result.
        M_temp = np.ones((len(parameter_names),
                          len(parameter_names)), 'float64')
        indices_final = np.zeros(len(parameter_names))
        indices_initial = np.zeros(len(covnames))
        # Remove names that are in parameter names but not in covnames, and
        # set to zero the corresponding columns of the final result.
        for k in range(len(parameter_names)):
            if parameter_names[k] in covnames:
                indices_final[k] = 1
        for zeros in np.where(indices_final == 0)[0]:
            M_temp[zeros, :] = 0
            M_temp[:, zeros] = 0
        # Remove names that are in covnames but not in param_names
        for h in range(len(covnames)):
            if covnames[h] in parameter_names:
                indices_initial[h] = 1
        # There, put a place holder number (we are using a pure imaginary
        # number: i, to avoid any problem) in the initial matrix, so that the
        # next step only copy the interesting part of the input to the final
        # matrix.
        for zeros in np.where(indices_initial == 0)[0]:
            M[zeros, :] = 1*j
            M[:, zeros] = 1*j
        # Now put in the temporary matrix, where the 1 were, the interesting
        # quantities from the input (the one that are not equal to i).
        M_temp[M_temp == 1] = M[M != 1*j]
        M = np.copy(M_temp)
        # on all other lines, that contain 0, just use sigma^2
        for zeros in np.where(indices_final == 0)[0]:
            M[zeros, zeros] = np.array(
                data.mcmc_parameters[parameter_names[zeros]]['initial'][3],
                'float64')**2

    # else, take sigmas^2.
    else:
        M = np.identity(len(parameter_names), 'float64')
        for elem in parameter_names:
            M[i][i] = np.array(
                data.mcmc_parameters[elem]['initial'][3], 'float64')**2
            i += 1

    # Final print out, the actually used covariance matrix
    sys.stdout.write('\nDeduced starting covariance matrix:\n')
    print(parameter_names)
    print(M)

    #inverse, and diagonalization
    eigv, eigV = np.linalg.eig(np.linalg.inv(M))
    return eigv, eigV, M


def get_new_position(data, eigv, U, k, Cholesky, Inverse_Cholesky, Rotation):
    """
    Obtain a new position in the parameter space from the eigen values of the
    inverse covariance matrix, or from the Cholesky decomposition (original
    idea by Anthony Lewis, in a yet unpublished paper, Efficient sampling of
    fast and slow cosmological parameters).

    .. note::

        U, eigv are not used anymore in v1.2.0, but might come back in v1.2.1.

    :Parameters:
        * **eigv** (`numpy array`) - eigenvalues previously computed *obsolete
                in v1.2.0*
        * **U** (`numpy_array`) - *obsolete in v1.2.0*, was the covariance
                matrix.
        * **k** (`int`) - Number of points so far in the chain, is used to
                rotate through parameters
        * **Cholesky** (`numpy_array`) - Cholesky decomposition of the
                covariance matrix, and its inverse
        * **Rotation** (`numpy_array`) - Not used yet

    """

    parameter_names = data.get_mcmc_parameters(['varying'])
    vector_new = np.zeros(len(parameter_names), 'float64')
    sigmas = np.zeros(len(parameter_names), 'float64')

    # Write the vector of last accepted points, or if it does not exist
    # (initialization routine), take the mean value
    vector = np.zeros(len(parameter_names), 'float64')
    try:
        for elem in parameter_names:
            vector[parameter_names.index(elem)] = \
                data.mcmc_parameters[elem]['last_accepted']
    except KeyError:
        for elem in parameter_names:
            vector[parameter_names.index(elem)] = \
                data.mcmc_parameters[elem]['initial'][0]

    # Initialize random seed
    rd.seed()

    # Choice here between sequential and global change of direction
    if data.jumping == 'global':
        for i in range(len(vector)):
            sigmas[i] = (math.sqrt(1/eigv[i]/len(vector))) * \
                rd.gauss(0, 1)*data.jumping_factor
    elif data.jumping == 'sequential':
        i = k % len(vector)
        sigmas[i] = (math.sqrt(1/eigv[i]))*rd.gauss(0, 1)*data.jumping_factor
    elif data.jumping == 'fast':
        i = k % len(vector)
        ###############
        # method fast+global
        for elem in data.blocks_parameters:
            if i < elem:
                index = data.blocks_parameters.index(elem)
                if index == 0:
                    Range = elem
                    Previous = 0
                else:
                    Range = elem-data.blocks_parameters[index-1]
                    Previous = data.blocks_parameters[index-1]
                for j in range(Range):
                    sigmas[j+Previous] = (math.sqrt(1./Range)) * \
                        rd.gauss(0, 1)*data.jumping_factor
                break
            else:
                continue
        ####################
        # method fast+sequential
        #sigmas[i] = rd.gauss(0,1)*data.jumping_factor
    else:
        print('\n\n Jumping method unknown (accepted : ')
        print('global (default), sequential, fast)')

    # Fill in the new vector
    if data.jumping in ['global', 'sequential']:
        vector_new = vector + np.dot(U, sigmas)
    else:
        vector_new = vector + np.dot(Cholesky, sigmas)

    # Check for boundaries problems
    flag = 0
    for elem in parameter_names:
        i = parameter_names.index(elem)
        value = data.mcmc_parameters[elem]['initial']
        if((str(value[1]) != str(-1) and value[1] is not None) and
                (vector_new[i] < value[1])):
            flag += 1  # if a boundary value is reached, increment
        elif((str(value[2]) != str(-1) and value[1] is not None) and
                vector_new[i] > value[2]):
            flag += 1  # same

    # At this point, if a boundary condition is not fullfilled, ie, if flag is
    # different from zero, return False
    if flag != 0:
        return False

    # Check for a slow step (only after the first time, so we put the test in a
    # try: statement: the first time, the exception KeyError will be raised)
    try:
        data.check_for_slow_step(vector_new)
    except KeyError:
        pass

    # If it is not the case, proceed with normal computation. The value of
    # new_vector is then put into the 'current' point in parameter space.
    for elem in parameter_names:
        i = parameter_names.index(elem)
        data.mcmc_parameters[elem]['current'] = vector_new[i]

    # Propagate the information towards the cosmo arguments
    data.update_cosmo_arguments()

    return True


def accept_step(data):
    """
    Transfer the 'current' point in the varying parameters to the last accepted
    one.

    """
    for elem in data.get_mcmc_parameters(['varying']):
        data.mcmc_parameters[elem]['last_accepted'] = \
            data.mcmc_parameters[elem]['current']
    for elem in data.get_mcmc_parameters(['derived']):
        data.mcmc_parameters[elem]['last_accepted'] = \
            data.mcmc_parameters[elem]['current']


######################
# MCMC CHAIN
######################
def chain(cosmo, data, command_line):
    """
    Run a Markov chain of fixed length.

    Main function of this module, this is the actual Markov chain procedure.
    After having selected a starting point in parameter space defining the
    first **last accepted** one, it will, for a given amount of steps :

    + choose randomnly a new point following the *proposal density*,
    + compute the cosmological *observables* through the cosmological module,
    + compute the value of the *likelihoods* of the desired experiments at this point,
    + *accept/reject* this point given its likelihood compared to the one of
      the last accepted one.

    Every time the code accepts :code:`data.write_step` number of points
    (quantity defined in the input parameter file), it will write the result to
    disk (flushing the buffer by forcing to exit the output file, and reopen it
    again.

    .. note::

        to use the code to set a fiducial file for certain fixed parameters,
        you can use two solutions. The first one is to put all input 1-sigma
        proposal density to zero (this method still works, but is not
        recommended anymore). The second one consist in using the flag "-f 0",
        to force a step of zero amplitude.

    """

    ## Initialisation
    loglike = 0

    # Recover the covariance matrix according to the input, if the varying set
    # of parameters is non-zero
    if (data.get_mcmc_parameters(['varying']) != []):
        sigma_eig, U, C = get_covariance_matrix(data, command_line)
        if data.jumping_factor == 0:
            io_mp.message(
                "The jumping factor has been set to 0. The above covariance \
                matrix will not be used.",
                "info")

    # In case of a fiducial run (all parameters fixed), simply run once and
    # print out the likelihood. This should not be used any more (one has to
    # modify the log.param, which is never a good idea. Instead, force the code
    # to use a jumping factor of 0 with the option "-f 0".
    else:
        io_mp.message(
            "You are running with no varying parameters... I will compute \
            only one point and exit",
            "info")
        data.update_cosmo_arguments()  # this fills in the fixed parameters
        loglike = compute_lkl(cosmo, data)
        io_mp.print_vector([data.out, sys.stdout], 1, loglike, data)
        return 1, loglike

    # In the fast-slow method, one need the Cholesky decomposition of the
    # covariance matrix. Return the Cholesky decomposition as a lower
    # triangular matrix
    Cholesky = None
    Inverse_Cholesky = None
    Rotation = None
    if command_line.jumping == 'fast':
        Cholesky = la.cholesky(C).T
        Inverse_Cholesky = np.linalg.inv(Cholesky)
        Rotation = np.identity(len(sigma_eig))

    # If restart wanted, pick initial value for arguments
    if command_line.restart is not None:
        read_args_from_chain(data, command_line.restart)

    # If restart from best fit file, read first point (overwrite settings of
    # read_args_from_chain)
    if command_line.bf is not None:
        read_args_from_bestfit(data, command_line.bf)

    # Pick a position (from last accepted point if restart, from the mean value
    # else), with a 100 tries.
    for i in range(100):
        if get_new_position(data, sigma_eig, U, i,
                            Cholesky, Inverse_Cholesky, Rotation) is True:
            break
        if i == 99:
            io_mp.message(
                "You should probably check your prior boundaries... because \
                no valid starting position was found after 100 tries",
                "error")

    # Compute the starting Likelihood
    loglike = compute_lkl(cosmo, data)

    # Choose this step as the last accepted value
    # (accept_step), and modify accordingly the max_loglike
    accept_step(data)
    max_loglike = loglike

    # If the jumping factor is 0, the likelihood associated with this point is
    # displayed, and the code exits.
    if data.jumping_factor == 0:
        io_mp.print_vector([data.out, sys.stdout], 1, loglike, data)
        return 1, loglike

    acc, rej = 0.0, 0.0  # acceptance and rejection number count
    N = 1   # number of time the system stayed in the current position

    # Print on screen the computed parameters
    io_mp.print_parameters(sys.stdout, data)

    k = 1
    # Main loop, that goes on while the maximum number of failure is not
    # reached, and while the expected amount of steps (N) is not taken.
    while k <= command_line.N:

        # Pick a new position ('current' flag in mcmc_parameters), and compute
        # its likelihood. If get_new_position returns True, it means it did not
        # encounter any boundary problem. Otherwise, just increase the
        # multiplicity of the point and start the loop again
        if get_new_position(
                data, sigma_eig, U, k, Cholesky,
                Inverse_Cholesky, Rotation) is True:
            newloglike = compute_lkl(cosmo, data)
        else:  # reject step
            rej += 1
            N += 1
            k += 1
            continue

        # Harmless trick to avoid exponentiating large numbers. This decides
        # whether or not the system should move.
        if (newloglike != data.boundary_loglike):
            if (newloglike >= loglike):
                alpha = 1.
            else:
                alpha = np.exp(newloglike-loglike)
        else:
            alpha = -1

        if ((alpha == 1.) or (rd.uniform(0, 1) < alpha)):  # accept step

            # Print out the last accepted step (WARNING: this is NOT the one we
            # just computed ('current' flag), but really the previous one.)
            # with its proper multiplicity (number of times the system stayed
            # there).
            io_mp.print_vector([data.out, sys.stdout], N, loglike, data)

            # Report the 'current' point to the 'last_accepted'
            accept_step(data)
            loglike = newloglike
            if loglike > max_loglike:
                max_loglike = loglike
            acc += 1.0
            N = 1  # Reset the multiplicity

        else:  # reject step
            rej += 1.0
            N += 1  # Increase multiplicity of last accepted point

        # Regularly (option to set in parameter file), close and reopen the
        # buffer to force to write on file.
        if acc % data.write_step == 0:
            io_mp.refresh_file(data)
        k += 1  # One iteration done
    # END OF WHILE LOOP

    # If at this moment, the multiplicity is higher than 1, it means the
    # current point is not yet accepted, but it also mean that we did not print
    # out the last_accepted one yet. So we do.
    if N > 1:
        io_mp.print_vector([data.out, sys.stdout], N-1, loglike, data)

    # Print out some information on the finished chain
    rate = acc / (acc + rej)
    sys.stdout.write('\n#  {0} steps done, acceptance rate: {1}\n'.
                     format(command_line.N, rate))

    # For a restart, erase the starting point to keep only the new, longer
    # chain.
    if command_line.restart is not None:
        os.remove(command_line.restart)
        sys.stdout.write('    deleting starting point of the chain {0}\n'.
                         format(command_line.restart))

    return
