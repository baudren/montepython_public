"""
.. module:: mcmc
   :synopsis: Monte Carlo procedure
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>

This module defines one key function, :func:`chain`, that handles the Markov
chain. So far, the code uses only one chain, as no parallelization is done.

The following routine is also defined in this module, which is called at
every step:

* :func:`get_new_position` returns a new point in the parameter space,
  depending on the proposal density.

The :func:`chain` in turn calls several helper routines, defined in
:mod:`sampler`. These are called just once:

* :func:`compute_lkl() <sampler.compute_lkl>` is called at every step in the Markov chain, returning
  the likelihood at the current point in the parameter space.
* :func:`get_covariance_matrix() <sampler.get_covariance_matrix>`
* :func:`read_args_from_chain() <sampler.read_args_from_chain>`
* :func:`read_args_from_bestfit() <sampler.read_args_from_bestfit>`
* :func:`accept_step() <sampler.accept_step>`

Their usage is described in :mod:`sampler`. On the contrary, the following
routines are called at every step:

The arguments of these functions will often contain **data** and/or **cosmo**.
They are both initialized instances of respectively :class:`data` and the
cosmological class. They will thus not be described for every function.
"""

import os
import sys
import math
import random as rd
import numpy as np
import warnings
import scipy.linalg as la
from pprint import pprint

import io_mp
import sampler


def get_new_position(data, eigv, U, k, Cholesky, Rotation):
    """
    Obtain a new position in the parameter space from the eigen values of the
    inverse covariance matrix, or from the Cholesky decomposition (original
    idea by Anthony Lewis, in `Efficient sampling of fast and slow
    cosmological parameters <http://arxiv.org/abs/1304.4473>`_ )

    The three different jumping options, decided when starting a run with the
    flag **-j**  are **global**, **sequential** and **fast** (by default) (see
    :mod:`parser_mp` for reference).

    .. warning::

        For running Planck data, the option **fast** is highly recommended, as
        it speeds up the convergence. Note that when using this option, the
        list of your likelihoods in your parameter file **must match** the
        ordering of your nuisance parameters (as always, they must come after
        the cosmological parameters, but they also must be ordered between
        likelihood, with, preferentially, the slowest likelihood to compute
        coming first).


    - **global**: varies all the parameters at the same time. Depending on the
      input covariance matrix, some degeneracy direction will be followed,
      otherwise every parameter will jump independently of each other.
    - **sequential**: varies every parameter sequentially. Works best when
      having no clue about the covariance matrix, or to understand which
      estimated sigma is wrong and slowing down the whole process.
    - **fast**: privileged method when running the Planck likelihood. Described
      in the aforementioned article, it separates slow (cosmological) and fast
      (nuisance) parameters.

    Parameters
    ----------
    eigv : numpy array
        Eigenvalues previously computed
    U : numpy_array
        Covariance matrix.
    k : int
        Number of points so far in the chain, is used to rotate through
        parameters
    Cholesky : numpy array
        Cholesky decomposition of the covariance matrix, and its inverse
    Rotation : numpy_array
        Not used yet

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
        #i = k % len(vector)
        j = k % len(data.over_sampling_indices)
        i = data.over_sampling_indices[j]
        ###############
        # method fast+global
        for index, elem in enumerate(data.block_parameters):
            # When the running index is below the maximum index of a block of
            # parameters, this block is varied, and **only this one** (note the
            # break at the end of the if clause, it is not a continue)
            if i < elem:
                if index == 0:
                    Range = elem
                    Previous = 0
                else:
                    Range = elem-data.block_parameters[index-1]
                    Previous = data.block_parameters[index-1]
                # All the varied parameters are given a random variation with a
                # sigma of 1. This will translate in a jump for all the
                # parameters (as long as the Cholesky matrix is non diagonal)
                for j in range(Range):
                    sigmas[j+Previous] = (math.sqrt(1./Range)) * \
                        rd.gauss(0, 1)*data.jumping_factor
                break
            else:
                continue
    else:
        print('\n\n Jumping method unknown (accepted : ')
        print('global, sequential, fast (default))')

    # Fill in the new vector
    if data.jumping in ['global', 'sequential']:
        vector_new = vector + np.dot(U, sigmas)
    else:
        vector_new = vector + np.dot(Cholesky, sigmas)

    # Check for boundaries problems
    flag = 0
    for i, elem in enumerate(parameter_names):
        value = data.mcmc_parameters[elem]['initial']
        if((str(value[1]) != str(-1) and value[1] is not None) and
                (vector_new[i] < value[1])):
            flag += 1  # if a boundary value is reached, increment
        elif((str(value[2]) != str(-1) and value[2] is not None) and
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
    for index, elem in enumerate(parameter_names):
        data.mcmc_parameters[elem]['current'] = vector_new[index]

    # Propagate the information towards the cosmo arguments
    data.update_cosmo_arguments()

    return True


######################
# MCMC CHAIN
######################
def chain(cosmo, data, command_line):
    """
    Run a Markov chain of fixed length with a Metropolis Hastings algorithm.

    Main function of this module, this is the actual Markov chain procedure.
    After having selected a starting point in parameter space defining the
    first **last accepted** one, it will, for a given amount of steps :

    + choose randomnly a new point following the *proposal density*,
    + compute the cosmological *observables* through the cosmological module,
    + compute the value of the *likelihoods* of the desired experiments at this
      point,
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

    # In case command_line.silent has been asked, outputs should only contain
    # data.out. Otherwise, it will also contain sys.stdout
    outputs = [data.out]
    if not command_line.silent:
        outputs.append(sys.stdout)

    # check for MPI
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        # suppress duplicate output from slaves
        if rank:
            command_line.quiet = True
    except ImportError:
        # set all chains to master if no MPI
        rank = 0

    # Recover the covariance matrix according to the input, if the varying set
    # of parameters is non-zero
    if (data.get_mcmc_parameters(['varying']) != []):
        sigma_eig, U, C = sampler.get_covariance_matrix(cosmo, data, command_line)
        if data.jumping_factor == 0:
            warnings.warn(
                "The jumping factor has been set to 0. The above covariance " +
                "matrix will not be used.")

    # In case of a fiducial run (all parameters fixed), simply run once and
    # print out the likelihood. This should not be used any more (one has to
    # modify the log.param, which is never a good idea. Instead, force the code
    # to use a jumping factor of 0 with the option "-f 0".
    else:
        warnings.warn(
            "You are running with no varying parameters... I will compute " +
            "only one point and exit")
        data.update_cosmo_arguments()  # this fills in the fixed parameters
        loglike = sampler.compute_lkl(cosmo, data)
        io_mp.print_vector(outputs, 1, loglike, data)
        return 1, loglike

    # In the fast-slow method, one need the Cholesky decomposition of the
    # covariance matrix. Return the Cholesky decomposition as a lower
    # triangular matrix
    Cholesky = None
    Rotation = None
    if command_line.jumping == 'fast':
        Cholesky = la.cholesky(C).T
        Rotation = np.identity(len(sigma_eig))

    # If the update mode was selected, the previous (or original) matrix should be stored
    if command_line.update:
        previous = (sigma_eig, U, C, Cholesky)

    # If restart wanted, pick initial value for arguments
    if command_line.restart is not None:
        sampler.read_args_from_chain(data, command_line.restart)

    # If restart from best fit file, read first point (overwrite settings of
    # read_args_from_chain)
    if command_line.bf is not None:
        sampler.read_args_from_bestfit(data, command_line.bf)

    # Pick a position (from last accepted point if restart, from the mean value
    # else), with a 100 tries.
    for i in range(100):
        if get_new_position(data, sigma_eig, U, i,
                            Cholesky, Rotation) is True:
            break
        if i == 99:
            raise io_mp.ConfigurationError(
                "You should probably check your prior boundaries... because " +
                "no valid starting position was found after 100 tries")

    # Compute the starting Likelihood
    loglike = sampler.compute_lkl(cosmo, data)

    # Choose this step as the last accepted value
    # (accept_step), and modify accordingly the max_loglike
    sampler.accept_step(data)
    max_loglike = loglike

    # If the jumping factor is 0, the likelihood associated with this point is
    # displayed, and the code exits.
    if data.jumping_factor == 0:
        io_mp.print_vector(outputs, 1, loglike, data)
        return 1, loglike

    acc, rej = 0.0, 0.0  # acceptance and rejection number count
    N = 1   # number of time the system stayed in the current position

    # define path and covmat
    input_covmat = command_line.cov
    base = os.path.basename(command_line.folder)
    # the previous line fails when "folder" is a string ending with a slash. This issue is cured by the next lines:
    if base == '':
        base = os.path.basename(command_line.folder[:-1])
    command_line.cov = os.path.join(
        command_line.folder, base+'.covmat')

    # Print on screen the computed parameters
    if not command_line.silent and not command_line.quiet:
        io_mp.print_parameters(sys.stdout, data)

    # Suppress non-informative output after initializing
    command_line.quiet = True

    k = 1
    # Main loop, that goes on while the maximum number of failure is not
    # reached, and while the expected amount of steps (N) is not taken.
    while k <= command_line.N:

        # If the number of steps reaches the number set in the update method,
        # then the proposal distribution should be adapted.
        if command_line.update:

            # master chain behavior
            if not rank:
                # Add the folder to the list of files to analyze, and switch on the
                # options for computing only the covmat
                from parser_mp import parse
                info_command_line = parse(
                    'info %s --minimal --noplot --keep-fraction 0.5 --keep-non-markovian --want-covmat' % command_line.folder)
                info_command_line.update = command_line.update
                # the +10 below is here to ensure that the first master update will take place before the first slave updates,
                # but this is a detail, the code is robust against situations where updating is not possible, so +10 could be omitted
                if not (k+10) % command_line.update and k > 10:
                    # Try to launch an analyze
                    try:
                        from analyze import analyze
                        R_minus_one = analyze(info_command_line)
                    except:
                        if not command_line.silent:
                            print 'Step ',k,' chain ', rank,': Failed to calculate covariant matrix'
                        pass

                if not (k-1) % command_line.update:
                    try:
                        # Read the covmat
                        sigma_eig, U, C = sampler.get_covariance_matrix(
                            cosmo, data, command_line)
                        if command_line.jumping == 'fast':
                            Cholesky = la.cholesky(C).T
                        # Test here whether the covariance matrix has really changed
                        # We should in principle test all terms, but testing the first one should suffice
                        if not C[0,0] == previous[2][0,0]:
                            previous = (sigma_eig, U, C, Cholesky)
                            if k == 1:
                                if not command_line.silent:
                                    if not input_covmat == None:
                                        warnings.warn(
                                            'Appending to an existing folder: using %s instead of %s. '
                                            'If new input covmat is desired, please delete previous covmat.'
                                            % (command_line.cov, input_covmat))
                                    else:
                                        warnings.warn(
                                            'Appending to an existing folder: using %s. '
                                            'If no starting covmat is desired, please delete previous covmat.'
                                            % command_line.cov)
                            else:
                                data.out.write('# After %d accepted steps: update proposal with max(R-1) = %f \n' % (int(acc), max(R_minus_one)))
                                if not command_line.silent:
                                    print 'After %d accepted steps: update proposal with max(R-1) = %f \n' % (int(acc), max(R_minus_one))
                                try:
                                    if stop-after-update:
                                        k = command_line.N
                                        print 'Covariant matrix updated - stopping run'
                                except:
                                    pass

                    except:
                        pass

                    command_line.quiet = True

            # slave chain behavior
            else:
                if not (k-1) % command_line.update:
                    try:
                        sigma_eig, U, C = sampler.get_covariance_matrix(
                            cosmo, data, command_line)
                        if command_line.jumping == 'fast':
                            Cholesky = la.cholesky(C).T
                        # Test here whether the covariance matrix has really changed
                        # We should in principle test all terms, but testing the first one should suffice
                        if not C[0,0] == previous[2][0,0] and not k == 1:
                            data.out.write('# After %d accepted steps: update proposal \n' % int(acc))
                            if not command_line.silent:
                                print 'After %d accepted steps: update proposal \n' % int(acc)
                            try:
                                if stop_after_update:
                                    k = command_line.N
                                    print 'Covariant matrix updated - stopping run'
                            except:
                                pass
                        previous = (sigma_eig, U, C, Cholesky)

                    except IOError:
                        pass

        # Pick a new position ('current' flag in mcmc_parameters), and compute
        # its likelihood. If get_new_position returns True, it means it did not
        # encounter any boundary problem. Otherwise, just increase the
        # multiplicity of the point and start the loop again
        if get_new_position(
                data, sigma_eig, U, k, Cholesky, Rotation) is True:
            newloglike = sampler.compute_lkl(cosmo, data)
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
            io_mp.print_vector(outputs, N, loglike, data)

            # Report the 'current' point to the 'last_accepted'
            sampler.accept_step(data)
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
            # Update the outputs list
            outputs[0] = data.out
        k += 1  # One iteration done
    # END OF WHILE LOOP

    # If at this moment, the multiplicity is higher than 1, it means the
    # current point is not yet accepted, but it also mean that we did not print
    # out the last_accepted one yet. So we do.
    if N > 1:
        io_mp.print_vector(outputs, N-1, loglike, data)

    # Print out some information on the finished chain
    rate = acc / (acc + rej)
    sys.stdout.write('\n#  {0} steps done, acceptance rate: {1}\n'.
                     format(command_line.N, rate))

    # In case the acceptance rate is too low, or too high, print a warning
    if rate < 0.05:
        warnings.warn("The acceptance rate is below 0.05. You might want to "
                      "set the jumping factor to a lower value than the "
                      "default (2.4), with the option `-f 1.5` for instance.")
    elif rate > 0.6:
        warnings.warn("The acceptance rate is above 0.6, which means you might"
                      " have difficulties exploring the entire parameter space"
                      ". Try analysing these chains, and use the output "
                      "covariance matrix to decrease the acceptance rate to a "
                      "value between 0.2 and 0.4 (roughly).")

    # For a restart, erase the starting point to keep only the new, longer
    # chain.
    if command_line.restart is not None:
        os.remove(command_line.restart)
        sys.stdout.write('    deleting starting point of the chain {0}\n'.
                         format(command_line.restart))

    return
