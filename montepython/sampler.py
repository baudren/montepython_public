"""
.. module:: sampler
    :synopsis: Generic sampler
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>
.. moduleauthor:: Surhudm More <>

This module defines one key function, :func:`run`, that distributes the work to
the desired actual sampler (Metropolis Hastings, or Nested Sampling so far).

It also defines a serie of helper functions, that aim to be generically used by
all different sampler methods:

* :func:`get_covariance_matrix`
* :func:`read_args_from_chain`
* :func:`read_args_from_bestfit`
* :func:`accept_step`
* :func:`compute_lkl`


"""
import numpy as np
import sys
import warnings

import io_mp
import os


def run(cosmo, data, command_line):
    """
    Depending on the choice of sampler, dispatch the appropriate information

    The :mod:`mcmc` module is used as previously, except the call to
    :func:`mcmc.chain`, or :func:`nested_sampling.run` is now within
    this function, instead of from within :mod:`MontePython`.

    In the long term, this function should contain any potential hybrid scheme.

    """

    if command_line.method == 'MH':
        import mcmc
        mcmc.chain(cosmo, data, command_line)
        data.out.close()
    elif command_line.method == 'NS':
        import nested_sampling as ns
        ns.run(cosmo, data, command_line)
    elif command_line.method == 'CH':
        import cosmo_hammer as hammer
        hammer.run(cosmo, data, command_line)
    elif command_line.method == 'IS':
        import importance_sampling as ims
        ims.run(cosmo, data, command_line)
    elif command_line.method == 'Der':
        import add_derived as der
        der.run(cosmo, data, command_line)
    else:
        raise io_mp.ConfigurationError(
            "Sampling method %s not understood" % command_line.method)


def read_args_from_chain(data, chain):
    """
    Pick up the last accepted values from an input chain as a starting point

    Function used only when the restart flag is set. It will simply read the
    last line of an input chain, using the tail command from the extended
    :class:`io_mp.File` class.

    .. warning::
        That method was not tested since the adding of derived parameters. The
        method :func:`read_args_from_bestfit` is the prefered one.

    .. warning::
        This method works because of the particular presentation of the chain,
        and the use of tabbings (not spaces). Please keep this in mind if you
        are having difficulties

    Parameters
    ----------
    chain : str
        Name of the input chain provided with the command line.

    """
    chain_file = io_mp.File(chain, 'r')
    parameter_names = data.get_mcmc_parameters(['varying'])

    i = 1
    for elem in parameter_names:
        data.mcmc_parameters[elem]['last_accepted'] = float(
            chain_file.tail(1)[0].split('\t')[i])
        i += 1


def read_args_from_bestfit(data, bestfit):
    """
    Deduce the starting point either from the input file, or from a best fit
    file.

    Parameters
    ----------
    bestfit : str
        Name of the bestfit file from the command line.

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


def get_covariance_matrix(cosmo, data, command_line):
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

    # Define quiet setting if not previously defined
    try:
        command_line.quiet
    except:
        command_line.quiet = False

    if command_line.fisher and not command_line.cov:
        # We will work out the fisher matrix for all the parameters and
        # write it to a file
        if not command_line.silent:
            warnings.warn("Fisher implementation is being tested")

        # Let us create a separate copy of data
        from copy import deepcopy
        # Do not modify data, instead copy
        temp_data = deepcopy(data)
        done = False

        # Create the center dictionary, which will hold the center point
        # information (or best-fit) TODO
        # This dictionary will be updated in case it was too far from the
        # best-fit, and found a non positive-definite symmetric fisher matrix.
        center = {}
        if not command_line.bf:
            for elem in parameter_names:
                temp_data.mcmc_parameters[elem]['current'] = (
                    data.mcmc_parameters[elem]['initial'][0])
                center[elem] = data.mcmc_parameters[elem]['initial'][0]
        else:
            read_args_from_bestfit(temp_data, command_line.bf)
            for elem in parameter_names:
                temp_data.mcmc_parameters[elem]['current'] = (
                    temp_data.mcmc_parameters[elem]['last_accepted'])
                center[elem] = temp_data.mcmc_parameters[elem]['last_accepted']

        # Have a security index that prevents looping indefinitely
        security = 0
        while not done and security < 10:
            security += 1
            # Compute the Fisher matrix and the gradient array at the center
            # point.
            fisher_matrix, gradient = compute_fisher(
                temp_data, cosmo, center, 0.01)

            # Compute inverse of the fisher matrix, catch LinAlgError exception
            fisher_invert_success = True
            try:
                if not command_line.silent:
                    print("Fisher matrix computed:")
                    print(fisher_matrix)
                cov_matrix = np.linalg.inv(fisher_matrix)
            except np.linalg.LinAlgError:
                raise io_mp.ConfigurationError(
                    "Could not find Fisher matrix, please remove the "
                    "option --fisher and run with Metropolis-Hastings "
                    "or another sampling method.")
                fisher_invert_success = False
                done = True

            # Write it to the file
            if fisher_invert_success:
                io_mp.write_covariance_matrix(
                    cov_matrix, parameter_names,
                    os.path.join(command_line.folder, 'covariance_fisher.mat'))

                command_line.cov = os.path.join(
                    command_line.folder, 'covariance_fisher.mat')

                done = True
                # Check if the diagonal elements are non-negative
                for h, elem in enumerate(parameter_names):
                    if cov_matrix[h][h] < 0:
                        warnings.warn(
                            "Covariance has negative values on diagonal, "
                            "moving to a better point and repeating "
                            "the Fisher computation")
                        done = False
                        break

                if not done:
                    # Solve for a step
                    step = np.dot(cov_matrix, gradient)
                    # Now modify data_parameters TODO HERE update center
                    for k, elem in enumerate(parameter_names):
                        data.mcmc_parameters[elem]['initial'][0] = data.mcmc_parameters[elem]['initial'][0]-step[k]
                        temp_data.mcmc_parameters[elem]['initial'][0] = temp_data.mcmc_parameters[elem]['initial'][0]-step[k]
                        print "Moved %s to:"%(elem),data.mcmc_parameters[elem]['initial'][0]

    # if the user provides a .covmat file or if user asks to compute a fisher matrix
    if command_line.cov is not None:

        cov = open('{0}'.format(command_line.cov), 'r')

        i = 0
        for line in cov:
            if line.find('#') != -1:
                # Extract the names from the first line
                covnames = line.strip('#').replace(' ', '').\
                    replace('\n', '').split(',')
                # Initialize the matrices
                matrix = np.zeros((len(covnames), len(covnames)), 'float64')
                rot = np.zeros((len(covnames), len(covnames)))
            else:
                line = line.split()
                for j in range(len(line)):
                    matrix[i][j] = np.array(line[j], 'float64')
                i += 1

        # First print out
        if not command_line.silent and not command_line.quiet:
            print('\nInput covariance matrix:')
            print(covnames)
            print(matrix)
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
        matrix = np.dot(invscales.T, np.dot(matrix, invscales))

        # Second print out, after having applied the scale factors
        if not command_line.silent and not command_line.quiet:
            print('\nFirst treatment (scaling)')
            print(covnames)
            print(matrix)

        # Rotate matrix for the parameters to be well ordered, even if some
        # names are missing or some are in extra.
        # First, store the parameter names in temp_names that also appear in
        # the covariance matrix, in the right ordering for the code (might be
        # different from the input matri)
        temp_names = [elem for elem in parameter_names if elem in covnames]

        # If parameter_names contains less things than covnames, we will do a
        # small trick. Create a second temporary array, temp_names_2, that will
        # have the same dimension as covnames, and containing:
        # - the elements of temp_names, in the order of parameter_names (h
        # index)
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
        matrix = np.dot(rot, np.dot(matrix, np.transpose(rot)))

        # Third print out
        if not command_line.silent and not command_line.quiet:
            print('\nSecond treatment (partial reordering and cleaning)')
            print(temp_names_2)
            print(matrix)

        # Final step, creating a temporary matrix, filled with 1, that will
        # eventually contain the result.
        matrix_temp = np.ones((len(parameter_names),
                               len(parameter_names)), 'float64')
        indices_final = np.zeros(len(parameter_names))
        indices_initial = np.zeros(len(covnames))
        # Remove names that are in parameter names but not in covnames, and
        # set to zero the corresponding columns of the final result.
        for k in range(len(parameter_names)):
            if parameter_names[k] in covnames:
                indices_final[k] = 1
        for zeros in np.where(indices_final == 0)[0]:
            matrix_temp[zeros, :] = 0
            matrix_temp[:, zeros] = 0
        # Remove names that are in covnames but not in param_names
        for h in range(len(covnames)):
            if covnames[h] in parameter_names:
                indices_initial[h] = 1
        # There, put a place holder number (we are using a pure imaginary
        # number: i, to avoid any problem) in the initial matrix, so that the
        # next step only copy the interesting part of the input to the final
        # matrix.
        max_value = np.finfo(np.float64).max
        for zeros in np.where(indices_initial == 0)[0]:
            matrix[zeros, :] = [max_value for _ in range(
                len(matrix[zeros, :]))]
            matrix[:, zeros] = [max_value for _ in range(
                len(matrix[:, zeros]))]
        # Now put in the temporary matrix, where the 1 were, the interesting
        # quantities from the input (the one that are not equal to i).
        matrix_temp[matrix_temp == 1] = matrix[matrix != max_value]
        matrix = np.copy(matrix_temp)
        # on all other lines, that contain 0, just use sigma^2
        for zeros in np.where(indices_final == 0)[0]:
            matrix[zeros, zeros] = np.array(
                data.mcmc_parameters[parameter_names[zeros]]['initial'][3],
                'float64')**2
    # else, take sigmas^2.
    else:
        matrix = np.identity(len(parameter_names), 'float64')
        for index, elem in enumerate(parameter_names):
            matrix[index][index] = np.array(
                data.mcmc_parameters[elem]['initial'][3], 'float64')**2


    # Final print out, the actually used covariance matrix
    if not command_line.silent and not command_line.quiet:
        sys.stdout.write('\nDeduced starting covariance matrix:\n')
        print(parameter_names)
        print(matrix)

    #inverse, and diagonalization
    eigv, eigV = np.linalg.eig(np.linalg.inv(matrix))
    return eigv, eigV, matrix


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


def check_flat_bound_priors(parameters, names):
    """
    Ensure that all varying parameters are bound and flat

    It is a necessary condition to use the code with Nested Sampling or the
    Cosmo Hammer.
    """
    is_flat = all(parameters[name]['prior'].prior_type == 'flat'
                  for name in names)
    is_bound = all(parameters[name]['prior'].is_bound()
                   for name in names)
    return is_flat, is_bound


def compute_lkl(cosmo, data):
    """
    Compute the likelihood, given the current point in parameter space.

    This function now performs a test before calling the cosmological model
    (**new in version 1.2**). If any cosmological parameter changed, the flag
    :code:`data.need_cosmo_update` will be set to :code:`True`, from the
    routine :func:`check_for_slow_step <data.Data.check_for_slow_step>`.

    Returns
    -------
    loglike : float
        The log of the likelihood (:math:`\\frac{-\chi^2}2`) computed from the
        sum of the likelihoods of the experiments specified in the input
        parameter file.

        This function returns :attr:`data.boundary_loglkie
        <data.data.boundary_loglike>`, defined in the module :mod:`data` if
        *i)* the current point in the parameter space has hit a prior edge, or
        *ii)* the cosmological module failed to compute the model. This value
        is chosen to be extremly small (large negative value), so that the step
        will always be rejected.


    """
    from classy import CosmoSevereError, CosmoComputationError

    # If the cosmological module has already been called once, and if the
    # cosmological parameters have changed, then clean up, and compute.
    if cosmo.state and data.need_cosmo_update is True:
        cosmo.struct_cleanup()

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
        # - CosmoSevereError is returned if a parameter was not properly set
        # during the initialisation (for instance, you entered Ommega_cdm
        # instead of Omega_cdm).  Then, the code exits, to prevent running with
        # imaginary parameters. This behaviour is also used in case you want to
        # kill the process.
        # - CosmoComputationError is returned if Class fails to compute the
        # output given the parameter values. This will be considered as a valid
        # point, but with minimum likelihood, so will be rejected, resulting in
        # the choice of a new point.
        try:
            cosmo.compute(["lensing"])
        except CosmoComputationError as failure_message:
            sys.stderr.write(str(failure_message)+'\n')
            sys.stderr.flush()
            return data.boundary_loglike
        except CosmoSevereError as critical_message:
            raise io_mp.CosmologicalModuleError(
                "Something went wrong when calling CLASS" +
                str(critical_message))
        except KeyboardInterrupt:
            raise io_mp.CosmologicalModuleError(
                "You interrupted execution")

    # For each desired likelihood, compute its value against the theoretical
    # model
    loglike = 0
    # This flag holds the information whether a fiducial model was written. In
    # this case, the log likelihood returned will be '1j', meaning the
    # imaginary number i.
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
        if value == 1j:
            flag_wrote_fiducial += 1

    # Compute the derived parameters if relevant
    if data.get_mcmc_parameters(['derived']) != []:
        try:
            derived = cosmo.get_current_derived_parameters(
                data.get_mcmc_parameters(['derived']))
            for name, value in derived.iteritems():
                data.mcmc_parameters[name]['current'] = value
        except AttributeError:
            # This happens if the classy wrapper is still using the old
            # convention, expecting data as the input parameter
            cosmo.get_current_derived_parameters(data)
        except CosmoSevereError:
            raise io_mp.CosmologicalModuleError(
                "Could not write the current derived parameters")
    for elem in data.get_mcmc_parameters(['derived']):
        data.mcmc_parameters[elem]['current'] /= \
            data.mcmc_parameters[elem]['scale']

    # If fiducial files were created, inform the user, and exit
    if flag_wrote_fiducial > 0:
        if flag_wrote_fiducial == len(data.lkl):
            raise io_mp.FiducialModelWritten(
                "Fiducial file(s) was(were) created, please start a new chain")
        else:
            raise io_mp.FiducialModelWritten(
                "Some previously non-existing fiducial files were created, " +
                "but potentially not all of them. Please check now manually" +
                " on the headers, of the corresponding that all parameters " +
                "are coherent for your tested models")

    return loglike


def compute_fisher(data, cosmo, center, step_size):

    parameter_names = data.get_mcmc_parameters(['varying'])
    fisher_matrix = np.zeros(
        (len(parameter_names), len(parameter_names)), 'float64')
    # Initialise the gradient field
    gradient = np.zeros(len(parameter_names), 'float64')

    for k, elem_k in enumerate(parameter_names):
        kdiff = center[elem_k]*step_size
        if kdiff == 0.0:
            kdiff = step_size
        for h, elem_h in enumerate(parameter_names):
            hdiff = center[elem_h]*step_size
            if hdiff == 0.0:
                hdiff = step_size
            # Since the matrix is symmetric, we only compute the
            # elements of one half of it plus the diagonal.
            if k > h:
                continue
            if k != h:
                fisher_matrix[k][h] = compute_fisher_element(
                    data, cosmo, center,
                    (elem_k, kdiff),
                    (elem_h, hdiff))
                fisher_matrix[h][k] = fisher_matrix[k][h]
            else:
                fisher_matrix[k][k], gradient[k] = compute_fisher_element(
                    data, cosmo, center,
                    (elem_k, kdiff))

    return fisher_matrix, gradient


def compute_fisher_element(data, cosmo, center, one, two=None):

    # Unwrap
    name_1, diff_1 = one
    if two:
        name_2, diff_2 = two
        data.mcmc_parameters[name_1]['current'] = (
            center[name_1]+diff_1)
        data.mcmc_parameters[name_2]['current'] = (
            center[name_2]+diff_2)
        data.update_cosmo_arguments()
        loglike_1 = compute_lkl(cosmo, data)

        data.mcmc_parameters[name_2]['current'] -= 2*diff_2
        data.update_cosmo_arguments()
        loglike_2 = compute_lkl(cosmo, data)

        data.mcmc_parameters[name_1]['current'] -= 2*diff_1
        data.mcmc_parameters[name_2]['current'] += 2*diff_2
        data.update_cosmo_arguments()
        loglike_3 = compute_lkl(cosmo, data)

        data.mcmc_parameters[name_2]['current'] -= 2*diff_2
        data.update_cosmo_arguments()
        loglike_4 = compute_lkl(cosmo, data)

        fisher_off_diagonal = -(
            loglike_1-loglike_2-loglike_3+loglike_4)/(4.*diff_1*diff_2)

        return fisher_off_diagonal
    # It is otherwise a diagonal component
    else:
        data.mcmc_parameters[name_1]['current'] = center[name_1]
        data.update_cosmo_arguments()
        loglike_1 = compute_lkl(cosmo, data)

        data.mcmc_parameters[name_1]['current'] += diff_1
        data.update_cosmo_arguments()
        loglike_2 = compute_lkl(cosmo, data)

        data.mcmc_parameters[name_1]['current'] -= 2*diff_1
        data.update_cosmo_arguments()
        loglike_3 = compute_lkl(cosmo, data)

        fisher_diagonal = -(
            loglike_2-2.*loglike_1+loglike_3)/(diff_1**2)
        gradient = -(loglike_2-loglike_3)/(2.*diff_1)

        return fisher_diagonal, gradient
