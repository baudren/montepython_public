"""
.. module:: sampler
    :synopsis: Generic sampler
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>

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

import io_mp
from classy import CosmoSevereError, CosmoComputationError


def run(cosmo, data, command_line):
    """
    First rudimentary implementation

    The :mod:`mcmc` module is used as previously, except the call to
    :func:`mcmc.chain`, or :func:`nested_sampling.run` is now within
    this function, instead of from within :mod:`MontePython`.

    In the long term, this function should contain any potential hybrid scheme,
    and any chain communication (defining different roles, etc)

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
                matrix = np.zeros((len(covnames), len(covnames)), 'float64')
                rot = np.zeros((len(covnames), len(covnames)))
            else:
                line = line.split()
                for j in range(len(line)):
                    matrix[i][j] = np.array(line[j], 'float64')
                i += 1

        # First print out
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
        for elem in parameter_names:
            matrix[i][i] = np.array(
                data.mcmc_parameters[elem]['initial'][3], 'float64')**2
            i += 1

    # Final print out, the actually used covariance matrix
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
            cosmo.get_current_derived_parameters(data)
        except NameError:
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
