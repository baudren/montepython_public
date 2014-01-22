"""
.. module:: nested_sampling
    :synopsis: Interface the MultiNest program with Monte Python

This implementation relies heavily on the existing Python wrapper for
MultiNest, called PyMultinest, written by Johannes Buchner, and available `at
this address <https://github.com/JohannesBuchner/PyMultiNest>`_ .

The main routine, :func:`run`, truly interfaces the two codes. It takes for
input the cosmological module, data and command line. It then defines
internally two functions, :func:`prior() <nested_sampling.prior>` and
:func:`loglike` that will serve as input for the run function of PyMultiNest.

.. moduleauthor:: Jesus Torrado Cacho <torradocacho@lorentz.leidenuniv.nl>
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>
"""
import pymultinest
import numpy as np
import os
import io_mp
import sampler


def from_ns_output_to_chains(command_line, basename):
    """
    Translate the output of MultiNest into readable output for Monte Python

    This routine will be called after the MultiNest run has been successfully
    completed.

    """
    # First, take care of post_equal_weights (accepted points)
    accepted_chain = os.path.join(command_line.folder,
                                  'chain_NS__accepted.txt')
    rejected_chain = os.path.join(command_line.folder,
                                  'chain_NS__rejected.txt')

    # creating chain of accepted points (straightforward reshuffling of
    # columns)
    with open(basename+'post_equal_weights.dat', 'r') as input_file:
        output_file = open(accepted_chain, 'w')
        array = np.loadtxt(input_file)
        output_array = np.ones((np.shape(array)[0], np.shape(array)[1]+1))
        output_array[:, 1] = -array[:, -1]
        output_array[:, 2:] = array[:, :-1]
        np.savetxt(
            output_file, output_array,
            fmt='%i '+' '.join(['%.6e' for _ in
                               range(np.shape(array)[1])]))
        output_file.close()

    # Extracting log evidence
    with open(basename+'stats.dat') as input_file:
        lines = [line for line in input_file if 'Global Log-Evidence' in line]
        if len(lines) > 1:
            lines = [line for line in lines if 'Importance' in line]
        log_evidence = float(lines[0].split(':')[1].split('+/-')[0])
        print 'Evidence is ', log_evidence

    # Creating chain from rejected points, with some interpretation of the
    # weight associated to each point arXiv:0809.3437 sec 3
    with open(basename+'ev.dat', 'r') as input_file:
        output = open(rejected_chain, 'w')
        array = np.loadtxt(input_file)
        output_array = np.zeros((np.shape(array)[0], np.shape(array)[1]-1))
        output_array[:, 0] = np.exp(array[:, -3]+array[:, -2]-log_evidence)
        output_array[:, 0] *= np.sum(output_array[:, 0])*np.shape(array)[0]
        output_array[:, 1] = -array[:, -3]
        output_array[:, 2:] = array[:, :-3]
        np.savetxt(
            output, output_array,
            fmt=' '.join(['%.6e' for _ in
                         range(np.shape(output_array)[1])]))
        output.close()


def run(cosmo, data, command_line):
    """
    Main call to prepare the information for the MultiNest run.

    Note the unusual set-up here, with the two following functions, `prior` and
    `loglike` having their docstrings written in the encompassing function.
    This trick was necessary as MultiNest required these two functions to be
    defined with a given number of parameters, so we could not add `data`. By
    defining them inside the run function, this problem was by-passed.

    .. function:: prior

        Generate the prior function for MultiNest

        It should transform the input unit cube into the parameter cube. This
        function actually wraps the method :func:`map_from_unit_interval()
        <prior.Prior.map_from_unit_interval>` of the class :class:`Prior
        <prior.Prior>`.

        :Parameters:
            **cube** (`array`) - Contains the current point in unit parameter
                space that has been selected within the MultiNest part.
            **ndim** (`int`) - Number of varying parameters
            **nparams** (`int`) - Total number of parameters, including the
                derived ones (not used, so hidden in `*args`)


    .. function:: loglike

        Generate the Likelihood function for MultiNest

        :Parameters:
            **cube** (`array`) - Contains the current point in the correct
                parameter space after transformation from :func:`prior`.
            **ndim** (`int`) - Number of varying parameters
            **nparams** (`int`) - Total number of parameters, including the
                derived ones (not used, so hidden in `*args`)

    """
    # Convenience variables
    varying_param_names = data.get_mcmc_parameters(['varying'])
    derived_param_names = data.get_mcmc_parameters(['derived'])

    # Check that all the priors are flat and that all the parameters are bound
    if not(all(data.mcmc_parameters[name]["prior"].prior_type == 'flat'
               for name in varying_param_names)):
        raise io_mp.ConfigurationError(
            "Nested Sampling with MultiNest is only possible with flat " +
            "priors. Sorry!")
    if not(all(data.mcmc_parameters[name]["prior"].is_bound()
               for name in varying_param_names)):
        raise io_mp.ConfigurationError(
            "Nested Sampling with MultiNest is only possible for bound " +
            "parameters. Set reasonable bounds for them in the '.param'" +
            "file.")

    def prior(cube, ndim, *args):
        """
        Please see the encompassing function docstring

        """
        for i, name in zip(range(ndim), varying_param_names):
            cube[i] = data.mcmc_parameters[name]["prior"]\
                .map_from_unit_interval(cube[i])

    def loglike(cube, ndim, *args):
        """
        Please see the encompassing function docstring

        """
        # Updates values: cube --> data
        for i, name in zip(range(ndim), varying_param_names):
            data.mcmc_parameters[name]['current'] = cube[i]
        # Propagate the information towards the cosmo arguments
        data.update_cosmo_arguments()
        lkl = sampler.compute_lkl(cosmo, data)
        for i, name in enumerate(derived_param_names):
            cube[ndim+i] = data.mcmc_parameters[name]["current"]
        return lkl

    # If absent, create the sub-folder NS
    ns_subfolder = os.path.join(command_line.folder, 'NS/')
    if not os.path.exists(ns_subfolder):
        os.makedirs(ns_subfolder)

    basename = os.path.join(
        ns_subfolder,
        command_line.folder.split(os.path.sep)[-2]+'-')

    # Prepare arguments for PyMultiNest
    # -- Automatic parameters
    data.ns_parameters['n_dims'] = len(varying_param_names)
    data.ns_parameters['n_params'] = (len(varying_param_names) +
                                      len(derived_param_names))
    data.ns_parameters['verbose'] = True
    data.ns_parameters['outputfiles_basename'] = basename
    # -- User-defined parameters
    parameters = ['n_live_points', 'sampling_efficiency', 'evidence_tolerance',
                  'importance_nested_sampling', 'const_efficiency_mode',
                  'log_zero', 'max_iter', 'seed', 'n_iter_before_update']
    prefix = 'NS_option_'
    for param in parameters:
        value = getattr(command_line, prefix+param)
        if value != -1:
            data.ns_parameters[param] = value
        # else: don't define them -> use PyMultiNest default value

    # Launch MultiNest, and recover the output code
    output = pymultinest.run(loglike, prior, **data.ns_parameters)

    # Assuming this worked, i.e. if output is `None`, translate the output
    # ev.txt into the same format as standard Monte Python chains for further
    # analysis.
    if output is None:
        from_ns_output_to_chains(command_line, basename)
