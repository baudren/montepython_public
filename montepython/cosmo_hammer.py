"""
.. module:: cosmo_hammer
    :synopsis: Invoke the Cosmo Hammer
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>

This module interface Monte Python with the CosmoHammer, available
`publicly <http://www.astro.ethz.ch/refregier/research/Software/cosmohammer>`_
and developped by Joel Akeret and Sebastian Seehars, who helped me a lot
creating this interface.

The link between the two codes is that some functions have been added to the
classes defined in Monte Python, so that they can become CoreModules and
LikelihoodModules of CosmoHammer. Mostly, all the important classes in Monte
Python (:class:`Data` <data.Data>, Class and :class:`Likelihood`
<likelihood_class.Likelihood>)

"""
import numpy as np
import os
import warnings
import logging

import io_mp
import sampler
from cosmoHammer.likelihood.chain.LikelihoodComputationChain import (
    LikelihoodComputationChain)
from cosmoHammer.sampler.CosmoHammerSampler import CosmoHammerSampler
from cosmoHammer.util.SampleFileUtil import SampleFileUtil

# Cosmo Hammer subfolder and name separator
CH_subfolder = 'CH'
CH_separator = '-'
# Cosmo Hammer file names ending, after the defined 'base_name'
name_arguments = '.arguments'
name_chain = 'chain_CH__sampling.txt'

# Cosmo Hammer option prefix
CH_prefix = 'CH_'
# User-defined arguments for the Hammer
CH_user_arguments = {
    'walkersRatio':
    {'help': 'Number of walkers',
     'type': int},
    'burninIterations':
    {'help': 'Number of burnin phase iterations',
     'type': int},
    'sampleIterations':
    {'help': 'Number of sample iterations',
     'type': int}}


def run(cosmo, data, command_line):
    """
    Sample with the CosmoHammer

    """
    # Store the parameters inside the format expected by CosmoHammer
    # TODO: about the derived params?
    parameter_names = data.get_mcmc_parameters(["varying"])

    # Ensure that their prior is bound and flat
    is_flat, is_bound = sampler.check_flat_bound_priors(
        data.mcmc_parameters, parameter_names)
    if not is_flat:
        raise io_mp.ConfigurationError(
            'The Cosmo Hammer is only available with flat ' +
            'priors. Sorry!')
    if not is_bound:
        raise io_mp.ConfigurationError(
            'The Cosmo Hammer is only available for bound ' +
            'parameters. Set reasonable bounds for them in the ".param"' +
            'file.')

    params = []
    for parameter in parameter_names:
        params.append(data.mcmc_parameters[parameter]['initial'])
    params = np.array(params)

    # If absent, create the sub-folder CH
    CH_folder = os.path.join(command_line.folder, CH_subfolder)
    if not os.path.exists(CH_folder):
        os.makedirs(CH_folder)

    # Initialize a chain object (Beware, the concept is quite different than
    # the chain of the module :mod:`mcmc`)
    chain = LikelihoodComputationChain(
        min=params[:, 1],
        max=params[:, 2])

    # Add data and cosmo as two core modules. Note that the order is important
    # here, since data must be called before cosmo.
    chain.addCoreModule(data)
    chain.addCoreModule(cosmo)

    # Add each likelihood class as a LikelihoodModule
    for likelihood in data.lkl.itervalues():
        chain.addLikelihoodModule(likelihood)

    # Define the file prefix
    chain_name = [a for a in command_line.folder.split(os.path.sep) if a][-1]
    file_prefix = os.path.join(command_line.folder, CH_subfolder, chain_name)

    # Recover the User options
    data.CH_arguments = {}
    for arg in CH_user_arguments:
        value = getattr(command_line, CH_prefix+arg)
        if value != -1:
            data.CH_arguments[arg] = value
        # else, do not define them, and leave the default Cosmo Hammer ones.

    # Write the CosmoHammer arguments
    with open(file_prefix+name_arguments, 'w') as arg_file:
        for arg in data.CH_arguments:
            arg_file.write(
                ' = '.join([str(arg), str(data.CH_arguments[arg])]) + '\n')

    # Create an extension to the SampleFileUtil from cosmoHammer
    derived_util = DerivedUtil(file_prefix)

    try:
        num_threads = int(os.environ['OMP_NUM_THREADS'])
    except KeyError:
        warnings.warn(
            "The environment variable OMP_NUM_THREADS is not set. "
            "To run the Cosmo Hammer meaningfully, you should better "
            "set it to something! Defaulting to 1 for now.")
        num_threads = 1

    # Create the Sampler object
    sampler_hammer = CosmoHammerSampler(
        params=params,
        likelihoodComputationChain=chain,
        filePrefix=file_prefix,
        walkersRatio=50,
        burninIterations=10,
        sampleIterations=30,
        storageUtil=derived_util,
        threadCount=num_threads,
        **data.CH_arguments)

    # create console handler and set level to debug (does not seem to appear)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    logging.getLogger().addHandler(console_handler)

    sampler_hammer.startSampling()

    #print sampler_hammer
    #for elem in dir(sampler_hammer):
        #print elem
    #_, data_from_context = chain()
    #print data_from_context


def from_CH_output_to_chains(folder):
    """
    Translate the output of the Cosmo Hammer into Monte Python chains

    This function will be called by the module :mod:`analyze`.
    """

    chain_name = [a for a in folder.split(os.path.sep) if a][-2]
    base_name = os.path.join(folder, chain_name)
    # Recover the points in parameter space
    with open(base_name+'.out', 'r') as param_values_file:
        chains = np.loadtxt(param_values_file)
    # Recover the associated likelihood (it is the -log likelihood)
    with open(base_name+'prob.out', 'r') as lkl_values_file:
        lkl = np.loadtxt(lkl_values_file)

    # Glue them together, with an additional column of ones, for the
    # multiplicity. This does not mean that the acceptance rate is one, but
    # points where the code stayed are duplicated in the file.

    ## First, reshape the lkl array
    lkl = np.array([[elem] for elem in lkl])

    ## Create the array of ones
    ones = np.array([[1] for _ in range(len(lkl))])

    ## Concatenate everything and save to file
    final = np.concatenate((ones, lkl, chains), axis=1)
    output_folder = os.path.join(folder, '..')
    output_chain_path = os.path.join(output_folder, name_chain)
    np.savetxt(output_chain_path, final)


class DerivedUtil(SampleFileUtil):
    """
    Extends the writing class from CosmoHammer to include derived parameters.
    """
    def persistValues(self, posFile, probFile, pos, prob, data):
        """
        Writes the walker positions and the likelihood to the disk
        """
        # extend the pos array to also contain the value of the derived
        # parameters
        derived = np.array(
            [[a for a in elem.itervalues()] for elem in data])
        final = np.concatenate((pos, derived), axis=1)

        posFile.write("\n".join(
            [" ".join([str(q) for q in p]) for p in final]))
        posFile.write("\n")
        posFile.flush()

        probFile.write("\n".join([str(p) for p in prob]))
        probFile.write("\n")
        probFile.flush()
