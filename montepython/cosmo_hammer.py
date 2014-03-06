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
LikelihoodModules of CosmoHammer.

"""
import numpy as np
import os
import logging

import io_mp
import sampler
from cosmoHammer.likelihood.chain.LikelihoodComputationChain import (
    LikelihoodComputationChain)
from cosmoHammer.sampler.CosmoHammerSampler import CosmoHammerSampler

# Cosmo Hammer subfolder and name separator
CH_subfolder = 'CH'
CH_separator = '-'
# Cosmo Hammer file names ending, after the defined 'base_name'
name_arguments = '.arguments'

# Cosmo Hammer option prefix
CH_prefix = 'CH_'
# User-defined arguments for the Hammer
CH_user_arguments = {
    'walkersRatio':
    {'metavar': 'Number of walkers',
     'type': int},
    'burninIterations':
    {'metavar': 'Number of burnin phase iterations',
     'type': int},
    'sampleIterations':
    {'metavar': 'Number of sample iterations',
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
    file_prefix = os.path.join(
        command_line.folder,
        os.path.join(CH_subfolder, chain_name))

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
            arg_file.write(' = '.join([str(arg), str(data.CH_arguments[arg])]))

    # Create the Sampler object
    sampler_hammer = CosmoHammerSampler(
        params=params,
        likelihoodComputationChain=chain,
        filePrefix=file_prefix,
        **data.CH_arguments)

    # create console handler and set level to debug (does not seem to appear)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    logging.getLogger().addHandler(console_handler)

    sampler_hammer.startSampling()
