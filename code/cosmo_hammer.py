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

from cosmoHammer.likelihood.chain.LikelihoodComputationChain import (
    LikelihoodComputationChain)
from cosmoHammer.sampler.CosmoHammerSampler import CosmoHammerSampler


def run(cosmo, data, command_line):
    """
    Sample with the CosmoHammer

    """
    # Store the parameters inside the format expected by CosmoHammer
    parameter_names = data.get_mcmc_parameters(["varying"])

    params = []
    for parameter in parameter_names:
        params.append(data.mcmc_parameters[parameter]['initial'])
    params = np.array(params)

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
    file_prefix = os.path.join(command_line.folder, 'cosmoHammer')

    # Create the Sampler object
    sampler = CosmoHammerSampler(
        params=params,
        likelihoodComputationChain=chain,
        filePrefix=file_prefix,
        walkersRatio=50,
        burninIterations=250,
        sampleIterations=250)

    # create console handler and set level to debug (does not seem to appear)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)
    logging.getLogger().addHandler(console_handler)

    sampler.startSampling()
