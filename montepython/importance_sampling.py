"""
.. module:: importance_sampling
    :synopsis: Perform an Importance Sampling from an existing folder

.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>
"""
try:
    from collections import OrderedDict as od
except ImportError:
    from ordereddict import OrderedDict as od
from copy import copy
from multiprocessing import Pool
import os
import warnings
import math

import io_mp
import sampler
from data import Data


def run(cosmo, data, command_line):
    """
    Performing the Importance Sampling

    The idea is to start from an existing run, constraining a certain model I,
    given a set of experiments. The new run will constrain the same model I,
    but adding one or several new experiments. In the case where it is expected
    that the final posterior distribution should not differ too greatly between
    the two parameter extractions, then using Importance Sampling can speed up
    significantly the second one.

    Instead of properly sampling randomly the parameter space, it instead reads
    the chains from the previous run, recompute the cosmology at this point,
    then adds the log-likelihood contributed by the new experiments to the
    previous ones. As an input of the method, with the flag
    `--IS-starting-folder`, you can thus specify either a folder containing a
    Monte Python run, or a set of chains that you want to be converted.

    The code will automatically compute the minimum amount of things. For
    instance, if the first run had all the Planck likelihoods, and the second,
    all the Planck likelihoods plus a prior on :math:`H_0`, it would be absurd
    to recompute also the cosmological perturbations: the only needed quantity
    is a background quantity.

    The new chains will hence store the same points in parameter space, but
    with a different value of the likelihood, and also of the multiplicity -
    that will become non-integer. Indeed, the multiplicity is also a probe of
    the posterior, and this new, higher likelihood should have had a higher
    multiplicity.
    """
    # Check that the command_line "--IS-starting-folder" points to an existing
    # Monte Python folder run, or a subset of files, and store in any case all
    # the chains to analyze in the chains.
    starting_folder = command_line.IS_starting_folder
    if not starting_folder:
        raise io_mp.ConfigurationError(
            "When running importance sampling, you should specify a folder or"
            " a set of chains with the option '--IS-starting-folder'")
    chains = []
    if len(starting_folder) == 1:
        starting_folder = starting_folder[0]
        if os.path.isdir(starting_folder):
            for elem in os.listdir(starting_folder):
                if elem.find("__") != -1:
                    chains.append(elem)
    else:
        chains = starting_folder
        starting_folder = os.path.sep.join(chains[0].split(os.path.sep)[:-1])
        chains = [elem.split(os.path.sep)[-1] for elem in chains]

    # Recovering only the extra likelihoods
    new_experiments = recover_new_experiments(
        data, command_line, starting_folder)
    if not new_experiments:
        raise io_mp.ConfigurationError(
            "To use Importance Sampling, your new folder should have at least"
            " one more likelihood than the previous one. Please ensure that "
            "it is the case.")

    # resetting the needed cosmo arguments, and deleting the dictionary of
    # likelihoods, only if new_experiments is smaller than the old ones.
    ignore_likelihood = True
    if set(new_experiments) < set(data.experiments):
        data.cosmo_arguments = {}
        data.initialise_likelihoods(new_experiments)
        ignore_likelihood = False
    elif set(new_experiments) == set(data.experiments):
        warnings.warn(
            "All likelihoods were found to be different than ones fron the "
            "starting folder. The previous value of the likelihood will be "
            "discarded.")
    else:
        io_mp.ConfigurationError(
            "You probably tried to run Importance Sampling with less "
            "experiments than what you started with. This is not supported.")

    # Multiprocessing part, to analyze all the chains in parralel. When not
    # specifying any 'processes' keyword argument to the Pool call, the system
    # uses as many as possible.
    pool = Pool()
    args = [(data, cosmo, command_line, starting_folder,
             elem, ignore_likelihood) for elem in chains]
    # Note the use of translate_chain_star, and not translate_chain, because of
    # the limitations of the `map` function (it only takes one argument). The
    # `_star` function simply unwraps the argument.
    print '\nStart extracting the chains:\n'
    pool.map(translate_chain_star, args)
    # Close the pool, and join everything (the join might not be needed)
    pool.close()
    pool.join()


def recover_new_experiments(data, command_line, starting_folder):
    """
    Given the input, extract the additional likelihoods

    """
    # Initialize the companion data structure, on a modified command line
    modified_command_line = copy(command_line)
    modified_command_line.folder = starting_folder
    modified_command_line.param = os.path.join(starting_folder, 'log.param')

    print 'Reading the starting folder'
    print '---------------------------'
    data2 = Data(modified_command_line, data.path)
    print '---------------------------'
    print 'Finished loading existing data'
    print
    print 'The likelihood will be computed for:'
    new_experiments = [
        elem for elem in data.experiments if elem not in data2.experiments]
    print ' ->',
    print ', '.join(new_experiments)

    return new_experiments


def translate_chain(data, cosmo, command_line,
                    starting_folder, chain_name, ignore_likelihood=False):
    """Translate the input to the output

    .. note::

        if the keyword argument `ignore_likelihood` is set to true, the
        previous value of the likelihood is discarded.

    """

    input_path = os.path.join(starting_folder, chain_name)
    output_path = os.path.join(command_line.folder, chain_name)
    print ' -> reading ', input_path
    parameter_names = data.get_mcmc_parameters(['varying'])
    with open(input_path, 'r') as input_chain:
        with open(output_path, 'w') as output_chain:
            for line in input_chain:
                params = line.split()
                # recover the likelihood of this point
                if not ignore_likelihood:
                    loglike = -float(params[1])
                else:
                    loglike = 0
                N = float(params[0])
                # Assign all the recovered values to the data structure
                for index, param in enumerate(parameter_names):
                    data.mcmc_parameters[param]['current'] = \
                        float(params[2+index])
                data.update_cosmo_arguments()

                newloglike = sampler.compute_lkl(cosmo, data)

                weight = math.exp(newloglike)
                newloglike += loglike
                # Accept the point
                sampler.accept_step(data)
                io_mp.print_vector([output_chain], N*weight, newloglike, data)
    print output_path, 'written'


def translate_chain_star(args):
    """Trick function for multiprocessing"""
    return translate_chain(*args)
