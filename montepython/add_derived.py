"""
.. module:: add_derived
   :synopsis: Read a chain and add new derived parameters from an existing run

.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>
"""
try:
    from collections import OrderedDict as od
except ImportError:
    from ordereddict import OrderedDict as od
from copy import copy
from multiprocessing import Pool
import os

import io_mp
import sampler
from data import Data
from classy import CosmoComputationError


def run(cosmo, data, command_line):
    """
    Rewrite chains with more derived parameters

    Starting from an existing folder, with some chains, constraining a certain
    model, and having some derived parameters, the idea is to recompute the
    cosmological code to follow additional derived parameters.
    """
    starting_folder = command_line.Der_starting_folder
    if not starting_folder:
        raise io_mp.ConfigurationError(
            "When running importance sampling, you should specify a folder or"
            " a set of chains with the option '--IS-starting-folder'")
    chains = []
    # If starting_folder is of length 1, it means it is either a whole folder,
    # or just one chain. If it is a folder, we recover all chains within.
    if len(starting_folder) == 1:
        starting_folder = starting_folder[0]
        if os.path.isdir(starting_folder):
            for elem in os.listdir(starting_folder):
                if elem.find("__") != -1:
                    chains.append(elem)
    # Else, it is a list of chains, of which we recover folder name, and store
    # all of them in chains.
    else:
        chains = starting_folder
        starting_folder = os.path.sep.join(chains[0].split(os.path.sep)[:-1])
        chains = [elem.split(os.path.sep)[-1] for elem in chains]

    # Read the additional derived parameter, remove the needs for output=mPk
    # except if sigma8 is there.
    new_derived = recover_new_derived_parameters(
        data, command_line, starting_folder)
    if not new_derived:
        raise io_mp.ConfigurationError(
            "You asked to add new derived parameters to a folder, but "
            "your new parameter file does not contain new derived "
            "parameters. Please reconsider your intent.")
    data.cosmo_arguments = {}
    if 'sigma8' in new_derived:
        data.cosmo_arguments.update({'output': 'mPk'})

    # Preparing the arguments for reading the files
    pool = Pool()
    args = [(data, cosmo, command_line, starting_folder,
             elem, new_derived) for elem in chains]
    # Note the use of translate_chain_star, and not translate_chain, because of
    # the limitations of the `map` function (it only takes one argument). The
    # `_star` function simply unwraps the argument.
    print '\nStart extracting the chains:\n'
    pool.map(extend_chain_star, args)
    # Close the pool, and join everything (the join might not be needed)
    pool.close()
    pool.join()


def recover_new_derived_parameters(data, command_line, starting_folder):
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
    new_derived = [elem for elem in data.get_mcmc_parameters(['derived'])
                   if elem not in data2.get_mcmc_parameters(['derived'])]
    print
    print 'The new folder will follow the additional derived parameters:'
    print ' ->' + ', '.join(new_derived)

    return new_derived


def extend_chain(data, cosmo, command_line, starting_folder, chain_name,
                 new_derived):
    """
    Reading the input point, and computing the new derived values

    """
    input_path = os.path.join(starting_folder, chain_name)
    output_path = os.path.join(command_line.folder, chain_name)
    print ' -> reading ', input_path
    # Put in parameter_names all the varying parameters, plus the derived ones
    # that are not part of new_derived
    parameter_names = data.get_mcmc_parameters(['varying'])
    parameter_names.extend([
        elem for elem in data.get_mcmc_parameters(['derived'])
        if elem not in new_derived])
    with open(input_path, 'r') as input_chain:
        with open(output_path, 'w') as output_chain:
            for line in input_chain:
                params = line.split()
                # recover the likelihood of this point
                loglike = -float(params[1])
                N = int(params[0])
                # Assign all the recovered values to the data structure
                for index, param in enumerate(parameter_names):
                    data.mcmc_parameters[param]['current'] = \
                        float(params[2+index])
                # Compute the cosmology
                data.update_cosmo_arguments()
                if cosmo.state:
                    cosmo.struct_cleanup()
                cosmo.set(data.cosmo_arguments)
                try:
                    cosmo.compute(["lensing"])
                except CosmoComputationError:
                    pass
                # Recover all the derived parameters
                derived = cosmo.get_current_derived_parameters(
                    new_derived)
                for name, value in derived.iteritems():
                    data.mcmc_parameters[name]['current'] = value
                for elem in new_derived:
                    data.mcmc_parameters[elem]['current'] /= \
                        data.mcmc_parameters[elem]['scale']
                # Accept the point
                sampler.accept_step(data)
                io_mp.print_vector([output_chain], N, loglike, data)
    print output_path, 'written'


def extend_chain_star(args):
    """Trick function for multiprocessing"""
    return extend_chain(*args)
