"""
.. module:: add_derived
   :synopsis: Read a chain and add new derived parameters from an existing run

.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>
"""
try:
    from collections import OrderedDict as od
except ImportError:
    from ordereddict import OrderedDict as od
from multiprocessing import Pool
import os

import io_mp
import sampler
from data import Data
from data import Parameter
from classy import CosmoComputationError


def run(cosmo, data, command_line):
    """
    Rewrite chains with more derived parameters

    Starting from an existing folder, with some chains, constraining a certain
    model, and having some derived parameters, the idea is to recompute the
    cosmological code to follow additional derived parameters.
    """
    target_folder = command_line.Der_target_folder
    # If it does not exist, create it
    if not os.path.isdir(target_folder):
        os.makedirs(target_folder)

    starting_folder = command_line.folder
    # Recover all chains in the starting folder
    chains = []
    #  If it exists, we recover all chains within.
    if os.path.isdir(starting_folder):
        for elem in os.listdir(starting_folder):
            if elem.find("__") != -1:
                chains.append(elem)

    # Read the additional derived parameter, remove the needs for output=mPk
    # except if sigma8 is there.
    new_derived = command_line.derived_parameters
    if not new_derived:
        raise io_mp.ConfigurationError(
            "You asked to add derived parameters, but did not specify a list "
            "of new ones to consider. Please use the flag `--Der-param-list`.")
    # Add them to the mcmc_parameters dict
    for param in new_derived:
        data.mcmc_parameters[param] = Parameter(
            [0, None, None, 0, 1, 'derived'], param)
    # Reset the cosmo_arguments dict output entry, and adapt it in case a
    # derived parameter requires a particular CLASS behaviour.
    data.cosmo_arguments.update({'output': ''})
    for key in ['lensing', 'l_max_scalars']:
        if key in data.cosmo_arguments.keys():
            data.cosmo_arguments.pop(key)
    if 'sigma8' in new_derived:
        data.cosmo_arguments.update({'output': 'mPk'})

    # Copy the log.param over from the starting folder, and add new lines
    # concerning the new derived parameters, for analysis.
    copy_log_file(starting_folder, target_folder, new_derived)

    # Preparing the arguments for reading the files
    pool = Pool()
    args = [(data, cosmo, command_line, target_folder,
             elem, new_derived) for elem in chains]
    # Note the use of translate_chain_star, and not translate_chain, because of
    # the limitations of the `map` function (it only takes one argument). The
    # `_star` function simply unwraps the argument.
    print '\nStart extracting the chains:\n'
    pool.map(extend_chain_star, args)
    # Close the pool, and join everything (the join might not be needed)
    pool.close()
    pool.join()


def extend_chain(data, cosmo, command_line, target_folder, chain_name,
                 new_derived):
    """
    Reading the input point, and computing the new derived values

    """
    input_path = os.path.join(command_line.folder, chain_name)
    output_path = os.path.join(target_folder, chain_name)
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
                if line[0] == '#':
                    output_chain.write(line)
                    continue
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
                    data.get_mcmc_parameters(['derived']))
                for name, value in derived.iteritems():
                    data.mcmc_parameters[name]['current'] = value
                for name in derived.iterkeys():
                    data.mcmc_parameters[elem]['current'] /= \
                        data.mcmc_parameters[elem]['scale']
                # Accept the point
                sampler.accept_step(data)
                io_mp.print_vector([output_chain], N, loglike, data)
    print output_path, 'written'


def extend_chain_star(args):
    """Trick function for multiprocessing"""
    return extend_chain(*args)


def copy_log_file(starting_folder, target_folder, new_derived):
    """
    Copy and extend the log.param from the starting to the target folder

    while adding new derived parameters
    """
    in_path = os.path.join(starting_folder, 'log.param')
    out_path = os.path.join(target_folder, 'log.param')
    with open(in_path, 'r') as input_log:
        with open(out_path, 'w') as output_log:
            # Read the whole file
            text = input_log.readlines()
            # Find the index of the last 'derived' line
            last_derived_line = 0
            for index, line in enumerate(text):
                if line.find('derived') != -1:
                    last_derived_line = index
            # Write everything up to this point to the output
            for index in range(last_derived_line+1):
                output_log.write(text[index])
            # Append lines for the new derived parameters
            for name in new_derived:
                output_log.write(
                    "data.parameters['%s'] = [0,-1,-1,0,1,'derived']\n" % (
                        name))
            # Write the rest
            for index in range(last_derived_line+1, len(text)):
                output_log.write(text[index])
