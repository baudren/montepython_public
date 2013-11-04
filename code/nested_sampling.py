import pymultinest
import numpy as np
import os
import io_mp
from pprint import pprint

def prior(parameter_vector):
    return [data.mcmc_parameters[name]["prior"].map_from_unit_interval(value)
            for value, name
            in zip(parameter_vector, data.get_mcmc_parameters(['varying']))]
        
def from_NS_output_to_chains(command_line, basename):
    # First, take care of post_equal_weights (accepted points)
    accepted_chain = os.path.join(command_line.folder, 'chain_NS__accepted.txt')
    rejected_chain = os.path.join(command_line.folder, 'chain_NS__rejected.txt')

    # creating chain of accepted points (straightforward reshuffling of
    # columns)
    with open(basename+'post_equal_weights.dat', 'r') as Input:
        output = open(accepted_chain, 'w')
        array = np.loadtxt(Input)
        output_array = np.ones((np.shape(array)[0], np.shape(array)[1]+1))
        output_array[:, 1] = -array[:, -1]
        output_array[:, 2:] = array[:, :-1]
        pprint(array)
        pprint(output_array)
        np.savetxt(output, output_array, fmt='%i '+' '.join(['%.6e' for i in
            range(np.shape(array)[1])]))
        output.close()

    # Extracting log evidence
    with open(basename+'stats.dat') as Input:
        lines = [line for line in Input if 'Global Log-Evidence' in line]
        if len(lines) > 1:
            lines = [line for line in lines if 'Importance' in line]
        log_evidence = float(lines[0].split(':')[1].split('+/-')[0])
        print 'Evidence is ', log_evidence

    # Creating chain from rejected points, with some interpretation of the
    # weight associated to each point arXiv:0809.3437 sec 3
    with open(basename+'ev.dat', 'r') as Input:
        output = open(rejected_chain, 'w')
        array = np.loadtxt(Input)
        output_array = np.zeros((np.shape(array)[0], np.shape(array)[1]-1))
        output_array[:, 0] = np.exp(array[:, -3]+array[:, -2]-log_evidence)
        output_array[:, 0] *= np.sum(output_array[:, 0])*np.shape(array)[0]
        output_array[:, 1] = -array[:, -3]
        output_array[:, 2:] = array[:, :-3]
        pprint(array)
        pprint(output_array)
        np.savetxt(output, output_array, fmt=' '.join(['%.6e' for i in
            range(np.shape(output_array)[1])]))
        output.close()


def run(cosmo, data, command_line):
    
    varying_param_names = data.get_mcmc_parameters(['varying'])
    derived_param_names = data.get_mcmc_parameters(['derived'])

    # Check that all the priors are flat and that all the parameters are bound
    if not(all(data.mcmc_parameters[name]["prior"].prior_type == 'flat'
               for name in varying_param_names)):
        io_mp.message(
            "Nested Sampling with MultiNest is only possible with flat priors. \
             Sorry!",
            "error")
    if not(all(data.mcmc_parameters[name]["prior"].is_bound()
               for name in varying_param_names)):
        io_mp.message(
            "Nested Sampling with MultiNest is only possible for bound parameters. \
             Set reasonable bounds for them in the '.param' file.",
            "error")

    # Generate the prior function for MultiNest
    def prior(cube, ndim, nparams):
        for i, name in zip(range(ndim), varying_param_names):
            cube[i] = data.mcmc_parameters[name]["prior"]\
                      .map_from_unit_interval(cube[i])

    # Generate the Likelihood function for MultiNest
    from mcmc import compute_lkl
    def loglik(cube, ndim, nparams):
        # Updates values: cube --> data
        for i, name in zip(range(ndim), varying_param_names):
            data.mcmc_parameters[name]['current'] = cube[i]
        # Propagate the information towards the cosmo arguments
        data.update_cosmo_arguments()
        lkl = compute_lkl(cosmo, data)
        for i, name in enumerate(derived_param_names):
            cube[ndim+i] = data.mcmc_parameters[name]["current"]
        return lkl
    
    # If absent, create the subfolder NS
    NS_subfolder = os.path.join(command_line.folder, 'NS/')
    if not os.path.exists(NS_subfolder):
        os.makedirs(NS_subfolder)

    basename = os.path.join(NS_subfolder,
        command_line.folder.split(os.path.sep)[-2]+'-')
    # Launch MultiNest
    pymultinest.run(loglik, prior, n_dims=len(varying_param_names),
                    n_params=len(varying_param_names)+len(derived_param_names),
                    outputfiles_basename=basename,
                    evidence_tolerance=0.5,
                    n_live_points=100,
                    verbose=True)
    
    # Assuming this worked, translate the output ev.txt into the same format as
    # standard Monte Python chains for further analysis.
    from_NS_output_to_chains(command_line, basename)

