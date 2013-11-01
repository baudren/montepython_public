import pymultinest
import io_mp

def prior(parameter_vector):
    return [data.mcmc_parameters[name]["prior"].map_from_unit_interval(value)
            for value, name
            in zip(parameter_vector, data.get_mcmc_parameters(['varying']))]
        

def run(cosmo, data, command_line):
    
    varying_param_names = data.get_mcmc_parameters(['varying'])

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
        return compute_lkl(cosmo, data)
            
    # Launch MultiNest
    pymultinest.run(loglik, prior, n_dims=len(varying_param_names),
                    outputfiles_basename=command_line.folder,
                    evidence_tolerance=0.5,
                    n_live_points=800,
                    verbose=True)

