import os
from montepython.likelihood_class import Likelihood_prior


class BK14priors(Likelihood_prior):

    # initialisation of the class is done within the parent Likelihood_prior. For
    # this case, it does not differ, actually, from the __init__ method in
    # Likelihood class.
    def loglkl(self, cosmo, data):
        BBbetadust = data.mcmc_parameters['BBbetadust']['current']*data.mcmc_parameters['BBbetadust']['scale']
        BBbetasync = data.mcmc_parameters['BBbetasync']['current']*data.mcmc_parameters['BBbetasync']['scale']
        loglkl = -0.5 * (BBbetadust - self.mean_BBbetadust) ** 2 / (self.sigma_BBbetadust ** 2) -0.5 * (BBbetasync - self.mean_BBbetasync) ** 2 / (self.sigma_BBbetasync ** 2) 
        return loglkl
