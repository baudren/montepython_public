from montepython.likelihood_class import likelihood


class test_nuisance2(likelihood):

    def loglkl(self, cosmo, data):

        h = cosmo.h()
        amplitude = (data.mcmc_parameters['amplitude']['current'] *
                     data.mcmc_parameters['amplitude']['scale'])
        other = (data.mcmc_parameters['other']['current'] *
                 data.mcmc_parameters['other']['scale'])

        loglkl = -0.5 * (other*h - amplitude*self.h) ** 2 / (self.sigma**2)
        return loglkl
