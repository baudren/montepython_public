from montepyhon.likelihood_class import likelihood
from numpy import matrix, dot


class test_gaussian(likelihood):

    def loglkl(self, cosmo, data):
        H0, ob, oc  = (
            data.mcmc_parameters[p]['current']*data.mcmc_parameters[p]['scale']
            for p in ['H0', 'omega_b', 'omega_cdm'])
        # Modes
        loglkl = 0
        for centre, covmat in zip([self.centre1, self.centre2, self.centre3],
                                  [self.covmat1, self.covmat2, self.covmat3]):
            diffvec = matrix([x-mu for x, mu in zip([H0, ob, oc], centre)])
            Hessian = matrix(covmat).I
            loglkl += -0.5 * (dot(diffvec, dot(Hessian, diffvec.T)))
        return loglkl

