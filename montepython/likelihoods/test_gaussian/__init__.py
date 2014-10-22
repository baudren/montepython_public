from montepython.likelihood_class import Likelihood
from numpy import matrix, dot, exp, log


class test_gaussian(Likelihood):

    def loglkl(self, cosmo, data):
        H0, ob, oc = (
            data.mcmc_parameters[p]['current']*data.mcmc_parameters[p]['scale']
            for p in ['H0', 'omega_b', 'omega_cdm'])
        # Modes
        lkl = 0
        for centre, covmat in zip([self.centre1, self.centre2, self.centre3],
                                  [self.covmat1, self.covmat2, self.covmat3]):
            diffvec = matrix([x-mu for x, mu in zip([H0, ob, oc], centre)])
            minusHessian = matrix(covmat).I
            lkl += exp(-0.5 * (dot(diffvec, dot(minusHessian, diffvec.T))))
        return log(lkl)
