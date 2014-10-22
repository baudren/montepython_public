"""
.. module:: JLA_simple
    :synopsis: Simplified JLA likelihood from Betoule et al. 2014

.. moduleauthor:: Benjamin Audren <benjamin.audren@gmail.com>

Copied from the original c++ code available at
`this address <http://supernovae.in2p3.fr/sdss_snls_jla/ReadMe.html>`_.

"""
import os
import numpy as np
import scipy.linalg as la
import montepython.io_mp as io_mp
from montepython.likelihood_class import Likelihood_sn


class JLA_simple(Likelihood_sn):

    def __init__(self, path, data, command_line):

        # This reads the configuration file as well
        try:
            Likelihood_sn.__init__(self, path, data, command_line)
        except IOError:
            raise io_mp.LikelihoodError(
                "The JLA data files were not found. Please download the "
                "following link "
                "http://supernovae.in2p3.fr/sdss_snls_jla/jla_likelihood_v4.tgz"
                ", extract it, and copy all files present in "
                "`jla_likelihood_v4/data` to `your_montepython/data/JLA`")

        # read the only matrix
        self.C00 = self.read_matrix(self.mu_covmat_file)

        # Read the simplified light-curve self.data_file
        self.light_curve_params = self.read_light_curve_parameters()

        # The covariance matrix can be already inverted, once and for all
        # (cholesky)
        self.C00 = la.cholesky(self.C00, lower=True, overwrite_a=True)

    def loglkl(self, cosmo, data):
        # Recover the distance moduli from CLASS (a size N vector of double
        # containing the predicted distance modulus for each SN in the JLA
        # sample, given the redshift of the supernova.)
        sn = self.light_curve_params
        size = sn.z.size

        moduli = np.empty((size, ))
        for index, z in enumerate(sn.z):
            moduli[index] = cosmo.luminosity_distance(z)
        moduli = 5 * np.log10(moduli) + 25

        # Convenience variables: store the nuisance parameters in short named
        # variables
        M = (data.mcmc_parameters['M']['current'] *
             data.mcmc_parameters['M']['scale'])

        residuals = sn.mu - (M+moduli)

        # Whiten the residuals
        residuals = la.solve_triangular(self.C00, residuals, lower=True,
                                        check_finite=False)

        chi2 = (residuals**2).sum()

        return -0.5 * chi2
