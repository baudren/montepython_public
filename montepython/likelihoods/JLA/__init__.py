"""
.. module:: JLA
    :synopsis: JLA likelihood from Betoule et al. 2014

.. moduleauthor:: Benjamin Audren <benjamin.audren@gmail.com>

Copied from the original c++ code available at
`this address <http://supernovae.in2p3.fr/sdss_snls_jla/ReadMe.html>`_.
The only modification was to keep the names of the covariance matrices to their
original names instead of going for C00, etc... Here is the conversion from
their original file names to theirs. Note that the convention is used for the
computations afterwards.

.. code::

    C00 = mag_covmat_file
    C11 = stretch_covmat_file
    C22 = colour_covmat_file
    C01 = mag_stretch_covmat_file
    C02 = mag_colour_covmat_file
    C12 = stretch_colour_covmat_file

.. note::

    Since there are a lot of file manipulation involved, the "pandas" library
    has to be installed -- it is an 8-fold improvement in speed over numpy, and
    a 2-fold improvement over a fast Python implementation. The "numexpr"
    library is also needed for doing the fast array manipulations, done with
    blas daxpy function in the original c++ code. Both can be installed with
    pip (Python package manager) easily.

"""
import numpy as np
import scipy.linalg as la
import montepython.io_mp as io_mp
try:
    import numexpr as ne
except ImportError:
    raise io_mp.MissingLibraryError(
        "This likelihood has intensive array manipulations. You "
        "have to install the numexpr Python package. Please type:\n"
        "(sudo) pip install numexpr --user")
from montepython.likelihood_class import Likelihood_sn


class JLA(Likelihood_sn):

    def __init__(self, path, data, command_line):

        # Unusual construction, since the data files are not distributed
        # alongside JLA (size problems)
        try:
            Likelihood_sn.__init__(self, path, data, command_line)
        except IOError:
            raise io_mp.LikelihoodError(
                "The JLA data files were not found. Please download the "
                "following link "
                "http://supernovae.in2p3.fr/sdss_snls_jla/jla_likelihood_v4.tgz"
                ", extract it, and copy all files present in "
                "`jla_likelihood_v4/data` to `your_montepython/data/JLA`")

        # Load matrices from text files, whose names were read in the
        # configuration file
        self.C00 = self.read_matrix(self.mag_covmat_file)
        self.C11 = self.read_matrix(self.stretch_covmat_file)
        self.C22 = self.read_matrix(self.colour_covmat_file)
        self.C01 = self.read_matrix(self.mag_stretch_covmat_file)
        self.C02 = self.read_matrix(self.mag_colour_covmat_file)
        self.C12 = self.read_matrix(self.stretch_colour_covmat_file)

        # Reading light-curve parameters from self.data_file (jla_lcparams.txt)
        self.light_curve_params = self.read_light_curve_parameters()

    def loglkl(self, cosmo, data):
        """
        Compute negative log-likelihood (eq.15 Betoule et al. 2014)

        """
        # Recover the distance moduli from CLASS (a size N vector of double
        # containing the predicted distance modulus for each SN in the JLA
        # sample, given the redshift of the supernova.)
        redshifts = self.light_curve_params.zcmb
        size = redshifts.size

        moduli = np.empty((size, ))
        for index, row in self.light_curve_params.iterrows():
            z_cmb = row['zcmb']
            moduli[index] = cosmo.luminosity_distance(z_cmb)
        moduli = 5 * np.log10(moduli) + 25

        # Convenience variables: store the nuisance parameters in short named
        # variables
        alpha = (data.mcmc_parameters['alpha']['current'] *
                 data.mcmc_parameters['alpha']['scale'])
        beta = (data.mcmc_parameters['beta']['current'] *
                data.mcmc_parameters['beta']['scale'])
        M = (data.mcmc_parameters['M']['current'] *
             data.mcmc_parameters['M']['scale'])
        Delta_M = (data.mcmc_parameters['Delta_M']['current'] *
                   data.mcmc_parameters['Delta_M']['scale'])

        # Compute the covariance matrix
        # The module numexpr is used for doing quickly the long multiplication
        # of arrays (factor of 3 improvements over numpy). It is used as a
        # replacement of blas routines cblas_dcopy and cblas_daxpy
        # For numexpr to work, we need (seems like a bug, but anyway) to create
        # local variables holding the arrays. This cost no time (it is a simple
        # pointer assignment)
        C00, C11, C22 = self.C00, self.C11, self.C22
        C01, C02, C12 = self.C01, self.C02, self.C12
        cov = ne.evaluate(
            "(C00 + alpha**2*C11 + beta**2*C22"
            "+2.*alpha*C01 -2.*beta*C02 -2.*alpha*beta*C12)")

        # Compute the residuals (estimate of distance moduli - exact moduli)
        residuals = np.empty((size,))
        sn = self.light_curve_params
        # This operation loops over all supernovae!
        # Compute the approximate moduli
        residuals = sn.mb - (
            M - alpha*sn.x1 + beta*sn.color + Delta_M*(
                sn.thirdvar > self.scriptmcut))
        # Remove from the approximate moduli the one computed from CLASS
        residuals -= moduli

        # Update the diagonal terms of the covariance matrix with the
        # statistical error
        cov += np.diag(sn.dmb**2 + (alpha*sn.dx1)**2 + (beta*sn.dcolor)**2
                       + 2.*alpha*sn.cov_m_s
                       - 2.*beta*sn.cov_m_c
                       - 2.*alpha*beta*sn.cov_s_c)

        # Whiten the residuals, in two steps
        # 1) Compute the Cholesky decomposition of the covariance matrix, in
        # place. This is a time expensive (0.015 seconds) part
        cov = la.cholesky(cov, lower=True, overwrite_a=True)

        # 2) Solve the triangular system, also time expensive (0.02 seconds)
        residuals = la.solve_triangular(cov, residuals, lower=True, check_finite=False)

        # Finally, compute the chi2 as the sum of the squared residuals
        chi2 = (residuals**2).sum()

        return -0.5 * chi2
