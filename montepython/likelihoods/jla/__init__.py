"""
.. module:: jla
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
    a 2-fold improvement over a fast Python implementation.

"""
import os
import numpy as np
import numexpr as ne
import scipy.linalg as la
import montepython.io_mp as io_mp
try:
    from pandas import read_table
except ImportError:
    raise io_mp.MissingLibraryError(
        "This likelihood has a lot of IO manipulation. You have "
        "to install the 'pandas' library to use it. Please type:\n"
        "`(sudo) pip install pandas --user`")
from montepython.likelihood_class import Likelihood


class jla(Likelihood):

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)

        # check if jla_simple is also in experiments, in which case, complain
        if 'jla_simple' in data.experiments:
            raise io_mp.LikelihoodError(
                'conflicting jla_simple measurements, you can only have'
                ' either "jla" or "jla_simple" as an experiment, not both')

        # Read the configuration file. Note that we unfortunately can not
        # immediatly execute the file, as it is not formatted as strings.
        self.read_configuration_file()

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

    def read_configuration_file(self):
        """
        Extract Python variables from the configuration file

        This routine performs the equivalent to the program "inih" used in the
        original c++ library.
        """
        settings_path = os.path.join(self.data_directory, self.settings)
        with open(settings_path, 'r') as config:
            for line in config:
                # Dismiss empty lines and commented lines
                if line and line.find('#') == -1:
                    lhs, rhs = [elem.strip() for elem in line.split('=')]
                    # lhs will always be a string, so set the attribute to this
                    # likelihood. The right hand side requires more work.
                    # First case, if set to T or F for True or False
                    if str(rhs) in ['T', 'F']:
                        rhs = True if str(rhs) == 'T' else False
                    # It can also be a path, starting with 'data/'. We remove
                    # this leading folder path
                    elif str(rhs).find('data/') != -1:
                        rhs = rhs.replace('data/', '')
                    else:
                        # Try  to convert it to a float
                        try:
                            rhs = float(rhs)
                        # If it fails, it is a string
                        except ValueError:
                            rhs = str(rhs)
                    # Set finally rhs to be a parameter of the class
                    setattr(self, lhs, rhs)

    def read_matrix(self, path):
        """
        extract the matrix from the path

        This routine uses the blazing fast pandas library (0.10 seconds to load
        a 740x740 matrix). If not installed, it uses a custom routine that is
        twice as slow (but still 4 times faster than the straightforward
        numpy.loadtxt method)

        .. note::

            the length of the matrix is stored on the first line... then it has
            to be unwrapped. The pandas routine read_table understands this
            immediatly, though.

        """
        path = os.path.join(self.data_directory, path)
        # The first line should contain the length.
        with open(path, 'r') as text:
            length = int(text.readline())

        # Note that this function does not require to skiprows, as it
        # understands the convention of writing the length in the first
        # line
        matrix = read_table(path).as_matrix().reshape((length, length))

        return matrix

    def read_light_curve_parameters(self):
        """
        Read the file jla_lcparams.txt containing the SN data

        .. note::

            the length of the resulting array should be equal to the length of
            the covariance matrices stored in C00, etc...

        """
        path = os.path.join(self.data_directory, self.data_file)

        # Recover the names of the columns. The names '3rdvar' and 'd3rdvar'
        # will be changed, because 3rdvar is not a valid variable name
        with open(path, 'r') as text:
            names = [e.strip().replace('3rd', 'third')
                     for e in text.readline()[1:].split()]

        lc_parameters = read_table(
            path, sep=' ', names=names, header=0)

        return lc_parameters

    def loglkl(self, cosmo, data):
        """
        Compute negative log-likelihood (eq.15 Betoule et al. 2014)

        """
        # Recover the distance moduli from CLASS (a size N vector of double
        # containing the predicted distance modulus for each SN in the JLA
        # sample, given the redshift of the supernova.)
        # FIXME what is this number + 25?
        # luminosity_distance returns (comoving_distance(z_cmb) * (1+z_cmb)),
        # whereas the quantity needed for the computation is
        # (comoving_distance(z_cmb) * (1+z_h))
        redshifts = self.light_curve_params.zcmb
        size = redshifts.size

        moduli = np.empty((size, ))
        for index, row in self.light_curve_params.iterrows():
            z_cmb, z_hel = row['zcmb'], row['zhel']
            moduli[index] = cosmo.luminosity_distance(z_cmb)*(
                1.+z_hel)/(1.+z_cmb)
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
            "+2.*alpha*C01 -2.*beta*C02"
            "-2.*alpha*beta*C12)")

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
        cov = la.cholesky(cov, lower=False, overwrite_a=True)
        # 2) Solve the triangular system, also time expensive (0.02 seconds)
        residuals = la.solve_triangular(cov, residuals, check_finite=False)

        # Finally, compute the chi2 as the sum of the squared residuals
        chi2 = (residuals**2).sum()

        return -0.5 * chi2
