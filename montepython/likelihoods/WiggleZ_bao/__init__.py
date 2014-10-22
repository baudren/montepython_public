from montepython.likelihood_class import Likelihood
import numpy as np
from math import log
import os


class WiggleZ_bao(Likelihood):
    """From 1401.0358v2"""

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)

        measurements = np.loadtxt(os.path.join(self.data_directory, self.file))
        self.z_eff, self.Dv = measurements[:, 0], measurements[:, 1]

        # Read the inverse covariance matrix, stored in data
        self.inverse_covmat = np.loadtxt(os.path.join(
            self.data_directory, self.inverse_covmat_file))
        # Multiply by 1e-4, as explained in Table 4
        self.inverse_covmat *= 1e-4

        # The matrix should be symmetric
        assert (self.inverse_covmat.T == self.inverse_covmat).all()

    def loglkl(self, cosmo, data):
        # Modes

        difference = []
        for z, Dv in zip(self.z_eff, self.Dv):
            # Recover the Dv and rs from CLASS
            da = cosmo.angular_distance(z)
            dr = z/cosmo.Hubble(z)
            dv = pow(da**2*(1.+z)**2*dr, 1./3)
            rs = cosmo.rs_drag()
            # To have the same units, we must multiply dv/rs from CLASS by
            # rs_fiducial
            difference.append(Dv - dv/rs*self.rs_fiducial)

        difference = np.array(difference)
        likelihood = np.exp(-0.5 * (
            np.dot(difference, np.dot(self.inverse_covmat, difference.T))))
        return log(likelihood)
