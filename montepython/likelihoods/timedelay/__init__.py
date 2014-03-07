import os
import numpy as np
from math import log, sqrt, pi
from montepython.likelihood_class import Likelihood


class timedelay(Likelihood):

    # initialization routine

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)

        # define array for values of z and data points
        self.zd = np.array([], 'float64')
        self.zs = np.array([], 'float64')
        self.lambda_d = np.array([], 'float64')
        self.mu_d = np.array([], 'float64')
        self.sigma_d = np.array([], 'float64')

        # read redshifts and data points
        for line in open(os.path.join(
                self.data_directory, self.file), 'r'):
            if (line.find("#") == -1):
                self.zd = np.append(self.zd, float(line.split()[0]))
                self.zs = np.append(self.zs, float(line.split()[1]))
                self.lambda_d = np.append(
                    self.lambda_d, float(line.split()[2]))
                self.mu_d = np.append(self.mu_d, float(line.split()[3]))
                self.sigma_d = np.append(self.sigma_d, float(line.split()[4]))

        # number of data points
        self.num_points = np.shape(self.zd)[0]

        # end of initialization

    # compute likelihood

    def loglkl(self, cosmo, data):

        lkl = 0.

        # for each point, compute angular distance da, radial distance dr,
        # volume distance dv, sound horizon at baryon drag rs_d,
        # theoretical prediction and chi2 contribution
        for i in range(self.num_points):
            Dd = cosmo.angular_distance(self.zd[i])
            Ds = cosmo.angular_distance(self.zs[i])
            Dds = ((1. + self.zs[i]) * Ds - (1 + self.zd[i]) * Dd) / (
                1. + self.zs[i])
            Dt = (1 + self.zd[i]) * Dd * Ds / Dds

            if (Dt > self.lambda_d[i]):
                lkl = lkl - (log(Dt - self.lambda_d[i]) - self.mu_d[i]) ** 2 / 2. / self.sigma_d[
                    i] ** 2 - log(sqrt(2. * pi) * (Dt - self.lambda_d[i]) * self.sigma_d[i])
            else:
                lkl = data.boundary_loglike

        return lkl
