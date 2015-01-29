import os
import numpy as np
import warnings
import montepython.io_mp as io_mp
from montepython.likelihood_class import Likelihood
import scipy.constants as conts

class bao_boss_aniso_gauss_approx(Likelihood):

    # initialization routine

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)

        # are there conflicting experiments?
        if 'bao_boss_aniso' in data.experiments:
            raise io_mp.LikelihoodError(
                'conflicting bao_boss_aniso measurments')

        # define array for values of z and data points
        self.z = np.array([], 'float64')
        self.DA_rdfid_by_rd_in_Mpc = np.array([], 'float64')
        self.DA_error = np.array([], 'float64')
        self.H_rd_by_rdfid_in_km_per_s_per_Mpc = np.array([], 'float64')
        self.H_error = np.array([], 'float64')
        self.cross_corr = np.array([], 'float64')
        self.rd_fid_in_Mpc = np.array([], 'float64')

        # read redshifts and data points
        i = 0
        with open(os.path.join(self.data_directory, self.file), 'r') as filein:
            for i, line in enumerate(filein):
                if line.strip() and line.find('#') == -1:
                    this_line = line.split()
                    # this_line[0] is some identifier
                    self.z = np.append(self.z, float(this_line[1]))
                    self.DA_rdfid_by_rd_in_Mpc = np.append(
                        self.DA_rdfid_by_rd_in_Mpc, float(this_line[2]))
                    self.DA_error = np.append(
                        self.DA_error, float(this_line[3]))
                    self.H_rd_by_rdfid_in_km_per_s_per_Mpc = np.append(
                        self.H_rd_by_rdfid_in_km_per_s_per_Mpc, float(this_line[4]))
                    self.H_error = np.append(
                        self.H_error, float(this_line[5]))
                    self.cross_corr = np.append(
                        self.cross_corr, float(this_line[6]))
                    self.rd_fid_in_Mpc = np.append(
                        self.rd_fid_in_Mpc, float(this_line[7]))

                    # is the cross correlation coefficient valid
                    if self.cross_corr[i] < -1.0 or self.cross_corr[i] > 1.0:
                        raise io_mp.LikelihoodError(
                            "invalid cross correlation coefficient in entry "
                            "%d: %f" % (i, self.cross_corr[i]))

        # number of data points
        self.num_points = np.shape(self.z)[0]

        # end of initialization

    # compute likelihood

    def loglkl(self, cosmo, data):

        chi2 = 0.

        # for each point, compute angular distance da, radial distance dr,
        # volume distance dv, sound horizon at baryon drag rs_d,
        # theoretical prediction and chi2 contribution
        for i in range(self.num_points):

            DA_at_z = cosmo.angular_distance(self.z[i])
            H_at_z = cosmo.Hubble(self.z[i]) * conts.c / 1000.0
            #dv = pow(da * da * (1 + self.z[i]) * (1 + self.z[i]) * dr, 1. / 3.)
            rd = cosmo.rs_drag() * self.rs_rescale

            theo_DA_rdfid_by_rd_in_Mpc = DA_at_z / rd * self.rd_fid_in_Mpc[i]
            theo_H_rd_by_rdfid = H_at_z * rd / self.rd_fid_in_Mpc[i]

            chi2 += ((theo_DA_rdfid_by_rd_in_Mpc - self.DA_rdfid_by_rd_in_Mpc[i]) / self.DA_error[i]) ** 2
            chi2 += ((theo_H_rd_by_rdfid - self.H_rd_by_rdfid_in_km_per_s_per_Mpc[i]) / self.H_error[i]) ** 2
            # account for cross correlation
            chi2 -= 2 * self.cross_corr[i] \
                * (theo_DA_rdfid_by_rd_in_Mpc - self.DA_rdfid_by_rd_in_Mpc[i]) \
                * (theo_H_rd_by_rdfid - self.H_rd_by_rdfid_in_km_per_s_per_Mpc[i]) \
                / self.DA_error[i] / self.H_error[i]

        # return ln(L)
        lkl = - 0.5 * chi2

        return lkl
