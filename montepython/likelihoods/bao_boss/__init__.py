import os
import numpy as np
from montepython.likelihood_class import Likelihood
import montepython.io_mp as io_mp
import warnings


class bao_boss(Likelihood):

    # initialization routine

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)

        # exclude the isotropic CMASS experiment when the anisotrpic
        # measurement is also used
        exclude_isotropic_CMASS = False

        conflicting_experiments = [
            'bao_boss_aniso', 'bao_boss_aniso_gauss_approx']
        for experiment in conflicting_experiments:
            if experiment in data.experiments:
                exclude_isotropic_CMASS = True

        if exclude_isotropic_CMASS:
            warnings.warn("excluding isotropic CMASS measurement")
            if not hasattr(self, 'exclude') or self.exclude == None:
                self.exclude = ['CMASS']
            else:
                self.exclude.append('CMASS')

        # define array for values of z and data points
        self.z = np.array([], 'float64')
        self.data = np.array([], 'float64')
        self.error = np.array([], 'float64')
        self.type = np.array([], 'int')

        # read redshifts and data points
        with open(os.path.join(self.data_directory, self.file), 'r') as filein:
            for line in filein:
                if line.strip() and line.find('#') == -1:
                    # the first entry of the line is the identifier
                    this_line = line.split()
                    # insert into array if this id is not manually excluded
                    if not this_line[0] in self.exclude:
                        self.z = np.append(self.z, float(this_line[1]))
                        self.data = np.append(self.data, float(this_line[2]))
                        self.error = np.append(self.error, float(this_line[3]))
                        self.type = np.append(self.type, int(this_line[4]))

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

            da = cosmo.angular_distance(self.z[i])
            dr = self.z[i] / cosmo.Hubble(self.z[i])
            dv = pow(da * da * (1 + self.z[i]) * (1 + self.z[i]) * dr, 1. / 3.)
            rs = cosmo.rs_drag()

            if self.type[i] == 3:
                theo = dv / rs

            elif self.type[i] == 4:
                theo = dv

            elif self.type[i] == 5:
                theo = da / rs

            elif self.type[i] == 6:
                theo = 1. / cosmo.Hubble(self.z[i]) / rs

            elif self.type[i] == 7:
                theo = rs / dv
            else:
                raise io_mp.LikelihoodError(
                    "In likelihood %s. " % self.name +
                    "BAO data type %s " % self.type[i] +
                    "in %d-th line not understood" % i)

            chi2 += ((theo - self.data[i]) / self.error[i]) ** 2

        # return ln(L)
        lkl = - 0.5 * chi2

        return lkl
