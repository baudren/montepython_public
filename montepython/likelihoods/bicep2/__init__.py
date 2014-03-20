import numpy as np
import math
from montepython.likelihood_class import Likelihood
import montepython.io_mp as io_mp
# import the python package of the BICEP2 collaboration
import bicep_util as bu


class bicep2(Likelihood):

    def __init__(self, path, data, command_line):

        # Define a default value, empty, for which field to use with the
        # likelihood. The proper value should be initialised in the .data file.
        self.field = ""
        Likelihood.__init__(self, path, data, command_line)

        # Recover all the relevant quantities for the likelihood computation
        # from the BICEP collaboration, which includes the band power window
        # functions (bpwf). It assumes that in the "root" directory, there is a
        # "windows" directory which contains the band power window functions
        # from BICEP2.
        C_l, C_l_hat, N_l, C_fl, M_inv, bpwf_l, bpwf_Cs_l = bu.init(
            "bicep2",
            self.field,
            self.data_directory)

        # Store the useful quantities
        self.C_l = C_l
        self.C_l_hat = C_l_hat
        self.N_l = N_l
        self.C_fl = C_fl
        self.M_inv = M_inv
        self.bpwf_l = bpwf_l
        self.bpwf_Cs_l = bpwf_Cs_l

        # Read the desired max ell from the band power window function.
        self.l_max = self.bpwf_l[-1]
        print self.l_max

        # Require tensor modes from Class
        arguments = {
            'output': 'tCl pCl lCl',
            'lensing': 'yes',
            'modes': 's, t',
            'l_max_scalars': self.l_max,
            'l_max_tensors': self.l_max,}
        self.need_cosmo_arguments(data, arguments)

    def loglkl(self, cosmo, data):

        # Recover Cl_s from CLASS, which is a dictionary, with the method
        # get_cl from the Likelihood class, because it already makes the
        # conversion to uK^2.
        dict_Cls = self.get_cl(cosmo, self.l_max)

        # Convert the dict to the same format expected by BICEP
        # that is:
        # 0: TT
        # 1: TE
        # 2: EE
        # 3: BB
        # 6: ET, and the rest to 0 (must be an array of width 9)
        # Warnings: BICEP2 expects l*(l+1)*Cl/(2*pi) in units of uK^2
        cosmo_Cls = np.zeros((self.l_max+1, 9))
        ell = np.arange(self.l_max+1)

        cosmo_Cls[:, 0] = ell*(ell+1)*dict_Cls['tt']/(2*math.pi)
        cosmo_Cls[:, 1] = ell*(ell+1)*dict_Cls['te']/(2*math.pi)
        cosmo_Cls[:, 2] = ell*(ell+1)*dict_Cls['ee']/(2*math.pi)
        cosmo_Cls[:, 3] = ell*(ell+1)*dict_Cls['bb']/(2*math.pi)
        cosmo_Cls[:, 6] = ell*(ell+1)*dict_Cls['te']/(2*math.pi)

        # Get the expectation value of the data considering this theoretical
        # model
        expectation_value = bu.calc_expvals(
            ell, cosmo_Cls,
            self.bpwf_l, self.bpwf_Cs_l)

        # Fill the C_l matrix
        if self.field == "T":
            self.C_l[:, 0, 0] = expectation_value[:, 0]
        elif self.field == 'E':
            self.C_l[:, 0, 0] = expectation_value[:, 2]
        elif self.field == 'B':
            self.C_l[:, 0, 0] = expectation_value[:, 3]
        elif self.field == "EB":
            self.C_l[:, 0, 0] = expectation_value[:, 2]
            self.C_l[:, 0, 1] = expectation_value[:, 5]
            self.C_l[:, 1, 0] = expectation_value[:, 5]
            self.C_l[:, 1, 1] = expectation_value[:, 3]
        elif self.field == "TB":
            self.C_l[:, 0, 0] = expectation_value[:, 0]
            self.C_l[:, 0, 1] = expectation_value[:, 4]
            self.C_l[:, 1, 0] = expectation_value[:, 5]
            self.C_l[:, 1, 1] = expectation_value[:, 3]
        elif self.field == "TE":
            self.C_l[:, 0, 0] = expectation_value[:, 0]
            self.C_l[:, 0, 1] = expectation_value[:, 1]
            self.C_l[:, 1, 0] = expectation_value[:, 1]
            self.C_l[:, 1, 1] = expectation_value[:, 2]
        else:
            raise io_mp.LikelihoodError(
                "BICEP2 requires a field to compute the likelihood"
                " which should so far be T, E or B, but was read to"
                " be '%s'" % self.field)

        # Add the noise
        self.C_l += self.N_l

        # Actually compute the likelihood
        loglkl = bu.evaluateLikelihood(
            self.C_l, self.C_l_hat, self.C_fl, self.M_inv)

        return loglkl
