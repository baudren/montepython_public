import os
import numpy as np
import montepython.io_mp as io_mp
from montepython.likelihood_class import Likelihood
import scipy.interpolate as interp
import scipy.constants as conts

class bao_boss_aniso(Likelihood):

    # initialization routine

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)

        # are there conflicting experiments?
        if 'bao_boss_aniso_gauss_approx' in data.experiments:
            raise io_mp.LikelihoodError(
                'conflicting bao_boss_aniso_gauss_approx measurments')

        # self.z, .hdif, .dafid, and .rsfid are read from the data file

        # load the ansio likelihood
        filepath = os.path.join(self.data_directory, self.file)
        # alpha_perp = D_A / rs (rs / DA)_fid
        # alpha_para = (H rs)_fid / (H rs)
        prob_dtype = [('alpha_perp', np.float64),
            ('alpha_para', np.float64),
            ('prob', np.float64)]
        prob = np.loadtxt(
            filepath, delimiter=None,
            comments='#', skiprows=0, dtype=prob_dtype)
        size = np.sqrt(len(prob))
        x = prob['alpha_perp'].reshape(size, size)[:,0]
        y = prob['alpha_para'].reshape(size, size)[0,:]
        Z = prob['prob'].reshape(size, size)
        normZ = np.max(Z)
        Z = Z/normZ
        # use the faster interp.RectBivariateSpline interpolation scheme
        self.prob_interp = interp.RectBivariateSpline(x, y, Z, kx=3, ky=3, s=0)

        # end of initialization

    # compute likelihood

    def loglkl(self, cosmo, data):

        Da = cosmo.angular_distance(self.z)
        H = cosmo.Hubble(self.z) * conts.c / 1000.0
        #dr = self.z / H
        #dv = pow(da * da * (1 + self.z) * (1 + self.z) * dr, 1. / 3.)
        rs = cosmo.rs_drag() * self.rs_rescale

        alpha_perp = Da / rs / (self.Dafid / self.rsfid)
        alpha_para = (self.Hfid * self.rsfid) / (H * rs)

        bao_aniso_like = self.prob_interp(alpha_perp, alpha_para)[0, 0]

        # return ln(L)
        if bao_aniso_like > 0.0:
            lkl = np.log(bao_aniso_like)
        else:
            lkl = data.boundary_loglike

        return lkl
