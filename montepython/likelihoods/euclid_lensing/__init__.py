########################################################
# Euclid_lensing likelihood
########################################################
# written by Benjamin Audren
# (adapted from J Lesgourgues's old COSMOS likelihood for CosmoMC)
#
# Modified by S. Clesse in March 2016 to add an optional form of n(z)
# motivated by ground based exp. (Van Waerbeke et al., 2013)
# See google doc document prepared by the Euclid IST - Splinter 2
#
# Modified by J. Lesgourgues in March 2016 to vectorise and speed up

from montepython.likelihood_class import Likelihood
import io_mp
#import time

import scipy.integrate
from scipy import interpolate as itp
import os
import numpy as np
import math


class euclid_lensing(Likelihood):

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)

        # Force the cosmological module to store Pk for redshifts up to
        # max(self.z) and for k up to k_max
        self.need_cosmo_arguments(data, {'output': 'mPk'})
        self.need_cosmo_arguments(data, {'z_max_pk': self.zmax})
        self.need_cosmo_arguments(data, {'P_k_max_h/Mpc': self.k_max_h_by_Mpc})

        # Compute non-linear power spectrum if requested
        if (self.use_halofit):
            self.need_cosmo_arguments(data, {'non linear':'halofit'})

        # Define array of l values, and initialize them
        # It is a logspace
        # find nlmax in order to reach lmax with logarithmic steps dlnl
        self.nlmax = np.int(np.log(self.lmax/self.lmin)/self.dlnl)+1
        # redefine slightly dlnl so that the last point is always exactly lmax
        self.dlnl = np.log(self.lmax/self.lmin)/(self.nlmax-1)
        self.l = self.lmin*np.exp(self.dlnl*np.arange(self.nlmax))

        ########################################################
        # Find distribution of dn_dz (not normalized) in each bin
        ########################################################
        # Assuming each bin contains the same number of galaxies, we find the
        # bin limits in z space
        # Compute the total number of galaxies until zmax (no normalization
        # yet), that is the integral of the galaxy distribution function from 0
        # to self.zmax
        n_tot, error = scipy.integrate.quad(
            self.galaxy_distribution, 0, self.zmax)
        assert error <= 1e-7,  (
            "The integration of the galaxy distribution is not as "
            "precise as expected.")

        # For each bin, compute the limit in z space

        # Create the array that will contain the z boundaries for each bin. The
        # first value is already correctly set to 0.
        self.z_bin_edge = np.zeros(self.nbin+1, 'float64')

        for Bin in xrange(self.nbin-1):
            bin_count = 0.
            z = self.z_bin_edge[Bin]
            while (bin_count <= n_tot/self.nbin):
                gd_1 = self.galaxy_distribution(z)
                gd_2 = self.galaxy_distribution(z+self.dz)
                bin_count += 0.5*(gd_1+gd_2)*self.dz
                z += self.dz
            self.z_bin_edge[Bin+1] = z
        self.z_bin_edge[self.nbin] = self.zmax

        # Fill array of discrete z values
        self.z = np.linspace(0, self.zmax, num=self.nzmax)

        # Fill distribution for each bin (convolving with photo_z distribution)
        self.eta_z = np.zeros((self.nzmax, self.nbin), 'float64')
        gal = self.galaxy_distribution(self.z, True)
        for Bin in xrange(self.nbin):
            low = self.z_bin_edge[Bin]
            hig = self.z_bin_edge[Bin+1]
            for nz in xrange(self.nzmax):
                z = self.z[nz]
                integrand = gal*self.photo_z_distribution(z, self.z, True)
                integrand = np.array([
                    elem if low <= self.z[index] <= hig else 0
                    for index, elem in enumerate(integrand)])
                self.eta_z[nz, Bin] = scipy.integrate.trapz(
                    integrand,
                    self.z)

        # integrate eta(z) over z (in view of normalizing it to one)
        self.eta_norm = np.zeros(self.nbin, 'float64')
        for Bin in xrange(self.nbin):
            self.eta_norm[Bin] = np.sum(0.5*(
                self.eta_z[1:, Bin]+self.eta_z[:-1, Bin])*(
                self.z[1:]-self.z[:-1]))

        ################
        # Noise spectrum
        ################

        # Number of galaxies per steradian
        self.noise = 3600.*self.gal_per_sqarcmn*(180./math.pi)**2

        # Number of galaxies per steradian per bin
        self.noise = self.noise/self.nbin

        # Noise spectrum (diagonal in bin*bin space, independent of l and Bin)
        self.noise = self.rms_shear**2/self.noise

        ###########
        # Read data
        ###########

        # If the file exists, initialize the fiducial values
        # It has been stored flat, so we use the reshape function to put it in
        # the right shape.
        self.Cl_fid = np.zeros((self.nlmax, self.nbin, self.nbin), 'float64')
        self.fid_values_exist = False
        fid_file_path = os.path.join(self.data_directory, self.fiducial_file)
        if os.path.exists(fid_file_path):
            self.fid_values_exist = True
            flat_Cl = np.loadtxt(fid_file_path)
            self.Cl_fid = flat_Cl.reshape((self.nlmax, self.nbin, self.nbin))

        return

    def galaxy_distribution(self, z, array=False):
        """
        Galaxy distribution returns the function D(z) from the notes

        If the array flag is set to True, z is then interpretated as an array,
        and not as a single value.

        Modified by S. Clesse in March 2016 to add an optional form of n(z) motivated by ground based exp. (Van Waerbeke et al., 2013)
        See google doc document prepared by the Euclid IST - Splinter 2
        """

        zmean = 0.9
        z0 = zmean/1.412

        if not array:
            galaxy_dist = z**2*math.exp(-(z/z0)**(1.5))
        elif self.nofz_method==1:
            return z**2*np.exp(-(z/z0)**(1.5))
        else:
            return self.a1*np.exp(-(z-0.7)**2/self.b1**2.)+self.c1*np.exp(-(z-1.2)**2/self.d1**2.)


        return galaxy_dist

    def photo_z_distribution(self, z, zph, array=True):
        """
        Photo z distribution

        If the array flag is set to True, z is then interpretated as an array,
        and not as a single value.
        """

        # Standard error on dz/(1+z)
        sigma_ph = 0.05

        # Note: you must normalize it yourself to one if you want to get nice
        # plots of the galaxy distribution function in each bin (otherwise, the
        # spectra will remain correct, but each D_i(x) will loot strangely
        # normalized when compared to the original D(z)
        if not array:
            photo_z_dist = math.exp(-0.5*(
                (z-zph)/sigma_ph/(1.+z))**2)/sigma_ph/(1.+z)/math.sqrt(
                2.*math.pi)
        else:
            photo_z_dist = np.exp(-0.5*(
                (z-zph)/sigma_ph/(1.+z))**2)/sigma_ph/(1.+z)/math.sqrt(
                2.*math.pi)

        return photo_z_dist

    def loglkl(self, cosmo, data):

        #start = time.time()

        # One wants to obtain here the relation between z and r, this is done
        # by asking the cosmological module with the function z_of_r
        self.r = np.zeros(self.nzmax, 'float64')
        self.dzdr = np.zeros(self.nzmax, 'float64')

        self.r, self.dzdr = cosmo.z_of_r(self.z)

        # Compute now the selection function eta(r) = eta(z) dz/dr normalized
        # to one. The np.newaxis helps to broadcast the one-dimensional array
        # dzdr to the proper shape. Note that eta_norm is also broadcasted as
        # an array of the same shape as eta_z
        self.eta_r = self.eta_z*(self.dzdr[:, np.newaxis]/self.eta_norm)

        # Compute function g_i(r), that depends on r and the bin
        # g_i(r) = 2r(1+z(r)) int_0^+\infty drs eta_r(rs) (rs-r)/rs
        # TODO is the integration from 0 or r ?
        g = np.zeros((self.nzmax, self.nbin), 'float64')
        for Bin in xrange(self.nbin):
            for nr in xrange(1, self.nzmax-1):
                fun = self.eta_r[nr:, Bin]*(self.r[nr:]-self.r[nr])/self.r[nr:]
                g[nr, Bin] = np.sum(0.5*(
                    fun[1:]+fun[:-1])*(self.r[nr+1:]-self.r[nr:-1]))
                g[nr, Bin] *= 2.*self.r[nr]*(1.+self.z[nr])

        # Get power spectrum P(k=l/r,z(r)) from cosmological module
        kmin_in_inv_Mpc = self.k_min_h_by_Mpc * cosmo.h()
        kmax_in_inv_Mpc = self.k_max_h_by_Mpc * cosmo.h()
        pk = np.zeros((self.nlmax, self.nzmax), 'float64')
        for index_l in xrange(self.nlmax):
            for index_z in xrange(1, self.nzmax):

        # These lines would return an error when you ask for P(k,z) out of computed range
        #        if (self.l[index_l]/self.r[index_z] > self.k_max):
        #            raise io_mp.LikelihoodError(
        #                "you should increase euclid_lensing.k_max up to at least %g" % (self.l[index_l]/self.r[index_z]))
        #        pk[index_l, index_z] = cosmo.pk(
        #            self.l[index_l]/self.r[index_z], self.z[index_z])

        # These lines set P(k,z) to zero out of [k_min, k_max] range
                k_in_inv_Mpc =  self.l[index_l]/self.r[index_z]
                if (k_in_inv_Mpc < kmin_in_inv_Mpc) or (k_in_inv_Mpc > kmax_in_inv_Mpc):
                    pk[index_l, index_z] = 0.
                else:
                    pk[index_l, index_z] = cosmo.pk(self.l[index_l]/self.r[index_z], self.z[index_z])

        # Recover the non_linear scale computed by halofit. If no scale was
        # affected, set the scale to one, and make sure that the nuisance
        # parameter epsilon is set to zero
        k_sigma = np.zeros(self.nzmax, 'float64')
        if (cosmo.nonlinear_method == 0):
            k_sigma[:] = 1.e6
        else:
            k_sigma = cosmo.nonlinear_scale(self.z, self.nzmax)

        # Define the alpha function, that will characterize the theoretical
        # uncertainty. Chosen to be 0.001 at low k, raise between 0.1 and 0.2
        # to self.theoretical_error
        alpha = np.zeros((self.nlmax, self.nzmax), 'float64')
        # self.theoretical_error = 0.1
        if self.theoretical_error != 0:
            for index_l in range(self.nlmax):
                k = self.l[index_l]/self.r[1:]
                alpha[index_l, 1:] = np.log(1.+k[:]/k_sigma[1:])/(
                    1.+np.log(1.+k[:]/k_sigma[1:]))*self.theoretical_error

        # recover the e_th_nu part of the error function
        e_th_nu = self.coefficient_f_nu*cosmo.Omega_nu/cosmo.Omega_m()

        # Compute the Error E_th_nu function
        if 'epsilon' in self.use_nuisance:
            E_th_nu = np.zeros((self.nlmax, self.nzmax), 'float64')
            for index_l in range(1, self.nlmax):
                E_th_nu[index_l, :] = np.log(
                    1.+self.l[index_l]/k_sigma[:]*self.r[:]) / (
                    1.+np.log(1.+self.l[index_l]/k_sigma[:]*self.r[:]))*e_th_nu

        # Add the error function, with the nuisance parameter, to P_nl_th, if
        # the nuisance parameter exists
                for index_l in range(self.nlmax):
                    epsilon = data.mcmc_parameters['epsilon']['current']*(
                        data.mcmc_parameters['epsilon']['scale'])
                    pk[index_l, :] *= (1.+epsilon*E_th_nu[index_l, :])

        # Start loop over l for computation of C_l^shear
        Cl_integrand = np.zeros((self.nzmax, self.nbin, self.nbin), 'float64')
        Cl = np.zeros((self.nlmax, self.nbin, self.nbin), 'float64')
        # Start loop over l for computation of E_l
        if self.theoretical_error != 0:
            El_integrand = np.zeros((self.nzmax, self.nbin, self.nbin),
                                    'float64')
            El = np.zeros((self.nlmax, self.nbin, self.nbin), 'float64')

        for nl in xrange(self.nlmax):

            # find Cl_integrand = (g(r) / r)**2 * P(l/r,z(r))
            for Bin1 in xrange(self.nbin):
                for Bin2 in xrange(Bin1,self.nbin):
                    Cl_integrand[1:, Bin1, Bin2] = g[1:, Bin1]*g[1:, Bin2]/(
                        self.r[1:]**2)*pk[nl, 1:]
                    if self.theoretical_error != 0:
                        El_integrand[1:, Bin1, Bin2] = g[1:, Bin1]*(
                            g[1:, Bin2])/(
                            self.r[1:]**2)*pk[nl, 1:]*alpha[nl, 1:]

            # Integrate over r to get C_l^shear_ij = P_ij(l)
            # C_l^shear_ij = 9/16 Omega0_m^2 H_0^4 \sum_0^rmax dr (g_i(r)
            # g_j(r) /r**2) P(k=l/r,z(r))
            # It it then multiplied by 9/16*Omega_m**2 to be in units of Mpc**4
            # and then by (h/2997.9)**4 to be dimensionless
            for Bin1 in xrange(self.nbin):
                for Bin2 in xrange(Bin1,self.nbin):
                    Cl[nl, Bin1, Bin2] = np.sum(0.5*(
                        Cl_integrand[1:, Bin1, Bin2] +
                        Cl_integrand[:-1, Bin1, Bin2])*(
                        self.r[1:]-self.r[:-1]))
                    Cl[nl, Bin1, Bin2] *= 9./16.*(cosmo.Omega_m())**2
                    Cl[nl, Bin1, Bin2] *= (cosmo.h()/2997.9)**4

                    if self.theoretical_error != 0:
                        El[nl, Bin1, Bin2] = np.sum(0.5*(
                            El_integrand[1:, Bin1, Bin2] +
                            El_integrand[:-1, Bin1, Bin2])*(
                            self.r[1:]-self.r[:-1]))
                        El[nl, Bin1, Bin2] *= 9./16.*(cosmo.Omega_m())**2
                        El[nl, Bin1, Bin2] *= (cosmo.h()/2997.9)**4
                    if Bin1 == Bin2:
                        Cl[nl, Bin1, Bin2] += self.noise

        # Write fiducial model spectra if needed (exit in that case)
        if self.fid_values_exist is False:
            # Store the values now, and exit.
            fid_file_path = os.path.join(
                self.data_directory, self.fiducial_file)
            with open(fid_file_path, 'w') as fid_file:
                fid_file.write('# Fiducial parameters')
                for key, value in data.mcmc_parameters.iteritems():
                    fid_file.write(
                        ', %s = %.5g' % (key, value['current']*value['scale']))
                fid_file.write('\n')
                for nl in range(self.nlmax):
                    for Bin1 in range(self.nbin):
                        for Bin2 in range(self.nbin):
                            fid_file.write("%.8g\n" % Cl[nl, Bin1, Bin2])
            print '\n\n /|\    Writing fiducial model in {0}'.format(
                fid_file_path)
            print '/_o_\ for {0} likelihood'.format(self.name)
            return 1j

        # Now that the fiducial model is stored, we add the El to both Cl and
        # Cl_fid (we create a new array, otherwise we would modify the
        # self.Cl_fid from one step to the other)

        # Spline Cl[nl,Bin1,Bin2] along l
        spline_Cl = np.empty((self.nbin, self.nbin), dtype=(list, 3))
        for Bin1 in xrange(self.nbin):
            for Bin2 in xrange(Bin1, self.nbin):
                spline_Cl[Bin1, Bin2] = list(itp.splrep(
                    self.l, Cl[:, Bin1, Bin2]))
                if Bin2 > Bin1:
                    spline_Cl[Bin2, Bin1] = spline_Cl[Bin1, Bin2]

        # Spline El[nl,Bin1,Bin2] along l
        if self.theoretical_error != 0:
            spline_El = np.empty((self.nbin, self.nbin), dtype=(list, 3))
            for Bin1 in xrange(self.nbin):
                for Bin2 in xrange(Bin1, self.nbin):
                    spline_El[Bin1, Bin2] = list(itp.splrep(
                        self.l, El[:, Bin1, Bin2]))
                    if Bin2 > Bin1:
                        spline_El[Bin2, Bin1] = spline_El[Bin1, Bin2]

        # Spline Cl_fid[nl,Bin1,Bin2]    along l
        spline_Cl_fid = np.empty((self.nbin, self.nbin), dtype=(list, 3))
        for Bin1 in xrange(self.nbin):
            for Bin2 in xrange(Bin1, self.nbin):
                spline_Cl_fid[Bin1, Bin2] = list(itp.splrep(
                    self.l, self.Cl_fid[:, Bin1, Bin2]))
                if Bin2 > Bin1:
                    spline_Cl_fid[Bin2, Bin1] = spline_Cl_fid[Bin1, Bin2]

        # Compute likelihood

        # Prepare interpolation for every integer value of l, from the array
        # self.l, to finally compute the likelihood (sum over all l's)
        dof = 1./(int(self.l[-1])-int(self.l[0])+1)

        ells = range(int(self.l[0]), int(self.l[-1])+1)

        # Define cov theory, observ and error on the whole integer range of ell
        # values
        Cov_theory = np.zeros((len(ells), self.nbin, self.nbin), 'float64')
        Cov_observ = np.zeros((len(ells), self.nbin, self.nbin), 'float64')
        Cov_error = np.zeros((len(ells), self.nbin, self.nbin), 'float64')

        for Bin1 in xrange(self.nbin):
            for Bin2 in xrange(Bin1, self.nbin):
                Cov_theory[:, Bin1, Bin2] = itp.splev(
                    ells, spline_Cl[Bin1, Bin2])
                Cov_observ[:, Bin1, Bin2] = itp.splev(
                    ells, spline_Cl_fid[Bin1, Bin2])
                if self.theoretical_error > 0:
                    Cov_error[:, Bin1, Bin2] = itp.splev(
                        ells, spline_El[Bin1, Bin2])
                if Bin2 > Bin1:
                    Cov_theory[:, Bin2, Bin1] = Cov_theory[:, Bin1, Bin2]
                    Cov_observ[:, Bin2, Bin1] = Cov_observ[:, Bin1, Bin2]
                    Cov_error[:, Bin2, Bin1] = Cov_error[:, Bin1, Bin2]

        chi2 = 0.

        # chi2 computation in presence of theoretical error
        # (in absence of it, computation more straightforward, see below)
        # TODO parallelize this
        if (self.theoretical_error > 0):

            for index, ell in enumerate(ells):

                det_theory = np.linalg.det(Cov_theory[index, :, :])
                det_observ = np.linalg.det(Cov_observ[index, :, :])

                det_cross_err = 0
                for i in range(self.nbin):
                    newCov = np.copy(Cov_theory)
                    newCov[:, i] = Cov_error[:, i]
                    det_cross_err += np.linalg.det(newCov)

                # Newton method to minimise chi2 over nuisance parameter epsilon_l
                # (only when using theoretical error scheme of 1210.2194)
                # Find starting point for the method:
                start = 0
                step = 0.001*det_theory/det_cross_err
                error = 1
                old_chi2 = -1.*data.boundary_loglike
                error_tol = 0.01
                epsilon_l = start
                while error > error_tol:
                    vector = np.array([epsilon_l-step,
                                       epsilon_l,
                                       epsilon_l+step])
                # Computing the function on three neighbouring points
                    function_vector = np.zeros(3, 'float64')
                    for k in range(3):
                        Cov_theory_plus_error = Cov_theory+vector[k]*Cov_error
                        det_theory_plus_error = np.linalg.det(
                            Cov_theory_plus_error)
                        det_theory_plus_error_cross_obs = 0
                        for i in range(self.nbin):
                            newCov = np.copy(Cov_theory_plus_error)
                            newCov[:, i] = Cov_observ[:, i]
                            det_theory_plus_error_cross_obs += np.linalg.det(
                                newCov)
                        function_vector[k] = (2.*ell+1.)*self.fsky*(det_theory_plus_error_cross_obs/det_theory_plus_error + math.log(det_theory_plus_error/det_observ) - self.nbin ) + dof*vector[k]**2

                    # Computing first
                    first_d    = (function_vector[2]-function_vector[0]) / (vector[2]-vector[0])
                    second_d = (function_vector[2]+function_vector[0]-2*function_vector[1]) / (vector[2]-vector[1])**2

                    # Updating point and error
                    epsilon_l = vector[1] - first_d/second_d
                    error = abs(function_vector[1] - old_chi2)
                    old_chi2 = function_vector[1]
                # End Newton

                Cov_theory_plus_error = Cov_theory + epsilon_l * Cov_error
                det_theory_plus_error = np.linalg.det(Cov_theory_plus_error)

                det_theory_plus_error_cross_obs = 0
                for i in range(self.nbin):
                    newCov = np.copy(Cov_theory_plus_error)
                    newCov[:, i] = Cov_observ[:, i]
                    det_theory_plus_error_cross_obs += np.linalg.det(newCov)

                chi2 += (2.*ell+1.)*self.fsky*(det_theory_plus_error_cross_obs/det_theory_plus_error + math.log(det_theory_plus_error/det_observ) - self.nbin ) + dof*epsilon_l**2


        # chi2 computation in absence of theoretical error (vectorized)
        else:

            det_theory = np.zeros(len(ells), 'float64')
            det_observ = np.zeros(len(ells), 'float64')
            det_cross_term = np.zeros((self.nbin, len(ells)), 'float64')
            det_cross = np.zeros(len(ells), 'float64')

            det_theory[:] = np.linalg.det(Cov_theory[:, :, :])
            det_observ[:] = np.linalg.det(Cov_observ[:, :, :])

            for i in xrange(self.nbin):
                newCov = np.copy(Cov_theory)
                newCov[:,:, i] = Cov_observ[:,:, i]
                det_cross_term[i,:] = np.linalg.det(newCov[:,:,:])

            det_cross = np.sum(det_cross_term,axis=0)

            for index, ell in enumerate(ells):
                chi2 += (2.*ell+1.)*self.fsky*(det_cross[index]/det_theory[index] + math.log(det_theory[index]/det_observ[index]) - self.nbin)

        # Finally adding a gaussian prior on the epsilon nuisance parameter, if
        # present
        if 'epsilon' in self.use_nuisance:
            epsilon = data.mcmc_parameters['epsilon']['current'] * \
                data.mcmc_parameters['epsilon']['scale']
            chi2 += epsilon**2

        #end = time.time()
        #print "Time in s:",end-start

        return -chi2/2.
