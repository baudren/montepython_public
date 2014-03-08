from montepython.likelihood_class import Likelihood
import montepython.io_mp as io_mp
import numpy as np


class wmap_9yr(Likelihood):

    def __init__(self, path, data, command_line):

        # Standard initialization, reads the .data
        Likelihood.__init__(self, path, data, command_line)

        # Extra needed cosmological paramters
        self.need_cosmo_arguments(
            data, {'output': 'tCl pCl lCl', 'lensing': 'yes'})

        try:
            import pywlik
        except ImportError:
            raise io_mp.MissingLibraryError(
                "You must first activate the binaries from the wmap wrapper." +
                " Please run : \n " +
                "]$ source /path/to/wrapper_wmap/bin/clik_profile.sh\n " +
                "and try again.")

        # try importing the wrapper_wmap
        self.wmaplike = pywlik.wlik(self.large_data_directory, self.ttmin,
                                    self.ttmax, self.temin, self.temax, self.use_gibbs, self.use_lowlpol)

        # self.cls = np.loadtxt(self.cl_test_file)

        # loglike = self.wmaplike(self.cls)
        # print "got %g expected %g"%(loglike,-845.483)

        self.l_max = max(self.ttmax, self.temax)
        self.need_cosmo_arguments(data, {'l_max_scalars': self.l_max})

        # deal with nuisance parameters
        try:
            self.use_nuisance
        except:
            self.use_nuisance = []
        self.read_contamination_spectra(data)

        pass

    def loglkl(self, cosmo, data):

        nuisance_parameter_names = data.get_mcmc_parameters(['nuisance'])

        # get Cl's from the cosmological code
        cl = self.get_cl(cosmo)

        # add contamination spectra multiplied by nuisance parameters
        cl = self.add_contamination_spectra(cl, data)

        # allocate array of Cl's and nuisance parameters
        tot = np.zeros(np.sum(self.wmaplike.get_lmax()) + 6)

        # fill with Cl's
        index = 0
        for i in range(np.shape(self.wmaplike.get_lmax())[0]):
            if (self.wmaplike.get_lmax()[i] > -1):
                for j in range(self.wmaplike.get_lmax()[i] + 1):
                    if (i == 0):
                        tot[index + j] = cl['tt'][j]
                    if (i == 1):
                        tot[index + j] = cl['ee'][j]
                    if (i == 2):
                        tot[index + j] = cl['bb'][j]
                    if (i == 3):
                        tot[index + j] = cl['te'][j]
                    if (i == 4):
                        tot[index + j] = cl['tb'][j]
                    if (i == 5):
                        tot[index + j] = cl['eb'][j]

                index += self.wmaplike.get_lmax()[i] + 1

        # fill with nuisance parameters
        # for nuisance in nuisance_parameter_names:
            # if nuisance == 'A_SZ':
                # nuisance_value = data.mcmc_parameters[nuisance]['current']*data.mcmc_parameters[nuisance]['scale']
                # tot[index]=nuisance_value
                # index += 1

        # compute likelihood
        lkl = self.wmaplike(tot)[0]

        # add prior on nuisance parameters
        lkl = self.add_nuisance_prior(lkl, data)

        return lkl
