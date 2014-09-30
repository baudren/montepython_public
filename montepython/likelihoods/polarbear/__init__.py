import os
from math import pi
import numpy as np
from montepython.likelihood_class import Likelihood


class polarbear(Likelihood):

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)

        # Read the four data points from the bandpower file
        self.bandpowers = np.loadtxt(os.path.join(
            self.data_directory, self.bandpower_file))

        # Read the band power window function (bpwf hereafter... yes, but
        # sometimes, explicit is too much)
        self.bpwf = self.load_bandpower_window_function(os.path.join(
            self.data_directory, self.bpwf_file))

        # l_max is now read from the bandpower window functions
        self.l_max = int(self.bpwf[0][-1, 0])

        # Require polarization from class
        arguments = {
            'output': 'tCl pCl lCl',
            'lensing': 'yes',
            'l_max_scalars': self.l_max}
        self.need_cosmo_arguments(data, arguments)

    def load_bandpower_window_function(self, path):
        """
        Read n^th blocks in the bpwf_file

        """
        size = len(self.bandpowers)

        empty_lines = 0
        blocks = []
        with open(path, 'r') as bpfw_file:
            for line in bpfw_file:
                # Check for comments
                if not line or line.startswith('#'):
                    # If it is the first one: new block
                    if empty_lines == 0:
                        blocks.append([])
                    empty_lines += 1
                # Non empty line: add line in current(last) block
                elif line.strip():
                        empty_lines = 0
                        clean_line = line.strip()
                        blocks[-1].append(
                            [float(e) for e in clean_line.split()])

        # Convert everything to numpy arrays
        blocks = [np.array(block) for block in blocks]

        # Check that sufficiently many blocks were read
        assert len(blocks) == size

        return blocks

    def loglkl(self, cosmo, data):
        # Recover the Cl_BB from CLASS
        cls = self.get_cl(cosmo, self.l_max)
        ell = cls['ell']
        cls_bb = cls['bb']*ell*(ell+1.)/(2.*pi)

        # Recover the predicted Cl_BB for each of the four bandpowers
        BB_th = []
        for block in self.bpwf:
            # each block contains a window function
            integrand = np.array(
                [block[index, 1]*cls_bb[e]
                 for index, e in enumerate(block[:, 0])])
            convolution = 0.5*((integrand[1:]+integrand[:-1])*(
                block[1:, 0]-block[:-1, 0])).sum()
            BB_th.append(convolution)

        BB_exp = self.bandpowers[:, 3]
        Delta_BB = self.bandpowers[:, 4]

        chi2 = ((BB_th - BB_exp)**2/(Delta_BB)**2).sum()

        return -0.5*chi2
