"""
.. module:: data
   :synopsis: Define the data class

.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>
"""
import os
import sys
import math
import random as rd
import numpy as np             # Numerical Python module

# A modified version of Python dictionary in order to keep track of the order
# in it (much the same as in an array). In case an older version of Python is
# used, this module does not belong to collections. Please remind to put that
# into your PYTHONPATH variable.
try:
    from collections import OrderedDict as od
except:
    from ordereddict import OrderedDict as od
from datetime import date

import io_mp  # Needs to talk to io_mp.py file for the logging of parameters


class data(object):
    """
    Store all relevant data to communicate between the different modules.

    """

    def __init__(self, command_line, path):
        """
        The data class holds the cosmological information, the parameters from
        the MCMC run, the information coming from the likelihoods. It is a wide
        collections of information, with in particular two main dictionaries:
        cosmo_arguments and mcmc_parameters.

        It defines several useful **methods**. The following ones are called
        just once, at initialization:

        * :func:`fill_mcmc_parameters`
        * :func:`from_input_to_mcmc_parameters`
        * :func:`read_file`
        * :func:`read_version`
        * :func:`group_parameters_in_blocks`

        On the other hand, these two following functions are called every step.

        * :func:`check_for_slow_step`
        * :func:`update_cosmo_arguments`

        Finally, the convenient method :func:`get_mcmc_parameters` will be
        called in many places, to return the proper list of desired parameters.

        It has a number of different **attributes**, and the more important
        ones are listed here:

        * :attr:`cosmo_arguments`
        * :attr:`mcmc_parameters`
        * :attr:`need_cosmo_update`
        * :attr:`log_flag`
        * :attr:`boundary_loglike`

        .. note::

            The `experiments` attribute is extracted from the parameter file,
            and contains the list of likelihoods to use

        To create an instance of this class, one must feed the following
        parameters and keyword arguments:

        :Parameters:
            - **command_line** (`dict`) - dictionary containing the input from the
              :mod:`parser_mp`. It stores the input parameter file, the
              jumping methods, the output folder, etc...
              Most of the information extracted from the command_file will
              be transformed into :class:`data` attributes, whenever it felt
              meaningful to do so.

            - **path** (`dict`) - contains a dictionary of important local paths.
              It is used here to find the cosmological module location.

        """

        # Initialisation of the random seed
        rd.seed()

        # Store the parameter file
        self.param = command_line.param

        # Recover jumping method from command_line
        self.jumping = command_line.jumping
        self.jumping_factor = command_line.jumping_factor
        self.path = path

        self.boundary_loglike = -1e30
        """
        Define the boundary loglike, the value used to defined a loglike that
        is out of bounds. If a point in the parameter space is affected to this
        value, it will be automatically rejected, hence increasing the
        multiplicity of the last accepted point.
        """

        # Creation of the two main dictionnaries:
        self.cosmo_arguments = {}
        """
        Simple dictionary that will serve as a communication interface with the
        cosmological code. It contains all the parameters for the code that
        will not be set to their default values.  It is updated from
        :attr:`mcmc_parameters`.

        :rtype:   dict
        """
        self.mcmc_parameters = od()
        """
        Ordered dictionary of dictionaries, it contains everything needed by
        the :mod:`mcmc` module for the MCMC procedure.  Every parameter name
        will be the key of a dictionary, containing the initial configuration,
        role, status, last accepted point and current point.

        :rtype: ordereddict
        """

        # Read from the parameter file to fill properly the mcmc_parameters
        # dictionary.
        self.fill_mcmc_parameters()

        # Determine which cosmological code is in use
        if path['cosmo'].find('class') != -1:
            self.cosmological_module_name = 'Class'
        else:
            self.cosmological_module_name = None

        # Recover the cosmological code version (and subversion if relevant).
        # To implement a new cosmological code, please add another case to the
        # test below.
        if self.cosmological_module_name == 'Class':
            svn_file = open(path['cosmo']+'/include/svnversion.h', 'r')
            self.subversion = svn_file.readline().split()[-1].\
                replace('"', '')
            svn_file.close()
            for line in open(path['cosmo']+'/include/common.h', 'r'):
                if line.find('_VERSION_') != -1:
                    self.version = line.split()[-1].replace('"', '')
                    break
        else:  # read in the existing parameter file
            self.read_version(self.param_file)

        # End of initialisation with the parameter file
        self.param_file.close()

        self.log_flag = False
        """
        Stores the information whether or not the likelihood data files need to
        be written down in the log.param file. Initially at False.

        :rtype: bool
        """

        self.need_cosmo_update = True
        """
        `added in version 1.1.1`. It stores the truth value of whether the
        cosmological block of parameters was changed from one step to another.
        See :meth:`group_parameters_in_blocks`

        :rtype: bool
        """

        sys.stdout.write('Testing likelihoods for:\n -> ')
        for i in range(len(self.experiments)):
            sys.stdout.write(self.experiments[i]+', ')
        sys.stdout.write('\n')

        # logging the parameter file (only if folder does not exist !)
        if command_line.folder[-1] != '/':
            command_line.folder += '/'
        if (os.path.exists(command_line.folder) and
                not os.path.exists(command_line.folder+'log.param')):
            if command_line.param is not None:
                io_mp.message(
                    "Detecting empty folder, logging the parameter file",
                    "warning")
                io_mp.log_parameters(self, command_line)
                self.log_flag = True
        if not os.path.exists(command_line.folder):
            os.mkdir(command_line.folder)
            # Logging of parameters
            io_mp.log_parameters(self, command_line)
            self.log_flag = True

        self.lkl = od()

        # adding the likelihood directory to the path, to import the module
        # then, for each library, calling an instance of the likelihood.
        # Beware, though, if you add new likelihoods, they should go to the
        # folder likelihoods/yourlike/yourlike.py, and contain a yourlike.data,
        # otherwise the following set of commands will not work anymore.

        # For the logging if log_flag is True, each likelihood will log its
        # parameters

        for elem in self.experiments:

            folder = os.path.abspath(
                path['MontePython'])+"/../likelihoods/%s" % elem
            # add the folder of the likelihood to the path of libraries to...
            if folder not in sys.path:
                sys.path.insert(0, folder)
            # ... import easily the likelihood.py program
            exec "import %s" % elem
            # Initialize the likelihoods. Depending on the values of
            # command_line and log_flag, the routine will call slightly different
            # things. If log_flag is True, the log.param will be appended.
            exec "self.lkl['%s'] = %s.%s('%s/%s.data',\
                self,command_line)" % (
                elem, elem, elem, folder, elem)

        # Storing parameters by blocks of speed
        self.group_parameters_in_blocks()

        # Finally, log the cosmo_arguments used. This comes in the end, because
        # it can be modified inside the likelihoods init functions
        if self.log_flag:
            io_mp.log_cosmo_arguments(self, command_line)
            io_mp.log_default_configuration(self, command_line)

    def fill_mcmc_parameters(self):
        """
        Initializes the ordered dictionary :attr:`mcmc_parameters` from
        the input parameter file.

        It uses :meth:`read_file`, and calls
        :meth:`from_input_to_mcmc_parameters` to actually fill in
        :attr:`mcmc_parameters`.

        """

        # Define temporary quantities, only to simplify the input in the
        # parameter file
        self.parameters = od()

        # Read from the parameter file everything
        try:
            self.param_file = open(self.param, 'r')
        except IOError:
            io_mp.message(
                "Error in initializing the data class, the parameter file \
                {0} does not point to a proper file".format(self.param),
                "error")
        self.read_file(self.param_file)

        # Transform from parameters dictionnary to mcmc_parameters dictionary
        # of dictionaries, method defined just below
        self.from_input_to_mcmc_parameters(self.parameters)

    def from_input_to_mcmc_parameters(self, dictionary):
        """
        Converts dictionary of raw quantities into a meaningful one.

        At the end of this initialization, every field but one is filled for
        every parameter, be it fixed or varying. The missing field is the
        'last_accepted' one, that will be filled in the module :mod:`mcmc`.

        The other fields are

        `initial`:
            initial array of input values defined in the parameter file.
            Contains (in this order) `mean`, `minimum`, `maximum`, `1-sigma`.
            If the min/max values (**TO CHECK** proposal density boundaries)
            are unimportant/unconstrained, use `None` or `-1` (without a period
            !)
        `scale`:
            5th entry of the initial array in the parameter file.
        `role`:
            6th entry of the initial array, can be `cosmo`, `nuisance` or
            `derived`. A `derived` parameter will not be considered as varying,
            but will be instead recovered from the cosmological code for each
            point in the parameter space.
        `tex_name`:
            A tentative tex version of the name, provided by the function
            :func:`io_mp.get_tex_name`.
        `status`:
            Depending on the `1-sigma` value in the initial array, it will be
            set to `fixed` or `varying` (resp. zero and non-zero)
        `current`:
            Stores the value at the current point in parameter space (`not
            allowed initially`)

        .. note::

            The syntax of the parameter files is defined here - if one
            wants to change it, one should report the changes in there.

        :Parameters:
            - **dictionary** (`dict`) - raw dictionary containing the input
              from the parameter file. Its content will be transformed and
              processed into the final :attr:`mcmc_parameters`

        """
        for key, value in dictionary.iteritems():
            self.mcmc_parameters[key] = od()
            self.mcmc_parameters[key]['initial'] = value[0:4]
            self.mcmc_parameters[key]['scale'] = value[4]
            self.mcmc_parameters[key]['role'] = value[-1]
            self.mcmc_parameters[key]['tex_name'] = io_mp.get_tex_name(key)
            if value[3] == 0:
                self.mcmc_parameters[key]['status'] = 'fixed'
                self.mcmc_parameters[key]['current'] = value[0]
            else:
                self.mcmc_parameters[key]['status'] = 'varying'

    def read_file(self, File):
        """
        Execute all lines concerning the data class from a parameter file

        All lines starting with `data.` will be replaced by `self.`, so the
        current instance of the class will contain all the information.

        .. note::

            A rstrip() was added at the end, because of an uncomprehensible bug
            on some systems that imagined some unexistant characters at the end
            of the line... Now should work

        """
        for line in File:
            if line.find('#') == -1:
                if line.split('=')[0].find('data.') != -1:
                    exec(line.replace('data.', 'self.').rstrip())
        File.seek(0)

    def group_parameters_in_blocks(self):
        """
        Regroup mcmc parameters by blocks of same speed

        This method divides all varying parameters from :attr:`mcmc_parameters`
        into as many categories as there are likelihoods, plus one (the slow
        block of cosmological parameters).

        It creates the attribute :attr:`blocks_parameters`, to be used in the
        module :mod:`mcmc`.

        .. note::

            It does not compute by any mean the real speed of each parameter,
            instead, every parameter belonging to the same likelihood will
            be considered as fast as its neighbour.

        .. warning::

            It assumes that the nuisance parameters are already written
            sequentially, and grouped together (not necessarilly in the order
            described in :attr:`experiments`). If you mix up the different
            nuisance parameters in the .param file, this routine will not
            method as intended. It also assumes that the cosmological
            parameters are written at the beginning of the file.

        """
        array = []
        # First obvious block is all cosmological parameters
        array.append(len(self.get_mcmc_parameters(['varying', 'cosmo'])))
        # Then, store all nuisance parameters
        nuisance = self.get_mcmc_parameters(['varying', 'nuisance'])

        for likelihood in self.lkl.itervalues():
            count = 0
            for elem in nuisance:
                if elem in likelihood.nuisance:
                    count += 1
            likelihood.varying_nuisance_parameters = count

        # Then circle through them
        index = 0
        while index < len(nuisance):
            elem = nuisance[index]
            flag = False
            # For each one, check if they belong to a likelihood
            for likelihood in self.lkl.itervalues():
                if elem in likelihood.nuisance:
                    # If yes, store the number of nuisance parameters needed
                    # for this likelihood.
                    flag = True
                    array.append(
                        likelihood.varying_nuisance_parameters+array[-1])
                    index += likelihood.varying_nuisance_parameters
                    continue
            if not flag:
                # If the loop reaches this part, it means this nuisance
                # parameter was associated with no likelihood: this should not
                # happen
                print '/!\ nuisance parameter {0}'.format(elem),
                print 'is associated to no likelihood'
                exit()
        # Store the result
        self.blocks_parameters = array

    def read_version(self, File):
        """
        Extract version and subversion from an existing log.param
        """
        # Read the first line (cosmological code version)
        first_line = File.readline()
        self.version = first_line.split()[1]
        self.subversion = first_line.split()[-1].replace(')', '').\
            replace('-', '')
        File.seek(0)

    def get_mcmc_parameters(self, table_of_strings):
        """
        Returns an ordered array of parameter names filtered by
        `table_of_strings`.

        :Parameters:
            - **table_of_strings** (`list`) - List of strings whose role and
              status must be matched by a parameter. For instance,

              >>> get_mcmc_parameters(['varying'])

              will return a list of all the varying parameters, both
              cosmological and nuisance ones (derived parameters being `fixed`,
              they wont be part of this list). Instead,

              >>> get_mcmc_parameters(['nuisance', 'varying'])

              will only return the nuisance parameters that are being varied.

        """
        table = []
        for key, value in self.mcmc_parameters.iteritems():
            number = 0
            for subkey, subvalue in value.iteritems():
                for string in table_of_strings:
                    if subvalue == string:
                        number += 1
            if number == len(table_of_strings):
                table.append(key)
        return table

    def check_for_slow_step(self, new_step):
        """
        Check whether the value of cosmological parameters were
        changed, and if no, skip computation of the cosmology.

        """
        parameter_names = self.get_mcmc_parameters(['varying'])
        cosmo_names = self.get_mcmc_parameters(['cosmo'])

        need_change = 0

        # For all elements in the varying parameters:
        for elem in parameter_names:
            i = parameter_names.index(elem)
            # If it is a cosmological parameter
            if elem in cosmo_names:
                if self.mcmc_parameters[elem]['current'] != new_step[i]:
                    need_change += 1

        # If any cosmological value was changed,
        if need_change > 0:
            self.need_cosmo_update = True
        else:
            self.need_cosmo_update = False

        for likelihood in self.lkl.itervalues():
            # If the cosmology changed, you need to recompute the likelihood
            # anyway
            if self.need_cosmo_update:
                likelihood.need_update = True
                continue
            # Otherwise, check if the nuisance parameters of this likelihood
            # were changed
            need_change = 0
            for elem in parameter_names:
                i = parameter_names.index(elem)
                if elem in likelihood.nuisance:
                    if self.mcmc_parameters[elem]['current'] != new_step[i]:
                        need_change += 1
            if need_change > 0:
                likelihood.need_update = True
            else:
                likelihood.need_update = False

    def update_cosmo_arguments(self):
        """
        Put in :attr:`cosmo_arguments` the current values of
        :attr:`mcmc_parameters`

        This method is called at every step in the Markov chain, to update the
        dictionary. In the Markov chain, the scale is not remembered, so one
        has to apply it before giving it to the cosmological code.

        .. note::

            When you want to define new parameters in the Markov chain that do
            not have a one to one correspondance to a cosmological name, you
            can redefine its behaviour here. You will find in the source
            several such examples.

        """
        # For all elements in any cosmological parameters
        for elem in self.get_mcmc_parameters(['cosmo']):
            # Fill in the dictionnary with the current value of parameters
            self.cosmo_arguments[elem] = \
                self.mcmc_parameters[elem]['current'] *\
                self.mcmc_parameters[elem]['scale']

        # For all elements in the cosmological parameters from the mcmc list,
        # translate any-one that is not directly a Class parameter into one.
        # The try: except: syntax ensures that the first call
        for elem in self.get_mcmc_parameters(['cosmo']):
            # infer h from Omega_Lambda and delete Omega_Lambda
            if elem == 'Omega_Lambda':
                omega_b = self.cosmo_arguments['omega_b']
                omega_cdm = self.cosmo_arguments['omega_cdm']
                Omega_Lambda = self.cosmo_arguments['Omega_Lambda']
                self.cosmo_arguments['h'] = math.sqrt(
                    (omega_b+omega_cdm) / (1.-Omega_Lambda))
                del self.cosmo_arguments[elem]
            # infer omega_cdm from Omega_L and delete Omega_L
            if elem == 'Omega_L':
                omega_b = self.cosmo_arguments['omega_b']
                h = self.cosmo_arguments['h']
                Omega_L = self.cosmo_arguments['Omega_L']
                self.cosmo_arguments['omega_cdm'] = (1.-Omega_L)*h*h-omega_b
                del self.cosmo_arguments[elem]
            if elem == 'ln10^{10}A_s':
                self.cosmo_arguments['A_s'] = math.exp(
                    self.cosmo_arguments[elem]) / 1.e10
                del self.cosmo_arguments[elem]
            if elem == 'exp_m_2_tau_As':
                tau_reio = self.cosmo_arguments['tau_reio']
                self.cosmo_arguments['A_s'] = self.cosmo_arguments[elem] * \
                    math.exp(2.*tau_reio)
                del self.cosmo_arguments[elem]
            if elem == 'f_cdi':
                self.cosmo_arguments['n_cdi'] = self.cosmo_arguments['n_s']
            if elem == 'beta':
                self.cosmo_arguments['alpha'] = 2.*self.cosmo_arguments['beta']
            # We only do that on xe_1, for there is at least one of them.
            if elem.find('xe_1') != -1:
                # To pass this option, you must have set a number of
                # cosmological settings reio_parametrization to reio_bins_tanh,
                # binned_reio_z set, and binned_reio_num First, you need to set
                # reio_parametrization to reio_bins_tanh
                    if (self.cosmo_arguments['reio_parametrization'] !=
                            'reio_bins_tanh'):
                        io_mp.message(
                            "You set binned_reio_xe to some values \
                            without setting reio_parametrization to \
                            reio_bins_tanh",
                            "error")
                    else:
                        try:
                            size = self.cosmo_arguments['binned_reio_num']
                        except (KeyError):
                            io_mp.message(
                                "You need to set reio_binnumber to the value \
                                corresponding to the one in binned_reio_xe",
                                "error")
                    string = ''
                    for i in range(1, size+1):
                        string += '%.4g' % self.cosmo_arguments['xe_%d' % i]
                        del self.cosmo_arguments['xe_%d' % i]
                        if i != size:
                            string += ','
                    self.cosmo_arguments['binned_reio_xe'] = string

    def __cmp__(self, other):
        """
        Redefinition of the 'compare' method for two instances of this class.

        It will decide which basic operations to perform when the code asked if
        two instances are the same (in case you want to launch a new chain in
        an existing folder, with your own parameter file) Comparing
        cosmological code versions (warning only, will not fail the comparison)

        """
        if self.version != other.version:
            io_mp.message(
                "You are running with a different version of your \
                cosmological code",
                "warning")

        # Defines unordered version of the dictionaries of parameters
        self.uo_parameters = {}
        other.uo_parameters = {}

        # Check if all the experiments are tested again,
        if len(list(set(other.experiments).symmetric_difference(
                set(self.experiments)))) == 0:
            # Check that they have been called with the same .data file, stored
            # in dictionary when initializing.
            for experiment in self.experiments:
                for elem in self.lkl[experiment].dictionary:
                    if self.lkl[experiment].dictionary[elem] != \
                            other.lkl[experiment].dictionary[elem]:
                        print 'in your parameter file: ',
                        print self.lkl[experiment].dictionary
                        print 'in log.param:           ',
                        print other.lkl[experiment].dictionary
                        return -1
            # Fill in the unordered version of dictionaries
            for key, elem in self.mcmc_parameters.iteritems():
                self.uo_parameters[key] = elem['initial']
            for key, elem in other.mcmc_parameters.iteritems():
                other.uo_parameters[key] = elem['initial']

            # And finally compare them (standard comparison between
            # dictionnaries, will return True if both have the same keys and
            # values associated to them.
            return cmp(self.uo_parameters, other.uo_parameters)
        else:
            return -1
