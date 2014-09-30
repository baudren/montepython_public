"""
.. module:: test_montepython

to run it, do
~] nosetest tests/test_montepython.py
"""
import unittest
import nose
import os
import datetime
import shutil
import re
import numpy as np
from itertools import count
import warnings
# This discards warning messages (maybe this should be tuned to discard only
# the ones specifically asked by this code..)
warnings.filterwarnings('ignore')

# ---- local imports -----#
from montepython import io_mp
from montepython import parser_mp
from montepython import sampler
from montepython.initialise import initialise
from montepython.run import run
from montepython.analyze import Information


class TestMontePython(unittest.TestCase):
    """
    Generic class modifying the display
    """
    def shortDescription(self):
        return (self.__class__.__name__ +
                ": " + self._testMethodDoc)


class Test01CommandLineInputBehaviour(TestMontePython):
    """
    Testing basic reaction to user input
    """

    def setUp(self):
        """set up the data used in the tests"""
        self.date = str(datetime.date.today())
        self.temp_folder_path = os.path.join(
            os.path.sep.join(os.path.realpath(__file__).split(
                os.path.sep)[:-1]),
            'test_config_'+self.date)

    def tearDown(self):
        """Remove the parser"""
        del self.date
        del self.temp_folder_path

    def test_no_output_folder(self):
        """
        If you ask the code to run without an output folder, it should stop
        """
        self.assertRaises(
            io_mp.ConfigurationError,
            parser_mp.parse,
            'run -N 10')

    def test_inexistant_output_folder(self):
        """
        Providing a non existant output folder, without a param file
        """
        # Check if the folder indeed does not exist
        self.assertFalse(os.path.exists(self.temp_folder_path))
        # Check for proper behaviour
        self.assertRaises(
            io_mp.ConfigurationError,
            parser_mp.parse,
            'run -N 10 -o tests/test_config_%s' % self.date)

    def test_no_log_param_in_output_folder(self):
        """
        Create a temp, empty directory without log.param in it
        """
        os.mkdir(self.temp_folder_path)
        self.assertRaises(
            io_mp.ConfigurationError,
            parser_mp.parse,
            'run -N 10 -o tests/test_config_%s' % self.date)
        os.rmdir(self.temp_folder_path)


class Test02Setup(TestMontePython):
    """Input from known cosmology on one single point"""
    def setUp(self):
        self.date = str(datetime.date.today())
        self.folder = os.path.join(
            'tests', 'test02_%s' % self.date)
        self.custom_command = (
            'run -N 1 -p test.param -o %s' % self.folder)
        try:
            self.cosmo, self.data, self.command_line, _ = initialise(
                self.custom_command)
        except io_mp.ConfigurationError:
            raise io_mp.ConfigurationError(
                "For the test suite to run, you need to have"
                " your default.conf properly configured")

    def tearDown(self):
        del self.custom_command
        del self.cosmo, self.data, self.command_line

        shutil.rmtree(self.folder)
        del self.date

    def test_folder_created(self):
        """
        Is the initialisation creating a folder?
        """
        assert os.path.exists(self.folder)

    def test_log_param_written(self):
        """
        Is the log.param properly written?
        """
        assert os.path.exists(
            os.path.join(self.folder, 'log.param'))

        # Check if the CLASS version is written properly in the log.param
        with open(os.path.join(self.folder, 'log.param'), 'r') as log_param:
            first_line = log_param.readline().strip()
            version = first_line.split()[1]
            assert re.match('v[0-9].[0-9].[0-9]', version) is not None

    def test_configuration_file(self):
        """
        Is the default.conf recovered and used properly?
        """
        # assert that the default.conf exists
        self.assertTrue(
            os.path.exists('default.conf'),
            "You need the file default.conf properly configured")
        # First, check that it has been written in the log.param
        conf_file_path = 'default.conf'
        path = {}
        ## Recover the path dictionnary
        with open(conf_file_path, 'r') as conf_file:
            for line in conf_file:
                exec(line)
        ## Compare with the one stored in data
        for key, value in path.iteritems():
            self.assertEqual(os.path.abspath(
                os.path.expanduser(value)), self.data.path[key])

    def test_likelihood_data_recovered(self):
        """
        TODO: Is the data from the likelihood folder properly handled?
        """
        # A rerun should read the data from the log.param, not from the
        # original likelihood folder
        pass


class Test03NoDefaultConf(TestMontePython):
    """
    Try removing the default.conf, and not specifying any other conf file

    .. warning::

        this Test prevents nosetests to be ran with --processes different from
        0
    """
    def setUp(self):
        self.date = str(datetime.date.today())
        self.folder = os.path.join(
            'tests', 'test03_%s' % self.date)
        self.custom_command = (
            'run -N 1 -p test.param -o %s' % self.folder)
        try:
            shutil.move("default.conf", "default_%s.conf" % self.date)
        except IOError:
            pass

    def tearDown(self):
        try:
            shutil.move("default_%s.conf" % self.date, "default.conf")
        except IOError:
            pass
        del self.date
        del self.custom_command

    def test_no_default_conf(self):
        """
        Is the execution stopped if no file called default.conf is found?
        """
        self.assertRaises(
            io_mp.ConfigurationError,
            initialise,
            self.custom_command)


class Test04CosmologicalCodeWrapper(TestMontePython):
    """
    Check the behaviour of the cosmological code through Monte Python
    """
    def setUp(self):
        self.date = str(datetime.date.today())
        self.folder = os.path.join(
            'tests', 'test04_%s' % self.date)
        self.custom_command = (
            'run -N 1 -p test.param -o %s' % self.folder)
        self.cosmo, self.data, self.command_line, _ = initialise(
            self.custom_command)

    def tearDown(self):
        shutil.rmtree(self.folder)
        self.cosmo.struct_cleanup()
        del self.cosmo, self.data, self.command_line

    def test_has_all_attributes(self):
        """Is the cosmological module containing all needed functions?"""
        # Core module
        methods = [
            'state', 'set', 'empty', 'compute', 'struct_cleanup',
            'get_current_derived_parameters', 'angular_distance',
            'Hubble', 'rs_drag', 'nonlinear_method', 'nonlinear_scale',
            'pk', 'get_pk', 'ionization_fraction', 'baryon_temperature',
            'lensed_cl', 'raw_cl', 'T_cmb', 'h']
        for method in methods:
            self.assertTrue(
                hasattr(self.cosmo, method),
                "%s need to be defined in the cosmological wrapper" % method)

    def test_default_parameters(self):
        """Has the cosmo code a good default behaviour?"""
        # The `state` flag should be initially set to False)
        self.assertFalse(self.cosmo.state)
        self.cosmo.set({})
        self.cosmo.compute()
        # After the computation call, the `state` flag should update to True
        self.assertTrue(self.cosmo.state)
        self.cosmo.struct_cleanup()
        # After a structure cleanup, the cosmological module should be
        # considered not ready
        self.assertFalse(self.cosmo.state)
        # Emptying the parameters arguments (should be empty here, but this
        # should not raise an error)
        self.cosmo.empty()

    def test_pk_behaviour(self):
        """Is the pk function well behaved?"""
        self.cosmo.set({'output': 'mPk'})
        self.cosmo.compute()
        self.assertTrue(self.cosmo.state)
        self.assertIsInstance(self.cosmo.pk(0.01, 0), float)
        self.cosmo.struct_cleanup()
        self.cosmo.empty()

    def test_raw_cl_behaviour(self):
        """Are the raw Cls well behaved?"""
        self.cosmo.set({'output': 'tCl'})
        self.cosmo.compute()
        self.assertTrue(self.cosmo.state)
        raw_cl = self.cosmo.raw_cl()
        self.assertIsInstance(raw_cl, dict)
        expected_keys = ['ell', 'tt', 'te', 'ee', 'bb', 'pp', 'tp']
        for key in raw_cl.iterkeys():
            self.assertIn(key, expected_keys)
        self.cosmo.struct_cleanup()
        self.cosmo.empty()

    def test_lensed_cl_behaviour(self):
        """Are the lensed Cls well behaved?"""
        self.cosmo.set({'output': 'lCl pCl', 'lensing': 'yes'})
        self.cosmo.compute()
        self.assertTrue(self.cosmo.state)
        lensed_cl = self.cosmo.lensed_cl()
        self.assertIsInstance(lensed_cl, dict)
        expected_keys = ['ell', 'tt', 'te', 'ee', 'bb', 'pp', 'tp']
        for key in lensed_cl.iterkeys():
            self.assertIn(key, expected_keys)
        self.cosmo.struct_cleanup()
        self.cosmo.empty()


class Test05DataModule(TestMontePython):
    """
    Check all functionnalities of the Data and Parameter class

    """
    def setUp(self):
        self.date = str(datetime.date.today())
        self.folder = os.path.join(
            'tests', 'test05_%s' % self.date)
        self.number = 30
        self.custom_command = (
            'run -N %d -p test.param -o %s -j fast' % (self.number, self.folder))
        self.cosmo, self.data, self.command_line, _ = initialise(
            self.custom_command)

    def tearDown(self):
        shutil.rmtree(self.folder)
        self.cosmo.struct_cleanup()
        self.cosmo.empty()
        del self.cosmo, self.data, self.command_line

    def test_cosmological_arguments(self):
        """
        Are the cosmological arguments well set?
        """
        # After the initialisation, the dictionnary mcmc_parameters should
        # contain five Parameter instances
        self.assertEqual(
            len(self.data.mcmc_parameters),
            5,
            "mcmc_parameters badly defined")

        # Run the sampler
        sampler.run(self.cosmo, self.data, self.command_line)
        self.assertEqual(
            self.data.mcmc_parameters['omega_b']['current'] *
            self.data.mcmc_parameters['omega_b']['scale'],
            self.data.cosmo_arguments['omega_b'],
            "the cosmo_arguments dict was not updated properly")

    def test_block_behaviour(self):
        """
        Are the mcmc arguments well grouped by block?
        """
        # The block separation must have selected three blocks, of size
        # respectively 2, 1, and 1 (the two test nuisance likelihoods are
        # sharing a nuisance parameter)
        self.assertEqual(
            self.data.block_parameters,
            [2, 3, 4],
            "The block selection went wrong")

        # By default, the over-sampling should be set to 1
        self.assertEqual(
            self.data.over_sampling,
            [1, 1, 1],
            "The default over sampling is messed up")

        # For the sake of the test, set them to something else
        self.data.over_sampling = [1, 3, 5]
        self.data.assign_over_sampling_indices()
        self.assertEqual(
            self.data.over_sampling_indices,
            [0, 1, 2, 2, 2, 3, 3, 3, 3, 3])

        # Run the sampler
        sampler.run(self.cosmo, self.data, self.command_line)

        # Verify that the output is well behaved
        output_path = os.path.join(
            self.folder, '%s_%d__1.txt' % (self.date, self.number))
        with open(output_path, 'r') as output_file:
            chain = np.loadtxt(output_file)
        for index in range(2, 4):
            self.assertLessEqual(
                len(np.unique(chain[:, index])),
                float(self.number)/10*2,
                "%s" % np.unique(chain[:, index]))
        self.assertLessEqual(
            len(np.unique(chain[:, 3])),
            float(self.number)/10*3)
        self.assertLessEqual(
            len(np.unique(chain[:, 4])),
            float(self.number)/10*5)


class Test06MetropolisHastingsImportanceSampling(TestMontePython):
    """
    Check that the default sampling method is working

    The call should create a chain, check that the number of points is properly
    used.
    """
    def setUp(self):
        self.date = str(datetime.date.today())
        self.folder = os.path.join(
            'tests', 'test06_%s' % self.date)
        self.number = 30
        self.custom_command = (
            'run -N %d -p test.param -o %s' % (self.number, self.folder))
        self.cosmo, self.data, self.command_line, _ = initialise(
            self.custom_command)
        sampler.run(self.cosmo, self.data, self.command_line)

    def tearDown(self):
        shutil.rmtree('%s' % self.folder)
        self.cosmo.struct_cleanup()
        self.cosmo.empty()
        del self.cosmo, self.data, self.command_line

    def test_short_run(self):
        """Are the MH and IS sampler basic functionalities working?"""
        # Check that a file with the proper name was created
        output_chain = self.date+'_%d__1.txt' % self.number
        self.assertTrue(os.path.exists(
            os.path.join(self.folder, output_chain)))
        # Recover the data in it
        data = np.loadtxt(os.path.join(self.folder, output_chain))
        # Check that the sum of the first column adds up to the desired number
        self.assertEqual(np.sum(data[:, 0]), self.number)
        # Check that the analyze module works for this
        custom_command = 'info %s' % self.folder
        run(custom_command)
        ## verify that this command created the appropriate files
        expected_files = [
            'test_' + self.date + '.' + elem
            for elem in ['bestfit', 'covmat', 'h_info', 'v_info',
                         'info', 'log', 'tex']]
        for output_file in os.listdir(self.folder):
            if os.path.isfile(output_file):
                if output_file.find('.txt') == -1:
                    self.assertIn(output_file, expected_files)
        ## check if the plot folder was created, with the two pdf files
        self.assertTrue(os.path.isdir(
            os.path.join(self.folder, 'plots')))
        for pdf_file in os.listdir(
                os.path.join(self.folder, 'plots')):
            self.assertTrue(
                pdf_file.find('_triangle') != -1 or pdf_file.find('_1d') != -1)
        # Reset the "Information._id"
        Information._ids = count(0)
        ## Run an Importance sampling run on this data
        is_folder = self.folder + '_is'
        custom_command = (
            'run -p test_is.param -o %s' % is_folder +
            ' -m IS --IS-starting-folder %s' % self.folder)
        run(custom_command)
        # Check that the file has been copied, and that both the likelihood and
        # the multiplicity has been changed
        self.assertIn(output_chain, os.listdir(is_folder),
                      "The IS did not create the right chain")
        # Use these two runs to test the several folders option
        custom_command = 'info %s %s' % (
            self.folder, is_folder)
        run(custom_command)
        # Check that all expected files were created
        short_folder = os.path.join(self.folder.split(os.path.sep)[-1])
        short_is_folder = os.path.join(is_folder.split(os.path.sep)[-1])
        expected_files = set([
            '%s-vs-%s_triangle.pdf' % (short_folder, short_is_folder),
            '%s-vs-%s_1d.pdf' % (short_folder, short_is_folder),
            '%s_triangle.pdf' % short_folder,
            '%s_1d.pdf' % short_folder])
        self.assertSetEqual(expected_files, set(
            os.listdir(os.path.join(self.folder, 'plots'))))
        shutil.rmtree('%s' % is_folder)


class Test07CosmoHammerBehaviour(TestMontePython):
    """
    Check if the modules are callable
    """
    def setUp(self):
        self.date = str(datetime.date.today())
        self.custom_command = (
            'run -N 1 -p test.param -o tests/test_%s' % self.date +
            ' -m CH')
        self.cosmo, self.data, self.command_line, _ = initialise(
            self.custom_command)

    def tearDown(self):
        shutil.rmtree('tests/test_%s' % self.date)
        self.cosmo.empty()
        del self.cosmo, self.data, self.command_line

    def test_callable_objects(self):
        """Are the cosmo and data objects callable?"""
        self.assertTrue(hasattr(self.cosmo, '__call__'))
        self.assertTrue(hasattr(self.data, '__call__'))
        for experiment in self.data.experiments:
            self.assertTrue(hasattr(
                self.data.lkl[experiment], 'computeLikelihood'))


class Test08NestedSamplingBehaviour(TestMontePython):
    """
    Check if Nested Sampling works
    """
    def setUp(self):
        self.date = str(datetime.date.today())
        self.folder = os.path.join('tests', 'test_%s' % self.date)
        self.custom_command = (
            'run -N 1 -p test_gaussian.param -o %s' % self.folder +
            ' -m NS --NS_n_live_points 30 --NS_max_iter 10')
        self.cosmo, self.data, self.command_line, _ = initialise(
            self.custom_command)

    def tearDown(self):
        shutil.rmtree(self.folder)
        self.cosmo.struct_cleanup()
        self.cosmo.empty()
        del self.cosmo, self.data, self.command_line

    def test_behaviour(self):
        """Check Nested Sampling global behaviour"""
        self.assertTrue(os.path.exists(
            os.path.join(self.folder, 'NS')))
        sampler.run(self.cosmo, self.data, self.command_line)


class Test09MPI(TestMontePython):
    """
    Check the MPI behaviour
    """
    pass


if __name__ == '__main__':
    nose.runmodule()
