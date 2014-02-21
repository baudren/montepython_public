"""
.. module:: test_montepython

to run it, do
~] nosetest code/test_montepython.py
"""
import unittest
import nose
import os
import datetime
import shutil
import re
import warnings
# This discards warning messages (maybe this should be tuned to discard only
# the ones specifically asked by this code..)
warnings.filterwarnings('ignore')

#---- local imports -----#
import io_mp
import parser_mp
from MontePython import initialise


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
            '-N 10')

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
            '-N 10 -o code/test_config_%s' % self.date)

    def test_no_log_param_in_output_folder(self):
        """
        Create a temp, empty directory without log.param in it
        """
        os.mkdir(self.temp_folder_path)
        self.assertRaises(
            io_mp.ConfigurationError,
            parser_mp.parse,
            '-N 10 -o code/test_config_%s' % self.date)
        os.rmdir(self.temp_folder_path)


class Test02Setup(TestMontePython):
    """Input from known cosmology on one single point"""
    def setUp(self):
        self.date = str(datetime.date.today())
        self.custom_command = (
            '-N 1 -p test.param -o code/test_%s' % self.date)
        self.cosmo, self.data, self.command_line = initialise(
            self.custom_command)

    def tearDown(self):
        del self.custom_command
        del self.cosmo, self.data, self.command_line
        shutil.rmtree('code/test_%s' % self.date)
        del self.date

    def test_folder_created(self):
        """
        Is the initialisation creating a folder?
        """
        assert os.path.exists(
            'code/test_%s' % self.date)

    def test_log_param_written(self):
        """
        Is the log.param properly written?
        """
        assert os.path.exists(
            'code/test_%s/log.param' % self.date)

        # Check if the CLASS version is written properly in the log.param
        with open('code/test_%s/log.param' % self.date, 'r') as log_param:
            first_line = log_param.readline().strip()
            version = first_line.split()[1]
            assert re.match('v[0-9].[0-9].[0-9]', version) is not None

    def test_configuration_file(self):
        """
        Is the .conf recovered and used properly?
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
            self.assertEqual(os.path.abspath(value), self.data.path[key])

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
    """
    def setUp(self):
        self.date = str(datetime.date.today())
        self.custom_command = (
            '-N 1 -p test.param -o code/test_%s' % self.date)
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
        self.custom_command = (
            '-N 1 -p test.param -o code/test_%s' % self.date)
        self.cosmo, self.data, self.command_line = initialise(
            self.custom_command)

    def tearDown(self):
        shutil.rmtree('code/test_%s' % self.date)
        self.cosmo.cleanup()
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
        expected_keys = ['tt', 'te', 'ee', 'bb', 'pp', 'tp']
        for key in raw_cl.iterkeys():
            self.assertIn(key, expected_keys)
        self.cosmo.struct_cleanup()
        self.cosmo.empty()

    def test_lensed_cl_behaviour(self):
        """Are the lensed Cls well behaved?"""
        self.cosmo.set({'output': 'lCl'})
        self.cosmo.compute()
        self.assertTrue(self.cosmo.state)
        lensed_cl = self.cosmo.lensed_cl()
        self.assertIsInstance(lensed_cl, dict)
        expected_keys = ['tt', 'te', 'ee', 'bb', 'pp', 'tp']
        for key in lensed_cl.iterkeys():
            self.assertIn(key, expected_keys)
        self.cosmo.struct_cleanup()
        self.cosmo.empty()


class Test05SamplingMethodBehaviour(TestMontePython):
    """
    Check that all existing sampling method are initializing
    """
    def setUp(self):
        pass

    def tearDown(self):
        pass


if __name__ == '__main__':
    nose.runmodule()
