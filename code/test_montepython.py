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
        shutil.rmtree('code/test_%s' % self.date)

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

    def test_likelihood_data_recovered(self):
        """
        Is the data from the likelihood folder properly handled?
        """
        # A rerun should read the data from the log.importparam, not from the
        # original likelihood folder
        pass

    def test_configuration_file(self):
        """
        Is the .conf recovered and used properly?
        """
        pass


class Test03CosmologicalCode(TestMontePython):
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


class Test04SamplingMethodBehaviour(TestMontePython):
    """
    Check that all existing sampling method are initializing
    """
    def setUp(self):
        pass

    def tearDown(self):
        pass


if __name__ == '__main__':
    nose.runmodule()
