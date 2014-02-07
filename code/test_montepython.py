import unittest
import os
import datetime
#---- local imports -----#
import io_mp
import parser_mp


class TestInputBehaviour(unittest.TestCase):
    """
    Testing basic reaction to user input

    to run it, do
    ~] nosetest code/test_montepython.py -v
    """

    def setUp(self):
        """set up the data used in the tests"""
        self.date = str(datetime.date.today())
        self.temp_folder_path = os.path.join(
            os.path.sep.join(os.path.realpath(__file__).split(
                os.path.sep)[:-1]),
            'test_'+self.date)

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
            '-N 10 -o code/test_%s' % self.date)

    def test_no_log_param_in_output_folder(self):
        """
        Create a temp, empty directory without log.param in it
        """
        os.mkdir(self.temp_folder_path)
        self.assertRaises(
            io_mp.ConfigurationError,
            parser_mp.parse,
            '-N 10 -o code/test_%s' % self.date)
        os.rmdir(self.temp_folder_path)


class TestReferenceCosmology(unittest.TestCase):
    """Input from know cosmology to expect results"""
    
    def setUp(self):
        pass

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
