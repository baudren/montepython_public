#!/usr/bin/python

##############################################################
# MontePython, a Monte Carlo Markov Chain code (with Class!)
# Version 1.0.0
# written by Benjamin Audren
# Additional parts by Julien Lesgourgues, Karim Benabed
##############################################################

#-------------------IMPORT-PACKAGES-------------------------------------

# Python modules
import os,sys

# Checking for python version
version = sys.version[:3]
if float(version) < 2.7:
  print '\n\n /|\  You must have Python >= 2.7,' 
  print '/_o_\ or install manually the following modules for your older distribution:'
  print '      argparse and OrderedDict (see MontePython.pdf)'

import parser	# parsing the input command line
import io	# all the input/output mechanisms
import mcmc	# the actual Monte Carlo chain procedure, along with the useful functions
import data	# data handling

#------------------MAIN-DEFINITION--------------------------------------
def main():
  # Parsing line argument
  command_line=parser.parse()

  # Default configuration
  path = {}
  # On execution, sys.path contains all the standard locations for the
  # libraries, plus, on the first position (index 0), the directory from where
  # the code is executed. By default, then, the data folder is located in the
  # same root directory. Any setting in the configuration file will overwrite
  # this one.
  path['MontePython'] = sys.path[0]+'/'
  path['data']        = path['MontePython'][:-5]+'data/'

  # Configuration file, defaulting to default.conf in your root directory. This
  # can be changed with the command line option -conf. All changes will be
  # stored into the log.param of your folder, and hence will be reused for an
  # ulterior run in the same directory
  conf_file = path['MontePython']+'/../'+command_line.config_file
  if os.path.isfile(conf_file):
    for line in open(conf_file):
      exec(line)
    for key,value in path.iteritems():
      if value[-1]!='/':
	path[key] = value+'/'
  else:
    print ' /|\  You must provide a .conf file (default.conf by default)'
    print '/_o_\ in your montepython directory that specifies'
    print '      the correct locations for your data folder, Class (, Clik), ...'

  sys.stdout.write('Running MontePython version 1.0\n')

  # If the info flag was used, read a potential chain (or set of chains) to be
  # analysed with default procedure. If the argument is a .info file, then it
  # will extract information from it (plots to compute, chains to analyse,
  # etc...)
  if command_line.files is not None:
    import analyze	# analysis module, only invoked if asked in the command line
    analyze.info(command_line)
    exit()

  # If the restart flag was used, load the cosmology directly from the
  # log.param file, and append to the existing chain. 
  if command_line.restart is not None:
    if command_line.restart[0] == '/':
      folder = ''
    else:
      folder = './'
    for elem in command_line.restart.split("/")[:-1]:
      folder += ''.join(elem+'/')
    command_line.param = folder+'log.param'
    command_line.folder = folder
    sys.stdout.write('Reading {0} file'.format(command_line.restart))
    Data = data.data(command_line,path)

  # Else, fill in data, starting from  parameter file, If output folder already
  # exists, first load a data instance with used param, and compare the two
  # instances. If different, exit: you are not able to run two different things
  # in one folder.
  else:
    Data=data.data(command_line,path)
    if command_line.param.find('log.param')==-1:
      Data_old=data.data(command_line,path,False)
      if Data!=Data_old:
        print '\n /|\  You are starting a chain in {0} with different parameters\n/_o_\ than used previously.\n      Exiting'.format(command_line.folder)
	exit()

  # Overwrite arguments from parameter file with the command line
  if command_line.N is None:
    try:
      command_line.N=Data.N
    except AttributeError:
      print '\n /|\  You did not provide a number of steps,\n/_o_\ neither via command line, neither in {0}'.format(command_line.param)
      exit()

  # Creating the file that will contain the chain
  io.create_output_files(command_line,Data)

  # Loading up the cosmological backbone. For the moment, only Class has been wrapped.

  #| Importing the python-wrapped Class from the correct folder, defined in the
  #| .conf file, or overwritten at this point by the log.param
  #| If there is a conflict between the log.param value and the .conf file, exiting.
  if Data.path != path:
    print ' /|\  Your log.param files is in contradiction with your .conf file,'
    print '/_o_\ This might mean that you are not computing from the correct Class version'
    print '      Exiting now'
    exit()

  # If the cosmological code is Class, do the following to import all relevant quantities
  if Data.cosmological_module_name == 'Class':
    for elem in os.listdir(Data.path['cosmo']+"python/build"):
      if elem.find("lib.") != -1:
        classy_path = path['cosmo']+"python/build/"+elem

    #| Inserting the previously found path into the list of folders to search for
    #| python modules.
    sys.path.insert(1,classy_path)
    try:
      from classy import Class
    except ImportError:
      print(" /|\  You must have compiled the classy.pyx")
      print("/_o_\ please go to /path/to/class/python and run")
      print("      python setup.py build, then python setup.py install --user")
      exit()

    _cosmo=Class()
  else:
    print(" /|\  Error: unrecognised cosmological module")
    print("/_o_\ Be sure to define the correct behaviour")
    print("      in MontePython.py, and data.py")
    exit()

  # MCMC chain 
  rate,min_LogLike=mcmc.chain(_cosmo,Data,command_line)
  
  # Closing up the file
  Data.out.close()


#-----------------MAIN-CALL---------------------------------------------
if __name__ == '__main__':
  main()
