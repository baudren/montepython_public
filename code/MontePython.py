#!/software/base/python-2.7.1-gnu/bin/python

####################
# Monte-Python, a Monte Carlo Markov Chain code (with Class!)
# Version 0.9
# written by Benjamin Audren
####################

#-------------------IMPORT-PACKAGES-------------------------------------

try:
  from classy import Class
except ImportError:
  print " /|\  You must have installed the classy.pyx"
  print "/_o_\ please go to /path/to/class/python and run"
  print "      python setup.py build, followed by python setup.py install --user"
  exit()

import parser	# parsing the input command line
import io	# all the input/output mechanisms
import mcmc	# the actual Monte Carlo chain procedure, along with the useful functions
import data	# data handling
#import analyze	# analysis module, only invoked if asked in the command line

import os,sys

#------------------DEFAULT-CONFIGURATION--------------------------------
path = {}
path['MontePython'] = sys.path[0]
conf_file = path['MontePython']+'/../default.conf'
if os.path.isfile(conf_file):
  for line in open(conf_file):
    exec(line)
else:
  print ' /|\  You must provide a default.conf file'
  print '/_o_\ in your montepython directory that specify'
  print '      the correct locations of MontePython, Class, Clik...'


#------------------MAIN-DEFINITION--------------------------------------
def main():
  # Parsing line argument
  command_line=parser.parse()

  # read the cleaning command: it will get rid of every file in the desired
  # folder that is not a completed MCM chain.
  if command_line.clean is not None:
    io.clean(command_line.clean)
    exit()

  # If info flag was used, read a potential chain (or set of chains) to be
  # analysed with default procedure. If the argument is a .info file, then it
  # will extract information from it (plots to compute, chains to analyse,
  # etc...)
  if command_line.files is not None:
    import analyze	# analysis module, only invoked if asked in the command line
    analyze.info(command_line.files,command_line.bins)
    exit()

  # If the restart flag was used, load the cosmology directly from the
  # log.param file, and append to the existing chain. The log.dat file will be
  # updated at the end of the run as well.
  if command_line.restart is not None:
    command_line.par = command_line.restart.split('/')[0]+'/log.param'
    Data = data.data(command_line,path)

  # Else, fill in data, starting from default_path, then custom parameter file,
  # then additionnal command line arguments.
  # If output folder already exists, first load a data instance with used
  # param, and compare the two instances. If different, displays a warning.
  else:
    Data=data.data(command_line,path)
    if command_line.par.find('log.param')==-1:
      if os.path.exists(command_line.folder+'log.dat'):
	Data_old=data.data(command_line,path,False)
	if Data!=Data_old:
	  print '\n /|\  You are starting a chain in {0} with different parameters\n/_o_\ than used previously.\n      Exiting'.format(command_line.folder)

  # Overwrite arguments from parameter file with the command line
  if command_line.N is None:
    try:
      command_line.N=Data.N
    except AttributeError:
      print '\n /|\  You did not provide a number of steps,\n/_o_\ neither via command line, neither in {0}'.format(command_line.par)
      exit()

  # Logging all configuration
  io.create_output_files(command_line,Data)

  # Load up the cosmological backbone. For the moment, only Class has been wrapped.
  _cosmo=Class()

  # Main chain
  rate,min_LogLike=mcmc.chain(_cosmo,Data,command_line)
  
  # Closing up and loggin the result of the chain stuff
  Data.out.close()
  io.write_log(Data,rate,min_LogLike)
  Data.log.close()

#-----------------MAIN-CALL---------------------------------------------

if __name__ == '__main__':
  main()
