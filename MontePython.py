#!/usr/bin/python

####################
# Monte-Python, a Monte Carlo Markov Chain code (with Class!)
# Version 0.9
# written by Benjamin Audren
####################


#------------------DEFAULT-CONFIGURATION--------------------------------
path={}
path['class']		= '/itp/baudren/Desktop/codes/class/'
path['clik']		= '/itp/baudren/Desktop/codes/clik_epfl/examples/'
path['clik_wmap']	= '/itp/baudren/Desktop/codes/clik_epfl/examples/wmap_full.clik'
path['MontePython']	= '/itp/baudren/Desktop/codes/montepython/'


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
import analyze	# analysis module, only invoked if asked in the command line

import os



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
    analyze.info(command_line.files,command_line.bins)
    exit()

  # If the restart flag was used, load the cosmology directly from the
  # log.param file, and append to the existing chain. The log.dat file will be
  # updated at the end of the run as well.
  if command_line.restart is not None:
    command_line.par = command_line.restart.split('/')[0]+'/log.param'
    Data = data.data(command_line.par,path)

  # Else, fill in data, starting from default_path, then custom parameter file,
  # then additionnal command line arguments.
  # If output folder already exists, first load a data instance with used
  # param, and compare the two instances. If different, displays a warning.
  else:
    Data=data.data(command_line.par,path)
    if os.path.isdir(command_line.folder):
      Data_old=data.data(command_line.folder+'/log.param',path,False)
      if Data!=Data_old:
	print '\n /|\  You are starting a chain in {0} with different parameters\n/_o_\ than used previously.\n      Proceeding with the computation'.format(command_line.folder)

  # Overwrite arguments from parameter file with the command line
  if command_line.N is None:
    try:
      command_line.N=Data.N
    except AttributeError:
      print '\n /|\  You did not provide a number of steps,\n/_o_\ neither via command line, neither in {0}'.format(command_line.par)
      exit()

  # Logging all configuration
  out,log,Data.out_name=io.create_output_files(command_line)
  if command_line.restart is None:
    io.log_parameters(Data,path,command_line)

  # Load up the cosmological backbone. For the moment, only Class has been wrapped.
  _cosmo=Class()

  # initialisation of output here
  io.class_output(Data)

  # Main chain
  rate,min_LogLike=mcmc.chain(_cosmo,Data,command_line,out)
  
  # Closing up and loggin the result of the chain stuff
  out.close()
  io.write_log(log,out,Data.param_names,rate,min_LogLike)
  log.close()

#-----------------MAIN-CALL---------------------------------------------

if __name__ == '__main__':
  main()
