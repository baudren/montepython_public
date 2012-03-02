import os,sys
import argparse # Python module to handle command line arguments

# Definition of the object, below will be added all the possible options. This
# will also create an automatic help command, available through the option -h
parser = argparse.ArgumentParser(description='Monte Python, a Monte Carlo code in Python')

###############
# MCMC basic
# -- number of steps	(OPTIONAL)
parser.add_argument('-N', metavar='steps',type=int,dest='N')
# -- output folder	(OBLIGATORY)
parser.add_argument('-o', metavar='output folder',type=str,dest='folder')
# -- parameter file	(OBLIGATORY)
parser.add_argument('-p', metavar='input param file',type=str,dest='param')
# -- covariance matrix	(OPTIONAL)
parser.add_argument('-c', metavar='input cov matrix',type=str,dest='cov')
# -- jumping method	(OPTIONAL)
parser.add_argument('-j', metavar='jumping method',type=str,dest='jumping',default='global')
# -- configuration file (OPTIONAL)
parser.add_argument('-conf', metavar='configuration file',type=str,dest='config_file',default='default.conf')

###############
# MCMC restart from chain
parser.add_argument('-r', metavar='restart from chain',type=str,dest='restart')

###############
# Information
# -- folder to analyze
parser.add_argument('-info', metavar='compute information of desired file',type=str,dest='files',nargs='*')
# -- number of bins (defaulting to 20)
parser.add_argument('-bins', metavar='desired number of bins, default is 20',type=int,dest='bins',default=20)
# -- possible comparison folder
parser.add_argument('-comp',metavar='comparison folder',type=str,dest='comp',nargs='*')
# -- if you just want the covariance matrix, use this option
parser.add_argument('-noplot',metavar='comparison folder',dest='plot',action='store_const',const=False,default=True)


def parse():
  # Recover all command line arguments in the args dictionnary
  args=parser.parse_args()

  # First of all, if the analyze module is invoked, there is no point in
  # checking for existing folder
  if args.files is None:
    
    # If the user wants to start over from an existing chain, the program will
    # use automatically the same folder, and the log.param in it
    if args.restart is not None:
      args.folder = args.restart.split('/')[0]+'/'
      args.param 	  = args.folder+'log.param'

    # Else, the user should provide an output folder
    else:
      if args.folder==None:
	print ' /|\   You must provide an output folder,\n/_o_\  because you do not want your main folder to look dirty, do you?'
	exit()

      # If he did so, 
      else:

	# check that the provided name is ending with a /,
	if args.folder[-1]!='/':
	  args.folder+='/'

	# and if the folder already exists, and that no parameter file was
	# provided, use the log.param
	if os.path.isdir(args.folder):
	  if args.param==None:
	    args.param=args.folder+'log.param'
	else:
	  if args.param==None:
	    print ' /|\   No log.param was found in your output folder,\n/_o_\  You must then provide a parameter file,\n       use the command line option -p any.param'
	    exit()

  return args
