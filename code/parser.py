import argparse
import os

parser = argparse.ArgumentParser(description='Monte Python, a Monte Carlo code in Python')

###############
# MCMC basic
# -- number of steps
parser.add_argument('-N', metavar='steps',type=int,dest='N')
# -- output folder
parser.add_argument('-o', metavar='output folder',type=str,dest='folder')
# -- parameter file
parser.add_argument('-p', metavar='input param file',type=str,dest='par')
# -- covariance matrix
parser.add_argument('-c', metavar='input cov matrix',type=str,dest='cov')
# -- jumping method
parser.add_argument('-j', metavar='jumping method',type=str,dest='jumping',default='global')

###############
# MCMC restart from chain
parser.add_argument('-r', metavar='restart from chain',type=str,dest='restart')

###############
# information
# -- folder to analyze
parser.add_argument('-info', metavar='compute information of desired file',type=str,dest='files',nargs='*')
# -- number of bins (defaulting to 20)
parser.add_argument('-bins', metavar='desired number of bins, default is 20',type=int,dest='bins',default=20)
# -- possible comparison folder
parser.add_argument('-comp',metavar='comparison folder',type=str,dest='comp',nargs='*')


def parse():
  args=parser.parse_args()
  if args.restart is None:
    if args.files is None:
      if args.folder==None:
	print ' /|\   You must provide an output folder,\n/_o_\  because you do not want your main folder to look dirty, do you?'
	exit()
      else:
	if args.folder[-1]!='/':
	  args.folder+='/'
	if os.path.isdir(args.folder):
	  if args.par==None:
	    args.par=args.folder+'log.param'
	else:
	  if args.par==None:
	    print ' /|\   No log.param was found in your output folder,\n/_o_\  You must then provide a parameter file,\n       use the command line option -p any.param'
	    exit()
  else:
    args.folder = args.restart.split('/')[0]+'/'
    args.par 	  = args.folder+'log.param'

  return args
