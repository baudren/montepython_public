import argparse
import os

parser = argparse.ArgumentParser(description='Monte Python, a Monte Carlo code in Python')
parser.add_argument('-N', metavar='steps',type=int,dest='N')
parser.add_argument('-o', metavar='output folder',type=str,dest='folder')
parser.add_argument('-p', metavar='input param file',type=str,dest='par')
parser.add_argument('-c', metavar='input cov matrix',type=str,dest='cov')
parser.add_argument('-r', metavar='restart from chain',type=str,dest='restart')
parser.add_argument('-j', metavar='jumping method',type=str,dest='jumping',default='global')

# cleaning
parser.add_argument('-clean', metavar='input folder to clean',type=str,dest='clean')

# for information part
parser.add_argument('-info', metavar='compute information of desired file',type=str,dest='files',nargs='*')
parser.add_argument('-bins', metavar='desired number of bins, default is 10',type=int,dest='bins',default=10)

parser.add_argument('-input',metavar='input file',type=str,dest='infile')


def parse():
  args=parser.parse_args()
  if args.clean is None:
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
	      print ' /|\   No runs were found in your output folder,\n/_o_\  You must provide a parameter file,\n       use the command line option -p any.param'
	      exit()
    else:
      args.folder = args.restart.split('/')[0]+'/'
      args.par 	  = args.folder+'log.param'

  return args
