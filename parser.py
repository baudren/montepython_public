import argparse

parser = argparse.ArgumentParser(description='Monte Python, a Monte Carlo code in Python')
parser.add_argument('-N', metavar='steps',type=int,dest='N')
parser.add_argument('-o', metavar='output folder',type=str,dest='folder')
parser.add_argument('-p', metavar='input param file',type=str,dest='par')
parser.add_argument('-c', metavar='input cov matrix',type=str,dest='cov')
parser.add_argument('-r', metavar='restart from chain',type=str,dest='restart')

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
	if args.par==None:
	  print ' /|\   You must provide a parameter file,\n/_o_\  use the command line option -p any.param'
	  exit()

  return args
