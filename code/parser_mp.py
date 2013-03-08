import os
import argparse  # Python module to handle command line arguments

# Definition of the object, below will be added all the possible options. This
# will also create an automatic help command, available through the option -h
parser = argparse.ArgumentParser(
    description='Monte Python, a Monte Carlo code in Python')

###############
# MCMC basic
# -- number of steps	(OPTIONAL)
parser.add_argument('-N', metavar='steps', type=int, dest='N')
# -- output folder	(OBLIGATORY)
parser.add_argument('-o', metavar='output folder', type=str, dest='folder')
# -- parameter file	(OBLIGATORY)
parser.add_argument('-p', metavar='input param file', type=str, dest='param')
# -- covariance matrix	(OPTIONAL)
parser.add_argument('-c', metavar='input cov matrix', type=str, dest='cov')
# -- jumping method	(OPTIONAL)
parser.add_argument('-j', metavar='jumping method', type=str,
                    dest='jumping', default='global')
# -- jumping factor	(OPTIONAL)
parser.add_argument('-f', metavar='jumping factor', type=float,
                    dest='jumping_factor', default=2.4)
# -- configuration file (OPTIONAL)
parser.add_argument('-conf', metavar='configuration file', type=str,
                    dest='config_file', default='default.conf')
# -- arbitraty numbering of an output chain (OPTIONAL)
parser.add_argument('-chain_number', metavar='chain number', type=str,
                    dest='chain_number', default=None)

###############
# MCMC restart from chain or best fit file
parser.add_argument('-r', metavar='restart from chain', type=str,
                    dest='restart')
parser.add_argument('-bf', metavar='restart from best fit file', type=str,
                    dest='bf')

###############
# Information
# -- folder to analyze
parser.add_argument('-info', metavar='compute information of desired file',
                    type=str, dest='files', nargs='*')
# -- number of bins (defaulting to 20)
parser.add_argument('-bins', metavar='desired number of bins, default is 20',
                    type=int, dest='bins', default=20)
# -- to remove the mean-likelihood line
parser.add_argument('-no_mean', metavar='remove the mean likelihood plot',
                    dest='mean_likelihood', action='store_const',
                    const=False, default=True)
# -- possible comparison folder
parser.add_argument('-comp', metavar='comparison folder', type=str,
                    dest='comp', nargs=1)
# -- possible plot file describing custom commands
parser.add_argument('-extra', metavar='plot file for custom needs',
                    type=str, dest='optional_plot_file', nargs=1)
# -- if you just want the covariance matrix, use this option
parser.add_argument('-noplot', metavar='ommit the plotting part',
                    dest='plot', action='store_const',
                    const=False, default=True)
# -- if you want to output every single subplots
parser.add_argument(
    '-all', metavar='plot every single subplot in a separate pdf file',
    dest='subplot', action='store_const', const=True, default=False)
# -- to change the extension used to output files (pdf is the default one, but
# takes long, valid options are png and eps)
parser.add_argument('-ext', metavar='change extension for the output file',
                    type=str, dest='extension', default='pdf')
# -- fontsize of plots (defaulting to 15)
parser.add_argument('-fontsize', metavar='desired fontsize, default is 15',
                    type=int, dest='fontsize', default=15)
# -- ticksize of plots (defaulting to 13)
parser.add_argument('-ticksize', metavar='desired ticksize, default is 13',
                    type=int, dest='ticksize', default=13)


def parse():
    # Recover all command line arguments in the args dictionnary
    args = parser.parse_args()

    # First of all, if the analyze module is invoked, there is no point in
    # checking for existing folder
    if args.files is None:

        # If the user wants to start over from an existing chain, the program
        # will use automatically the same folder, and the log.param in it
        if args.restart is not None:
            args.folder = args.restart.split('/')[0]+'/'
            args.param = args.folder+'log.param'

        # Else, the user should provide an output folder
        else:
            if args.folder is None:
                print ' /|\   You must provide an output folder,'
                print '/_o_\  because you do not want your main folder '
                print '       to look dirty, do you?'
                exit()

            # If he did so,
            else:
                # check that the provided name is ending with a /,
                if args.folder[-1] != '/':
                    args.folder += '/'

            # and if the folder already exists, and that no parameter file was
            # provided, use the log.param
            if os.path.isdir(args.folder):
                if os.path.exists(args.folder+'log.param'):
                    old_param = args.param
                    args.param = args.folder+'log.param'
                    if args.param is not None:
                        print ' /!\ Appending to an existing folder: '
                        print '     using the log.param instead of %s' % (
                            old_param)
                else:
                    if args.param is None:
                        print ' /|\  The requested output folder appears to be'
                        print '/_o_\ empty. You must then provide a parameter '
                        print '      file (command line option -p any.param)'
                        exit()
            else:
                if args.param is None:
                    print ' /|\  The requested output folder appears to be '
                    print '/_o_\ non-existent, you must then provide a param'
                    print '      file (command line option -p any.param)'
                    exit()

    return args
