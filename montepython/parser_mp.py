"""
.. module:: parser_mp
    :synopsis: Definition of the command line options
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>

Defines the command line options and their help messages in
:func:`create_parser` and read the input command line in :func:`parse`, dealing
with different possible configurations.

"""
import os
import sys
import textwrap
import argparse  # Python module to handle command line arguments
import warnings

import io_mp


# -- custom Argument Parser that throws an io_mp.ConfigurationError
# -- for unified look within montepython
class MpArgumentParser(argparse.ArgumentParser):
    """Extension of the default ArgumentParser"""

    def error(self, message):
        """Override method to raise error
        Parameters
        ----------
        message: string
            error message
        """
        raise io_mp.ConfigurationError(message)

    def safe_parse_args(self, args=None):
        """
        Allows to set a default subparser

        This trick is there to maintain the previous way of calling
        MontePython.py
        """
        args = self.set_default_subparser('run', args)
        return self.parse_args(args)

    def set_default_subparser(self, default, args=None):
        """
        If no subparser option is found, add the default one

        .. note::

            This function relies on the fact that all calls to MontePython will
            start with a `-`. If this came to change, this function should be
            revisited

        """
        if not args:
            args = sys.argv[1:]
        if args[0] not in ['-h', '--help', '--version', '-info']:
            if args[0].find('-') != -1:
                args.insert(0, default)
        elif args[0] == '-info':
            args[0] = 'info'
        return args


# -- custom argparse types
# -- check that the argument is a positive integer
def positive_int(string):
    """
    Check if the input is integer positive
    Parameters
    ----------
    string: string
        string to parse

    output: int
        return the integer
    """
    try:
        value = int(string)
        if value <= 0:
            raise ValueError
        return value
    except ValueError:
        raise argparse.ArgumentTypeError(
            "You asked for a non-positive number of steps. "
            "I am not sure what to do, so I will exit. Sorry.")


# -- check that the argument is an existing file
def existing_file(fname):
    """
    Check if the file exists. If not raise an error
    Parameters
    ----------
    fname: string
        file name to parse

    output: fname
    """
    if os.path.isfile(fname):
        return fname
    else:
        msg = "The file '{}' does not exist".format(fname)
        raise argparse.ArgumentTypeError(msg)


def create_parser():
    """
    Definition of the parser command line options

    The main parser has so far two subparsers, corresponding to the two main
    modes of operating the code, namely `run` and `info`. If you simply call
    :code:`python montepython/MontePython.py -h`, you will find only this piece
    of information. To go further, and find the command line options specific
    to these two submodes, one should then do: :code:`python
    montepython/MontePython.py run -h`, or :code:`info -h`.

    All command line arguments are defined below, for each of the two
    subparsers. This function create the automatic help command.

    Each flag outputs the following argument to a destination variable,
    specified by the `dest` keyword argument in the source code. Please check
    there to understand the variable names associated with each option.

    **Options**:

        **run**:

            - **-N** (`int`) - number of steps in the chain (**OBL**). Note
              that when running on a cluster, your run might be stopped before
              reaching this number.
            - **-o** (`str`) - output folder (**OBL**). For instance :code:`-o
              chains/myexperiments/mymodel`. Note that in this example, the
              folder :code:`chains/myexperiments` must already exist.
            - **-p** (`str`) - input parameter file (**OBL**). For example
              :code:`-p input/exoticmodel.param`.
            - **-c** (`str`) - input covariance matrix (*OPT*). A covariance
              matrix is created when analyzing previous runs.

              .. note::
                    The list of parameters in the input covariance matrix and
                    in the run do not necessarily coincide.

            - **-j** (`str`) - jumping method (`global` (default),
              `sequential` or `fast`) (*OPT*).

              With the `global` method the code generates a new random
              direction at each step, with the `sequential` one it cycles over
              the eigenvectors of the proposal density (= input covariance
              matrix).

              The `global` method the acceptance rate is usually lower but
              the points in the chains are less correlated. We recommend using
              the sequential method to get started in difficult cases, when the
              proposal density is very bad, in order to accumulate points and
              generate a covariance matrix to be used later with the `default`
              jumping method.

              The `fast` method implements the Cholesky decomposition presented
              in http://arxiv.org/abs/1304.4473 by Antony Lewis.
            - **-m** (`string`) - sampling method used, by default 'MH' for
              Metropolis-Hastings, can be set to 'NS' for Nested Sampling
              (using Multinest wrapper PyMultiNest) or 'CH' for Cosmo Hammer
              (using the Cosmo Hammer wrapper to emcee algorithm).
            - **-f** (`float`) - jumping factor (>= 0, default to 2.4)
              (*OPT*).

              the proposal density is given by the input covariance matrix (or
              a diagonal matrix with elements given by the square of the input
              sigma's) multiplied by the square of this factor. In other words,
              a typical jump will have an amplitude given by sigma times this
              factor.

              The default is the famous factor 2.4, found by **TO CHECK**
              Dunkley et al. to be an optimal trade-off between high acceptance
              rate and high correlation of chain elements, at least for
              multivariate gaussian posterior probabilities. It can be a good
              idea to reduce this factor for very non-gaussian posteriors.

              Using :code:`-f 0 -N 1` is a convenient way to get the likelihood
              exactly at the starting point passed in input.
            - **--conf** (`str`) - configuration file (default to
              `default.conf`) (*OPT*). This file contains the path to your
              cosmological module directory.
            - **--chain_number** (`str`) - arbitrary numbering of the output
              chain, to overcome the automatic one (*OPT*).

              By default, the chains are named :code:`yyyy-mm-dd_N__i.txt` with
              year, month and day being extracted, :code:`N` being the number
              of steps, and :code:`i` an automatically updated index.

              This means that running several times the code with the same
              command will create different chains automatically.

              This option is a way to enforce a particular number :code:`i`.
              This can be useful when running on a cluster: for instance you
              may ask your script to use the job number as :code:`i`.
            - **-r** (`str`) - start a new chain from the last point of the
              given one, to avoid the burn-in stage (*OPT*).

              At the beginning of the run, the previous chain will be deleted,
              and its content transfered to the beginning of the new chain.
            - **-b** (`str`) - start a new chain from the bestfit computed in
              the given file (*OPT*)

        **info**:

              Replaces the old **-info** command, which is deprecated but still
              available.

              You can specify either single files, or a complete folder, for
              example :code:`-info chains/my-run/2012-10-26*`, or :code:`-info
              chains/my-run`
            - **--bins** (`int`) - number of bins in the histograms used to
              derive posterior probabilities and credible intervals (default to
              20). Decrease this number for smoother plots at the expense of
              masking details.
            - **--no_mean** (`None`) - by default, when plotting marginalised
              1D posteriors, the code also shows the mean likelihood per bin
              with dashed lines; this flag switches off the dashed lines
            - **--comp** (`str`) - pass the name of another folder (or another
              set of chains, same syntax as -info) if you want to compare 1D
              posteriors on the same plot.

              The lists of parameters in the two folders to compare do not need
              to coincide. It is limited so far to two folders to compare in
              total.
            - **--extra** (`str`) - extra file to customize the output plots.

            .. code::

                info.to_change={'oldname1':'newname1','oldname2':'newname2',...}
                info.to_plot=['name1','name2','newname3',...]
                info.new_scales={'name1':number1,'name2':number2,...}

            - **--noplot** (`None`) - do not produce plot, and compute only the
              covariance matrix (flag)
            - **--plot-2d** (`str`) - output triangle plot of 2d contours
              (`no`, `not_if_comp` (default), `always`, `overplot_comp`)
            - **--alpha** (`float`) - transparency of the second 2d posterior
              distribution in case of a comparison. Values accepted between 0
              and 1 (default to 0.8).
            - **--all** (`None`) - output every subplot in a separate file
              (flag)
            - **--ext** (`str`) - specify the extension of the figures (`pdf`
              (default), `png` (faster))
            - **--fontsize** (`int`) - adjust fontsize (default to 15)
            - **--ticksize** (`int`) - adjust ticksize (default to 13)

    """
    # Customized usage, for more verbosity concerning these subparsers options.
    usage = """%(prog)s [-h] [--version] {run,info} ... """
    usage += textwrap.dedent("""\n
        From more help, type:
        %(prog)s run -h
        %(prog)s info -h""")

    # parser = argparse.ArgumentParser(
    parser = MpArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Monte Python, a Monte Carlo code in Python',
        usage=usage)

    # -- version
    path_file = os.path.sep.join(
        os.path.abspath(__file__).split(os.path.sep)[:-2])
    with open(os.path.join(path_file, 'VERSION'), 'r') as version_file:
        version = version_file.readline()
        parser.add_argument('--version', action='version', version=version)

    # -- add the subparsers
    subparser = parser.add_subparsers(dest='subparser_name')

    ###############
    # run the MCMC
    runparser = subparser.add_parser('run', help="run the MCMC chains")

    # -- number of steps (OPTIONAL)
    runparser.add_argument('-N', help='number of steps',
                           type=positive_int, dest='N')
    # -- output folder (OBLIGATORY)
    runparser.add_argument('-o', '--output', help='output folder',
                           type=str, dest='folder')
    # -- parameter file (OBLIGATORY)
    runparser.add_argument('-p', '--param', help='input param file',
                           type=existing_file, dest='param')
    # -- covariance matrix (OPTIONAL)
    runparser.add_argument('-c', '--covmat', help='input cov matrix',
                           type=existing_file, dest='cov')
    # -- jumping method (OPTIONAL)
    runparser.add_argument('-j', '--jumping', help='jumping method',
                           dest='jumping', default='global',
                           choices=['global', 'sequential', 'fast'])
    # -- sampling method (OPTIONAL)
    runparser.add_argument('-m', '--method', help='sampling method',
                           dest='method', default='MH',
                           choices=['MH', 'NS', 'CH'])
    # -- jumping factor (OPTIONAL)
    runparser.add_argument('-f', help='jumping factor', type=float,
                           dest='jumping_factor', default=2.4)
    # -- configuration file (OPTIONAL)
    runparser.add_argument('--conf', help='configuration file',
                           type=str, dest='config_file',
                           default='default.conf')
    # -- arbitrary numbering of an output chain (OPTIONAL)
    runparser.add_argument('--chain-number', help='chain number')

    ###############
    # MCMC restart from chain or best fit file
    runparser.add_argument('-r', help='restart from chain',
                           type=existing_file, dest='restart')
    runparser.add_argument('-b', '--bestfit', dest='bf',
                           help='restart from best fit file',
                           type=existing_file)

    ###############
    # MultiNest arguments (all OPTIONAL and ignored if not "-m=NS")
    # The default values of -1 mean to take the PyMultiNest default values
    try:
        from nested_sampling import NS_prefix, NS_user_arguments
        NSparser = runparser.add_argument_group(
            title="MultiNest",
            description="Run the MCMC chains using MultiNest"
            )
        for arg in NS_user_arguments:
            NSparser.add_argument('--'+NS_prefix+arg,
                                  default=-1,
                                  **NS_user_arguments[arg])
    except ImportError:
        # Not defined if not installed
        pass

    ###############
    # CosmoHammer arguments (all OPTIONAL and ignored if not "-m=CH")
    # The default values of -1 mean to take the CosmoHammer default values
    try:
        from cosmo_hammer import CH_prefix, CH_user_arguments
        CHparser = runparser.add_argument_group(
            title="CosmoHammer",
            description="Run the MCMC chains using the CosmoHammer framework")
        for arg in CH_user_arguments:
            CHparser.add_argument('--'+CH_prefix+arg,
                                  default=-1,
                                  **CH_user_arguments[arg])
    except ImportError:
        # Not defined if not installed
        pass

    ###############
    # Information
    infoparser = subparser.add_parser('info',
                                      help="analyze the MCMC chains")

    # -- folder to analyze
    infoparser.add_argument('files',
                            help='compute information of desired file',
                            nargs='+')
    # -- number of bins (defaulting to 20)
    infoparser.add_argument('--bins',
                            help='desired number of bins, default is 20',
                            type=int, default=20)
    # -- to remove the mean-likelihood line
    infoparser.add_argument('--no-mean',
                            help='remove the mean likelihood plot',
                            dest='mean_likelihood', action='store_false',)
    # -- possible comparison folder
    infoparser.add_argument('--comp', help='comparison folder')
    # -- possible plot file describing custom commands
    infoparser.add_argument('--extra', help='plot file for custom needs',
                            dest='optional_plot_file')
    # -- if you just want the covariance matrix, use this option
    infoparser.add_argument('--noplot', help='omit the plotting part',
                            dest='plot', action='store_false')
    # -- if you want to output 2d contours plots (the 'triangle' plot)
    # default: only as long as -comp is not specified
    infoparser.add_argument('--plot-2d',
                            help='plot the triangle plot of 2d contours',
                            dest='plot_2d', type=str,
                            choices=['no', 'not_if_comp',
                                     'always', 'overplot_comp'],
                            default='not_if_comp')
    # -- when comparing two folders, decide on the alpha setting of the second
    # plot. Defaults to 0.8
    infoparser.add_argument('--alpha',
                            help='choose the transparency of the compared run',
                            dest='alpha', type=float,
                            default=0.8)
    # -- if you want to output every single subplots
    infoparser.add_argument(
        '--all', help='plot every single subplot in a separate pdf file',
        dest='subplot', action='store_true')
    # -- to change the extension used to output files (pdf is the default one,
    # but takes long, valid options are png and eps)
    infoparser.add_argument(
        '--ext', help='''change extension for the output file.
        Any extension handled by `matplotlib` can be used''',
        type=str, dest='extension', default='pdf')
    # -- fontsize of plots (defaulting to 15)
    infoparser.add_argument('--fontsize', help='desired font size',
                            type=int, default=15)
    # -- ticksize of plots (defaulting to 13)
    infoparser.add_argument('--ticksize', help='desired tick size',
                            type=int, default=13)

    return parser


def parse(custom_command=''):
    """
    Check some basic organization of the folder, and exit the program in case
    something goes wrong.

    Keyword Arguments
    -----------------
    custom_command : str
        For testing purposes, instead of reading the command line argument,
        read instead the given string. It should ommit the start of the
        command, so e.g.: '-N 10 -o toto/'

    """
    # Create the parser
    parser = create_parser()

    # Recover all command line arguments in the args dictionary, except for a
    # test, where the custom_command string is read.
    # Note that the function safe_parse_args is read instead of parse_args. It
    # is a function defined in this file to allow for a default subparser.
    if not custom_command:
        args = parser.safe_parse_args()
    else:
        args = parser.safe_parse_args(custom_command.split(' '))

    # Some check to perform when running the MCMC chains is requested
    if args.subparser_name == "run":

        # If the user wants to start over from an existing chain, the program
        # will use automatically the same folder, and the log.param in it
        if args.restart is not None:
            args.folder = os.path.sep.join(
                args.restart.split(os.path.sep)[:-1])
            args.param = os.path.join(args.folder, 'log.param')
            warnings.warn(
                "Restarting from %s." % args.restart +
                " Using associated log.param.")

        # Else, the user should provide an output folder
        else:
            if args.folder is None:
                raise io_mp.ConfigurationError(
                    "You must provide an output folder, because you do not " +
                    "want your main folder to look dirty, do you ?")

            # and if the folder already exists, and that no parameter file was
            # provided, use the log.param
            if os.path.isdir(args.folder):
                if os.path.exists(
                        os.path.join(args.folder, 'log.param')):
                    # if the log.param exists, and that a parameter file was
                    # provided, take instead the log.param, and notify the
                    # user.
                    old_param = args.param
                    args.param = os.path.join(
                        args.folder, 'log.param')
                    if old_param is not None:
                        warnings.warn(
                            "Appending to an existing folder: using the "
                            "log.param instead of %s" % old_param)
                else:
                    if args.param is None:
                        raise io_mp.ConfigurationError(
                            "The requested output folder seems empty. "
                            "You must then provide a parameter file (command"
                            " line option -p any.param)")
            else:
                if args.param is None:
                    raise io_mp.ConfigurationError(
                        "The requested output folder appears to be non "
                        "existent. You must then provide a parameter file "
                        "(command line option -p any.param)")

    return args
