"""
.. module:: parser_np
    :synopsis: Definition of the command line options
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>

Defines the command line options and their help messages in
:func:`create_parser` and read the input command line in :func:`parse`, dealing
with different possible configurations.

"""
import os
import argparse  # Python module to handle command line arguments
import io_mp


def create_parser():
    """
    Definition of the parser command line options

    All command line arguments are defined below. This will also create an
    automatic help command, available through the call :code:`python
    code/Montepython -h`, listing the parameters and their function.

    Each flag outputs the following argument to a destination variable,
    specified by the `dest` keyword argument in the source code. Please check
    there to understand the variable names associated with each option.

    **Options**:

        **MCMC**:

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
                    The list of parameters in the input covariance matrix and in
                    the run do not necessarily coincide.

            - **-j** (`str`) - jumping method (`global` (default),
              `sequential`) (*OPT*).

              With the `global` method the code generates a new random direction
              at each step, with the `sequential` one it cycles over the
              eigenvectors of the proposal density (= input covariance matrix).

              The `global` method the acceptance rate is usually lower but
              the points in the chains are less correlated. We recommend using
              the sequential method to get started in difficult cases, when the
              proposal density is very bad, in order to accumulate points and
              generate a covariance matrix to be used later with the `default`
              jumping method.
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
            - **-conf** (`str`) - configuration file (default to
              `default.conf`) (*OPT*). This file contains the path to your
              cosmological module directory.
            - **-chain_number** (`str`) - arbitrary numbering of the output
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
            - **-bf** (`str`) - start a new chain from the bestfit computed in
              the given file (*OPT*)

        **Information**:

            - **-info** (`str`) - compute the information of given file/folder.

              You can specify either single files, or a complete folder, for
              example :code:`-info chains/my-run/2012-10-26*`, or :code:`-info
              chains/my-run`
            - **-bins** (`int`) - number of bins in the histograms used to
              derive posterior probabilities and credible intervals (default to
              20). Decrease this number for smoother plots at the expense of
              masking details.
            - **-no_mean** (`None`) - by default, when plotting marginalised 1D
              posteriors, the code also shows the mean likelihood per bin with
              dashed lines; this flag switches off the dashed lines
            - **-comp** (`str`) - pass the name of another folder (or another
              set of chains, same syntax as -info) if you want to compare 1D
              posteriors on the same plot.

              The lists of parameters in the two folders to compare do not need
              to coincide. It is limited so far to two folders to compare in
              total.
            - **-extra** (`str`) - extra file to customize the output plots.

            .. code::

                info.to_change={'oldname1':'newname1','oldname2':'newname2',...}
                info.to_plot=['name1','name2','newname3',...]
                info.new_scales={'name1':number1,'name2':number2,...}

            - **-noplot** (`None`) - do not produce plot, and compute only the
              covariance matrix (flag)
            - **-all** (`None`) - output every subplot in a separate file
              (flag)
            - **-ext** (`str`) - specify the extension of the figures (`pdf`
              (default), `png` (faster))
            - **-fontsize** (`int`) - adjust fontsize (default to 15)
            - **-ticksize** (`int`) - adjust ticksize (default to 13)

    """
    parser = argparse.ArgumentParser(
        description='Monte Python, a Monte Carlo code in Python')

    # -- number of steps (OPTIONAL)
    parser.add_argument('-N', metavar='number of steps', type=int, dest='N')
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

    return parser


def parse():
    """
    Check some basic organization of the folder, and exit the program in case
    something goes wrong.

    """
    # Create the parser
    parser = create_parser()

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
                io_mp.message(
                    "You must provide an output folder, because you do not \
                    want your main folder to look dirty, do you ?",
                    "error")

            # If he did so,
            else:
                # check that the provided name is ending with a /,
                if args.folder[-1] != '/':
                    args.folder += '/'

            # and if the folder already exists, and that no parameter file was
            # provided, use the log.param
            if os.path.isdir(args.folder):
                if os.path.exists(args.folder+'log.param'):
                    # if the log.param exists, and that a parameter file was
                    # provided, take instead the log.param, and notify the
                    # user.
                    old_param = args.param
                    args.param = args.folder+'log.param'
                    if args.param is not None:
                        io_mp.message(
                            "Appending to an existing folder: using the \
                            log.param instead of %s" % old_param,
                            "info")
                else:
                    if args.param is None:
                        io_mp.message(
                            "The requested output folder appears to be empty. \
                            You must then provide a parameter file (command \
                            line option -p any.param)",
                            "error")
            else:
                if args.param is None:
                    io_mp.message(
                        "The requested output folder appears to be non \
                        existent. You must then provide a parameter file \
                        (command line option -p any.param)",
                        "error")

    return args
