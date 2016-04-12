"""
.. module:: parser_mp
    :synopsis: Definition of the command line options
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>
.. moduleauthor:: Francesco Montesano <franz.bergesund@gmail.com>

Defines the command line options and their help messages in
:func:`create_parser` and read the input command line in :func:`parse`, dealing
with different possible configurations.

The fancy short/long help formatting, as well as the automatic help creation
from docstrings is entirely due to Francesco Montesano.

"""
import os
import sys
import textwrap as tw
import re
import argparse as ap  # Python module to handle command line arguments
import warnings

import io_mp


# -- custom Argument Parser that throws an io_mp.ConfigurationError
# -- for unified look within montepython
class MpArgumentParser(ap.ArgumentParser):
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
                msg = "Defaulting to the 'run' command. Please update the"
                msg += " call of MontePython. For more info, see the help"
                msg += " string and/or the documentation "
                warnings.warn(msg)
                args.insert(0, default)
        elif args[0] == '-info':
            msg = "The info option has been turned into a command. "
            msg += "Please substitute '-info' with 'info' when running "
            msg += "MontePython"
            warnings.warn(msg)
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
        raise ap.ArgumentTypeError(
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

    Returns
    -------
    fname : string
    """
    if os.path.isfile(fname):
        return fname
    else:
        msg = "The file '{}' does not exist".format(fname)
        raise ap.ArgumentTypeError(msg)


def parse_docstring(docstring, key_symbol="<**>", description_symbol="<++>"):
    """
    Extract from the docstring the keys and description, return it as a dict

    Parameters
    ----------
    docstring : str
    key_symbol : str
        identifies the key of an argument/option
    description_symbol : str
        identify the description of an argument/option

    output
    ------
    helpdic : dict
        help strings for the parser
    """
    # remove new lines and multiple whitespaces
    whitespaces = re.compile(r"\s+")
    docstring = whitespaces.sub(" ", docstring)

    # escape special characters
    key_symbol = re.escape(key_symbol)
    description_symbol = re.escape(description_symbol)

    # define the regular expressions to match the key and the description
    key_match = r'{0}-{{0,2}}(.+?){0}'
    re_key = re.compile(key_match.format(key_symbol))

    desc_match = r'({0}.+?{0}.+?{0})'
    re_desc = re.compile(desc_match.format(description_symbol))

    # get all and check that the keys and descriptions have the same lenghts
    keys = re_key.findall(docstring)
    descriptions = re_desc.findall(docstring)
    if len(keys) != len(descriptions):
        msg = "The option keys and their descriptions have different lenghts.\n"
        msg += "Make sure that there are as many string surrounded by '{0}'"
        msg += " as there are surrounded by '{1}"
        raise ValueError(msg.format(key_symbol, description_symbol))

    helpdict = dict(zip(keys, descriptions))
    return helpdict


def custom_help(split_string="<++>"):
    """
    Create a custom help action.

    It expects *split_string* to appear in groups of three.
    If the option string is '-h', then uses the short description
    between the first two *split_string*.
    If the option string is '-h', then uses all that is between
    the first and the third *split_string*, stripping the first one.

    Parameters
    ----------
    split_string: str
        string to use to select the help string and how to select them.
        They must appear in groups of *3*

    output
    ------
    CustomHelp: class definition
    """
    class CustomHelp(ap._HelpAction):
        def __call__(self, parser, namespace, values, option_string=None):

            # create the help string and store it into a string
            from StringIO import StringIO
            fstr = StringIO()
            try:
                parser.print_help(file=fstr)
                help_str = fstr.getvalue()
            finally:
                fstr.close()

            # create the regular expression to match the desciption
            descmatch = r'{0}(.+?){0}(.+?){0}'
            # escape possible dangerous characters
            esplit_string = re.escape(split_string)
            re_desc = re.compile(descmatch.format(esplit_string),
                                 flags=re.DOTALL)

            # select the case according to which option_string is selected
            if option_string == '-h':
                to_sub = r'\1'
            elif option_string == '--help':
                to_sub = r'\1\2'

            print(re_desc.sub(to_sub, help_str))
            parser.exit()

    return CustomHelp


def add_subparser(sp, name, **kwargs):
    """
    Add a parser to the subparser *sp* with *name*.

    All the logic common to all subparsers should go here

    Parameters
    ----------
    sp: subparser instance
    name: str
        name of the subparser
    kwargs: dict
        keywords to pass to the subparser

    output
    ------
    sparser: Argparse instance
        new subparser
    """
    kwargs["add_help"] = False
    kwargs['formatter_class'] = ap.ArgumentDefaultsHelpFormatter
    sparser = sp.add_parser(name, **kwargs)

    sparser.add_argument("-h", "--help", action=custom_help(),
                         help="print the short or long help")

    return sparser


def get_dict_from_docstring(key_symbol="<**>", description_symbol="<++>"):
    """
    Create the decorator

    Parameters
    ----------
    key_symbol : str
        identifies the key of a argument/option
    description_symbol: str
        identify the description of a argument/option

    Returns
    ------
    wrapper: function
    """
    def wrapper(func):
        """
        Decorator that wraps the function that implement the parser, parses the
        `__doc__` and construct a dictionary with the help strings.  The
        dictionary is added as an attribute of `func` and can be accessed in
        the function

        Parameters
        ----------
        func: function
            function with the docs to be parsed

        Returns
        ------
        func: function
            function with the dictionary added. *key_symbol* and
            *description_symbol* strings are removed
        """
        docstring = func.__doc__
        helpdict = parse_docstring(
            docstring, key_symbol=key_symbol,
            description_symbol=description_symbol)
        func.helpdict = helpdict
        # remove markers
        docstring = docstring.replace(key_symbol, '')
        func.__doc__ = docstring.replace(description_symbol, '')
        return func
    return wrapper


def initialise_parser(**kwargs):
    """
    Create the argument parser and returns it
    Parameters
    ----------
    kwargs: dictionary
        keyword to pass to the parser
    output
    ------
    p: MpArgumentParser instance
        parser with some keyword added
    """
    kwargs['formatter_class'] = ap.ArgumentDefaultsHelpFormatter
    p = MpArgumentParser(**kwargs)

    # -- version
    path_file = os.path.sep.join(
        os.path.abspath(__file__).split(os.path.sep)[:-2])
    with open(os.path.join(path_file, 'VERSION'), 'r') as version_file:
        version = version_file.readline()
        p.add_argument('--version', action='version', version=version)

    p.add_argument("-v", "--verbose", action="store_true", help="Verbose mode")

    return p


@get_dict_from_docstring()
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

    Options
    -------

    **run**

        <**>-N<**> : int
            <++>number of steps in the chain<++> (**OBL**). Note that when
            running on a cluster, your run might be stopped before reaching
            this number.<++>
        <**>-o<**> : str
            <++>output folder<++> (**OBL**). For instance :code:`-o
            chains/myexperiments/mymodel`. Note that in this example, the
            folder :code:`chains/myexperiments` must already exist.<++>
        <**>-p<**> : str
            <++>input parameter file<++> (**OBL**). For example :code:`-p
            input/exoticmodel.param`.<++>
        <**>-c<**> : str
            <++>input covariance matrix<++> (*OPT*). A covariance matrix is
            created when analyzing previous runs.

            Note that the list of parameters in the input covariance matrix and
            in the run do not necessarily coincide.<++>
        <**>-j<**> : str
            <++>jumping method<++> (`global` (default), `sequential` or `fast`)
            (*OPT*).

            With the `global` method the code generates a new random direction
            at each step, with the `sequential` one it cycles over the
            eigenvectors of the proposal density (= input covariance matrix).

            The `global` method the acceptance rate is usually lower but the
            points in the chains are less correlated. We recommend using the
            sequential method to get started in difficult cases, when the
            proposal density is very bad, in order to accumulate points and
            generate a covariance matrix to be used later with the `default`
            jumping method.

            The `fast` method implements the Cholesky decomposition presented
            in http://arxiv.org/abs/1304.4473 by Antony Lewis.<++>
        <**>-m<**> : str
            <++>sampling method<++>, by default 'MH' for Metropolis-Hastings,
            can be set to 'NS' for Nested Sampling (using Multinest wrapper
            PyMultiNest), 'CH' for Cosmo Hammer (using the Cosmo Hammer wrapper
            to emcee algorithm), and finally 'IS' for importance sampling.

            Note that when running with Importance sampling, you need to
            specify a folder to start from.<++>
        <**>--update<**> : int
            <++>update frequency for Metropolis Hastings.<++>
            If greater than zero, number of steps after which the proposal covariance
            matrix is updated automatically (*OPT*).<++>
        <**>-f<**> : float
            <++>jumping factor<++> (>= 0, default to 2.4) (*OPT*).

            The proposal density is given by the input covariance matrix (or a
            diagonal matrix with elements given by the square of the input
            sigma's) multiplied by the square of this factor. In other words, a
            typical jump will have an amplitude given by sigma times this
            factor.

            The default is the famous factor 2.4, advertised by Dunkley
            et al. to be an optimal trade-off between high acceptance rate and
            high correlation of chain elements, at least for multivariate
            gaussian posterior probabilities. It can be a good idea to reduce
            this factor for very non-gaussian posteriors.

            Using :code:`-f 0 -N 1` is a convenient way to get the likelihood
            exactly at the starting point passed in input.<++>
        <**>--conf<**> : str
            <++>configuration file<++> (default to `default.conf`) (*OPT*).
            This file contains the path to your cosmological module
            directory.<++>
        <**>--chain-number<**> : str
            arbitrary <++>number of the output chain<++>, to overcome the
            automatic one (*OPT*).

            By default, the chains are named :code:`yyyy-mm-dd_N__i.txt` with
            year, month and day being extracted, :code:`N` being the number of
            steps, and :code:`i` an automatically updated index.

            This means that running several times the code with the same
            command will create different chains automatically.

            This option is a way to enforce a particular number :code:`i`.
            This can be useful when running on a cluster: for instance you may
            ask your script to use the job number as :code:`i`.<++>
        <**>-r<**> : str
            <++>restart from last point in chain<++>, to avoid the burn-in
            stage (*OPT*).

            At the beginning of the run, the previous chain will be deleted,
            and its content transfered to the beginning of the new chain.<++>
        <**>-b<**> : str
            <++>start a new chain from the bestfit file<++> computed with
            analyze.  (*OPT*)<++>
        <**>--fisher<**> : None
            <++>Calculates the inverse of the fisher matrix<++> to use as
            proposal distribution<++>
        <**>--silent<**> : None
            <++>silence the standard output<++> (useful when running on
            clusters)<++>
        <**>--Der-target-folder<**> : str
            <++>Add additional derived params to this folder<++>. It has to be
            used in conjunction with `Der-param-list`, and the method set to
            Der: :code:`-m Der`. (*OPT*)<++>
        <**>--Der-param-list<**> : str
            <++>Specify a number of derived parameters to be added<++>. A
            complete example would be to add Omega_Lambda as a derived
            parameter:
            :code:`python montepython/MontePython.py run -o existing_folder
            -m Der --Der-target-folder non_existing_folder --Der-param-list
            Omega_Lambda`<++>
        <**>--IS-starting-folder<**> : str
            <++>Perform Importance Sampling from this folder or set of
            chains<++> (*OPT*)<++>
        <**>--stop-after-update<**> : bool
            <++>When using update mode, stop run after updating the covariant matrix.<++>
            Useful if you want to change settings after the first guess (*OPT*) (flag)<++>

        For Nested Sampling and Cosmo Hammer arguments, see
        :mod:`nested_sampling` and :mod:`cosmo_hammer`.

    **info**

              Replaces the old **-info** command, which is deprecated but still
              available.

        <**>files<**> : string/list of strings
            <++>you can specify either single files, or a complete folder<++>,
            for example :code:`info chains/my-run/2012-10-26*`, or :code:`info
            chains/my-run`.

            If you specify several folders (or set of files), a comparison
            will be performed.<++>
        <**>--minimal<**> : None
            <++>use this flag to avoid computing the posterior
            distribution.<++> This will decrease the time needed for the
            analysis, especially when analyzing big folders.<++>
        <**>--bins<**> : int
            <++>number of bins in the histograms<++> used to derive posterior
            probabilities and credible intervals (default to 20). Decrease this
            number for smoother plots at the expense of masking details.<++>
        <**>--no-mean<**> : None
            <++>remove the mean likelihood from the plot<++>. By default, when
            plotting marginalised 1D posteriors, the code also shows the mean
            likelihood per bin with dashed lines; this flag switches off the
            dashed lines.<++>
        <**>--extra<**> : str
            <++>extra file to customize the output plots<++>. You can actually
            set all the possible options in this file, including line-width,
            ticknumber, ticksize, etc... You can specify four fields,
            `info.redefine` (dict with keys set to the previous variable, and
            the value set to a numerical computation that should replace this
            variable), `info.to_change` (dict with keys set to the old variable
            name, and value set to the new variable name), `info.to_plot` (list
            of variables with new names to plot), and `info.new_scales` (dict
            with keys set to the new variable names, and values set to the
            number by which it should be multiplied in the graph).<++> For
            instance,

            .. code::

                info.to_change={'oldname1':'newname1','oldname2':'newname2',...}
                info.to_plot=['name1','name2','newname3',...]
                info.new_scales={'name1':number1,'name2':number2,...}

        <**>--noplot<**> : bool
            <++>do not produce any plot, simply compute the posterior<++>
            (*OPT*) (flag)<++>
        <**>--noplot-2d<**> : bool
            <++>produce only the 1d posterior plot<++> (*OPT*) (flag)<++>
        <**>--contours-only<**> : bool
            <++>do not fill the contours on the 2d plots<++> (*OPT*) (flag)<++>
        <**>--all<**> : None
            <++>output every subplot and data in separate files<++> (*OPT*)
            (flag)<++>
        <**>--ext<**> : str
            <++>change the extension for the output file. Any extension handled
            by :code:`matplotlib` can be used<++>. (`pdf` (default), `png`
            (faster))<++>
        <**>--fontsize<**> : int
            <++>desired fontsize<++> (default to 16)<++>
        <**>--ticksize<**> : int
            <++>desired ticksize<++> (default to 14)<++>
        <**>--line-width<**> : int
            <++>set line width<++> (default to 4)<++>
        <**>--decimal<**> : int
            <++>number of decimal places on ticks<++> (default to 3)<++>
        <**>--ticknumber<**> : int
            <++>number of ticks on each axis<++> (default to 3)<++>
        <**>--legend-style<**> : str
            <++>specify the style of the legend<++>, to choose from `sides` or
            `top`.<++>
        <**>--keep-non-markovian<**> : bool
            <++>Use this flag to keep the non-markovian part of the chains produced
            at the beginning of runs with --update mode<++>
            This option is only relevant when the chains were produced with --update (*OPT*) (flag)<++>
        <**>--keep-fraction<**> : float
            <++>after burn-in removal, analyze only last fraction of each chain.<++>
            (between 0 and 1). Normally one would not use this for runs with --update mode,
            unless --keep-non-markovian is switched on (*OPT*)<++>
        <**>--want-covmat<**> : bool
            <++>calculate the covariant matrix when analyzing the chains.<++>
            Warning: this will interfere with ongoing runs utilizing update mode (*OPT*) (flag)<++>

    Returns
    -------
    args : NameSpace
        parsed input arguments

    """
    helpdict = create_parser.helpdict
    # Customized usage, for more verbosity concerning these subparsers options.
    usage = """%(prog)s [-h] [--version] {run,info} ... """
    usage += tw.dedent("""\n
        From more help on each of the subcommands, type:
        %(prog)s run -h
        %(prog)s info -h\n\n""")

    # parser = ap.ArgumentParser(
    #parser = MpArgumentParser(
        #formatter_class=ap.ArgumentDefaultsHelpFormatter,
        #description='Monte Python, a Monte Carlo code in Python',
        #usage=usage)
    parser = initialise_parser(
        description='Monte Python, a Monte Carlo code in Python', usage=usage)

    # -- add the subparsers
    subparser = parser.add_subparsers(dest='subparser_name')

    ###############
    # run the MCMC
    runparser = add_subparser(subparser, 'run', help="run the MCMC chains")

    # -- number of steps (OPTIONAL)
    runparser.add_argument('-N', help=helpdict['N'], type=positive_int,
                           dest='N')
    # -- output folder (OBLIGATORY)
    runparser.add_argument('-o', '--output', help=helpdict['o'], type=str,
                           dest='folder')
    # -- parameter file (OBLIGATORY)
    runparser.add_argument('-p', '--param', help=helpdict['p'],
                           type=existing_file, dest='param')
    # -- covariance matrix (OPTIONAL)
    runparser.add_argument('-c', '--covmat', help=helpdict['c'],
                           type=existing_file, dest='cov')
    # -- jumping method (OPTIONAL)
    runparser.add_argument('-j', '--jumping', help=helpdict['j'],
                           dest='jumping', default='fast',
                           choices=['global', 'sequential', 'fast'])
    # -- sampling method (OPTIONAL)
    runparser.add_argument('-m', '--method', help=helpdict['m'],
                           dest='method', default='MH',
                           choices=['MH', 'NS', 'CH', 'IS', 'Der'])
    # -- update Metropolis Hastings (OPTIONAL)
    runparser.add_argument('--update', help=helpdict['update'], type=int,
                           default=0)
    # -- jumping factor (OPTIONAL)
    runparser.add_argument('-f', help=helpdict['f'], type=float,
                           dest='jumping_factor', default=2.4)
    # -- fisher (EXPERIMENTAL)
    runparser.add_argument('--fisher', help=helpdict['fisher'],
                           action='store_true')
    # -- configuration file (OPTIONAL)
    runparser.add_argument('--conf', help=helpdict['conf'],
                           type=str, dest='config_file',
                           default='default.conf')
    # -- arbitrary numbering of an output chain (OPTIONAL)
    runparser.add_argument('--chain-number', help=helpdict['chain-number'])
    # -- stop run after first successful update using --update (EXPERIMENTAL)
    runparser.add_argument('--stop-after-update', help=helpdict['stop-after-update'],
                           dest='stop_after_update', action='store_true')

    ###############
    # MCMC restart from chain or best fit file
    runparser.add_argument('-r', help=helpdict['r'],
                           type=existing_file, dest='restart')
    runparser.add_argument('-b', '--bestfit', dest='bf', help=helpdict['b'],
                           type=existing_file)

    ###############
    # Silence the output (no print on the console)
    runparser.add_argument('--silent', help=helpdict['silent'],
                           action='store_true')
    ###############
    # Adding new derived parameters to a run
    runparser.add_argument(
        '--Der-target-folder', dest="Der_target_folder",
        help=helpdict['Der-target-folder'], type=str, default='')
    runparser.add_argument(
        '--Der-param-list', dest='derived_parameters',
        help=helpdict['Der-param-list'], type=str, default='', nargs='+')

    ###############
    # Importance Sampling Arguments
    runparser.add_argument(
        '--IS-starting-folder', dest='IS_starting_folder',
        help=helpdict['IS-starting-folder'], type=str, default='', nargs='+')

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
    infoparser = add_subparser(subparser, 'info',
                               help="analyze the MCMC chains")

    # -- folder to analyze
    infoparser.add_argument('files', help=helpdict['files'],
                            nargs='+')
    # Silence the output (no print on the console)
    infoparser.add_argument('--silent', help=helpdict['silent'],
                            action='store_true')
    # -- to only write the covmat and bestfit, without computing the posterior
    infoparser.add_argument('--minimal', help=helpdict['minimal'],
                            action='store_true')
    # -- number of bins (defaulting to 20)
    infoparser.add_argument('--bins', help=helpdict['bins'],
                            type=int, default=20)
    # -- to remove the mean-likelihood line
    infoparser.add_argument('--no-mean', help=helpdict['no-mean'],
                            dest='mean_likelihood', action='store_false')
    # -- possible plot file describing custom commands
    infoparser.add_argument('--extra', help=helpdict['extra'],
                            dest='optional_plot_file', default='')
    # -- if you just want the covariance matrix, use this option
    infoparser.add_argument('--noplot', help=helpdict['noplot'],
                            dest='plot', action='store_false')
    # -- if you just want to output 1d posterior distributions (faster)
    infoparser.add_argument('--noplot-2d', help=helpdict['noplot-2d'],
                            dest='plot_2d', action='store_false')
    # -- when plotting 2d posterior distribution, use contours and not contours
    # filled (might be useful when comparing several folders)
    infoparser.add_argument('--contours-only', help=helpdict['contours-only'],
                            dest='contours_only', action='store_true')
    # -- if you want to output every single subplots
    infoparser.add_argument('--all', help=helpdict['all'], dest='subplot',
                            action='store_true')
    # -- to change the extension used to output files (pdf is the default one,
    # but takes long, valid options are png and eps)
    infoparser.add_argument('--ext', help=helpdict['ext'],
                            type=str, dest='extension', default='pdf')
    # -- only analyze the markovian part of the chains
    infoparser.add_argument('--keep-non-markovian', help=helpdict['keep-non-markovian'],
                            dest='markovian', action='store_false')
    # -- fraction of chains to be analyzed after burn-in removal (defaulting to 1.0)
    infoparser.add_argument('--keep-fraction', help=helpdict['keep-fraction'],
                            type=float, dest='keep_fraction', default=1.0)
    # -- calculate the covariant matrix when analyzing the chains
    infoparser.add_argument('--want-covmat', help=helpdict['want-covmat'],
                            dest='want_covmat', action='store_true')
    # -------------------------------------
    # Further customization
    # -- fontsize of plots (defaulting to 16)
    infoparser.add_argument('--fontsize', help=helpdict['fontsize'],
                            type=int, default=16)
    # -- ticksize of plots (defaulting to 14)
    infoparser.add_argument('--ticksize', help=helpdict['ticksize'],
                            type=int, default=14)
    # -- linewidth of 1d plots (defaulting to 4, 2 being a bare minimum for
    # legible graphs
    infoparser.add_argument('--line-width', help=helpdict['line-width'],
                            type=int, default=4)
    # -- number of decimal places that appear on the tick legend. If you want
    # to increase the number of ticks, you should reduce this number
    infoparser.add_argument('--decimal', help=helpdict['decimal'], type=int,
                            default=3)
    # -- number of ticks that appear on the graph.
    infoparser.add_argument('--ticknumber', help=helpdict['ticknumber'],
                            type=int, default=3)
    # -- legend type, to choose between top (previous style) to sides (new
    # style). It modifies the place where the name of the variable appear.
    infoparser.add_argument('--legend-style', help=helpdict['legend-style'],
                            type=str, choices=['sides', 'top'],
                            default='sides')

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

    # check for MPI
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
    except ImportError:
        # set all chains to master if no MPI
        rank = 0

    # Some check to perform when running the MCMC chains is requested
    if args.subparser_name == "run":

        # If the user wants to start over from an existing chain, the program
        # will use automatically the same folder, and the log.param in it
        if args.restart is not None:
            args.folder = os.path.sep.join(
                args.restart.split(os.path.sep)[:-1])
            args.param = os.path.join(args.folder, 'log.param')
            if not args.silent:
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
                        if not args.silent and not rank:
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
