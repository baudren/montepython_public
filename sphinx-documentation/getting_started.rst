Getting Started
===============


Foreword
--------

Python has a very nice way of handling errors in the execution.
Instead of a segmentation fault as in C, when the code breaks, you
have access to the whole stack of actions that lead to the error. This
helps you pin-point which function was called, which line was
responsible for the error.

It can however be lengthy, and to help everyone reading it, a
messaging system was implemented in Monte Python. After a blank line,
a summary of the actual error will be displayed. When reporting for an
error, please attach the entire output, as this is priceless for
debugging.


Input parameter file
--------------------

An example of input parameter file is provided with the download
package, under the name :code:`example.param`. Input files are
organised as follows:

.. code::

    data.experiments = ['experiment1', 'experiment2', ...]

    data.parameters['cosmo_name']       = [mean, min, max, sigma, scale, 'cosmo']
    ...

    data.parameters['nuisance_name']    = [mean, min, max, sigma, scale, 'nuisance']
    ...

    data.parameters['cosmo_name']       = [mean, min, max, sigma, scale, 'derived']
    ...

    data.cosmo_arguments['cosmo_name']           = value

    data.N = 10
    data.write_step = 5

The first command is rather explicit. You will list there all the
experiments you want to take into account. Their name should coincide
with the name of one of the several sub-directories in the
:code:`montepython/likelihoods/` directory. Likelihoods will be explained in the
:doc:`likelihood_class`

In :attr:`data.parameters`, you can list all the cosmo and nuisance
parameter that you want to vary in the Markov chains. For each of them
you must give an array with six elements, in this order:

    * **mean value** (your guess for the best fitting value, from
      which the first jump will start)
    * **minimum value** (set to `-1` or `None` for unbounded prior edge),
    * **maximum value** (set to `-1` or `None` for unbounded prior edge),
    * **sigma** (your guess for the standard deviation of the
      posterior of this parameter, its square will be used as the
      variance of the proposal density when there is no covariance
      matrix including this parameter passed as an input),
    * **scale** (most of the time, it will be 1, but occasionnaly you
      can use a rescaling factor for convenience, for instance {\tt
      1.e-9} if you are dealing with :math:`A_s` or :code:`0.01` if
      you are dealing with :math:`\omega_b`)
    * **role** (:code:`cosmo` for MCMC parameters used by the Boltzmann
      code, :code:`nuisance` for MCMC parameters used only by the
      likelihoods, and :code:`derived` for parameters not directly varied by
      the MCMC algorithm, but to be kept in the chains for memory).

In  :attr:`data.cosmo_arguments`, you can pass to the Boltzmann code
any parameter that you want to fix to a non-default value
(cosmological parameter, precision parameter, flag, name of input file
needed by the Bolztmann code, etc.). The names and values should be
the same as in a |CLASS| input file, so the values can be numbers or a
strings, e.g:

.. code::

    data.cosmo_arguments['Y_He']           = 0.25

or

.. code::

    data.cosmo_arguments['Y_He']      = 'BBN'
    data.cosmo_arguments['sBBN file'] = data.path['cosmo']+'/bbn/sBBN.dat'

All elements you input with a :code:`cosmo`, :code:`derived` or
:code:`cosmo_arguments` role will be interpreted by the cosmological
code (only |CLASS| so far). They are not coded anywhere inside |MP|.
|MP| takes parameter names, assigns values, and passes all of these to
|CLASS| as if they were written in a |CLASS| input file. The
advantages of this scheme are obvious. If you need to fix or vary
whatever parameter known by |CLASS|, you don't need to edit |MP|, you
only need to write these parameters in the input parameter file. Also,
|CLASS| is able to interpret input parameters from a |CLASS| input
file with a layer of simple logic, allowing to specify different
parameter combinations.  Parameters passed from the parameter file of
|MP| go through the same layer of logic.

If a :code:`cosmo`, :code:`derived` or :code:`cosmo_arguments`
parameter is not understood by the Boltzmann code, |MP| will stop
and return an explicit error message. A similar error will occur if
one of the likelihoods requires a :code:`nuisance` parameter that is
not passed in the list.

You may wish occasionally to use in the MCMC runs a new parameter
that is not a |CLASS|  parameter, but can be mapped to one or
several |CLASS| parameters (e.g. you may wish to use in your chains
:math:`\log(10^{10}A_s)` instead of :math:`A_s`). There is a function,
in the module :mod:`data`, that you can edit to define such
mappings: it is called  :func:`update_cosmo_arguments
<data.data.update_cosmo_arguments>`. Before calling \CLASS, this
function will simply substitute in the list of arguments your
customized parameters by some |CLASS| parameters.  Several exemple of
such mappings are already implemented, allowing you for instance to
use :code:`'Omega_Lambda'`, :code:`'ln10^{10}A_s'` or
:code:`'exp_m_2_tau_As'` in your chains. Looking at these examples,
the user can easily write new ones even without knowing python.

The last two lines of the input parameter file are the number of steps
you want your chain to contain (:code:`data.N`) and the number of
accepted steps the system should wait before writing it down to a file
(:code:`data.write_step`). Typically, you will need a rather low
number here, e.g. :code:`data.write_step = 5` or :code:`10`. The
reason for not setting this parameter to one is just to save a bit of
time in writing on the disk.

In general, you will want to specify the number of steps in the
command line, with the option :code:`-N` (see section~\ref{commands}).
This will overwrite the value passed in the input parameter file. The
value by default in the parameter file, :code:`data.N = 10`, is
intentionnaly low, simply to prevent doing any mistake while testing
the program on a cluster.


Output directory
----------------

You are assumed to use the code in the following way: for every set of
experiments and parameters you want to test, including different
priors, some parameters fixed, etc\ldots you should use one output
folder. This way, the folder will keep track of the exact calling of
the code, allowing you to reproduce the data at later times, or to
complete the existing chains. All important data are stored in your
:code:`folder/log.param` file.

Incidentaly, if you are starting the program in an existing folder,
already containing a :code:`log.param` file, then you do not even have
to specify a parameter file: the code will use it automatically. This
will avoid mixing things up. If you are using one anyway, the code
will warn you that it did not read it: it will always only use the
:code:`log.param` file.

In the folder :code:`montepython`, you can create a folder
:code:`chains` where you will organize your runs e.g. in the
following way:

.. code::

    montepython/chains/set_of_experiments1/model1
    montepython/chains/set_of_experiments1/model2
    ...
    montepython/chains/set_of_experiments2/model1
    montepython/chains/set_of_experiments2/model2
    ...

The minimum amount of command lines for running |MP| is an input file,
an output directory and a configuration file: if you have already
edited :code:`defaut.conf` or copied it to your own
:code:`my-machine.conf`, you may already try a mini-run with the
command

.. code::

    montepython]$ montepython/MontePython.py --conf my-machine.conf -p example.param -o test

If your configuration file is called :code:`defaut.conf`, you may even omit it (it is the default) and write only

.. code::

    montepython]$ montepython/MontePython.py -p example.param -o test

Analyzing chains and plotting
-----------------------------


Once you have accumulated a few chains, you can analyse the run to get
convergence estimates, best-fit values, minimum credible intervals, a
covariance matrix  and some plots of the marginalised posterior
probability. You can run again |MP| with the :code:`info` prefix
followed by the name of a directory or of several chains, e.g.
:code:`info chains/myrun/` or :code:`info chains/myrun/2012-10-26*
chains/myrun/2012-10-27*`. There is no need to pass an input file
with parameter names since they have all been stored in the
:code:`log.param`.

Information on the acceptance rate and minimum :math:`-\log{\cal
L}=\chi^2_{\rm eff}/2` is written in :code:`chains/myrun/myrun.log`.
Information on the convergence (Gelman-Rubin test for each chain
parameter), on the best fit, mean and minimum credible interval for
each parameter at the 68.26\%, 95.4\%, 99.7\% level are written in
horizontal presentation in :code:`chains/myrun/myrun.h_info`, and in
vertical presentation in :code:`chains/myrun/myrun.v_info` (without
99.7\% in the vertical one). A latex file to produce a table with
parameter names, means and 68\% errors in written in
:code:`chains/myrun/myrun.tex`.

The covariance matrix of the run is written in
:code:`chains/myrun/myrun.covmat`. It can be used as an input for the
proposal density in a future run. The first line, containing the
parameter name, will be read when the covariance matrix will be passed
in input. This means that the list of parameters in the input
covariance matrix and in the run don't need to coincide: the code will
automatically eliminate, add and reorder parameters (see
:func:`mcmc.get_covariance_matrix`). Note that the rescaling factors
passed in the input file are used internally during the run and also
in the presentation of results in the :code:`.h_info`,
:code:`.v_info`, :code:`.tex` files, but not in the covariance matrix
file, which refers to the true parameters.

The 1D posteriors and 2D posterior contours are plotted in
:code:`chains/myrun/plots/myrun_1D.pdf` and
:code:`chains/myrun/plots/myrun_triangle.pdf`. You will find in the
:doc:`parser_mp` documentation a list of commands to customize the
plots.

When the chains are not very converged and the posterior probability
has local maxima, the code will fail to compute minimum credible
intervals and say it in a warning. The two solutions are either to
re-run and increase the number of samples, or maybe just to decrease
the number of bins with the :code:`--bins` option.


Global running strategy
-----------------------

In the current version of |MP|, we deliberately  choose not to use MPI
communication between instances of the code. Indeed the use of MPI
usually makes the installation step more complicated, and the gain is,
in our opinion, not worth it. Several chains are launched as
individual serial runs (if each instance of |MP| is launched on
several cores, |CLASS| and the WMAP likelihood will parallelize since
they use OpenMP). They can be run with the same command since chain
names  are created automatically with different numbers for each
chain: the chain names are in  the form :code:`yyyy-mm-dd_N__i.txt`
where :code:`yyyy` is the year, :code:`mm` the month, :code:`dd` the
day, :code:`N` the requested number of steps and :code:`i` the
smallest available integer at the time of starting a new run.

However the absence of communication between chains implies that the
proposal density cannot be updated automatically during the initial
stage of a run. Hence the usual strategy consists in launching a first
run with a poor (or no) covariance matrix, and a low acceptance rate;
then to analyze this run and produce a better covariance matrix; and
then to launch a new run with high acceptance rate, leading to nice
plots. Remember that in order to respect strictly markovianity and the
Metropolis Hastings algorithm, one should not mix up chains produced
with different covariance matrices: this is easy if one takes
advantage of the :code:`info` syntax, for example :code:`info
chains/myrun/2012-10-26_10000*`. However mixing runs that started from
very similar covariance matrices is harmless.

It is also possible to run on several desktops instead of a single
cluster. Each desktop should have a copy of the output folder and with
the same :code:`log.param` file, and after running the chains can be
grouped on a single machine and analyse. In this case, take care of
avoiding that chains are produced with the same name (easy to ensure
with either the :code:`-N` or :code:`--chain-number` options). This is
a good occasion to keep the desktops of your department finally busy.


.. |CLASS| replace:: *Class*
.. |MP| replace:: *Monte Python*
