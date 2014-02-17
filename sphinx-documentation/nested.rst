Using MultiNest with MontePython
================================

Monte Python uses the implementation of `MultiNest <http://ccpforge.cse.rl.ac.uk/gf/project/multinest/>`__ by F. Feroz M. Hobson and M. Bridges, through the Python wrapper `PyMultiNest <http://github.com/JohannesBuchner/PyMultiNest>`__ by J. Buchner.

.. NOTE::
   By using MultiNest and PyMultiNest, you agree to their respective licenses, that can be found in the :code:`LICENCE` files into the respective installation folders.

Installation
------------

1. MultiNest
~~~~~~~~~~~~

Download MultiNest from `here <http://github.com/JohannesBuchner/MultiNest>`__, either using the `releases page <http://github.com/JohannesBuchner/MultiNest/releases>`__ or cloning with git

.. code::

    $ git clone http://github.com/JohannesBuchner/MultiNest

.. NOTE:: MultiNest requires the libraries :code:`lapack` and :code:`mpi` (optional), and the compilation tool :code:`cmake`.

We now follow the instructions in the :code:`README` file to compile:

.. code::

    $ cd /path/to/MultiNest  # Folder where MultiNest was downloaded
    $ cd build
    $ cmake ..  # -DCMAKE_Fortran_COMPILER=gfortran
    $ make

It is not necessary to install (:code:`make install`), but it is so to add the folder in which the library was installed to the list of paths in which compilers will looll for libraries to link. In Linux (and other systems using a bash shell), this consists simply of adding at the end of the file :code:`\home\<your user>\.bashrc` the line

.. code::

    export LD_LIBRARY_PATH=/path/to/MultiNest/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}

By default, MultiNest will try to use the Intel Fortran Compiler if installed. If you want to use :code:`gfortran` instead, add the flag :code:`-DCMAKE_Fortran_COMPILER=gfortran` to :code:`cmake`.

2. PyMultiNest
~~~~~~~~~~~~~~

Download PyMultiNest from `here <http://github.com/JohannesBuchner/PyMultiNest>`__, either using the `releases
page <http://github.com/JohannesBuchner/PyMultiNest/releases>`__ or cloning with git

.. code::

    $ git clone http://github.com/JohannesBuchner/PyMultiNest

Now go to the installation directory and install:

.. code::

    $ cd PyMultiNest
    $ python setup.py install  # --user

You may need the flag :code:`--user` in the last command if you do not have admin privileges.

If everything went ok, you should be able to run :code:`import pymultinest` in a python console without any output. In that case, you are good to go!


Basic usage and parameters
--------------------------

**TODO**
   Why would anyone like to use MultiNest?

The MultiNest sampling is invoked with the command line option :code:`-m NS`. As in the MCMC case, the parameter file is read and the sampling is launched. The output files are created inside a subfolder :code:`NS` inside the chain folder.

**TODO**
   What to do after finished.

**TODO**
   Derived parameters can be used as in MCMC

**TODO**
   Caveat about column ordering when usering clustering parameters

Thorough descriptions of the parameters can be found in the README file on MultiNest and the documentation of PyMultiNest. Here is how Monte Python manages them:

Automatic parameters
~~~~~~~~~~~~~~~~~~~~

(Technical section, you can skip)

The following parameters are defined automatically by the content of the :code:`.param` file, and you should not care about them:

-  :code:`ndims | n_dims` : number of varying parameters.
-  :code:`nPar | n_params` : number of varying parameters.
-  :code:`root | outputfiles_basename` : prefix of the MultiNest output files: name of the chain plus a hyphen.
-  :code:`outfile | write_output` : whether to write output files (yes, of course).
-  :code:`resume | resume` : whether to allow for resuming a previously killed run, enabled by default.
-  :code:`initMPI | init_MPI` : try to use MPI (only if MultiNest was compiled with MPI on), enabled by default.
-  :code:`feedback | verbose (True)` : print information periodically.

Manually set parameters
~~~~~~~~~~~~~~~~~~~~~~~

The following parameters can be changed by hand to adjust the sampling to one's needs. In the following, they are presented as

.. code::

    [MultiNest name] | [PyMultiNest name] (default value)

and are set in every run by command line options as

.. code::

    --NS_option_[PyMultiNest name] [value]

E.g. to set the number of "live points" to 100, one should add to the command :code:`python MontePython.py [...] -m NS` the option

.. code::

    --NS_option_n_live_points 100

.. NOTE::
   The default values are those defined in PyMultiNest (at least most of them), and are not hard-coded in Monte Python.

.. NOTE::
   The parameters not appearing in the following lists are not managed in the current implementation.

General sampling options
^^^^^^^^^^^^^^^^^^^^^^^^

-  :code:`nlive | n_live_points (400)` : number of points used in every iteration.
-  :code:`IS | importance_nested_sampling (True)` : whether to use Importance Nested Samplin (see `arXiv:1306.2144 <http://arxiv.org/abs/1306.2144>`__).
-  :code:`efr | sampling_efficiency (0.8)` : defines the sampling efficiency (see 'Use cases' below).
-  :code:`ceff | const_efficiency_mode (True)` : constant efficiency mode -- slower, but more accurate evidence estimation.
-  :code:`seed | seed (-1)`: seed of the random number generator (if negative, uses system clock).
-  :code:`logZero | log_zero (-1e90)` : if the log-likelihood of a sample is smaller than this value, the sample is ignored.
-  :code:`updInt | n_iter_before_update (100)` : number of iteration after which the output files are updated.

Ending conditions
^^^^^^^^^^^^^^^^^

-  :code:`tol | evidence_tolerance (0.5)`
-  :code:`maxiter | max_iter (0)`

The sampling ends after :code:`maxiter` iterations, or when the tolerance condition on the evidence defined by :code:`tol` is fulfilled, whatever happens first.

Multi-modal sampling
^^^^^^^^^^^^^^^^^^^^

-  :code:`mmodal | multimodal (False)` : whether to try to find separate modes in the posterior.
-  :code:`maxModes | max_modes (100)` : maximum number of separate modes to consider.
-  :code:`Ztol | mode_tolerance (-1e90)` : if the local log-evidence is greater than this value, a mode is created.

.. NOTE::
   Here, multi-modal sampling is disabled by default. If enabled, Imporance Nested Sampling will be automatically disabled, since both modes are not compatible.

We left out the option concerning the *clustering parameters*, i.e. on which parameters's subspace is MultiNest to look for posterior mode separation:

.. code::

   nCdims | n_clustering_params

In (Py)MultiNest, clustering parameters are specified as the :code:`n` first ones, which **must** be at the beginning of the parameters list. Here, instead, we override that limitation, and the clustering parameters are specified as

.. code::

   --NS_option_clustering_params param1 param2 ...

The reason for doing it this way is giving more flexibility to the user, being able to change the clustering parameters without having to modify the ordering of the parameters in the :code:`param` file to put the clustering parameters at the beginnig. But this comes at a price: the raw MultiNest chain files have the parameters ordered with the clustering parameters at the beginning, and then the rest as they appear in the :code:`.param` file. The ordering of the parameters is save to a file :code:`[chain name].paramnames` in the :code:`NS` subfolder. If you intend to use MustiNest's raw output files, you must take this into account! If, instead, you use nested sampling simply as a means to get a covariance matrix and some sample points (saved in :code:`chain_NS__[accepted/rejected].txt`), you do not need to care about this.


Usage cases (and suggested values of the options)
-------------------------------------------------

**TODO**

**I want a good covariance matrix**

(low sampling)

**I want to sample the posterior thoroughly, and I know there is only one mode**

**I want to map a posterior to find the different modes**

(multimodal, clustering, medium sampling)

**I want to sample a multi-modal posterior thoroughly**

(multimodal, clustering, high sampling)

**I want to evaluate the evidence of a model**

**(Other cases...)**


