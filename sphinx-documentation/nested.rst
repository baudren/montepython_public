Using MultiNest with Monte Python
=================================

Monte Python can easily use the implementation of `MultiNest <http://ccpforge.cse.rl.ac.uk/gf/project/multinest/>`__ by F. Feroz and M. Hobson [1]_, through the Python wrapper `PyMultiNest <http://github.com/JohannesBuchner/PyMultiNest>`__ by J. Buchner [2]_.

Some hints about why and how to use MultiNest can be found in `A Basic usage and parameters`_. A more thorough description of the MultiNest sampler can be found in the MultiNest papers [1]_. The `PyMultiNest tutorial <http://johannesbuchner.github.io/pymultinest-tutorial/>`_ is also worth checking out, as well as the respective :code:`README` files of both MultiNest and PyMultiNest.

.. NOTE::
   By using MultiNest and PyMultiNest, you agree to their respective licenses, that can be found in the :code:`LICENSE` (or :code:`LICENCE`) files into the respective installation folders.


Installation
------------

This basically follows the installation procedure in the `PyMultiNest documentation <http://johannesbuchner.github.io/PyMultiNest/pymultinest.html>`_.

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

By default, MultiNest will try to use the Intel Fortran Compiler if installed. If you want to use :code:`gfortran` instead, add the flag :code:`-DCMAKE_Fortran_COMPILER=gfortran` to :code:`cmake`. If you want to force compilation with MPI support (though :code:`gfortran` should autodetect it and implement it in most systems), use :code:`mpif90` as compiler name.

2. PyMultiNest
~~~~~~~~~~~~~~

.. NOTE::
   In a future implementation (hopefully soon), PyMultiNest will be installed automatically. Right now, the installations has to be done manually as it follows.

Download PyMultiNest from `here <http://github.com/JohannesBuchner/PyMultiNest>`__, either using the `releases
page <http://github.com/JohannesBuchner/PyMultiNest/releases>`__ or cloning with git

.. code::

    $ git clone http://github.com/JohannesBuchner/PyMultiNest

Now go to the installation directory and install:

.. code::

    $ cd PyMultiNest
    $ python setup.py install  # --user

You may need the flag :code:`--user` in the last command if you do not have admin privileges.

If everything went ok, you should be able to run :code:`import pymultinest` in a python console without getting any output. In that case, you are good to go!


Basic usage and parameters
--------------------------

The MultiNest sampling is invoked with the command line option :code:`-m NS`.
As in the MCMC case, the parameter file is read and the sampling is launched.
The output files are created inside a subfolder :code:`NS` inside the chain
folder. This will create the expected :code:`log.param` file inside the chain's
root folder, and the expected raw MultiNest files in the :code:`NS` subfolder
(see MultiNest's :code:`README`), along with two more files: :code:`[chain
name].paramnames`, which contains the ordering of the parameters in the nested
sampling chain files (not necessarily the ordering in which they appear in the
:code:`log.param`, since clustering parameters must go first), and
:code:`[chain name].arguments`, which contains the user defined MultiNest
arguments and their values (see below).

.. NOTE::
   If the sampling has been interrupted, simply run it again and MultiNest
   should be able to restart where it finished. If you intend to start a new
   sampling with different parameters for MultiNest, it is safer to delete the
   :code:`NS` subfolder (otherwise, the behaviour is not well defined).

.. NOTE::
   MultiNest can benefit greatly from being run in parallel with MPI. If it has
   been correctly compiled with MPI (try to run the examples distributed with
   the MultiNest code with MPI), it is possible to take advantage of it using
   Monte Python: simply run the sampler :code:`python MontePython.py` preceded
   by the appropriate MPI runner (:code:`mpirun` for Open MPI, :code:`mpiexec`
   for MPICH, etc.).

Once the sampling has finished, the output of it can be analised as in the MCMC
case with :code:`MontePython.py -info [chain_folder]/NS` (notice that one must
specify the :code:`NS` subfolder). This will create a chain file in the chain
root folder containing the (accepted) points of the nested sampling, and it
will be automatically analysed as a MCMC chain, producing the expected files
and plots.

The MultiNest parameters are added after the :code:`-m NS` flag in the command
line. They are described in the next section (more thorough descriptions are to
be looked for within the MultiNest documentation).

Automatic parameters
~~~~~~~~~~~~~~~~~~~~

(Technical section, you can skip)

The following parameters are defined automatically by the content of the
:code:`.param` file, and you should not care about them:

-  :code:`ndims | n_dims` : number of varying parameters.
-  :code:`nPar | n_params` : number of varying parameters.
-  :code:`root | outputfiles_basename` : prefix of the MultiNest output files: name of the chain plus a hyphen.
-  :code:`outfile | write_output` : whether to write output files (yes, of course).
-  :code:`resume | resume` : whether to allow for resuming a previously killed run, enabled by default.
-  :code:`initMPI | init_MPI` : initialise MPI within MultiNest (disabled: MPI, if requested, is initialised by Monte Python).
-  :code:`feedback | verbose (True)` : print information periodically.

Manually set parameters
~~~~~~~~~~~~~~~~~~~~~~~

The following parameters can be changed by hand to adjust the sampling to one's needs. In the following, they are presented as

.. code::

    [MultiNest name] | [PyMultiNest name] (default value)

and are set in every run by command line options as

.. code::

    --NS_[PyMultiNest name] [value]

E.g. to set the number of "live points" to 100, one should add to the command :code:`python MontePython.py [...] -m NS` the option

.. code::

    --NS_n_live_points 100

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

   --NS_clustering_params param1 param2 ...

The reason for doing it this way is giving more flexibility to the user, being
able to change the clustering parameters without having to modify the ordering
of the parameters in the :code:`param` file to put the clustering parameters at
the beginnig. But this comes at a price: the raw MultiNest chain files have the
parameters ordered with the clustering parameters at the beginning, and then
the rest as they appear in the :code:`.param` file. The ordering of the
parameters is save to a file :code:`[chain name].paramnames` in the :code:`NS`
subfolder. If you intend to use MustiNest's raw output files, you must take
this into account! If, instead, you use nested sampling simply as a means to
get a covariance matrix and some sample points (saved in
:code:`chain_NS__[accepted/rejected].txt`), you do not need to care about this.

References
----------

.. [1] `arXiv:0704.3704 <http://arxiv.org/abs/0704.3704>`_,
       `arXiv:0809.3437 <http://arxiv.org/abs/0809.3437>`_ and
       `arXiv:1306.2144 <http://arxiv.org/abs/1306.2144>`_.

.. [2] `arXiv:1402.0004 <http://arxiv.org/abs/1402.0004>`_.
       
