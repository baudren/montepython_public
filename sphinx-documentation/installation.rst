Installation Guide
==================

Prerequisites
-------------

Python
^^^^^^

First of all, you need a clean installation of Python_ (version 2.7 is
better, though it works also with 2.6 and 2.5 (see below). version 3.0
is not supported), with at least the numpy_ module (version :math:`\geq 1.4.1`) and
the cython_ module. This last one is to convert the C code CLASS_ into
a Python class.

If you also want the output plot to have cubic interpolation for
analyzing chains, you should also have the scipy_ module (at least
version 0.9.0). In case this one is badly installed, you will have an
error message when running the analyze module of |MP|, and obtain only
linear interpolation. Though not fatal, this problem produces ugly
plots.

To test for the presence of the modules **numpy**,  **scipy**,
**cython** on your machine, you can type

.. code::

   $ python
   $ >>> import numpy
   $ >>> import scipy
   $ >>> import cython
   $ >>> exit()

If one of these steps fails, go to the corresponding websites, and
follow the instructions (if you have the privilege to have the root
password on your machine, an `apt-get install python-numpy`,
`python-scipy` and `cython` will do the trick. Otherwise, all these
packages can also be downloaded and installed locally, with the
command :code:`python setup.py install --user`).

Note that you can use the code with Python 2.6 also, even though you
need to download two packages separately ordereddict_ and argparse_.
For this, it is just a matter of downloading the two files
`ordereddict.py` and `argparse.py`), and placing them in
your code directory without installation steps.

Class
^^^^^

Next in line, you must compile the python wrapper of CLASS_. Download
the latest version (:math:`\geq 1.5.0`), and follow the basic instruction.
Instead of  :code:`make class`, type :code:`make -j`. If you are using a |CLASS|
version :math:`\geq 2.3.0`, the wrapper will be installed by doing this step,
so skip ahead.

In case you are using an older version of |CLASS|, the compilation only created
an archiv `.ar` of the code, useful in the next step. After this, do:

.. code::

   class]$ cd python/
   python]$ python setup.py build
   python]$ python setup.py install --user

If you have correctly installed cython, this should add Classy as a new python
module. You can check the success of this operation by running the following
command:

.. code::

  ~]$ python
  >>> from classy import Class

If the installation was successful, this should work within any
directory. If you get no error message from this line, you know
everything is fine.

.. note::

  If the step :code:`python setup.py install --user` does not succeed, but that
  the :code:`build` is successful, then as far as |MP| is concerned, there are
  no issues. The code will be found nonetheless.

If at some point you have several different coexisting versions of
|CLASS| on the system, and you are worried that |MP| is not using the
good one, rest reassured. As long as you run |MP| with the proper
path to the proper |CLASS| in your configuration file (see
Installation_) then it will use this one.




Installation
------------

Main code
^^^^^^^^^

Move the latest release of |MP| to one of your folders, called e.g.
:code:`code/` (for instance, this could be the folder containing also
:code:`class/`), and untar its content:

.. code::

  code]$ bunzip montepython-v1.0.0.tar.bz2 
  code]$ tar -xvf montepython-v1.0.0.tar
  code]$ cd montepython

You will have to create one file holding the path of the codes you want to use.
There is a predefined template, :code:`default.conf.template`, in the root
directory of the code. You should copy it to a new file called
:code:`default.conf`, which will tell |MP|, where your other programs (in
particular |CLASS|) are installed, and where you are storing the data for the
likelihoods. It will be interpreted as a python file, so be careful to
reproduce the syntax exactly. At minimum, **default.conf** should contain one
line, filled with the path of your :code:`class/` directory:

.. code::

   path['cosmo']   = 'path/to/your/class/'
 
To check that |MP| is ready to work, simply type :code:`python
montepython/MontePython.py --help` (or just :code:`montepython/MontePython.py
--help`). This will provide you with a short description of the
available command line arguments, explained in :doc:`parser_mp`. 


Planck likelihood
^^^^^^^^^^^^^^^^^

With the release of Planck data comes the release of its likelihood.
It is distributed from this `ESA website
<http://www.sciops.esa.int/index.php?project=planck&page=Planck_Legacy_Archive>`_,
along with the data. Download all `tar.gz` files, extract them to the
place of your convenience.

The Planck Likelihood Code (**plc**) is based on a library called
`clik`. It will be extracted, alongside several `.clik` folders that
contain the likelihoods. The installation of the code is described in
the archive, and it uses an auto installer device, called `waf`.

.. warning::

  Note that you **are strongly advised** to configure `clik` with the
  Intel mkl library, and not with lapack. There is a massive gain in
  execution time: without it, the code is dominated by the execution
  of the low-l polarisation data from WMAP.


In your |MP| configuration file, to use this
code, you should add the following line

.. code:: python

  path['clik'] = 'path/to/your/plc/folder/'

The four likelihoods defined in |MP| for Planck are `Planck_highl`,
`Planck_lowl`, `Planck_lensing`, `lowlike` (the polarization data from
WMAP). In each of the respective data files for these likelihood,
please make sure that the line, for instance,

.. code:: python

  Planck_highl.path_clik = data.path['clik']+'../something.clik'

points to the correct clik file. Now, before trying to run this
likelihood, you will need to source the code to your system, by
typing:

.. code::

   ~]$ source /path/to/your/plc/folder/bin/clik_profile.sh
    
Once you made sure of this, you can then use the base.param file
distributed with MontePython, that defines all the needed nuisance
parameters, the covariance matrix as well as the bestfit file, in this
command:

.. code::

  python montepython/MontePython.py -o planck/ -p base.param -c covmat/base.covmat \
  -bf bestfit/base.bestfit --conf default.conf -f 1.5

.. note::

  The use of the factor 1.5 is to increase the acceptance rate, due to
  the non gaussianity of the nuisance parameters posterior.


WMAP likelihood
^^^^^^^^^^^^^^^

.. warning::

  As of version 1.2.5, with Planck data being available, installing
  this likelihood might not be so important anymore. You might prefer
  to skip this, at it is an **optional** part of the installation
  process.

.. warning::

  So far, the use of the WMAP wrapper is separated from the Planck
  wrapper, but it might be merged in the future, as it is based on the
  same code `clik` developped internally for Planck by Karim Benabed.

To use the likelihood of WMAP, we propose a python wrapper, located in
the :code:`wrapper_wmap` directory. Just like with the |CLASS|
wrapper, you need to install it, although the procedure differs. Go to
the wrapper directory, and enter:

.. code::

  wrapper_wmap]$ ./waf configure install_all_deps

This should read the configuration of your distribution, and install
the WMAP likelihood code and its dependencies (cfitsio) automatically
on your machine. For our purpose, though, we prefer using the intel
mkl libraries, which are much faster. To tell the code about your
local installation of mkl libraries, please add to the line above some
options:

.. code::

   --lapack_mkl=/path/to/intel/mkl/10.3.8 --lapack_mkl_version=10.3

Once the configuration is done properly, finalize the installation by typing:

.. code::

  wrapper_wmap]$ ./waf install

The code will generate a configuration file, that you will need to
source before using the WMAP likelihood with |MP|. The file is
:code:`clik_profile.sh`, and is located in :code:`wrapper_wmap/bin/`.
So if you want to use the likelihood :code:`'wmap'`, before any call
to |MP| (or inside your scripts), you should execute

.. code::

  ~]$ source /path/to/MontePython/wrapper_wmap/bin/clik_profile.sh

The wrapper will use the original version of the WMAP likelihood codes
downloaded and placed in the folder
:code:`wrapper_wmap/src/likelihood_v4p1/` during the installation
process. This likelihood will be compiled later, when you will call it
for the first time from the |MP| code. Before calling it for the first
time, you could eventually download the WMAP patch from Wayne Hu's web
site, for a faster likelihood.

You should finally download the WMAP data files by yourself, place
them anywhere on your system, and specify the path to these data files
in the file :code:`likelihoods/wmap/wmap.data`.


.. _Python: http://www.python.org/
.. _numpy: http://www.numpy.org/
.. _cython: http://www.cython.org/
.. _scipy: http://www.scipy.org/
.. _argparse: https://pypi.python.org/pypi/argparse
.. _ordereddict: http://code.activestate.com/recipes/576693/
.. _CLASS: http://www.class-code.net/
.. |CLASS| replace:: *Class*
.. |MP| replace:: *Monte Python*
