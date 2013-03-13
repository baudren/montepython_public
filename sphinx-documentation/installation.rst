Installation Guide
==================

Prerequisites
-------------

Python
^^^^^^

First of all, you need a clean installation of Python_ (version 2.7 is
better, though it works also with 2.6 and 2.5 (see below). version 3.0
is not supported), with at least the numpy_ module (version |geq| 1.4.1) and
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
`python-scipy` and `cython` will do the trick.).

Note that you can use the code with Python 2.6 also, even though you
need to download two packages separately ordereddict_ and argparse_.
For this, it is just a matter of downloading the two files
`ordereddict.py` and `argparse.py`), and placing them in
your code directory without installation steps.

Class
^^^^^

Next in line, you must compile the python wrapper of CLASS_. Download
the latest version (|geq| 1.5.0), and follow the basic instruction.
Instead of  :code:`make class`, type :code:`make`. This will also
create an archiv `.ar` of the code, useful in the next step. After
this, do:

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

If the installation was successfull, this should work within any
directory. If you get no error message from this line, you know
everything is fine.

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

You will have to edit two files (the first, once for every new distribution of |MP|, and
the second, once and for all). The first to edit is
:code:`code/MontePython.py`. Its first line reads:

.. code::

  #!/usr/bin/python

You should eventually replace this path with the one of your python 2.7 executable, if different.
This modification is not crucial, it simply allows to run the code by simply typing :code:`code/Montepython.py`.
If, instead, you run it through python (\emph{i.e.}: :code:`python`
:code:`code/MontePython.py`), then this line will be disregarded.

The second file to change, and this one is crucial, is
:code:`default.conf`, in the root directory of the code. This file will
tell |MP|, where your other programs (in particular |CLASS|) are
installed, and where you are storing the data for the likelihoods. It
will be interpreted as a python file, so be careful to reproduce the
syntax exactly. At minimum, {\bf default.conf} should contain one
line, filled with the path of your :code:`class/` directory:

.. code::

   path['cosmo']   = 'path/to/your/class/'
 
For members of the Planck collaboration only: if you have installed Clik, then you should also add:

.. code::

  path['clik']    = 'path/to/your/clik/examples/'

To check that |MP| is ready to work, simply type :code:`python
code/MontePython.py --help` (or just :code:`code/MontePython.py
--help`). This will provide you with a short description of the
available command line arguments, explained in :doc:`parser_mp`. 


Planck likelihood
^^^^^^^^^^^^^^^^^

TODO

WMAP likelihood
^^^^^^^^^^^^^^^

To use the likelihood of WMAP, we propose a python wrapper, located in the
:code:`wrapper_wmap` directory. Just like with the |CLASS| wrapper, you need to
install it, although the procedure differs. Go to the wrapper directory, and enter:

.. code::

  wrapper_wmap]$ ./waf configure install_all_deps

This should read the configuration of your distribution, and install the WMAP likelihood code and its
dependencies (cfitsio) automatically on your machine. For our purpose,
though, we prefer using the intel mkl libraries, which are much faster. To
tell the code about your local installation of mkl libraries, please add to the line above some options:

.. code::

   --lapack_mkl=/path/to/intel/mkl/10.3.8 --lapack_mkl_version=10.3

Once the configuration is done properly, finalize the installation by typing:

.. code::

  wrapper_wmap]$ ./waf install

The code will generate a configuration file, that you will need to source
before using the WMAP likelihood with |MP|. The file is :code:`clik_profile.sh`, and is located
in :code:`wrapper_wmap/bin/`. So if you want to use the likelihood :code:`'wmap'`, before any call to |MP| (or inside your
scripts), you should execute

.. code::

  ~]$ source /path/to/MontePython/wrapper_wmap/bin/clik_profile.sh

The wrapper will use the original version of the WMAP likelihood codes downloaded and placed in the folder :code:`wrapper_wmap/src/likelihood_v4p1/` during the installation process. This likelihood will be compiled later, when you will call it for the first time from the |MP| code. Before calling it for the first time, you could eventually download the WMAP patch from Wayne Hu's web site, for a faster likelihood.

You should finally download the WMAP data files by yourself, place them anywhere on your system, and specify the path to these data files in the file :code:`likelihoods/wmap/wmap.data`.


.. _Python: http://www.python.org/
.. _numpy: http://www.numpy.org/
.. _cython: http://www.cython.org/
.. _scipy: http://www.scipy.org/
.. _argparse: https://pypi.python.org/pypi/argparse
.. _ordereddict: http://code.activestate.com/recipes/576693/
.. _CLASS: http://www.class-code.net/
.. |CLASS| replace:: *Class*
.. |MP| replace:: *Monte Python*

.. |geq| unicode:: U+2265
.. |leq| unicode:: U+2264
