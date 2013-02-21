===========================================================
Monte Python, a Monte Carlo Markov Chain code (with Class!)
===========================================================

:Author: Benjamin Audren <benjamin.audren@epfl.ch>
:Version: Version 1.1


If you are searching for specific examples of a work session, please refer to
the pdf documentation.


Prerequisites
-------------

* You need the python program **version 2.7** or above, but less than 3.0.
  Note that lower versions of python will work, down to 2.6 (tested), if you
  add manually two extra packages (
  `ordereddict <http://code.activestate.com/recipes/576693/>`_ and 
  `argparse <https://pypi.python.org/pypi/argparse/1.2.1>`_).

* Your python of choice must have numpy (version >= 1.4.1), and cython. The
  later is used to wrap CLASS in python.

* *[optional]* If you want to use fully the plotting capabilities of Monte Python,
  you also need the scipy module, with interpolate.


The MontePython part
--------------------

Move the .tar.bz2 file to the place of your convenience, untar its content

.. code::

    $ bunzip2 montepython-vx.y.tar.bz2
    $ tar -xvf montepython-vx.y.tar

This will create a directory named montepython into your current directory.
You can add the following line to your .bashrc file:

.. code::

    export PATH=/path/to/MontePython/code/:$PATH

to be able to call the program from anywhere.

You will need to adapt only two files to your local configuration. The first
is the main file of the code `code/MontePython.py`, and it will be the only
time you will have to edit it, and it is simply to accomodate different
possible configurations of your computer.

Its first line reads

.. code::

    #!/usr/bin/python

This should be changed to wherever is your prefered python distribution
installed. For standard distribution, this should already be working. Now,
you should be able to execute directly the file, i.e. instead of calling:

The second file to modify is located in the root directory of Monte Python :
`default.conf`. This file will be read (and stored) whenever you execute the
program, and will search for your cosmological code path, your data path, and
your wmap wrapper path. You can alternatively create a second one, `my.conf`,
containing your setup, and then run the code providing this file (with the flag
`-conf`)


The Class part
--------------

Go to your class directory, and do **make clean**, then **make**. This builds the
`libclass.a`, needed for the next step. From there, 

.. code::

    $ cd python/
    $ python setup.py build
    $ python setup.py install --user

This will compile the file `classy.pyx`, which is the python wrapper for CLASS,
into a library, `classy.so`, located in the `build/` subdirectory. This is the
library called in Monte Python afterwards.

If this step fails, check that you have cython installed, numpy (a numerical
package for python), python (well... did I say this code was in python ?) with
a version > 2.6.  If this step fails again, kindly ask your sys.admin, (s)he
is there for this, after all. Note that the installation (last command) is
not stricly speaking mandatory.

Remember that if you modify CLASS to implement some new physics, you will need to
perform this part again for the new CLASS.


The wmap wrapper part
---------------------

Go to your `wrapper_wmap` sub-folder, and execute:

.. code::

    $ ./waf configure --install_all_deps

This will automatically install the wrapper. Please refer to the
MontePython.pdf documentation for further details, and more options concerning
this wrapper.

Do not forget to source your wrapper everytime you want to use it:

.. code::

    $ source YourWlikPath/bin/clik_profile.sh

You can put this line in your .bashrc file, and you should put it in your
scripts for cluster computing.


For Planck collaborators (for future Planck users)
--------------------------------------------------

Replace the upper section by Clik. You also have to source the file everytime
you want Monte Python to use it.


Enjoying the difference
-----------------------

Now the code is installed. Go anywhere, and just call

.. code::

    $ ./MontePython.py --help

To see a list of all commands. There are two essential ones, without which
the program will not start. At minimum, you should precise an output folder
('-o') and a parameter file ('-p'). An example of parameter file is found in
the main directory of MontePython (planck.param, for instance).

A typical call would then be:

.. code::

    $ ./MontePython.py -o planck -p planck.param

If non existant, the `planck/` folder will be created, and a run with the
number of steps described in `planck.param` will be started. To run a chain with
less steps, one can type:

.. code::

    $ ./MontePython.py -o planck -p planck.param -N 100

If you want to analyse the run, then just type

.. code::

    $ ./MontePython.py -info planck/


Details and Examples
--------------------

Please refer to the pdf documentation for further details.
