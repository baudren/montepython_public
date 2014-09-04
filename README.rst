===========================================================
Monte Python, a Monte Carlo Markov Chain code (with Class!)
===========================================================

:Author: Benjamin Audren <benjamin.audren@epfl.ch>
:License: MIT


If you are searching for specific examples of a work session, please refer to
the pdf documentation. The code is under the MIT license. As an additional
clause, you are also required to cite the original release paper when using it
in a scientific publication: `Conservative Constraints on Early Cosmology` (see
the tail of this document)


Prerequisites
-------------

* You need the python program **version 2.7** or above, but less than 3.0.
  Note that lower versions of python will work, down to 2.6 (tested), if you
  add manually two extra packages (
  `ordereddict <http://code.activestate.com/recipes/576693/>`_ and 
  `argparse <https://pypi.python.org/pypi/argparse/1.2.1>`_).

* Your python of choice must have `numpy` (version >= 1.4.1) and `cython`. The
  later is used to wrap CLASS in python.

* *[optional]* If you want to use fully the plotting capabilities of Monte Python,
  you also need the `scipy`, with interpolate, and `matplotlib` modules.

* *[optional]* You can now use Multi Nest and the CosmoHammer with Monte
  Python, though you need to install them. Please refer to the documentation.


The MontePython part
--------------------

Move the `.tar.bz2` file to the place of your convenience, untar its content

.. code::

    $ bunzip2 montepython-vx.y.tar.bz2
    $ tar -xvf montepython-vx.y.tar

This will create a directory named `montepython` into your current directory.
You can add the following line to your `.bashrc` file:

.. code::

    export PATH=/path/to/MontePython/montepython/:$PATH

to be able to call the program from anywhere.

You will need to adapt only two files to your local configuration. The first
is the main file of the code `montepython/MontePython.py`, and it will be the only
time you will have to edit it, and it is simply to accommodate different
possible configurations of your computer.

Its first line reads

.. code::

    #!/usr/bin/python

This should be changed to wherever is your preferred python distribution
installed. For standard distribution, this should already be working. Now,
you should be able to execute directly the file, i.e. instead of calling:

The second file to modify is located in the root directory of Monte Python :
`default.conf`. This file will be read (and stored) whenever you execute the
program, and will search for your cosmological code path, your data path, and
your wmap wrapper path. You can alternatively create a second one, `my.conf`,
containing your setup, and then run the code providing this file (with the flag
`--conf`)


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

If this step fails, check that you have `cython` installed, `numpy` (a numerical
package for python), python (well... did I say this code was in python ?) with
a version > 2.6.  If this step fails again, kindly ask your system admin, (s)he
is there for this, after all. Note that the installation (last command) is
not strictly speaking mandatory.

Remember that if you modify `CLASS` to implement some new physics, you will need to
perform this part again for the new `CLASS`.


The Planck likelihood part
---------------------------

The release of the Planck data comes with a likelihood program, called
Clik, that one can recover from the `ESA website
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

Go to your plc folder, and execute the following line, taking into
account the mkl installation

.. code::

    $ ./waf configure --install_all_deps --mkl=...

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

points to the correct clik file. Do not forget to source your Planck
likelihood every time you want to use it:

.. code::

    $ source Your/Plc/bin/clik_profile.sh

You can put this line in your .bashrc file, and you should put it in your
scripts for cluster computing.



Enjoying the difference
-----------------------

Now the code is installed. Go anywhere, and just call

.. code::

    $ python montepython/MontePython.py --help
    $ python montepython/MontePython.py run --help
    $ python montepython/MontePython.py info --help

To see a list of all commands. For the `run` subcommand, there are two
essential ones, without which the program will not start. At minimum, you
should precise an output folder (`-o`) and a parameter file (`-p`). An example
of parameter file is found in the main directory of MontePython (`test.param`,
for instance).

A typical call would then be:

.. code::

    $ python montepython/MontePython.py run -o test -p example.param

If non existent, the `test/` folder will be created, and a run with the number
of steps described in `example.param` will be started. To run a chain with more
steps, one can type:

.. code::

    $ python montepython/MontePython.py run -o test -p example.param -N 100

If you want to analyse the run, then just type

.. code::

    $ python montepython/MontePython.py info test/

Note that you probably want more than a hundred points before analyzing a
folder.

Details and Examples
--------------------

Please refer to the pdf or online documentation for further details. The `wiki
<https://github.com/baudren/montepython_public/wiki>`_ contains additional
details on installation. The `forum
<https://github.com/baudren/montepython_public/issues>`_ also contains a
collection of already answered questions, and can be used to discuss the code.


Bibtex entry
------------

When using *Monte Python* in a publication, please acknowledge the code by citing
the following paper. If you used *Class*, *Nested Sampling* or *Cosmo Hammer*,
you should also cite the original works.

.. code::

    @article{Audren:2012wb,
          author         = "Audren, Benjamin and Lesgourgues, Julien and Benabed,
                            Karim and Prunet, Simon",
          title          = "{Conservative Constraints on Early Cosmology: an
                            illustration of the Monte Python cosmological parameter
                            inference code}",
          journal        = "JCAP",
          volume         = "1302",
          pages          = "001",
          doi            = "10.1088/1475-7516/2013/02/001",
          year           = "2013",
          eprint         = "1210.7183",
          archivePrefix  = "arXiv",
          primaryClass   = "astro-ph.CO",
          reportNumber   = "CERN-PH-TH-2012-290, LAPTH-048-12",
          SLACcitation   = "%%CITATION = ARXIV:1210.7183;%%",
    }

