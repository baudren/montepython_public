.. Monte Python documentation master file, created by
   sphinx-quickstart on Thu Mar  7 14:13:29 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Monte Python's documentation!
========================================

The main page lives `here <http://baudren.github.io/montepython.html>`_, from which
you can download the code, see the changelog. The Github page is available
`there <https://github.com/baudren/montepython_public/>`_.

All useful information concerning the installation, some tips on how
to organize the folder, and the complete description of the code
source is found below.

For the list of command line arguments, please see the documentation of the
:func:`create parser function <parser_mp.create_parser>`. You can also ask this
same information interactively by asking:

.. code::

  python montepython/MontePython.py -h / --help
  python montepython/MontePython.py run -h
  python montepython/MontePython.py info -h

The first one gives you all the possible modes for running (:code:`run`, or
:code:`info`), while the other two give you the information specific for each
modes. Note that asking for :code:`-h` or :code:`--help` will result in using a
short or long format for the help.

Contents:

.. toctree::
   :maxdepth: 2

   installation
   getting_started
   example
   nested
   hammer
   likelihoods
   documentation


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
