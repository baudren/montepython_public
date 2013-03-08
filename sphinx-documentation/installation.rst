Installation
============

Prerequisites
-------------

Python
^^^^^^

First of all, you need a clean installation of Python_ (version 2.x,
x<=7, but <3.0), with at least the numpy_ module (version >=1.4.1) and
the cython_ module. This last one is to convert the C code CLASS_
into a Python class.

If you also want the output plot to have cubic interpolation for
analyzing chains, you should also have the scipy_ module (at least
version 0.9.0). In case this one is badly installed, you will have an
error message when running the analyze module of |MP|, and
obtain only linear interpolation. Though not fatal, this problem
produces ugly plots.

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


.. _Python: http://www.python.org/
.. _numpy: http://www.numpy.org/
.. _cython: http://www.cython.org/
.. _scipy: http://www.scipy.org/
.. _argparse: https://pypi.python.org/pypi/argparse
.. _ordereddict: http://code.activestate.com/recipes/576693/
.. _CLASS: http://www.class-code.net/
.. |MP| replace:: **Monte Python**
