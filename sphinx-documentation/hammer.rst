Using the Cosmo Hammer with Monte Python
========================================

Monte Python can now use the software `Cosmo Hammer
<http://www.astro.ethz.ch/refregier/research/Software/cosmohammer/>`_ written
by J.  Akeret and S. Seehars, which is based on the `emcee sampler
<http://arxiv.org/abs/1202.3665>`_, itself based on the `Affine Invariant
Markov Chain Monte Carlo <http://msp.berkeley.edu/camcos/2010/5-1/p04.xhtml>`_

.. note::

    By using the Cosmo Hammer, you agree to abide by the GNU General Public
    License v3.0 or higher (see their website)


Please look at their website for specifics about the installation.
    

Using with Monte Python
-----------------------

you can choose to use the Cosmo Hammer by specifying the argument: `-m CH`

You should probably always set the environment variable OMP_NUM_THREADS to your
maximum number of cores:

.. code::

    $] export OMP_NUM_THREADS=4

before running.
