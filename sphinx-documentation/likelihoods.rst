Existing likelihoods, and how to create new ones
================================================


This page is intended to explain in more concrete terms the
information contained in the :doc:`likelihood_class` documentation.
More specifically, you should be able to write new likelihood files
and understand the structure of existing ones.

One likelihood is one directory, one :code:`.py` and one :code:`.data` file
----------------------------------------------------------------------------

We have seen already that cosmological parameters are passed directly from the
input file to *Class*, and do not appear anywhere in the code itself, i.e. in
the files located in the :code:`montepython/` directory. The situation is the same for
likelihoods. You can write the name of a likelihood in the input file, and
|MP| will directly call one of the external likelihood codes implemented in the
:code:`montepython/likelihoods/` directory. This means that when you add some new
likelihoods, you don't need to declare them in the code. You implement them in
the :code:`likelihoods` directory, and they are ready to be used if mentioned in
the input file.

For This to work, a precise syntax must be respected. Each likelihood
is associated to a name, e.g. :code:`hst`,  :code:`wmap`,
:code:`WiggleZ` (the name is case-sensitive). This name is used:

* for calling the likelihood in the input file, e.g. :code:`data.experiments = ['hst', ...]`,
* for naming the directory of the likelihood, e.g. :code:`montepython/likelihoods/hst/`,
* for naming the input data file describing the characteristics of the
  experiment,  :code:`montepython/likelihoods/hst/hst.data` (this file can point
  to raw data files located in the :code:`data` directory)
* for naming the class declared in :code:`montepython/likelihoods/hst/__init__.py` and used also in :code:`montepython/likelihoods/hst/hst.data`

.. warning::

    Note that since release 2.0.0, the likelihood python source is not called
    any longer :code:`hst.py`, but :code:`__init__.py`. The reason was for
    packaging and ease of use when calling from a Python console.

When  implementing new likelihoods, you will have to follow this rule. You
could already wish to have two Hubble priors/likelihoods in your folder. For
instance, the distributed version of :code:`hst` corresponds to a gaussian prior
with standard deviation :math:`h=0.738\pm0.024`. If you want to change these numbers,
you can simply edit :code:`montepython/likelihoods/hst/hst.data`. But you could
also keep :code:`hst` unchanged and create a new likelihood called
e.g. :code:`spitzer`. We will come back to the creation of likelihoods
later, but just to illustrate the structure of likelihoods, let us see
how to create such a prior/likelihood:

.. code::

  $ mkdir likelihoods/spitzer
  $ cp likelihoods/hst/hst.data likelihoods/spitzer/spitzer.data
  $ cp likelihoods/hst/__init__.py likelihoods/spitzer/__init__.py

Then edit :code:`montepython/likelihoods/spitzer/__init__.py` and replace in the initial declaration the class name :code:`hst` by :code:`spitzer`:

.. code::

  class spitzer(Likelihood_prior):

Edit also :code:`montepython/likelihoods/spitzer/spitzer.data`, replace the class name :code:`hst` by :code:`spitzer`, and the numbers by your constraint:

.. code::

  spitzer.h = 0.743
  spitzer.sigma = 0.021

You are done. You can simply add :code:`data.experiments = [...,'spitzer', ...]` to the list of experiments in the input parameter file and the likelihood will be used.



Existing likelihoods
--------------------

We release the first version of |MP| with the likelihoods:

* :code:`spt`, :code:`bicep`, :code:`cbi`, :code:`acbar`,
  :code:`bicep`, :code:`quad`, the latest public versions of CMB data
  from SPT, Bicep, CBI, ACBAR, BICEP and Quad; for the SPT likelihoods
  we include three nuisance parameters obeying to gaussian priors,
  like in the original SPT paper, and for ACBAR one nuisance parameter
  with top-hat prior. These experiments are described by the very same
  files as in a ComsoMC implementation. They are located in the
  :code:`data/` directory. For each experiment, there is a master file
  :code:`xxx.dataset` containing several variables and the names of
  other files with the raw data. In the files
  :code:`likelihoods/xxx/xxx.data`, we just give the name of the
  different :code:`xxx.dataset` files, that |MP| is able to read just
  like CosmoMC.
* :code:`wmap`, original likelihood file accessed through the wmap
  wrapper. The file :code:`likelihoods/wmap/wmap.data` allows you to
  call this likelihood with a few different options (e.g. switching
  on/off Gibbs sampling, choosing the minimum and maximum multipoles
  to include, etc.) As usual, we implemented the nuisance parameter
  :code:`A_SZ` with a flat prior. In the input parameter file, you can
  decide to vary this parameter in the range 0-2, or to fix it to some
  value.
* :code:`hst` is the HST Key Project gaussian prior on :math:`h`,
* :code:`sn` constains the luminosity distance-redhsift relation using
  the Union 2 data compilation,
* :code:`WiggleZ` constraints the matter power spectrum :math:`P(k)`
  in four different redshift bins using recent WiggleZ data,

plus a few other likelihoods referring to future experiments,
described in the next subsection. All these likelihoods are strictly
equivalent to those in the CosmoMC patches released by the various
experimental collaborations.


Mock data likelihoods
---------------------

We also release simplified likelihoods :code:`fake_planck_bluebook`,
:code:`euclid_lensing` and :code:`euclid_pk` for doing forecasts for
Planck, Euclid (cosmic shear survey) and Euclid (redshift survey).

In the case of Planck, we use a simple gaussian likelihood for TT, TE,
EE (like in `astro-ph/0606227
<http://arxiv.org/abs/astro-ph/0606227>`_ with no lensing extraction)
with sensitivity parameters matching the numbers published in the
Planck bluebook. In the case of Euclid, our likelihoods and
sensitivity parameters are specified in `the Euclid Red Book
<http://arxiv.org/abs/1110.3193>`_. The sensitivity parameters can
always be modified by the user, by simply editing the :code:`.data`
files.

These likelihoods compare theoretical spectra to a fiducial spectrum
(and **not** to random data generated given the fiducial model: this
approach is simpler and leads to the same forecast error bars, see
`this paper again <http://arxiv.org/abs/astro-ph/0606227>`_).

Let us illustrate the way in which this works with
:code:`fake_planck_bluebook`, although the two Euclid likelihoods obey
exactly to the same logic.


When you download the code, the file
:code:`montepython/likelihoods/fake_planck_bluebook/fake_planck_bluebook.data` has
a field :code:`fake_planck_bluebook.fiducial_file` pointing to the
file :code:`'fake_planck_bluebook_fiducial.dat'`. You downloaded this
file together with the code: it is located in :code:`data` and it
contains the TT/TE/EE spectrum of a particular fiducial model (with
parameter values logged in the first line of the file). If you launch
a run with this likelihood, it will work immediately and fit the
various models to this fiducial spectrum.

But you probably wish to choose your own fiducial model. This is extremely
simple with |MP|. You can delete the provided fiducial
file\\:code:`'fake_planck_bluebook_fiducial.dat'`, or alternatively,
you can change the name of the fiducial file in
:code:`likelihoods/fake_planck_bluebook/fake_planck_bluebook.data`. When you
start the next run, the code will notice that there is no input fiducial
spectrum. It will then generate one automatically, write it in the correct file
with the correct location, and stop after this single step. Then, you can
launch new chains, they will fit this fiducial spectrum. 

When you generate the fiducial model, you probably want to control
exactly fiducial parameter values. If you start from an ordinary input
file with no particular options, |MP| will perform one random jump and
generate the fiducial model. Fiducial parameter values will be logged
in the first line of the fiducial file. But you did not choose them
yourself. However, when you call |MP| with the intention of generating
a fiducial spectrum, you can pass the command line option :code:`-f
0`. This sets the variance of the proposal density to zero. Hence the
fiducial model will have precisely the parameter values specified in
the input parameter file. The fiducial file is even logged in the
:code:`log.param` of all the runs that have been using it.


Creating new likelihoods belonging to pre-defined category
----------------------------------------------------------

A likelihood is a class (let's call it generically :code:`xxx`), declared and
defined in :code:`montepython/likelihoods/xxx/__init__.py`, using input numbers and input files
names specified in :code:`montepython/likelihoods/xxx/xxx.data`. The actual data files
should usually be placed in the :code:`data/` folder (with the exception of WMAP
data). Such a class will always inherit from the properties of the most generic
class defined inside :code:`montepython/likelihoods_class.py`. But it may fall in the
category of some pre-defined likelihoods and inherit more properties. In this
case the coding will be extremely simple, you won't need to write a specific
likelihood code.

In the current version, pre-defined classes are:

:class:`Likelihood_newdat <likelihood_class.Likelihood_newdat>` 
  suited for all CMB experiments described by a file in the
  :code:`.newdat` format (same files as in CosmoMC).
:class:`Likelihood_mock_cmb <likelihood_class.Likelihood_mock_cmb>`
  suited for all CMB experiments dexcribed with a simplified gaussian
  likelihood, like our :code:`fake_planck_bluebook` likelihood.
:class:`Likelihood_mpk <likelihood_class.Likelihood_mpk>`
  suited for matter power spectrum data that would be described with a
  :code:`.dataset` file in CosmoMC. This generic likelihood contains a
  piece of code following closely the routine :code:`mpk` developped
  for CosmoMC. In the released version of |MP|, this likelihood type
  is only used by each of the four redshift bins of the WiggleZ data,
  but it is almost ready for being used with other data set in this
  format.

Suppose, for instance, that a new CMB dataset :code:`nextcmb` is
released in the :code:`.newdat` format. You will then copy the
:code:`.newdat` file and other related files (with window functions,
etc.) in the folder :code:`data/`. You will then create a new
likelihood, starting from an existing one, e.g cbi:

.. code::

  $ mkdir likelihoods/nextcmb
  $ cp likelihoods/cbi/cbi.data likelihoods/nextcmb/nextcmb.data
  $ cp likelihoods/cbi/__init__.py likelihoods/nextcmb/__init__.py

The python file should only be there to tell the code that nextcmb is
in the :code:`.newdat` format. Hence it should only contain:

.. code::

  from montepython.likelihood_class import Likelihood_newdat 
  class nextcmb(Likelihood_newdat):
      pass

This is enough: the likelihood is fully defined. The data file should
only contain the name of the :code:`.newdat` file:

.. code::

  nextcmb.data_directory  = data.path['data']
  nextcmb.file            = 'next-cmb-file.newdat'

Once you have edited these few lines, you are done! No need to tell
|MP| that there is a new likelihood! Just call it in your next run by
adding :code:`data.experiments = [...,'nextcmb', ...]` to the list of
experiments in the input parameter file, and the likelihood will be
used.

You can also define nuisance parameters, contamination spectra and
nuisance priors for this likelihood, as explained in the next section.


Creating new likelihoods from scratch
-------------------------------------

The likelihood :code:`sn` is an example of individual likelihood code:
the actual code is explicitly written in :code:`sn.py`. To create your
own likelihood files, the best to is look at such examples and follow
them. We do not provide a full tutorial here, and encourage you to ask
for help if needed. Here are however some general indications.

Your customised likelihood should inherit from generic likelihood
properties through:

.. code::

  from montepython.likelihood_class import Likelihood
  class my-likelihood(Likelihood):


Implementing the likelihood amounts in developing in the python file
:code:`my-likelihood.py` the properties of two essential functions,
:code:`__init__` and :code:`loglkl`. But you don't need to code
everything from scratch, because the generic :class:`likelihood
<likelihood_class.Likelihood>` already knows the most generic steps.
The previous link will give you all the functions defined from this
base class, that your daughter class will inherit from. Here follows a
detailled explanation about how to use these.

One thing is that you don't need to write from scratch the parser
reading the :code:`.data` file: this will be done automatically at the
beginning of the initialization of your likelihood. Consider that any
field defined with a line in the :code:`.data` file, e.g.
:code:`my-likelihood.variance = 5`, are known in the likelihood code:
in this example you could write in the python code something like
:code:`chi2+=result**2/self.variance`.

You don't need either to write from scratch an interface with *Class*.
You just need to write somewhere in the initialization function  some
specific parameters that should be passed to *Class*. For instance, if
you need the matter power spectrum, write

.. code::

  self.need_cosmo_arguments(data,{'output':'mPk'})

that uses the method :func:`need_cosmo_arguments
<likelihood_class.Likelihood.need_cosmo_arguments>`. If this
likelihood is used, the field :code:`mPk` will be appended to the list
of output fields (e.g. :code:`output=tCl,pCl,mPk`), unless it was
already there. If you write

.. code::

  self.need_cosmo_arguments(data,{'l_max_scalars':3300})

the code will check if :code:`l_max_scalars` was already set at least
to 3300, and if not, it will increase it to 3300. But if another
likelihood needs more it will be more.
   
You don't need to redefine functions like for instance those defining
the role of nuisance parameters (especially for CMB experiments).  If
you write in the :code:`.data` file

.. code::

  my-likelihood.use_nuisance           = ['N1','N2'] 

the code will know that this likelihood cannot work if these two
nuisance parameters are not specified  in the parameter input file
(they can be varying or fixed; fix them by writing a 0 in the sigma
entry). If you try to run without them, the code will stop with an
explicit error message.  If the parameter :code:`N1` has a top-hat
prior, no need to write it: just specify prior edges in the input
parameter file. If :code:`N2` has a gaussian prior, specify it in the
:code:`.data` file, e.g.:

.. code::

  my-likelihood.N2_prior_center  = 1
  my-likelihood.N2_prior_variance = 2  

Since these fields refer to pre-defined properties of the likelihood,
you don't need to write explicitly in the code something like
:code:`chi2 += (N2-center)**2/variance`, adding the prior is done
automatically. Finally, if these nuisance parameters are associated to
a CMB dataset, they may stand for a multiplicative factor in front of
a contamination spectrum to be added to the theoretical
:math:`C_{\ell}`'s. This is the case for the nuisance parameters of
the :code:`acbar`, :code:`spt` and :code:`wmap` likelihoods delivered
with the code, so you can look there for concrete examples. To assign
this role to these nuisance parameters, you just need to write

.. code::

  my-likelihood.N1_file = 'contamination_corresponding_to_N1.data'

and the code will understand what it should do with the parameter
:code:`N1` and the file
:code:`data/contamination_corresponding_to_N1.data`. Optionally, the
factor in front of the contamination spectrum can be rescaled by a
constant number using the syntax:

.. code::

  my-likelihood.N1_scale = 0.5

Creating new likelihoods requires a basic knowledge of python. If you
are new in python, once you know the basics, you will realise how
concise a code can be. You can compare the length of the likelihood
codes that we provide with their equivalent in Fortran in the CosmoMC
package.

.. |MP| replace:: *Monte Python*
