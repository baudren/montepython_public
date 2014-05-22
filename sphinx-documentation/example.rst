Example of a complete work session
==================================


I just downloaded and installed |MP|, read the previous pages, and I
wish to launch and analyse my first run.

I can first create a few folders in order to keep my :code:`montepython` directory tidy in the future. I do a

:code:`$ mkdir chains` 
  for storing all my chains
:code:`$ mkdir chains/planck` 
  if the first run I want to launch is based on the fake planck likelihood proposed in the :code:`example.param` file
:code:`$ mkdir input`  
  for storing all my input files
:code:`$ mkdir scripts` 
  for storing all my scripts for running the code in batch mode

I then copy :code:`example.param` in my input folder, with a name of my choice, e.g. :code:`lcdm.param`, and edit it if needed:

.. code::

  $ cp example.param input/lcdm.param

I then launch a short chain with

.. code::

  $ montepython/Montepython.py run -p input/lcdm.param -o chains/planck/lcdm -N 5

I can see on the screen the evolution of the initialization of the
code. At the end I check that I have a chain and a :code:`log.param`
written in my :code:`chains/planck/lcdm/log.param` directory. I can
immediately repeat the experience with the same command. The second
chain is automatically created with number 2 instead of 1. I can also
run again without the input file:

.. code::

  $ montepython/Montepython.py run -o chains/planck/lcdm -N 5

This works equally well because all information is taken from the :code:`log.param` file.

In some cases, initally, I don't have a covariance matrix to pass in
input [#f1]_ . But in this particular example I can try
the one delivered with the |MP| package, in the :code:`covmat/` directory:

.. code::

  $ montepython/Montepython.py run -p input/lcdm.param\
          -o chains/planck/lcdm -c covmat/fake_planck_lcdm.covmat -N 5


I don't have yet a covariance matrix to pass in input, otherwise I
would have run with 

.. code::

  $ montepython/Montepython.py run -p input/lcdm.param -o chains/planck/lcdm -c mycovmat.covmat -N 5


I now wish to launch longer runs on my cluster or on a powerful desktop.
The syntax of the script depends on the cluster. In the simplest case it
will only contain some general commands concerning the job name, wall time
limit etc., and the command line above (I can use the one without input
file, provided that I made already one short interactive run, and that the
:code:`log.param` already exists; but I can now increase the number of
steps, e.g. to 5000 or 10000). On some cluster, the chain file is created
immediately in the output directory at start up. In this case, the
automatic numbering of chains proposed by |MP| will be satisfactory.

.. warning::

  On some clusters, the automatic numbering will conflict when the chains
  are created too fast. Please look at the section on how to use
  :code:`mpi_run` for guidance

In other clusters, the chains are created on a temporary file, and
then copied at the end to the output file. In this case, if I do
nothing, there is a risk that chain names are identical and clash. I
should then relate the chain name to the job number, with an
additional command line :code:`--chain_number $JOBID`. Some clusters,
:code:`$JOBID` is a string, but the job number can be extracted with a
line like :code:`export JOBNUM="$(echo $PBS_JOBID|cut -d'.' -f1)"`,
and passed to |MP| as  :code:`--chain_number $JOBNUM`.

If I use in a future run the Planck likelihood, I should not forget to
add in the script (before calling |MP|) the line

.. code::

  source /path/to/my/plc/bin/clik_profile.sh


I then launch a chain by submitting the script, with e.g. :code:`qsub
scripts/lcdm.sh`. I can launch many chains in one command with

.. code::

  $ for i in {1..10}; do qsub scripts/lcdm.sh;done

If you cluster creates the chains too fast, there might be conflicts in the
chain names. One way to go around this issue is to run with :code:`mpi`, which is
a parallelization process. The chains will be initialised one after the
other, each one sending a **go** signal to the next in line.

To launch a job with :code:`mpi`, the syntax is exactly the same than
without, except that you will start the whole command with, depending on
your installation, :code:`mpirun` or :code:`mpiexec`:

.. code::

  mpirun -np 4 python montepython/MontePython.py run -o chains/...

will simply launch 4 chains, each using the environment variable
:code:`$OMP_NUM_THREADS` for the number of cores to compute *Class*.

When the runs have stopped, I can analyse them with

.. code::

  $ montepython/Montepython.py info chains/planck/lcdm 

If I had been running without a covariance matrix, the results would probably
be bad, with a very low acceptance rate and few points. It would have however
created a covariance matrix :code:`chains/planck/lcdm/lcdm.covmat`. I can decide
to copy it in order to keep track of it even after analysing future runs, 

.. code::

  cp chains/planck/lcdm/lcdm.covmat chains/planck/lcdm/lcdm_run1.covmat

I now add to my script, in the line starting with :code:`montepython/Montepyhton.py`, the option 

.. code::

  -c chains/planck/lcdm/lcdm_run1.covmat

run on the same day as the previous one, it might be smart to change also a bit
the number of steps (e.g. from 5000 to 5001) in order to immediately identify
chains belonging to the same run.

When this second run is finished, I analyse it with e.g.

.. code::

  montepython/Montepython.py info chains/planck/lcdm/2012-10-27_5001*

If all R-1 numbers are small (typically :math:`<0.05`) and plots look nice, I am
done. If not, there can be two reasons: the covariance matrix is still bad, or
I just did not get enough samples.

I can check the acceptance rate of this last run by looking at the
:code:`chains/planck/lcdm/lcdm.log` file. If I am in a case with nearly gaussian
posterior (i.e. nearly ellipsoidal contours), an acceptance rate :math:`<0.2` or
:math:`>0.3` can be considered as bad. In other cases, even 0.1 might be the best
that I can expect. If the acceptance rate is bad, I must re-run with an
improved covariance matrix in order to converge quicker. I copy the last
covariance matrix to :code:`lcdm_run2.covmat` and use this one for the next run.
If the acceptance rate is good but the chains are not well converged because
they are simply too short, then I should better rerun with the same covariance
matrix :code:`lcdm_run1.covmat`: in this way, I know that the proposal density
is frozen since the second run, and I can safely analyse the second and third
runs altogether.

If I do two or three runs in that way, I always loose running time, because
each new chain will have a new burn-in phase (i.e. a phase when the log
likelihood is very bad and slowly decreasing towards values close to the
minimum). If this is a concern, I can avoid it in three ways:

* before launching the new run, I set the input mean value of each
  parameter in the input file to the best-fit value found in the previous run.
  The runs will then start from the best-fit value plus or minus the size of
  the first jump drown from the covariance matrix, and avoid burn-in. Since I
  have changed the input file, I must rerun with a new output directory, e.g.
  :code:`chain/lcdm2`. This is a clean method.
* I might prefer a less clean but slightly quicker variant: I modify the
  mean values, like in the previous item, but directly in the :code:`log.param`
  file, and I rerun in the same directory without an input file. This will
  work, but it is advisable not to edit the :code:`log.param` manually, since it
  is supposed to keep all the information from previous runs.
* I may restart the new chains from the previous chains using the :code:`-r`
  command line option. The name of previous chains can be written after
  :code:`-r` manually or through a script.
* I can also restart from the best-fit found previously, using the
  :code:`-bf` command line option, specifying the :code:`.bestfit`
  file to use.

When I am pleased with the final plots and result, I can customize the plot
content and labels by writing a short file :code:`plot_files/lcdm.plot` passed
through the :code:`-extra` command line option, and paste the latex file
produced by |MP| in my paper.

.. |MP| replace:: *Monte Python*

.. rubric:: Footnotes

.. [#f1] If I am also a CosmoMC user, I might have an adequate covmat
  to start with, before using the covmat that |MP| will produce. Fot
  this I just need to edit the first line, add comas between paramater
  names, and for parameter that are identical to those in my run,
  replace CosmoMC parameter names with equivalent *Class* parameter
  names.}
