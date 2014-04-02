"""
.. module:: run_mp
    :synopsis: Run the code, initialising everything and calling the sampler
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>

"""
from initialise import initialise
import sys
import io_mp
import re


def run(custom_command=''):
    """
    Main call of the function

    It recovers the initialised instances of cosmo Class, :class:`Data` and the
    NameSpace containing the command line arguments, feeding into the sampler.

    .. note::
        A possible parallelization would take place here.

    Parameters
    ----------
        custom_command: str
            allows for testing the code
    """
    # MPI stuff
    mpi_asked = False
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        nprocs = comm.Get_size()
        if nprocs > 1:
            mpi_asked = True
    except ImportError:
        pass

    if mpi_asked:
        mpi_run(custom_command)
        return

    # Initialisation routine for a normal openMP run
    else:
        cosmo, data, command_line, success = safe_initialisation(custom_command)

        # If success is False, it means either that the initialisation was not
        # successful, or that it was simply an analysis call. The run should
        # stop
        if not success:
            return

        # Once that the initialisation phase is done, one can import the sampler
        import sampler

        # Generic sampler call
        sampler.run(cosmo, data, command_line)

        return


def mpi_run(custom_command=""):
    """TODO"""

    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    nprocs = comm.Get_size()
    rank = comm.Get_rank()
    # this will be the master cpu. This guy will create - or append - a folder,
    # being sure to be the first to do so.
    if rank == 0:
        # First initialisation
        cosmo, data, command_line, success = safe_initialisation(
            custom_command)

        if not success:
            return

        regexp = re.match(".*__(\w*).txt", data.out_name)
        suffix = regexp.groups()[0]
        # Send an "OK" signal to all the other processes, actually giving the
        # suffix of this master chain. All the other will add 1 to this number
        for index in range(1, nprocs):
            comm.send(suffix, dest=index, tag=1)
    else:
        # If the rank is not 0, it is a slave process. It waits to receive the
        # "OK" message, which is immediatly discarded.
        suffix = comm.recv(source=0, tag=1)

        if not custom_command:
            custom_command = " ".join(sys.argv[1:])
        custom_command += " -chain_number %s" % (int(suffix)+rank)
        cosmo, data, command_line, success = safe_initialisation(
            custom_command)

    import sampler
    sampler.run(cosmo, data, command_line)

    return

def safe_initialisation(custom_command=""):
    """Wrapper around the init function to handle errors"""
    try:
        cosmo, data, command_line, success = initialise(custom_command)
    except io_mp.ConfigurationError as message:
        print str(message)
        raise io_mp.ConfigurationError(
            "The initialisation was not successful, resulting in a "
            "potentially half created `log.param`. Please see the "
            "above error message. If you run the exact same command, it"
            " will not work. You should solve the problem, and try again.")
    except KeyError:
        raise io_mp.ConfigurationError(
            "You are running in a folder that was created following "
            "a non-successful initialisation (wrong parameter name, "
            "wrong likelihood, etc...). If you have solved the issue, you "
            "should remove completely the output folder, and try again.")
    return cosmo, data, command_line, success
