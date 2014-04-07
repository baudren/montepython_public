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

    Parameters
    ----------
        custom_command: str
            allows for testing the code
    """
    # Create all the instances of the needed classes to run the code. The safe
    # initialisation handles the errors.
    cosmo, data, command_line, success = safe_initialisation(
        custom_command)

    # If success is False, it means either that the initialisation was not
    # successful, or that it was simply an analysis call. The run should
    # stop
    if not success:
        return

    # Once that the initialisation phase is done, one can import the
    # sampler
    import sampler

    # Generic sampler call
    sampler.run(cosmo, data, command_line)

    return


def mpi_run(custom_command=""):
    """
    Launch a simple MPI run, with no communication of covariance matrix

    It simply allows the first process to create the folder - so that the
    log.param is properly written. A signal is then send to the other
    processes, that contains the chain number of the parent run.

    In order to be sure to have different chain numbers, it adds the rank of
    the process and the initial job number - this should avoid conflict, but
    can be subject of future improvements
    """

    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    nprocs = comm.Get_size()
    rank = comm.Get_rank()
    # this will be the master cpu. This guy will create - or append - a folder,
    # being sure to be the first to do so.
    if rank == 0:
        # First initialisation
        cosmo, data, command_line, success = safe_initialisation(
            custom_command, comm, nprocs)

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

        # If a failed message was passed, exit the process
        if suffix == 'failed':
            return

        # Concatenate the rank to the suffix, and not the opposite, this should
        # avoid any conflicting name
        if not custom_command:
            custom_command = " ".join(sys.argv[1:])
        custom_command += " --chain-number %s" % str(rank)+suffix
        cosmo, data, command_line, success = initialise(custom_command)

    import sampler
    sampler.run(cosmo, data, command_line)

    return


def safe_initialisation(custom_command="", comm=None, nprocs=1):
    """
    Wrapper around the init function to handle errors

    KeyWord Arguments
    -----------------
    custom_command : str
        testing purposes
    comm : MPI.Intracomm
        object that helps communicating between the processes
    nprocs : int
        number of processes
    """
    try:
        cosmo, data, command_line, success = initialise(custom_command)
    except io_mp.ConfigurationError as message:
        if comm:
            for index in range(1, nprocs):
                comm.send('failed', dest=index, tag=1)
        print str(message)
        raise io_mp.ConfigurationError(
            "The initialisation was not successful, resulting in a "
            "potentially half created `log.param`. Please see the "
            "above error message. If you run the exact same command, it"
            " will not work. You should solve the problem, and try again.")
    except KeyError:
        if comm:
            for index in range(1, nprocs):
                comm.send('failed', dest=index, tag=1)
        raise io_mp.ConfigurationError(
            "You are running in a folder that was created following "
            "a non-successful initialisation (wrong parameter name, "
            "wrong likelihood, etc...). If you have solved the issue, you "
            "should remove completely the output folder, and try again.")
    return cosmo, data, command_line, success
