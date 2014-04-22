"""
.. module:: run_mp
    :synopsis: Run the code, initialising everything and calling the sampler
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>

"""
from initialise import initialise
import io_mp
import sys
import warnings
import os
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

        # Check that the run asked is compatible with mpirun.
        failed = False
        if command_line.subparser_name == 'info':
            warnings.warn(
                "Analyzing the chains is not supported in mpirun"
                " so this will run on one core only.")
            failed = True
        elif command_line.method in ['NS', 'CH', 'IS']:
            warnings.warn(
                "The methods NS, CH and IS are not supported in mpirun"
                " so this will run on one core only.")
            failed = True

        if failed:
            suffix = 'failed'
        else:
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
            success = False
        else:
            # Concatenate the rank to the suffix, and not the opposite, this
            # should avoid any conflicting name
            if not custom_command:
                custom_command = " ".join(sys.argv[1:])
            custom_command += " --chain-number %s" % str(int(suffix)+rank)
            cosmo, data, command_line, success = initialise(custom_command)

    if success:
        import sampler
        sampler.run(cosmo, data, command_line)

    success = comm.gather(success, root=0)
    return


def mock_update_run(custom_command=""):
    """
    Tentative covmat update run

    Not reachable yet by any option.
    """
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # store the command_line space
    if not custom_command:
        custom_command = " ".join(sys.argv[1:])

    # Do a first run
    mpi_run(custom_command)

    # Compute the covariance matrix
    if rank == 0:
        info_command = from_run_to_info(custom_command)
        initialise(info_command)

    # Make sure that the covariance matrix used is the one just computed
    with_covmat = add_covariance_matrix(custom_command)

    mpi_run(with_covmat)
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


def from_run_to_info(command):
    """
    Translate a command corresponding to a run into one for analysis

    """
    original = command.split()
    new = ['info']

    # Recover the folder
    index = original.index('-o')
    folder = original[index+1]

    # Recover chain number
    index = original.index('-N')
    number = original[index+1]

    # Recover all relevant chains
    chains = [os.path.join(folder, elem) for elem in os.listdir(folder)
              if elem.find(number) != -1]
    new.extend(chains)

    # Do not plot the pdf
    new.append("--noplot")
    return " ".join(new)


def add_covariance_matrix(command):
    """
    Make sure that the command uses the covariance matrix from the folder

    """
    original = command.split()

    # recover the folder
    index = original.index('-o')
    folder = original[index+1]
    # get the name of the folder
    name = folder.split(os.path.sep)[-1]
    covname = os.path.join(
        folder, name+'.covmat')
    bfname = os.path.join(
        folder, name+'.bestfit')

    if '-c' in original:
        index = original.index('-c')
        original[index+1] = covname
    else:
        original.extend(['-c', covname])

    if '-b' in original:
        index = original.index('-b')
        original[index+1] = bfname
    elif '--bestfit' in original:
        index = original.index('--bestfit')
        original[index+1] = bfname
    else:
        original.extend(['-b', bfname])

    # Change the number of steps asked
    index = original.index('-N')
    number = list(str(original[index+1]))
    number[0] = str(int(number[0])+1)
    newnumber = "".join(number)
    original[index+1] = str(newnumber)
    return " ".join(original)
