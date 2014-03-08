"""
.. module:: run_mp
    :synopsis: Run the code, initialising everything and calling the sampler
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>

"""
import sampler     # generic sampler that calls different
                                    # sampling algorithms
from initialise import initialise


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
    # Initialisation routine
    cosmo, data, command_line, success = initialise(custom_command)

    # If success is False, it means either that the initialisation was not
    # successful, or that it was simply an analysis call. The run should stop
    if not success:
        return

    # Generic sampler call
    sampler.run(cosmo, data, command_line)

    return
