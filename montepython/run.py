"""
.. module:: run_mp
    :synopsis: Run the code, initialising everything and calling the sampler
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>

"""
from initialise import initialise
import io_mp


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
    try:
        cosmo, data, command_line, success = initialise(custom_command)
    except io_mp.ConfigurationError as e:
        print str(e)
        raise io_mp.ConfigurationError(
            "The initialisation was not successful, resulting in a "
            "potentially half created `log.param`. Please see the "
            "above error message. If you run the exact same command, it will "
            "not work. You should solve the problem, and try again.")
    except KeyError:
        raise io_mp.ConfigurationError(
            "You are running in a folder that was created following "
            "a non-successful initialisation (wrong parameter name, "
            "wrong likelihood, etc...). If you have solved the issue, you "
            "should remove completely the output folder, and try again.")

    # If success is False, it means either that the initialisation was not
    # successful, or that it was simply an analysis call. The run should stop
    if not success:
        return

    # Once that the initialisation phase is done, one can import the sampler
    import sampler

    # Generic sampler call
    sampler.run(cosmo, data, command_line)

    return
