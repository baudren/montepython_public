"""
.. module:: sampler
    :synopsis: Generic sampler
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>

This module defines one key function, :func:`run`, that distributes the work to
the desired actual sampler (Metropolis Hastings, or Nested Sampling so far)

"""

import io_mp

def run(cosmo, data, command_line):
    """
    First rudimentary implementation

    The old mcmc module is used as previously, except the call is now within
    this function, instead of from within :mod:MontePython

    """

    if command_line.method == 'MH':
        import mcmc
        mcmc.chain(cosmo, data, command_line)
    elif command_line.method == 'NS':
        import nested_sampling as ns
        ns.run(cosmo, data, command_line)
    else:
        io_mp.message(
            "Sampling method %s not understood" % command_line.method,
            "error")
