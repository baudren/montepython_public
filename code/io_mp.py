"""
Input-Output handling

Handles all the input/output of the code (at least most of it). If something is
printed that does not satisfy you (number of decimals, for instance, in the
output files), you only have to find the called function and change a number.

Whenever the arguments of the functions are :code:`command_line` or
:code:`data`, no mention of them will be done - as it is now clear. On the
contrary, if there are more arguments, they will be detailled.

This module also defines a new class :class:`File`, that extends
:py:class:`file`, which provides a tail function. It is used in
:func:`mcmc.read_args_from_chain`.

Finally, the way the error messages are displayed is set there, along with
ascii-art for the exclamation mark sign.
"""

import os
import sys
import re  # Module to handle regular expressions
import random as rd
import numpy as np
try:
    from collections import OrderedDict as od
except:
    from ordereddict import OrderedDict as od
from datetime import date
import fcntl
import textwrap  # used to format the error messages

# Ascii art for error display
start_line = {}
start_line['error'] = [' /|\   ',
                       '/_o_\  ',
                       '       ']
start_line['warning'] = [' /!\ ',
                         '     ']
start_line['info'] = [' /!\ ',
                      '     ']

standard_length = 80  # standard, increase if you have a big screen

def log_parameters(data, command_line):
    """
    Write the first lines of the log.param

    Writes the beginning of log.param, starting with the header with the
    cosmological code version and subversion, and then recopies entirely the
    input parameter file.

    """
    log = open(command_line.folder+'/log.param', 'w')
    param_file = open(command_line.param, 'r')
    log.write("#-----{0} {1} (subversion {2})-----\n\n".format(
        data.cosmological_module_name, data.version, data.subversion))
    for line in param_file:
        log.write(line)
    param_file.close()
    log.close()


def log_likelihood_parameters(likelihood, command_line):
    """Write down the .data file of the input likelihood to log.param"""
    log = open(command_line.folder+'log.param', 'a')
    tolog = open(likelihood.path, 'r')
    log.write("\n\n#-----Likelihood-{0}-----\n".format(likelihood.name))
    for line in tolog:
        log.write(line)
    tolog.seek(0)
    tolog.close()
    log.close()


def log_cosmo_arguments(data, command_line):
    """
    Write down the `cosmo_arguments` used to log.param

    Third function called when writing log.param. It is understood here that
    all the other parameters for the cosmological modules are set to their
    default value directly in the program.

    It is written as an update for the dictionary cosmo_arguments (i.e. as
    :code:`dict.update()` and not as :code:`dict =`) in order not to erase
    previously initialized data.

    """
    if len(data.cosmo_arguments) >= 1:
        log = open(command_line.folder+'/log.param', 'a')
        log.write('\n\n#-----------Cosmological-arguments---------\n')
        log.write('data.cosmo_arguments.update({0})\n'.format(
            data.cosmo_arguments))
        log.close()


def log_default_configuration(data, command_line):
    """
    Log the .conf file to log.param

    Fourth and last function called when writing log.param. Only useful if you
    have several versions of your cosmological code installed in different
    locations, or different versions of Clik. But, as you never know what might
    go wrong, it is logged everytime !

    """
    log = open(command_line.folder+'/log.param', 'a')
    log.write('\n\n#--------Default-Configuration------\n')
    for key, value in data.path.iteritems():
        log.write("data.path['{0}']\t= '{1}'\n".format(key, value))
    log.close()


def print_parameters(out, data):
    """
    Will print the parameter names. In the code, :code:`out` is simply the
    standard output, as this information will not be printed on the output
    file.

    Indeed, you will be able to recover these information from the log.param.

    .. warning::

        Please pay attention to the fact that, once launched, the order of the
        parameters in log.param is crucial, as is it the only place where it is
        stored.

    """
    param = data.get_mcmc_parameters(['varying'])
    for elem in data.get_mcmc_parameters(['derived']):
        param.append(elem)
    out.write('\n#  -LogLkl\t')
    for i in range(len(param)):
        if data.mcmc_parameters[param[i]]['scale'] != 1:
            number = data.mcmc_parameters[param[i]]['scale']
            if (number > 100. or number < 0.01):
                string = '%0.e%s' % (1./number, param[i])
            else:
                string = '%0.2g%s' % (1./number, param[i])
        else:
            string = '%s' % param[i]
        out.write("%-16s" % string)
    out.write('\n')


def print_vector(out, N, loglkl, data):
    """
    Print the last accepted values to :code:`out`

    :Parameters:
        - **out** (`list`) - array containing both standard output and the output file.

          This way, if you run in interactive mode, you will be able to monitor
          the progress of the chain.
        - **N** (`int`) - multiplicity of the point, `i.e.` number of times the
          code stayed at this particular place.
        - **loglkl** (`float`) - value of the log likelihood at this point

    .. note::

        It is the `last_accepted` point that is printed, and **not** the
        `current` one (obviously, as one does not know yet the multiplicity of
        the current one !)

    """

    for j in range(len(out)):
        out[j].write('%d  %.6g\t' % (N, -loglkl))
        for elem in data.get_mcmc_parameters(['varying']):
            out[j].write('%.6e\t' %
                         data.mcmc_parameters[elem]['last_accepted'])
        for elem in data.get_mcmc_parameters(['derived']):
            out[j].write('%.6e\t' %
                         data.mcmc_parameters[elem]['last_accepted'])
        out[j].write('\n')


def refresh_file(data):
    """
    Closes and reopen the output file to write any buffered quantities

    """
    data.out.close()
    data.out = open(data.out_name, 'a')


def create_output_files(command_line, data):
    """
    Automatically create a new name for the chain.

    This routine takes care of organising the folder for you. It will
    automatically generate names for the new chains according to the date, number
    of points chosen.

    .. warning::

        The way these names are generated (with the proper number of _, __, -, and
        their placement) is exploited in the rest of the code in various places.
        Please keep that in mind if ever you are in the mood of changing things here.

    """
    if command_line.restart is None:
        number = command_line.N
    else:
        number = int(command_line.restart.split('/')[-1].split('__')[0].
                     split('_')[1]) + command_line.N

    # output file
    outname_base = '{0}_{1}__'.format(date.today(), number)
    suffix = 0
    Try = True
    if command_line.chain_number is None:
        for files in os.listdir(command_line.folder):
            if files.find(outname_base) != -1:
                if int(files.split('__')[-1].split('.')[0]) > suffix:
                    suffix = int(files.split('__')[-1].split('.')[0])
        suffix += 1
        while Try:
            data.out = open(command_line.folder+outname_base +
                            str(suffix)+'.txt', 'w')
            try:
                lock(data.out, fcntl.LOCK_EX | fcntl.LOCK_NB)
                Try = False
            except LockException:
                suffix += 1
        sys.stdout.write('Creating {0}{1}{2}.txt\n'.format(
            command_line.folder, outname_base, suffix))
        data.out_name = '{0}{1}{2}.txt'.format(
            command_line.folder, outname_base, suffix)
    else:
        data.out = open(command_line.folder+outname_base +
                        command_line.chain_number+'.txt', 'w')
        sys.stdout.write('Creating {0}{1}{2}.txt\n'.format(
            command_line.folder, outname_base, command_line.chain_number))
        data.out_name = '{0}{1}{2}.txt'.format(
            command_line.folder, outname_base, command_line.chain_number)
    # in case of a restart, copying the whole thing in the new file
    if command_line.restart is not None:
        for line in open(command_line.restart, 'r'):
            data.out.write(line)


def get_tex_name(name, number=1):
    """
    Simplistic tex name transformer.

    Essentially tries to add a backslash in front of known possible greek
    letters, and insert curly brackets { } around statement following an
    _ or a ^. It will also try to include the scale into the name in a nice
    way.

    .. note::

        This might easily fail on simple names, like `beta_plus_lambda`. In
        this case, please use an extra plot file with the command line option
        :code:`-extra plot_file`, or come up with a better function !

    :Parameters:
        - **name** (`str`) - input name

    :Keywords:
        - **number** (`float`) - scale
    """
    tex_greek = ['omega', 'tau', 'alpha', 'beta', 'delta', 'nu',
                 'Omega', 'Lambda', 'lambda']
    for elem in tex_greek:
        if elem in name:
            position = name.find(elem)
            name = name[:position]+"""\\"""+name[position:]
    if name.find('_') != -1:
        temp_name = name.split('_')[0]+'_{'
        for i in range(1, len(name.split('_'))):
            temp_name += name.split('_')[i]+' '
        name = temp_name + '}'
    if number == 1:
        name = "${0}$".format(name)
    elif (number < 1000 and number > 1):
        name = "$%0.d~%s$" % (number, name)
    else:
        temp_name = "$%0.e%s$" % (number, name)
        m = re.search(r'(?:\$[0-9]*e\+[0]*)([0-9]*)(.*)', temp_name)
        sign = '+'
        if m is None:
            m = re.search(r'(?:\$[0-9]*e\-[0]*)([0-9]*)(.*)', temp_name)
            sign = '-'
        name = '$10^{'+sign+m.groups()[0]+'}'+m.groups()[1]
    return name

def message(string, status):
    """
    Outputs the string formatted according to its status

    The input is a potentially long message, describing the problem.
    According to the severity of its status (so far, 'error' will exit the
    program, whereas 'warning' and 'info' will go through anyway).

    Standard length has been defined globally, as well as the ascii-art
    dictionary of arrays start_line.

    """

    length = standard_length-len(start_line[status][0])
    # Remove unwanted spaces (coming from carriage returns in the input string)
    # and handle voluntary carriage returns specified with \n
    first_cleanup = [' '.join(elem.lstrip(' ').split()) for elem in string.split('\n')]
    splitted = []
    # Recover the lines splitted at correct length
    for elem in first_cleanup:
        splitted.extend(textwrap.wrap(elem, length))
    
    if status == 'error':
        # Add a blank line so that the error displays better
        print

    # Add in front the appropriate fancy display
    index = 0
    for line in splitted:
        # If the number of needed lines is bigger than the ascii-art, the last
        # line of it (empty) will be used.
        if index < len(start_line[status]): 
            start_index = index
        else: 
            start_index = len(start_line[status])-1
        print start_line[status][start_index]+line 
        index += 1
    if status == 'error':
        # In case of a severe error, the program should stop the execution
        exit()


class File(file):
    """
    New class of file, to provide an equivalent of the tail command (on linux).
    It will be used when starting from an existing chain, and avoids circling
    through an immense file.
    """

    def tail(self, lines_2find=1):
        self.seek(0, 2)     # go to end of file
        bytes_in_file = self.tell()
        lines_found, total_bytes_scanned = 0, 0
        while (lines_2find+1 > lines_found and
                bytes_in_file > total_bytes_scanned):
            byte_block = min(1024, bytes_in_file-total_bytes_scanned)
            self.seek(-(byte_block+total_bytes_scanned), 2)
            total_bytes_scanned += byte_block
            lines_found += self.read(1024).count('\n')
        self.seek(-total_bytes_scanned, 2)
        line_list = list(self.readlines())
        return line_list[-lines_2find:]


class LockException(Exception):
    """
    .. warning::

        in the process of being tested
    """
    # Error codes:
    LOCK_FAILED = 1


def lock(file, flags):
    """
    Lock a given file to prevent other instances of the code to write to the
    same file.

    .. warning::

        in the process of being tested

    """
    import fcntl
    try:
        fcntl.flock(file.fileno(), flags)
    except IOError, exc_value:
        # The exception code varies on different systems so we'll catch
        # every IO error
        raise LockException(*exc_value)


def unlock(file):
    """
    Unlock a previously locked file.

    .. warning::

        in the process of being tested

    """
    import fcntl
    fcntl.flock(file.fileno(), fcntl.LOCK_UN)
