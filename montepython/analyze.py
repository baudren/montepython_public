"""
.. module:: analyze
   :synopsis: Extract data from chains and produce plots

.. moduleauthor:: Karim Benabed <benabed@iap.fr>
.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>

Collection of functions needed to analyze the Markov chains.

This module defines as well a class :class:`Information`, that stores useful
quantities, and shortens the argument passing between the functions.

.. note::
    Some of the methods used in this module are directly adapted from the
    `CosmoPmc <http://www.cosmopmc.info>`_ code from Kilbinger et. al.

"""
import os
import math
import numpy as np
# The root plotting module, to change options like font sizes, etc...
import matplotlib
# The following line suppresses the need for an X server
matplotlib.use("Agg")
# Module for handling display
import matplotlib.pyplot as plt
# Module to handle warnings from matplotlib
import warnings
import importlib
import io_mp
from debug import timeit

# Defined to remove the burnin for all the points that were produced before the
# first time where log-likelihood >= max-log-likelihood-LOG_LKL_CUTOFF
LOG_LKL_CUTOFF = 3

def analyze(command_line):
    """
    Main function, does the entire analysis.

    It calls in turn all the other routines from this module. To limit the
    arguments of each function to a reasonnable size, a :class:`Information`
    instance is used. This instance is initialized in this function, then
    appended by the other routines.

    """
    # Create an instance of the Information class, that will hold all relevant
    # information, and be used as a compact way of exchanging information
    # between functions
    info = Information()

    # Check if the scipy module has the interpolate method correctly
    # installed (should be the case on every linux distribution with
    # standard numpy)
    try:
        from scipy.interpolate import interp1d
        info.has_interpolate_module = True
    except ImportError:
        info.has_interpolate_module = False
        warnings.warn(
            'No cubic interpolation done (no interpolate method found ' +
            'in scipy), only linear')

    # Save the extension of output files and the bin number
    info.extension = command_line.extension
    info.binnumber = command_line.bins

    # Read a potential file describing changes to be done for the parameter
    # names, and number of paramaters plotted (can be let empty, all will
    # then be plotted).
    if command_line.optional_plot_file is not None:
        execfile(command_line.optional_plot_file)

    # Prepare the files, according to the case, load the log.param, and
    # prepare the output (plots folder, .covmat, .info and .log files).
    # After this step, info.files will contain all chains.
    status = prepare(command_line.files, info)
    # If the preparation step generated new files (for instance, translating
    # from NS or CH to Markov Chains) this routine should stop now.
    if not status:
        return

    # Compute the mean, maximum of likelihood, 1-sigma variance for this
    # main folder. This will create the info.spam chain
    convergence(info)

    # Create the main chain, which consists in all elements of info.spam
    # put together. This will serve for the plotting.
    chain = np.vstack(info.spam)

    # In case of comparison, launch the prepare and convergence methods, for
    # the other folders
    # TODO continue fixing by doing this over an arbitrary number of compared
    # folders
    if command_line.comp is not None:
        comp_info = Information()
        prepare(command_line.comp, comp_info)
        convergence(comp_info)
        # Create comp_chain
        comp_chain = np.vstack(comp_info.spam)

    # Total number of steps, after burnin
    weight = chain[:, 0].sum()

    # Covariance matrix computation (for the whole chain)
    info.mean = info.mean[0]
    info.covar = np.zeros((len(info.ref_names), len(info.ref_names)))

    print '--> Computing covariance matrix'
    for i in xrange(len(info.ref_names)):
        for j in xrange(i, len(info.ref_names)):
            info.covar[i, j] = (
                chain[:, 0]*(
                    (chain[:, i+2]-info.mean[i]) *
                    (chain[:, j+2]-info.mean[j]))).sum()
            if i != j:
                info.covar[j, i] = info.covar[i, j]
    info.covar /= weight

    # Removing scale factors in order to store true parameter covariance
    info.covar = np.dot(info.scales.T, np.dot(info.covar, info.scales))

    # Writing it out in name_of_folder.covmat
    with open(info.cov_path, 'w') as cov:
        cov.write('# %s\n' % ', '.join(
            ['%16s' % name for name in info.backup_names]))

        for i in range(len(info.ref_names)):
            for j in range(len(info.ref_names)):
                if info.covar[i][j] > 0:
                    cov.write(' %.5e\t' % info.covar[i][j])
                else:
                    cov.write('%.5e\t' % info.covar[i][j])
            cov.write('\n')

    # Sorting by likelihood: a will hold the list of indices where the
    # points are sorted with increasing likelihood.
    a = chain[:, 1].argsort(0)

    # Writing the best-fit model in name_of_folder.bestfit
    info.best_fit.write('# ')
    for i in range(len(info.ref_names)):
        string = info.backup_names[i]
        if i != len(info.ref_names)-1:
            string += ','
        info.best_fit.write('%-16s' % string)
    info.best_fit.write('\n')
    # Removing scale factors in order to store true parameter values
    for i in range(len(info.ref_names)):
        bfvalue = chain[a[0], 2+i]*info.scales[i, i]
        if bfvalue > 0:
            info.best_fit.write(' %.5e\t' % bfvalue)
        else:
            info.best_fit.write('%.5e\t' % bfvalue)
    info.best_fit.write('\n')

    # Defining the sigma contours (1, 2 and 3-sigma)
    info.lvls = (68.26, 95.4, 99.7)

    # Computing 1,2 and 3-sigma errors, and plot. This will create the
    # triangle and 1d plot by default.  If you also specified a comparison
    # folder, it will create a versus plot with the 1d comparison of all
    # the common parameters, plus the 1d distibutions for the others.
    info.bounds = np.zeros((len(info.ref_names), len(info.lvls), 2))
    if command_line.plot is True:
        if command_line.comp is None:
            plot_triangle(
                info, chain, command_line,
                bin_number=info.binnumber, levels=info.lvls)
        else:
            plot_triangle(
                info, chain, command_line, bin_number=info.binnumber,
                levels=info.lvls,
                comp_chain=comp_chain,
                comp_ref_names=comp_ref_names,
                comp_tex_names=comp_tex_names,
                comp_backup_names=comp_backup_names,
                comp_plotted_parameters=comp_plotted_parameters,
                comp_folder=comp_folder, comp_boundaries=comp_boundaries,
                comp_mean=comp_mean)

    # Creating the array indices to hold the proper ordered list of plotted
    # parameters
    indices = []
    for i in range(len(info.plotted_parameters)):
        indices.append(info.ref_names.index(info.plotted_parameters[i]))

    print '--> Writing .info and .tex files'
    # Write down to the .h_info file all necessary information
    info.h_info.write(' param names\t:\t')
    info.v_info_names = []
    for i in indices:
        if info.scales[i, i] != 1:
            if (float(info.scales[i, i]) > 100. or
                    (info.scales[i, i]) < 0.01):
                string = ' %0.e%s' % (
                    1./info.scales[i, i], info.ref_names[i])
            elif float(info.scales[i, i]) < 1:
                string = ' %2d%s' % (
                    1./info.scales[i, i], info.ref_names[i])
            else:
                string = ' %2g%s' % (
                    1./info.scales[i, i], info.ref_names[i])
        else:
            string = ' %s' % info.ref_names[i]
        info.v_info_names.append(string)
        info.h_info.write("%-16s" % string)

    write_h(info.h_info, indices, 'R-1 values', '%.6f', info.R)
    write_h(info.h_info, indices, 'Best Fit  ', '%.6e',
            chain[a[0], 2:])
    write_h(info.h_info, indices, 'mean      ', '%.6e', info.mean)
    write_h(info.h_info, indices, 'sigma     ', '%.6e',
            (info.bounds[:, 0, 1]-info.bounds[:, 0, 0])/2.)
    info.h_info.write('\n')
    write_h(info.h_info, indices, '1-sigma - ', '%.6e',
            info.bounds[:, 0, 0])
    write_h(info.h_info, indices, '1-sigma + ', '%.6e',
            info.bounds[:, 0, 1])
    write_h(info.h_info, indices, '2-sigma - ', '%.6e',
            info.bounds[:, 1, 0])
    write_h(info.h_info, indices, '2-sigma + ', '%.6e',
            info.bounds[:, 1, 1])
    write_h(info.h_info, indices, '3-sigma - ', '%.6e',
            info.bounds[:, 2, 0])
    write_h(info.h_info, indices, '3-sigma + ', '%.6e',
            info.bounds[:, 2, 1])
    # bounds
    info.h_info.write('\n')
    write_h(info.h_info, indices, '1-sigma > ', '%.6e',
            info.mean+info.bounds[:, 0, 0])
    write_h(info.h_info, indices, '1-sigma < ', '%.6e',
            info.mean+info.bounds[:, 0, 1])
    write_h(info.h_info, indices, '2-sigma > ', '%.6e',
            info.mean+info.bounds[:, 1, 0])
    write_h(info.h_info, indices, '2-sigma < ', '%.6e',
            info.mean+info.bounds[:, 1, 1])
    write_h(info.h_info, indices, '3-sigma > ', '%.6e',
            info.mean+info.bounds[:, 2, 0])
    write_h(info.h_info, indices, '3-sigma < ', '%.6e',
            info.mean+info.bounds[:, 2, 1])

    info.bestfit = np.zeros(len(info.ref_names))
    for i in xrange(len(info.ref_names)):
        info.bestfit[i] = chain[a[0], :][2+i]

    # Write vertical info file
    info.v_info.write('%-15s\t: %-6s %-10s %-10s %-10s %-11s %-10s %-11s %-10s %-10s %-10s %-10s %-10s' % (
        'param names', 'R-1', 'Best fit', 'mean', 'sigma', '1-sigma -',
        '1-sigma +', '2-sigma -', '2-sigma +', '1-sigma >', '1-sigma <',
        '2-sigma >', '2-sigma <'))
    for i in xrange(len(info.v_info_names)):
        name = info.v_info_names[i]
        #index = info.v_info_names.index(name)
        index = indices[i]
        info.v_info.write('\n%-15s\t: %.4f %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e' % (
            name, info.R[index], chain[a[0], 2:][index],
            info.mean[index],
            (info.bounds[:, 0, 1][index]-info.bounds[:, 0, 0][index])/2.,
            info.bounds[:, 0, 0][index], info.bounds[:, 0, 1][index],
            info.bounds[:, 1, 0][index], info.bounds[:, 1, 1][index],
            info.mean[index]+info.bounds[:, 0, 0][index],
            info.mean[index]+info.bounds[:, 0, 1][index],
            info.mean[index]+info.bounds[:, 1, 0][index],
            info.mean[index]+info.bounds[:, 1, 1][index]))

    # Writing the .tex file that will have this table prepared to be
    # imported in a tex document.
    write_tex(info, indices)


@timeit
def prepare(files, info):
    """
    Scan the whole input folder, and include all chains in it.

    Since you can decide to analyze some file(s), or a complete folder, this
    function first needs to separate between the two cases.

    .. warning::
        If someday you change the way the chains are named, remember to change
        here too, because this routine assumes the chains have a double
        underscore in their names.

    .. note::
        Only files ending with .txt will be selected, to keep compatibility
        with CosmoMC format

    .. note::
        New in version 2.0.0: if you ask to analyze a Nested Sampling
        sub-folder (i.e. something that ends in `NS` with capital letters), the
        analyze module will translate the output from Nested Sampling to
        standard chains for Monte Python, and stops. You can then run the
        `-- info` flag on the whole folder. **This procedure is not necessary
        if the run was complete, but only if the Nested Sampling run was killed
        before completion**.

    Parameters
    ----------
    files : list
        list of potentially only one element, containing the files to analyze.
        This can be only one file, or the encompassing folder, files
    info : Information instance
        Used to store the result

    """
    # First test if the folder is a Nested Sampling or CosmoHammer folder. If
    # so, call the module's own routine through the clean conversion function,
    # which will translate the output of this other sampling into MCMC chains
    # that can then be analyzed.
    modules = ['nested_sampling', 'cosmo_hammer']
    tags = ['NS', 'CH']
    for module_name, tag in zip(modules, tags):
        action_done = clean_conversion(module_name, tag, files[0])
        if action_done:
            return False

    # If the input command was an entire folder, then grab everything in it.
    # Too small files (below 600 octets) and subfolders are automatically
    # removed.
    folder, files = recover_folder_and_files(files)

    info.files = files
    info.folder = folder

    # Check if the log.param file exists
    parameter_file_path = os.path.join(folder, 'log.param')
    if os.path.isfile(parameter_file_path):
        if os.path.getsize(parameter_file_path) == 0:
            raise io_mp.AnalyzeError(
                "The log param file %s " % os.path.join(folder, 'log.param') +
                "seems empty")
    else:
        raise io_mp.AnalyzeError(
            "The log param file %s " % os.path.join(folder, 'log.param') +
            "is missing in the analyzed folder?")

    # If the folder has no subdirectory, then go for a simple infoname,
    # otherwise, call it with the last name
    basename = (os.path.basename(folder) if os.path.basename(folder) != '.'
                else os.path.basename(os.path.abspath(
                    os.path.join(folder, '..'))))

    info.v_info_path = os.path.join(folder, basename+'.v_info')
    info.h_info_path = os.path.join(folder, basename+'.h_info')
    info.tex_path = os.path.join(folder, basename+'.tex')
    info.cov_path = os.path.join(folder, basename+'.covmat')
    info.log_path = os.path.join(folder, basename+'.log')
    info.best_fit_path = os.path.join(folder, basename+'.bestfit')
    info.param_path = parameter_file_path

    return True


@timeit
def convergence(info):
    """
    Compute convergence for the desired chains, using Gelman-Rubin diagnostic

    Chains have been stored in the info instance of :class:`Information`. Note
    that the G-R diagnostic can be computed for a single chain, albeit it will
    most probably give absurd results. To do so, it separates the chain into
    three subchains.

    """
    # Recovering parameter names and scales, creating tex names,
    extract_parameter_names(info)

    # Circle through all files to find the global maximum of likelihood
    print '--> Finding global maximum of likelihood'
    find_maximum_of_likelihood(info)

    # Restarting the circling through files, this time removing the burnin,
    # given the maximum of likelihood previously found and the global variable
    # LOG_LKL_CUTOFF. spam now contains all the accepted points that were
    # explored once the chain moved within max_lkl - LOG_LKL_CUTOFF
    print '--> Removing burn-in'
    spam = remove_burnin(info)

    # Now that the list spam contains all the different chains removed of
    # their respective burn-in, proceed to the convergence computation

    # 2D arrays for mean and var, one column will contain the total (over
    # all chains) mean (resp. variance), and each other column the
    # respective chain mean (resp. chain variance). R only contains the
    # values for each parameter. Therefore, mean and var will have len(spam)+1
    # as a first dimension
    mean = np.zeros((len(spam)+1, info.number_parameters))
    var = np.zeros((len(spam)+1, info.number_parameters))
    R = np.zeros(info.number_parameters)

    # Store the total number of points, and the total in each chain
    total = np.zeros(len(spam)+1)
    for j in xrange(len(spam)):
        total[j+1] = spam[j][:, 0].sum()
    total[0] = total[1:].sum()

    # Compute mean and variance for each chain
    print '--> Computing mean values'
    compute_mean(mean, spam, total)

    print '--> Computing variance'
    compute_variance(var, mean, spam, total)

    print '--> Computing convergence criterium (Gelman-Rubin)'
    # Gelman Rubin Diagnostic:
    # Computes a quantity linked to the ratio of the mean of the variances of
    # the different chains (within), and the variance of the means (between)
    # Note: This is not strictly speaking the Gelman Rubin test, defined for
    # same-length MC chains. Our quantity is defined without the square root,
    # which should not change much the result: a small sqrt(R) will still be a
    # small R. The same convention is used in CosmoMC, except for the weighted
    # average: we decided to do the average taking into account that longer
    # chains should count more
    within = 0
    between = 0

    for i in xrange(np.shape(mean)[1]):
        for j in xrange(len(spam)):
            within += total[j+1]*var[j+1, i]
            between += total[j+1]*(mean[j+1, i]-mean[0, i])**2
        within /= total[0]
        between /= (total[0]-1)

        R[i] = between/within
        if i == 0:
            print ' -> R is %.6f' % R[i], '\tfor ', info.ref_names[i]
        else:
            print '         %.6f' % R[i], '\tfor ', info.ref_names[i]

    # Log finally the total number of steps, and absolute loglikelihood
    with open(info.log_path, 'a') as log:
        log.write("--> Total    number    of    steps: %d\n" % (
            info.steps))
        log.write("--> Total number of accepted steps: %d\n" % (
            info.accepted_steps))
        log.write("--> Minimum of -logLike           : %.2f" % info.max_lkl)

    # Store the remaining members in the info instance, for further writing to
    # files.
    info.spam = spam
    info.mean = mean
    info.R = R
    info.total = total


def plot_triangle(
        info, chain, command_line, bin_number=20, scales=(),
        levels=(68.26, 95.4, 99.7), aspect=(16, 16), fig=None,
        tick_at_peak=False, comp_chain=None, comp_ref_names=None,
        comp_tex_names=None, comp_backup_names=None,
        comp_plotted_parameters=None, comp_folder=None,
        comp_boundaries=None, comp_mean=None):
    """
    Plotting routine, also computes the sigma errors.

    Partly imported from Karim Benabed in CosmoPmc.

    """

    # If comparison is asked, don't plot 2d levels
    # unless explicitely wanted by specifying `-plot-2d <mode>`
    if command_line.comp is not None:
        plot_2d = command_line.plot_2d == 'always' or command_line.plot_2d == 'overplot_comp'
        overplot_comp_contour = command_line.plot_2d == 'overplot_comp'
        comp = True
        comp_done = False
    else:
        plot_2d = command_line.plot_2d != 'no'
        comp = False
        comp_done = False

    # Pre configuration of the output, note that changes to the font size
    # will occur later on as well, to obtain a nice scaling.
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=11)
    matplotlib.rc('xtick', labelsize='8')
    matplotlib.rc('ytick', labelsize='8')
    lvls = np.array(levels)/100.

    # Create the figures
    if plot_2d:
        if fig:
            fig2d = plt.figure(fig, aspect)
        else:
            fig2d = plt.figure(1, figsize=aspect)

    #exps = ', '.join([elem.replace('_', ' ') for elem in info.experiments])
    #plt.figtext(0.4,0.95,'Experiments: '+exps,fontsize=40,alpha=0.6)
    #plt.figtext(0.9, 0.7,'Monte Python',fontsize=70,\
    #rotation=90,alpha=0.15)
    fig1d = plt.figure(2, figsize=aspect)

    # clear figure
    plt.clf()

    # Recover the total number of parameters to potentially plot
    n = np.shape(chain)[1]-2
    if not scales:
        scales = np.ones(n)
    scales = np.array(scales)

    # 1D plot
    max_values = np.max(chain[:, 2:], axis=0)*scales
    min_values = np.min(chain[:, 2:], axis=0)*scales
    span = (max_values-min_values)

    best_minus_lkl = np.min(chain[:, 1], axis=0)

    if comp:
        comp_max_values = np.max(comp_chain[:, 2:], axis=0)
        comp_min_values = np.min(comp_chain[:, 2:], axis=0)
        comp_span = (comp_max_values-comp_min_values)

    # Define the place of ticks
    if tick_at_peak:
        pass
    else:
        ticks = np.array(
            (min_values+span*0.1, (max_values+min_values)/2.,
             max_values-span*0.1)).T
        x_range = np.array((min_values, max_values)).T
        if comp:
            comp_ticks = np.array(
                (comp_min_values+comp_span*0.1,
                 (comp_max_values+comp_min_values)/2.,
                 comp_max_values-comp_span*0.1)).T
            comp_x_range = np.array((comp_min_values, comp_max_values)).T
            for i in range(np.shape(comp_ticks)[0]):
                if abs(comp_x_range[i][0]-comp_boundaries[i][0]) < \
                        comp_span[i]/bin_number:
                    comp_ticks[i][0] = comp_boundaries[i][0]
                    comp_x_range[i][0] = comp_boundaries[i][0]
                if abs(comp_x_range[i][1]-comp_boundaries[i][1]) < \
                        comp_span[i]/bin_number:
                    comp_ticks[i][2] = comp_boundaries[i][1]
                    comp_x_range[i][1] = comp_boundaries[i][1]

    for i in range(np.shape(ticks)[0]):
        if abs(x_range[i][0]-info.boundaries[i][0]) < span[i]/bin_number:
            ticks[i][0] = info.boundaries[i][0]
            x_range[i][0] = info.boundaries[i][0]
        if abs(x_range[i][1]-info.boundaries[i][1]) < span[i]/bin_number:
            ticks[i][2] = info.boundaries[i][1]
            x_range[i][1] = info.boundaries[i][1]

    # Borders stuff, might need adjustement for printing on paper.
    fig1d.subplots_adjust(
        bottom=0.03, left=.07, right=0.98, top=0.93, hspace=.35)
    if plot_2d:
        fig2d.subplots_adjust(
            bottom=0.03, left=.07, right=0.98, top=0.93, hspace=.35)

    # In case of a comparison, figure out which names are shared, which are
    # unique and thus require a simple treatment.
    if comp:
        backup_comp_names = np.copy(comp_plotted_parameters)

        for i in xrange(len(info.plotted_parameters)):
            if info.plotted_parameters[i] in comp_plotted_parameters:
                comp_plotted_parameters.remove(info.plotted_parameters[i])

        num_columns = int(round(math.sqrt(
            len(info.plotted_parameters) + len(comp_plotted_parameters))))
        num_lines = int(math.ceil(
            (len(info.plotted_parameters)+len(comp_plotted_parameters)) *
            1.0/num_columns))
    else:
        num_columns = int(round(math.sqrt(len(info.plotted_parameters))))
        num_lines = int(math.ceil(
            len(info.plotted_parameters)*1.0/num_columns))

    # If plots/ folder in output folder does not exist, create it
    if os.path.isdir(info.folder+'plots') is False:
        os.mkdir(info.folder+'plots')

    # Actual plotting
    print '-----------------------------------------------'
    for i in xrange(len(info.plotted_parameters)):

        print ' -> Computing histograms for ',\
            info.plotted_parameters[i]

        index = info.ref_names.index(info.plotted_parameters[i])
        # Adding the subplots to the respective figures, this will be the
        # diagonal for the triangle plot.
        if plot_2d:
            ax2d = fig2d.add_subplot(
                len(info.plotted_parameters),
                len(info.plotted_parameters),
                i*(len(info.plotted_parameters)+1)+1,
                yticks=[])
        ax1d = fig1d.add_subplot(
            num_lines, num_columns, i+1, yticks=[])

        # normalized histogram
        hist, bin_edges = np.histogram(
            chain[:, index+2], bins=bin_number,
            weights=chain[:, 0], normed=False)
        bincenters = 0.5*(bin_edges[1:]+bin_edges[:-1])

        # interpolated histogram (if available)
        interp_hist, interp_grid = cubic_interpolation(
            info, hist, bincenters)
        interp_hist /= np.max(interp_hist)

        if comp:
            try:
                # For the names in common, the following line will not
                # output an error. Then compute the comparative
                # histogram
                comp_index = comp_ref_names.index(
                    info.plotted_parameters[i])
                comp_hist, comp_bin_edges = np.histogram(
                    comp_chain[:, comp_index+2], bins=bin_number,
                    weights=comp_chain[:, 0], normed=False)
                comp_bincenters = 0.5*(
                    comp_bin_edges[1:]+comp_bin_edges[:-1])
                interp_comp_hist, interp_comp_grid = \
                    cubic_interpolation(info, comp_hist, comp_bincenters)
                interp_comp_hist /= interp_comp_hist.max()
                comp_done = True
            except ValueError:
                # If the name was not found, return the error. This will be
                # then plotted at the end
                comp_done = False
        if comp:
            if not comp_done:
                print('{0} was not found in the second folder'.format(
                    info.plotted_parameters[i]))

        # minimum credible interval (method by Jan Haman). Fails for
        # multimodal histograms
        bounds = minimum_credible_intervals(hist, bincenters, lvls)
        if bounds is False:
            # print out the faulty histogram (try reducing the binnumber to
            # avoir this)
            print hist
        else:
            for elem in bounds:
                for j in (0, 1):
                    elem[j] -= info.mean[index]
            info.bounds[index] = bounds

        if comp_done:
            comp_bounds = minimum_credible_intervals(
                comp_hist, comp_bincenters, lvls)
            if comp_bounds is False:
                print comp_hist
            else:
                for elem in comp_bounds:
                    for j in (0, 1):
                        elem[j] -= comp_mean[comp_index]

        # plotting
        if plot_2d:
            ax2d.set_xticks(ticks[index])
            # First, assign to it the default value
            fontsize2d, ticksize2d = get_fontsize(info, '2d')
            # Then potentially overwrite with provided command line arguments
            if command_line.fontsize != -1:
                fontsize2d = command_line.fontsize
            if command_line.ticksize != -1:
                ticksize2d = command_line.ticksize
            ax2d.set_xticklabels(['%.3g' % s for s in ticks[index]],
                                 fontsize=ticksize2d)
            ax2d.set_title('%s= $%.3g^{+%.3g}_{%.3g}$' % (
                info.tex_names[index], info.mean[index],
                bounds[0][1], bounds[0][0]), fontsize=fontsize2d)
            ax2d.plot(interp_grid, interp_hist, color='red',
                      linewidth=2, ls='-')
            ax2d.axis([x_range[index][0], x_range[index][1], 0, 1.05])

        fontsize1d, ticksize1d = get_fontsize(info, '1d', num_columns)

        if command_line.fontsize != -1:
            fontsize1d = command_line.fontsize
        if command_line.ticksize != -1:
            ticksize1d = command_line.ticksize
        ax1d.set_title('%s= $%.3g^{+%.3g}_{%.3g}$' % (
            info.tex_names[index], info.mean[index],
            bounds[0][1], bounds[0][0]), fontsize=fontsize1d)
        ax1d.set_xticks(ticks[index])
        ax1d.set_xticklabels(['%.3g' % s for s in ticks[index]],
                             fontsize=ticksize1d)
        ax1d.axis([x_range[index][0], x_range[index][1], 0, 1.05])

        if comp_done:
            # complex variation of intervals
            comp_index = comp_ref_names.index(info.plotted_parameters[i])
            if comp_x_range[comp_index][0] > x_range[index][0]:
                comp_ticks[comp_index][0] = ticks[index][0]
                comp_x_range[comp_index][0] = x_range[index][0]
            if comp_x_range[comp_index][1] < x_range[index][1]:
                comp_ticks[comp_index][2] = ticks[index][2]
                comp_x_range[comp_index][1] = x_range[index][1]
            comp_ticks[comp_index][1] = (
                comp_x_range[comp_index][1]+comp_x_range[comp_index][0])/2.
            ax1d.set_xticks(comp_ticks[comp_index])
            ax1d.set_xticklabels(['%.3g' % s for s in comp_ticks[comp_index]],
                                 fontsize=ticksize1d)
            ax1d.axis([comp_x_range[comp_index][0], comp_x_range[comp_index][1], 0, 1.05])

        ax1d.plot(
            interp_grid, interp_hist, color='black', linewidth=2, ls='-')
        if comp_done:
            ax1d.plot(
                interp_comp_grid, interp_comp_hist, color='red',
                linewidth=2, ls='-')
            if plot_2d:
                ax2d.plot(
                    interp_comp_grid, interp_comp_hist, color='blue',
                    linewidth=2, ls='-')

        # mean likelihood (optional, if comparison, it will not be printed)
        if not comp and command_line.mean_likelihood:
            try:
                lkl_mean, _ = np.histogram(
                    chain[:, index+2],
                    bins=bin_edges,
                    normed=False,
                    weights=np.exp(best_minus_lkl-chain[:, 1])*chain[:, 0])
                lkl_mean /= lkl_mean.max()
                interp_lkl_mean, interp_grid = cubic_interpolation(
                    info, lkl_mean, bincenters)
                if plot_2d:
                    ax2d.plot(interp_grid, interp_lkl_mean, color='red',
                              ls='--', lw=2)
                ax1d.plot(interp_grid, interp_lkl_mean, color='red',
                          ls='--', lw=4)
            except:
                print 'could not find likelihood contour for ',
                print info.plotted_parameters[i]

        if command_line.subplot is True:
            if not comp:
                extent2d = ax2d.get_window_extent().transformed(
                    fig2d.dpi_scale_trans.inverted())
                fig2d.savefig(info.folder+'plots/{0}_{1}.{2}'.format(
                    info.folder.split('/')[-2],
                    info.plotted_parameters[i], info.extension),
                    bbox_inches=extent2d.expanded(1.1, 1.4))
            else:
                extent1d = ax1d.get_window_extent().transformed(
                    fig1d.dpi_scale_trans.inverted())
                fig1d.savefig(info.folder+'plots/{0}_{1}.{2}'.format(
                    info.folder.split('/')[-2],
                    info.plotted_parameters[i], info.extension),
                    bbox_inches=extent1d.expanded(1.1, 1.4))
            # Store the function in a file
            hist_file_name = os.path.join(
                info.folder, 'plots/{0}_{1}.hist'.format(
                    info.folder.split(os.path.sep)[-2],
                    info.ref_names[index]))
            print hist_file_name
            write_histogram(hist_file_name, interp_grid, interp_hist)

        # Now do the rest of the triangle plot
        if plot_2d:
            for j in xrange(i):
                second_index = info.ref_names.index(
                    info.plotted_parameters[j])
                ax2dsub = fig2d.add_subplot(
                    len(info.plotted_parameters),
                    len(info.plotted_parameters),
                    (i)*len(info.plotted_parameters)+j+1)
                n, xedges, yedges = np.histogram2d(
                    chain[:, index+2], chain[:, second_index+2],
                    weights=chain[:, 0], bins=(bin_number, bin_number),
                    normed=False)
                extent = [x_range[second_index][0],
                          x_range[second_index][1],
                          x_range[index][0], x_range[index][1]]
                x_centers = 0.5*(xedges[1:]+xedges[:-1])
                y_centers = 0.5*(yedges[1:]+yedges[:-1])

                ax2dsub.set_xticks(ticks[second_index])
                if i == len(info.plotted_parameters)-1:
                    ax2dsub.set_xticklabels(
                        ['%.3g' % s for s in ticks[second_index]],
                        fontsize=ticksize2d)
                else:
                    ax2dsub.set_xticklabels([''])

                ax2dsub.set_yticks(ticks[index])
                if j == 0:
                    ax2dsub.set_yticklabels(
                        ['%.3g' % s for s in ticks[index]],
                        fontsize=ticksize2d)
                else:
                    ax2dsub.set_yticklabels([''])
                #ax2dsub.imshow(n, extent=extent,
                               #aspect='auto', interpolation='gaussian',
                               #origin='lower', cmap=matplotlib.cm.Reds)

                # plotting contours, using the ctr_level method (from Karim
                # Benabed). Note that only the 1 and 2 sigma contours are
                # displayed (due to the line with lvls[:2])
                try:
                    contours = ax2dsub.contourf(
                        y_centers, x_centers, n,
                        extent=extent, levels=ctr_level(n, lvls[:2]),
                        zorder=4, cmap=plt.cm.autumn_r,)
                except Warning:
                    warnings.warn(
                        "The routine could not find the contour of the " +
                        "'%s-%s' 2d-plot" % (
                            info.plotted_parameters[i],
                            info.plotted_parameters[j]))
                    pass

                if comp_done:
                    # i.e. comp_index is valid
                    try:
                        # For the names in common, the following line will not
                        # output an error. Then compute the comparative
                        # histogram
                        comp_second_index = comp_ref_names.index(
                            info.plotted_parameters[j])
                        comp_n, comp_xedges, comp_yedges = np.histogram2d(
                            comp_chain[:, comp_index+2], comp_chain[:, comp_second_index+2],
                            weights=comp_chain[:, 0], bins=(bin_number, bin_number),
                            normed=False)
                        comp_extent = [comp_x_range[comp_second_index][0],
                            comp_x_range[comp_second_index][1],
                            comp_x_range[comp_index][0], comp_x_range[comp_index][1]]
                        comp_x_centers = 0.5*(comp_xedges[1:]+comp_xedges[:-1])
                        comp_y_centers = 0.5*(comp_yedges[1:]+comp_yedges[:-1])
                        comp_done_other = True
                    except ValueError:
                        # If the name was not found, return the error. This will be
                        # then plotted at the end
                        comp_done_other = False

                    if comp_done_other:
                        try:
                            if overplot_comp_contour:
                                contours = ax2dsub.contour(
                                    comp_y_centers, comp_x_centers, comp_n,
                                    extent=comp_extent, levels=ctr_level(comp_n, lvls[:2]),
                                    zorder=5, cmap=plt.cm.Blues)
                            else:
                                contours = ax2dsub.contourf(
                                    comp_y_centers, comp_x_centers, comp_n,
                                    extent=comp_extent, levels=ctr_level(comp_n, lvls[:2]),
                                    zorder=5, cmap=plt.cm.Blues,
                                    alpha=command_line.alpha)
                        except Warning:
                            warnings.warn(
                                "The routine could not find the contour of the " +
                                "'%s-%s' 2d-plot (using the comparison data)" % (
                                info.plotted_parameters[i],
                                info.plotted_parameters[j]))
                            pass
                        ax2dsub.axis([x_range[second_index][0], x_range[second_index][1],
                            x_range[index][0], x_range[index][1]])

                else:
                    comp_done_other = False

                if command_line.subplot is True:
                    # Store the individual 2d plots. Note that the tick and
                    # fontsize are hardcoded here since they will always keep
                    # the same size.
                    fig_temp = plt.figure(3, figsize=(6, 6))
                    fig_temp.clf()
                    ax_temp = fig_temp.add_subplot(111)
                    ax_temp.set_xticks(ticks[second_index])
                    ax_temp.set_yticks(ticks[index])
                    ax_temp.set_xticklabels(
                        ['%.3g' % s for s in ticks[second_index]],
                        fontsize=16)
                    ax_temp.set_yticklabels(
                        ['%.3g' % s for s in ticks[index]],
                        fontsize=16)
                    ax_temp.set_title(
                        '%s vs %s' % (
                            info.tex_names[index],
                            info.tex_names[second_index]),
                        fontsize=16)
                    try:
                        contours = ax_temp.contourf(
                            y_centers, x_centers, n, extent=extent,
                            levels=ctr_level(n, lvls[:2]),  # colors="k",
                            zorder=4, cmap=plt.cm.autumn_r,)
                    except Warning:
                        warnings.warn(
                            "The routine could not find the contour of the " +
                            "'%s-%s' 2d-plot" % (
                                info.plotted_parameters[i],
                                info.plotted_parameters[j]))
                        pass

                    fig_temp.savefig(
                        info.folder+'plots/{0}_2d_{1}-{2}.{3}'.format(
                            info.folder.split('/')[-2],
                            info.ref_names[index],
                            info.ref_names[second_index], info.extension))

                    if comp_done_other:
                        try:
                            if overplot_comp_contour:
                                contours = ax_temp.contour(
                                    comp_y_centers, comp_x_centers, comp_n,
                                    extent=comp_extent, levels=ctr_level(comp_n, lvls[:2]),
                                    zorder=5, cmap=plt.cm.Blues)
                            else:
                                contours = ax_temp.contourf(
                                    comp_y_centers, comp_x_centers, comp_n,
                                    extent=comp_extent, levels=ctr_level(comp_n, lvls[:2]),
                                    zorder=5, cmap=plt.cm.Blues,
                                    alpha=command_line.alpha)
                        except Warning:
                            warnings.warn(
                                "The routine could not find the contour of the " +
                                "'%s-%s' 2d-plot (using the comparison data)" % (
                                info.plotted_parameters[i],
                                info.plotted_parameters[j]))
                            pass
                        ax_temp.axis([x_range[second_index][0],
                                      x_range[second_index][1],
                                      x_range[index][0], x_range[index][1]])

                        fig_temp.savefig(
                            info.folder+'plots/{0}-vs-{1}_2d_{2}-{3}.{4}'.format(
                            info.folder.split('/')[-2],
                            comp_folder.split('/')[-2],
                            info.ref_names[index],
                            info.ref_names[second_index], info.extension))

                    # store the coordinates of the points for further
                    # plotting.
                    plot_file = open(
                        info.folder+'plots/{0}_2d_{1}-{2}.dat'.format(
                            info.folder.split('/')[-2],
                            info.ref_names[index],
                            info.ref_names[second_index]), 'w')
                    plot_file.write(
                        '# contour for confidence level {0}\n'.format(
                            levels[1]))
                    for elem in contours.collections[0].get_paths():
                        points = elem.vertices
                        for k in range(np.shape(points)[0]):
                            plot_file.write("%.8g\t %.8g\n" % (
                                points[k, 0], points[k, 1]))
                            # stop to not include the inner contours
                            if k != 0:
                                if all(points[k] == points[0]):
                                    plot_file.write("\n")
                                    break

                    plot_file.write("\n\n")

                    plot_file.write(
                        '# contour for confidence level {0}\n'.format(
                            levels[0]))
                    for elem in contours.collections[1].get_paths():
                        points = elem.vertices
                        for k in range(np.shape(points)[0]):
                            plot_file.write("%.8g\t %.8g\n" % (
                                points[k, 0], points[k, 1]))
                            if k != 0:
                                if all(points[k] == points[0]):
                                    plot_file.write("\n")
                                    break
                    plot_file.write("\n\n")
                    plot_file.close()

                    # Store also directly the histogram
                    #contours = ax2dsub.contourf(
                        #y_centers, x_centers, n,
                        #extent=extent, levels=ctr_level(n, lvls[:2]),
                        #zorder=5, cmap=plt.cm.autumn_r,
                        #alpha=command_line.alpha)
                    hist_file_name = os.path.join(
                        info.folder, 'plots/{0}_2d_{1}-{2}.hist'.format(
                            info.folder.split('/')[-2],
                            info.ref_names[index],
                            info.ref_names[second_index]))
                    write_histogram_2d(
                        hist_file_name, x_centers, y_centers, extent, n)

    # Plot the remaining 1d diagram for the parameters only in the comp
    # folder
    if comp:
        #if len(info.plotted_parameters) == len(info.ref_names):
        for i in xrange(
                len(info.plotted_parameters),
                len(info.plotted_parameters)+len(comp_plotted_parameters)):

            ax1d = fig1d.add_subplot(
                num_lines, num_columns, i+1, yticks=[])
            comp_index = comp_ref_names.index(
                comp_plotted_parameters[i-len(info.plotted_parameters)])

            comp_hist, comp_bin_edges = np.histogram(
                comp_chain[:, comp_index+2], bins=bin_number,
                weights=comp_chain[:, 0], normed=False)
            comp_bincenters = 0.5*(comp_bin_edges[1:]+comp_bin_edges[:-1])
            interp_comp_hist, interp_comp_grid = cubic_interpolation(
                info, comp_hist, comp_bincenters)
            interp_comp_hist /= interp_comp_hist.max()

            comp_bounds = minimum_credible_intervals(
                comp_hist, comp_bincenters, lvls)
            if comp_bounds is False:
                print comp_hist
            else:
                for elem in comp_bounds:
                    for j in (0, 1):
                        elem[j] -= comp_mean[comp_index]
            ax1d.set_xticks(comp_ticks[comp_index])
            ax1d.set_xticklabels(['%.3g' % s for s in comp_ticks[comp_index]],
                                 fontsize=ticksize1d)
            ax1d.axis([comp_x_range[comp_index][0], comp_x_range[comp_index][1], 0, 1.05])
            ax1d.set_title(
                '%s= $%.3g^{+%.3g}_{%.3g}$' % (
                    comp_tex_names[comp_index], comp_mean[comp_index], comp_bounds[0][1],
                    comp_bounds[0][0]), fontsize=fontsize1d)
            ax1d.plot(interp_comp_grid, interp_comp_hist, color='red',
                      linewidth=2, ls='-')
            if command_line.subplot is True:
                extent1d = ax1d.get_window_extent().transformed(
                    fig1d.dpi_scale_trans.inverted())
                fig1d.savefig(info.folder+'plots/{0}_{1}.{2}'.format(
                    info.folder.split('/')[-2],
                    info.ref_names[i], info.extension),
                    bbox_inches=extent1d.expanded(1.1, 1.4))

    print '-----------------------------------------------'
    print '--> Saving figures to .{0} files'.format(info.extension)
    if comp:
        if plot_2d:
            fig2d.savefig(
                info.folder+'plots/{0}-vs-{1}_triangle.{2}'.format(
                info.folder.split('/')[-2],
                comp_folder.split('/')[-2], info.extension),
                bbox_inches=0, )
        fig1d.savefig(
            info.folder+'plots/{0}-vs-{1}_1d.{2}'.format(
                info.folder.split('/')[-2],
                comp_folder.split('/')[-2], info.extension),
            bbox_inches=0)
    else:
        if plot_2d:
            fig2d.savefig(
                info.folder+'plots/{0}_triangle.{1}'.format(
                info.folder.split('/')[-2], info.extension),
                bbox_inches=0, )
        fig1d.savefig(
            info.folder+'plots/{0}_1d.{1}'.format(
                info.folder.split('/')[-2], info.extension),
            bbox_inches=0)


def ctr_level(histogram2d, lvl, infinite=False):
    """
    Extract the contours for the 2d plots (Karim Benabed)

    """

    hist = histogram2d.flatten()*1.
    hist.sort()
    cum_hist = np.cumsum(hist[::-1])
    cum_hist /= cum_hist[-1]

    alvl = np.searchsorted(cum_hist, lvl)[::-1]
    clist = [0]+[hist[-i] for i in alvl]+[hist.max()]
    if not infinite:
        return clist[1:]
    return clist


def minimum_credible_intervals(histogram, bincenters, levels):
    """
    Extract minimum credible intervals (method from Jan Haman)
    """
    bounds = np.zeros((len(levels), 2))
    j = 0
    delta = bincenters[1]-bincenters[0]
    left_edge = np.max(histogram[0] - 0.5*(histogram[1]-histogram[0]), 0.)
    right_edge = np.max(histogram[-1] + 0.5*(histogram[-1]-histogram[-2]), 0.)
    failed = False
    for level in levels:
        norm = float(
            (np.sum(histogram)-0.5*(histogram[0]+histogram[-1]))*delta)
        norm += 0.25*(left_edge+histogram[0])*delta
        norm += 0.25*(right_edge+histogram[-1])*delta
        water_level_up = np.max(histogram)*1.0
        water_level_down = np.min(histogram)*1.0
        top = 0.

        iterations = 0
        while (abs((top/norm)-level) > 0.0001) and not failed:
            top = 0.
            water_level = (water_level_up + water_level_down)/2.
            #ontop = [elem for elem in histogram if elem > water_level]
            indices = [i for i in range(len(histogram))
                       if histogram[i] > water_level]
            # check for multimodal posteriors
            if ((indices[-1]-indices[0]+1) != len(indices)):
                warnings.warn(
                    "could not derive minimum credible intervals " +
                    "for this multimodal posterior")
                failed = True
                break
            top = (np.sum(histogram[indices]) -
                   0.5*(histogram[indices[0]]+histogram[indices[-1]]))*(delta)

            # left
            if indices[0] > 0:
                top += (0.5*(water_level+histogram[indices[0]]) *
                        delta*(histogram[indices[0]]-water_level) /
                        (histogram[indices[0]]-histogram[indices[0]-1]))
            else:
                if (left_edge > water_level):
                    top += 0.25*(left_edge+histogram[indices[0]])*delta
                else:
                    top += (0.25*(water_level + histogram[indices[0]]) *
                            delta*(histogram[indices[0]]-water_level) /
                            (histogram[indices[0]]-left_edge))

            # right
            if indices[-1] < (len(histogram)-1):
                top += (0.5*(water_level + histogram[indices[-1]]) *
                        delta*(histogram[indices[-1]]-water_level) /
                        (histogram[indices[-1]]-histogram[indices[-1]+1]))
            else:
                if (right_edge > water_level):
                    top += 0.25*(right_edge+histogram[indices[-1]])*delta
                else:
                    top += (0.25*(water_level + histogram[indices[-1]]) *
                            delta * (histogram[indices[-1]]-water_level) /
                            (histogram[indices[-1]]-right_edge))

            if top/norm >= level:
                water_level_down = water_level
            else:
                water_level_up = water_level
            # safeguard, just in case
            iterations += 1
            if (iterations > 1000):
                warnings.warn(
                    "the loop to check for sigma deviations was " +
                    "taking too long to converge")
                break

        # min
        if indices[0] > 0:
            bounds[j][0] = bincenters[indices[0]] - delta*(histogram[indices[0]]-water_level)/(histogram[indices[0]]-histogram[indices[0]-1])
        else:
            if (left_edge > water_level):
                bounds[j][0] = bincenters[0]-0.5*delta
            else:
                bounds[j][0] = bincenters[indices[0]] - 0.5*delta*(histogram[indices[0]]-water_level)/(histogram[indices[0]]-left_edge)

        # max
        if indices[-1] < (len(histogram)-1):
            bounds[j][1] = bincenters[indices[-1]] + delta*(histogram[indices[-1]]-water_level)/(histogram[indices[-1]]-histogram[indices[-1]+1])
        else:
            if (right_edge > water_level):
                bounds[j][1] = bincenters[-1]+0.5*delta
            else:
                bounds[j][1] = bincenters[indices[-1]] + \
                    0.5*delta*(histogram[indices[-1]]-water_level) / \
                    (histogram[indices[-1]]-right_edge)

        j += 1

    return bounds


def write_h(info_file, indices, name, string, quantity, modifiers=None):
    """
    Write one horizontal line of output

    """
    info_file.write('\n '+name+'\t:\t')
    for i in indices:
        if quantity[i] >= 0:
            space_string = ' '
        else:
            space_string = ''
        info_file.write(space_string+string % quantity[i]+'\t')


def write_tex(info, indices):
    """
    Write a tex table containing the main results

    """
    info.tex.write("\\begin{tabular}{|l|c|c|c|c|} \n \\hline \n")
    info.tex.write("Param & best-fit & mean$\pm\sigma$ & 95\% lower & 95\% upper \\\\ \\hline \n")
    for i in indices:
        info.tex.write("%s &" % info.tex_names[i])
        info.tex.write("$%.4g$ & $%.4g_{%.2g}^{+%.2g}$ & $%.4g$ & $%.4g$ \\\\ \n" % (
            info.bestfit[i], info.mean[i], info.bounds[:, 0, 0][i],
            info.bounds[:, 0, 1][i], info.mean[i]+info.bounds[:, 1, 0][i],
            info.mean[i]+info.bounds[:, 1, 1][i]))

    info.tex.write("\\hline \n \\end{tabular} \\\\ \n")
    info.tex.write(
        "$-\ln{\cal L}_\mathrm{min} =%.6g$, minimum $\chi^2=%.4g$ \\\\ \n" %
        (info.max_lkl, info.max_lkl*2.))


def cubic_interpolation(info, hist, bincenters):
    """
    Small routine to accomodate the absence of the interpolate module

    """
    if info.has_interpolate_module:
        interp_grid = np.linspace(
            bincenters[0], bincenters[-1], len(bincenters)*10)
        from scipy.interpolate import interp1d
        f = interp1d(bincenters, hist, kind='cubic')
        interp_hist = f(interp_grid)
        return interp_hist, interp_grid
    else:
        return hist, bincenters


def get_fontsize(info, geometry, width=None):
    """
    Empirical method to adjust font size on the plots to fit the number of
    parameters. Feel free to modify to your needs.

    """
    # Approximate values to have roughly nice displays font size
    if geometry == '2d':
        size = 75./len(info.plotted_parameters)
    if geometry == '1d':
        if width:
            size = 80./width
        else:
            raise io_mp.AnalyzeError(
                "This routine expects a number of columns "
                "for the 1d plot")
    return size, size


def write_histogram(hist_file_name, x_centers, hist):
    """
    Store the posterior distribution to a file

    """
    with open(hist_file_name, 'w') as hist_file:
        hist_file.write("# 1d posterior distribution\n")
        hist_file.write("\n# x_centers\n")
        hist_file.write(", ".join(
            [str(elem) for elem in x_centers])+"\n")
        hist_file.write("\n# Histogram\n")
        hist_file.write(", ".join(
            [str(elem) for elem in hist])+"\n")
    print 'wrote ', hist_file_name


def read_histogram(histogram_path):
    """
    Recover a stored 1d posterior

    """
    with open(histogram_path, 'r') as hist_file:
        for line in hist_file:
            if line:
                if line.find("# x_centers") != -1:
                    x_centers = [float(elem) for elem in
                                 hist_file.next().split(",")]
                elif line.find("# Histogram") != -1:
                    hist = [float(elem) for elem in
                            hist_file.next().split(",")]
    x_centers = np.array(x_centers)
    hist = np.array(hist)

    return x_centers, hist


def write_histogram_2d(hist_file_name, x_centers, y_centers, extent, hist):
    """
    Store the histogram information to a file, to plot it later

    """
    with open(hist_file_name, 'w') as hist_file:
        hist_file.write("# Interpolated histogram\n")
        hist_file.write("\n# x_centers\n")
        hist_file.write(", ".join(
            [str(elem) for elem in x_centers])+"\n")

        hist_file.write("\n# y_centers\n")
        hist_file.write(", ".join(
            [str(elem) for elem in y_centers])+"\n")

        hist_file.write("\n# Extent\n")
        hist_file.write(", ".join(
            [str(elem) for elem in extent])+"\n")

        hist_file.write("\n# Histogram\n")
        for line in hist:
            hist_file.write(", ".join(
                [str(elem) for elem in line])+"\n")


def read_histogram_2d(histogram_path):
    """
    Read the histogram information that was stored in a file.

    To use it, call something like this:

    .. code::

        x_centers, y_centers, extent, hist = read_histogram_2d_from_file(path)
        fig, ax = plt.subplots()
        ax.contourf(
            y_centers, x_centers, hist, extent=extent,
            levels=ctr_level(hist, [0.68, 0.95]),
            zorder=5, cma=plt.cm.autumn_r)
        plt.show()

    """
    with open(histogram_path, 'r') as hist_file:
        length = 0
        for line in hist_file:
            if line:
                if line.find("# x_centers") != -1:
                    x_centers = [float(elem) for elem in
                                 hist_file.next().split(",")]
                    length = len(x_centers)
                elif line.find("# y_centers") != -1:
                    y_centers = [float(elem) for elem in
                                 hist_file.next().split(",")]
                elif line.find("# Extent") != -1:
                    extent = [float(elem) for elem in
                              hist_file.next().split(",")]
                elif line.find("# Histogram") != -1:
                    hist = []
                    for index in range(length):
                        hist.append([float(elem) for elem in
                                     hist_file.next().split(",")])
    x_centers = np.array(x_centers)
    y_centers = np.array(y_centers)
    extent = np.array(extent)
    hist = np.array(hist)

    return x_centers, y_centers, extent, hist


def clean_conversion(module_name, tag, folder):
    """
    Execute the methods "convert" from the different sampling algorithms

    Returns True if something was made, False otherwise
    """
    has_module = False
    subfolder_name = tag+"_subfolder"
    try:
        module = importlib.import_module(module_name)
        subfolder = getattr(module, subfolder_name)
        has_module = True
    except ImportError:
        # The module is not installed, the conversion can not take place
        pass

    if has_module and os.path.isdir(folder):
        # Remove any potential trailing slash
        folder = os.path.join(
            *[elem for elem in folder.split(os.path.sep) if elem])
        if folder.split(os.path.sep)[-1] == subfolder:
            try:
                getattr(module, 'from_%s_output_to_chains' % tag)(folder)
            except IOError:
                raise io_mp.AnalyzeError(
                    "You asked to analyze a %s folder which " % tag +
                    "seems to come from an unfinished run, or to be empty " +
                    "or corrupt. Please make sure the run went smoothly " +
                    "enough.")
            warnings.warn(
                "The content of the %s subfolder has been " % tag +
                "translated for Monte Python. Please run an "
                "analysis of the entire folder now.")
            return True
    else:
        return False


def recover_folder_and_files(files):
    """Distinguish the cases when analyze is called with files or folder"""
    # The following list defines the substring that a chain should contain for
    # the code to recognise it as a proper chain.
    substrings = ['.txt', '__']
    limit = 600
    # If the first element is a folder, grab all chain files inside
    if os.path.isdir(files[0]):
        folder = os.path.normpath(files[0])
        files = [os.path.join(folder, elem) for elem in os.listdir(folder)
                 if not os.path.isdir(os.path.join(folder, elem))
                 and not os.path.getsize(os.path.join(folder, elem)) < limit
                 and all([x in elem for x in substrings])]
    # Otherwise, extract the folder from the chain file-name.
    else:
        folder = os.path.relpath(
            os.path.dirname(os.path.realpath(files[0])), os.path.curdir)
        files = [os.path.join(folder, elem) for elem in os.listdir(folder)
                 if os.path.join(folder, elem) in np.copy(files)
                 and not os.path.isdir(os.path.join(folder, elem))
                 and not os.path.getsize(os.path.join(folder, elem)) < limit
                 and all([x in elem for x in substrings])]
    return folder, files


def extract_array(line):
    """
    Return the array on the RHS of the line

    >>> extract_array("toto = ['one', 'two']\n")
    ['one', 'two']
    >>> extract_array('toto = ["one", 0.2]\n')
    ['one', 0.2]

    """
    # Recover RHS of the equal sign, and remove surrounding spaces
    rhs = line.split('=')[-1].strip()
    # Remove array signs
    rhs = rhs.strip(']').lstrip('[')
    # Recover each element of the list
    sequence = [e.strip().strip('"').strip("'") for e in rhs.split(',')]
    for index, elem in enumerate(sequence):
        try:
            sequence[index] = int(elem)
        except ValueError:
            try:
                sequence[index] = float(elem)
            except ValueError:
                pass
    return sequence


def extract_dict(line):
    """
    Return the key and value of the dictionary element contained in line

    >>> extract_dict("something['toto'] = [0, 1, 2, -2, 'cosmo']")
    'toto', [0, 1, 2, -2, 'cosmo']

    """
    # recovering the array
    sequence = extract_array(line)
    # Recovering only the LHS
    lhs = line.split('=')[0].strip()
    # Recovering the name from the LHS
    name = lhs.split('[')[-1].strip(']')
    name = name.strip('"').strip("'")

    return name, sequence


def extract_parameter_names(info):
    """
    Reading the log.param, store in the Information instance the names
    """
    backup_names = []
    plotted_parameters = []
    boundaries = []
    ref_names = []
    tex_names = []
    scales = []
    with open(info.param_path, 'r') as param:
        for line in param:
            if line.find('#') == -1:
                if line.find('data.experiments') != -1:
                    info.experiments = extract_array(line)
                if line.find('data.parameters') != -1:
                    name, array = extract_dict(line)
                    original = name
                    # Rename the names according the .extra file (opt)
                    if name in info.to_change.iterkeys():
                        name = info.to_change[name]
                    # If the name corresponds to a varying parameter (fourth
                    # entry in the initial array being non-zero, or a derived
                    # parameter (could be designed as fixed, it does not make
                    # any difference)), then continue the process of analyzing.
                    if array[3] != 0 or array[5] == 'derived':
                        # The real name is always kept, to have still the class
                        # names in the covmat
                        backup_names.append(original)
                        # With the list "to_plot", we can potentially restrict
                        # the variables plotted. If it is empty, though, simply
                        # all parameters will be plotted.
                        if info.to_plot == []:
                            plotted_parameters.append(name)
                        else:
                            if name in info.to_plot:
                                plotted_parameters.append(name)

                        # Append to the boundaries array
                        boundaries.append([elem if elem != 'None' else -1
                                           for elem in array[1:3]])

                        ref_names.append(name)
                        # Take care of the scales
                        scale = array[4]
                        if name in info.new_scales.iterkeys():
                            scale = info.new_scales[name]
                        scales.append(scale)

                        # Given the scale, decide for the pretty tex name
                        number = 1./scale
                        tex_names.append(
                            io_mp.get_tex_name(name, number=number))
    scales = np.diag(scales)

    info.ref_names = ref_names
    info.tex_names = tex_names
    info.boundaries = boundaries
    info.plotted_parameters = plotted_parameters
    info.number_parameters = len(plotted_parameters)
    info.backup_names = backup_names
    info.scales = scales


@timeit
def find_maximum_of_likelihood(info):
    """
    Finding the global maximum of likelihood

    max_lkl will be appended with all the maximum likelihoods of files,
    then will be replaced by its own maximum. This way, the global
    maximum likelihood will be used as a reference, and not each chain's
    maximum.
    """
    max_lkl = []
    for chain_file in info.files:
        # cheese will brutally contain everything (- log likelihood) in the
        # file chain_file being scanned.
        # This could potentially be faster with pandas, but is already quite
        # fast TODO
        cheese = (np.array([float(line.split()[1].strip())
                            for line in open(chain_file, 'r')]))

        max_lkl.append(cheese[:].min())

    # beware, it is the min because we are talking about
    # '- log likelihood'
    # Selecting only the true maximum.
    try:
        max_lkl = min(max_lkl)
    except ValueError:
        raise io_mp.AnalyzeError(
            "No decently sized chain was found in the desired folder. " +
            "Please wait to have more accepted point before trying " +
            "to analyze it.")

    info.max_lkl = max_lkl


@timeit
def remove_burnin(info):
    """
    Create an array with all the points from the chains, after burnin

    """
    # spam will brutally contain all the chains with sufficient number of
    # points, after the burn-in was removed.
    spam = list()

    # Recover the longest file name, for pleasing display
    max_name_length = max([len(e) for e in info.files])

    # Total number of steps done:
    steps = 0
    accepted_steps = 0

    for index, chain_file in enumerate(info.files):
        # To improve presentation, and print only once the full path of the
        # analyzed folder, we recover the length of the path name, and
        # TODO this could go in IO
        # create an empty complementary string of this length
        total_length = 18+max_name_length
        empty_length = 18+len(os.path.dirname(chain_file))+1

        basename = os.path.basename(chain_file)
        if index == 0:
            exec "print '--> Scanning file %-{0}s' % chain_file,".format(
                max_name_length)
        else:
            exec "print '%{0}s%-{1}s' % ('', basename),".format(
                empty_length, total_length-empty_length)
        # cheese will brutally contain everything in the chain chain_file being
        # scanned
        cheese = (np.array([[float(elem) for elem in line.split()]
                            for line in open(chain_file, 'r')]))
        # If the file contains a broken line with a different number of
        # elements, the previous array generation might fail, and will not have
        # the correct shape. Hence the following command will fail. To avoid
        # that, the error is caught.
        try:
            local_max_lkl = cheese[:, 1].min()
        except IndexError:
            raise io_mp.AnalyzeError(
                "Error while scanning %s." % chain_file +
                " This file most probably contains "
                "an incomplete line, rendering the analysis impossible. "
                "I think that the following line(s) is(are) wrong:\n %s" % (
                    '\n '.join(
                        ['-> %s' % line for line in
                         open(chain_file, 'r') if
                         len(line.split()) != len(info.backup_names)+2])))
        line_count = float(sum(1 for line in open(chain_file, 'r')))

        # Logging the information obtained until now.
        with open(info.log_path, 'a') as log:
            number_of_steps = cheese[:, 0].sum()
            log.write("%s\t " % os.path.basename(chain_file))
            log.write(" Number of steps:%d\t" % number_of_steps)
            log.write(" Steps accepted:%d\t" % line_count)
            log.write(" acc = %.2g\t" % (float(line_count)/number_of_steps))
            log.write("min(-loglike) = %.2f\n" % local_max_lkl)
            steps += number_of_steps
            accepted_steps += line_count

        # Removing burn-in
        start = 0
        try:
            while cheese[start, 1] > info.max_lkl+LOG_LKL_CUTOFF:
                start += 1
            print ': Removed {0}\t points of burn-in'.format(start)
        except IndexError:
            print ': Removed everything: chain not converged'

        # ham contains cheese without the burn-in, if there are any points
        # left (more than 5)
        if np.shape(cheese)[0] > start+5:
            ham = np.copy(cheese[start::])

            # Deal with single file case
            if len(info.files) == 1:
                warnings.warn("Convergence computed for a single file")
                bacon = np.copy(cheese[::3, :])
                egg = np.copy(cheese[1::3, :])
                sausage = np.copy(cheese[2::3, :])

                spam.append(bacon)
                spam.append(egg)
                spam.append(sausage)
                continue

            # Adding resulting table to spam
            spam.append(ham)

    # Test the length of the list
    if len(spam) == 0:
        raise io_mp.AnalyzeError(
            "No decently sized chain was found. " +
            "Please wait a bit to analyze this folder")

    # Applying now new rules for scales
    for name in info.new_scales.iterkeys():
        index = info.ref_names.index(name)
        for i in xrange(len(spam)):
            spam[i][:, index+2] *= 1./info.scales[index, index]

    info.steps = steps
    info.accepted_steps = accepted_steps

    return spam


@timeit
def compute_mean(mean, spam, total):
    """
    """
    for i in xrange(np.shape(mean)[1]):
        for j in xrange(len(spam)):
            submean = np.sum(spam[j][:, 0]*spam[j][:, i+2])
            mean[j+1, i] = submean / total[j+1]
            mean[0, i] += submean
        mean[0, i] /= total[0]


@timeit
def compute_variance(var, mean, spam, total):
    """
    """
    for i in xrange(np.shape(var)[1]):
        for j in xrange(len(spam)):
            var[0, i] += np.sum(
                spam[j][:, 0]*(spam[j][:, i+2]-mean[0, i])**2)
            var[j+1, i] = np.sum(
                spam[j][:, 0]*(spam[j][:, i+2]-mean[j+1, i])**2) / \
                (total[j+1]-1)
        var[0, i] /= (total[0]-1)


class Information(object):
    """
    Hold all information for analyzing runs

    """
    def __init__(self):
        """
        The following initialization creates the three tables that can be
        customized in an extra plot_file (see :mod:`parser_mp`).

        """
        self.to_change = {}
        """
        Dictionary whose keys are the old parameter names, and values are the
        new ones. For instance :code:`{'beta_plus_lambda':'beta+lambda'}`

        """
        self.to_plot = []
        """
        Array of names of parameters to plot. If left empty, all will be
        plotted.

        .. warning::
            If you changed a parameter name with :attr:`to_change`, you need to
            give the new name to this array

        """
        self.new_scales = {}
        """
        Dictionary that redefines some scales. The keys will be the parameter
        name, and the value its scale.

        """

        # Follows a bunch of initialisation to provide default members
        self.ref_names, self.backup_names = [], []
        self.scales, self.plotted_parameters = [], []
        self.spam = []
        self.cov, self.v_info, self.h_info, self.R = None, None, None, None
        self.best_fit = None
