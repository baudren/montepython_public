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
from itertools import count
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
from itertools import ifilterfalse
from itertools import ifilter

# Defined to remove the burnin for all the points that were produced before the
# first time where -log-likelihood <= min-minus-log-likelihood+LOG_LKL_CUTOFF
LOG_LKL_CUTOFF = 3

NUM_COLORS = 4


def analyze(command_line):
    """
    Main function, does the entire analysis.

    It calls in turn all the other routines from this module. To limit the
    arguments of each function to a reasonnable size, a :class:`Information`
    instance is used. This instance is initialized in this function, then
    appended by the other routines.

    """
    # Check if the scipy module has the interpolate method correctly
    # installed (should be the case on every linux distribution with
    # standard numpy)
    try:
        from scipy.interpolate import interp1d
        Information.has_interpolate_module = True
    except ImportError:
        Information.has_interpolate_module = False
        warnings.warn(
            'No cubic interpolation done (no interpolate method found ' +
            'in scipy), only linear')

    # Determine how many different folders are asked through the 'info'
    # command, and create as many Information instances
    files = separate_files(command_line.files)

    # Create an instance of the Information class for each subgroup found in
    # the previous function. They will each hold all relevant information, and
    # be used as a compact way of exchanging information between functions
    information_instances = []
    for item in files:
        info = Information(command_line)
        information_instances.append(info)

        # Prepare the files, according to the case, load the log.param, and
        # prepare the output (plots folder, .covmat, .info and .log files).
        # After this step, info.files will contain all chains.
        status = prepare(item, info)
        # If the preparation step generated new files (for instance,
        # translating from NS or CH to Markov Chains) this routine should stop
        # now.
        if not status:
            return

        # Compute the mean, maximum of likelihood, 1-sigma variance for this
        # main folder. This will create the info.chain object, which contains
        # all the points computed stacked in one big array.
        convergence(info)

        # check if analyze() is called directly by the user, or by the mcmc loop during an updating phase
        try:
            # command_line.update is defined when called by the mcmc loop
            command_line.update
        except:
            # in case it was not defined (i.e. when analyze() is called directly by user), set it to False
            command_line.update = 0

        # compute covariance matrix, excepted when we are in update mode and convergence is too bad or too good
        if command_line.update and (np.amax(info.R) > 3. or np.amax(info.R) < 0.4):
            print '--> Not computing covariance matrix'
        else:
            try:
                if command_line.want_covmat:
                    print '--> Computing covariance matrix'
                    info.covar = compute_covariance_matrix(info)
                    # Writing it out in name_of_folder.covmat
                    io_mp.write_covariance_matrix(
                        info.covar, info.backup_names, info.cov_path)
            except:
                print '--> Computing covariance matrix failed'
                pass

        # Store an array, sorted_indices, containing the list of indices
        # corresponding to the line with the highest likelihood as the first
        # element, and then as decreasing likelihood
        info.sorted_indices = info.chain[:, 1].argsort(0)

        # Writing the best-fit model in name_of_folder.bestfit
        bestfit_line = [elem*info.scales[i, i] for i, elem in
                        enumerate(info.chain[info.sorted_indices[0], 2:])]
        io_mp.write_bestfit_file(bestfit_line, info.backup_names,
                                 info.best_fit_path)

    if not command_line.minimal:
        # Computing 1,2 and 3-sigma errors, and plot. This will create the
        # triangle and 1d plot by default.
        compute_posterior(information_instances)

        print '--> Writing .info and .tex files'
        for info in information_instances:
            info.write_information_files()

    # when called by MCMC in update mode, return R values so that they can be written for information in the chains
    if command_line.update:
        return info.R

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
    folder, files, basename = recover_folder_and_files(files)

    info.files = files
    info.folder = folder
    info.basename = basename

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
    # Now that the number of parameters is known, the array containing bounds
    # can be initialised
    info.bounds = np.zeros((len(info.ref_names), len(info.levels), 2))

    # Circle through all files to find the global maximum of likelihood
    print '--> Finding global maximum of likelihood'
    find_maximum_of_likelihood(info)

    # Restarting the circling through files, this time removing the burnin,
    # given the maximum of likelihood previously found and the global variable
    # LOG_LKL_CUTOFF. spam now contains all the accepted points that were
    # explored once the chain moved within min_minus_lkl - LOG_LKL_CUTOFF.
    # If the user asks for a keep_fraction <1, this is also the place where
    # a fraction (1-keep_fraction) is removed at the beginning of each chain.
    print '--> Removing burn-in'
    spam = remove_bad_points(info)

    info.remap_parameters(spam)
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
            print ' -> R-1 is %.6f' % R[i], '\tfor ', info.ref_names[i]
        else:
            print '           %.6f' % R[i], '\tfor ', info.ref_names[i]

    # Log finally the total number of steps, and absolute loglikelihood
    with open(info.log_path, 'a') as log:
        log.write("--> Total    number    of    steps: %d\n" % (
            info.steps))
        log.write("--> Total number of accepted steps: %d\n" % (
            info.accepted_steps))
        log.write("--> Minimum of -logLike           : %.2f" % (
            info.min_minus_lkl))

    # Store the remaining members in the info instance, for further writing to
    # files, storing only the mean and total of all the chains taken together
    info.mean = mean[0]
    info.R = R
    info.total = total[0]

    # Create the main chain, which consists in all elements of spam
    # put together. This will serve for the plotting.
    info.chain = np.vstack(spam)


def compute_posterior(information_instances):
    """
    computes the marginalized posterior distributions, and optionnally plots
    them

    Parameters
    ----------
    information_instances : list
        list of information objects, initialised on the given folders, or list
        of file, in input. For each of these instance, plot the 1d and 2d
        posterior distribution, depending on the flags stored in the instances,
        comming from command line arguments or read from a file.
    """
    # For convenience, store as `conf` the first element of the list
    # information_instances, since it will be called often to check for
    # configuration parameters
    conf = information_instances[0]

    # Pre configuration of the output, note that changes to the font size
    # will occur later on as well, to obtain a nice scaling.
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=11)
    matplotlib.rc('xtick', labelsize='8')
    matplotlib.rc('ytick', labelsize='8')

    # Recover max and min values for each instance, defining the a priori place
    # of ticks (in case of a comparison, this should change)
    for info in information_instances:
        info.define_ticks()
        # If plots/ folder in output folder does not exist, create it
        if os.path.isdir(os.path.join(info.folder, 'plots')) is False:
            os.mkdir(os.path.join(info.folder, 'plots'))

    # Determine the total number of parameters to plot, based on the list
    # without duplicates of the plotted parameters of all information instances
    plotted_parameters = []
    # For printing not in latex
    ref_names = []
    for info in information_instances:
        for index, name in enumerate(info.plotted_parameters):
            if name not in plotted_parameters:
                plotted_parameters.append(name)
                ref_names.append(info.ref_names[index])

    if len(plotted_parameters) == 0:
        raise io_mp.AnalyzeError(
            "You provided no parameters to analyze, probably by selecting"
            " wrong parameters names in the '--extra' file.")
    # Find the appropriate number of columns and lines for the 1d posterior
    # plot
    num_columns = int(round(math.sqrt(len(plotted_parameters))))
    num_lines = int(math.ceil(len(plotted_parameters)*1.0/num_columns))

    # For special needs, you can impose here a different number of columns and lines in the 1d plot
    # Here is a commented example:
    # if (len(plotted_parameters) == 10):
    #     num_columns = 5
    #     num_lines = 2

    # Create the figures
    # which will be 3*3 inches per subplot, quickly growing!
    if conf.plot:
        fig1d = plt.figure(num=1, figsize=(
            3*num_columns,
            3*num_lines), dpi=80)
    if conf.plot_2d:
        fig2d = plt.figure(num=2, figsize=(
            3*len(plotted_parameters),
            3*len(plotted_parameters)), dpi=80)

    # Create the name of the files, concatenating the basenames with
    # underscores.
    file_name = "_".join(
        [info.basename for info in information_instances])

    # Loop over all the plotted parameters
    # There will be two indices at all time, the one running over the plotted
    # parameters, `index`, and the one corresponding to the actual column in
    # the actual file, `native_index`. For instance, if you try to plot only
    # two columns of a several columns file, index will vary from 0 to 1, but
    # the corresponding native indices might be anything.
    # Obviously, since plotted parameters contain potentially names not
    # contained in some files (in case of a comparison), native index might be
    # undefined.
    # Defined the legends object, which will store the plot style, to display
    # at the level of the figure
    legends = [None for _ in range(len(information_instances))]
    if not conf.legendnames:
        legend_names = [info.basename.replace('_', ' ')
                        for info in information_instances]
    else:
        legend_names = conf.legendnames
    print '-----------------------------------------------'
    for index, name in enumerate(plotted_parameters):

        # Adding the subplots to the respective figures, this will correspond
        # to the diagonal on the triangle plot.
        if conf.plot_2d:
            ax2d = fig2d.add_subplot(
                len(plotted_parameters),
                len(plotted_parameters),
                index*(len(plotted_parameters)+1)+1,
                yticks=[])
            ax2d.set_color_cycle(conf.cm)
        if conf.plot:
            ax1d = fig1d.add_subplot(
                num_lines, num_columns, index+1, yticks=[])
            ax1d.set_color_cycle(conf.cm)

        # check for each instance if the name is part of the list of plotted
        # parameters, and if yes, store the native_index. If not, store a flag
        # to ignore any further plotting or computing issues concerning this
        # particular instance.
        for info in information_instances:
            try:
                info.native_index = info.ref_names.index(name)
                info.ignore_param = False
                standard_name = info.backup_names[info.native_index]
            except ValueError:
                info.ignore_param = True

        adjust_ticks(name, information_instances)

        # normalized histogram
        print ' -> Computing histograms for ', name
        for info in information_instances:
            if not info.ignore_param:
                info.hist, info.bin_edges = np.histogram(
                    info.chain[:, info.native_index+2], bins=info.bins,
                    weights=info.chain[:, 0], normed=False)
                info.bincenters = 0.5*(info.bin_edges[1:]+info.bin_edges[:-1])

                # interpolated histogram (if available)
                info.interp_hist, info.interp_grid = cubic_interpolation(
                    info, info.hist, info.bincenters)
                info.interp_hist /= np.max(info.interp_hist)

                # minimum credible interval (method by Jan Haman). Fails for
                # multimodal histograms #FIXME
                bounds = minimum_credible_intervals(info)
                info.bounds[info.native_index] = bounds

        # plotting
        for info in information_instances:
            if not info.ignore_param:
                if conf.plot_2d:
                    plot = ax2d.plot(
                        info.interp_grid, info.interp_hist,
                        linewidth=info.line_width, ls='-')
                    legends[info.id] = plot[0]
                    ax2d.set_xticks(info.ticks[info.native_index])
                    if conf.legend_style == 'top':
                        ax2d.set_title(
                            '%s=$%.{0}g^{{+%.{0}g}}_{{%.{0}g}}$'.format(
                                info.decimal) % (
                                info.tex_names[info.native_index],
                                info.mean[info.native_index],
                                info.bounds[info.native_index, 0, -1],
                                info.bounds[info.native_index, 0, 0]),
                            fontsize=info.fontsize)
                        ax2d.set_xticklabels(
                            ['%.{0}g'.format(info.decimal) % s
                             for s in info.ticks[info.native_index]],
                            fontsize=info.ticksize)
                    elif conf.legend_style == 'sides':
                        # Except for the last 1d plot (bottom line), don't
                        # print ticks
                        if index == len(plotted_parameters)-1:
                            ax2d.set_xticklabels(
                                ['%.{0}g'.format(info.decimal) % s
                                 for s in info.ticks[info.native_index]],
                                fontsize=info.ticksize)
                            ax2d.set_xlabel(
                                info.tex_names[info.native_index],
                                fontsize=info.fontsize)
                        else:
                            ax2d.set_xticklabels([])
                    ax2d.axis([info.x_range[info.native_index][0],
                               info.x_range[info.native_index][1],
                               0, 1.05])

                if conf.plot:
                    # Note the use of double curly brackets {{ }} to produce
                    # the desired LaTeX output. This is necessary because the
                    # format function would otherwise understand single
                    # brackets as fields.
                    ax1d.set_title(
                        '%s=$%.{0}g^{{+%.{0}g}}_{{%.{0}g}}$'.format(
                            info.decimal) % (
                            info.tex_names[info.native_index],
                            info.mean[info.native_index],
                            info.bounds[info.native_index, 0, -1],
                            info.bounds[info.native_index, 0, 0]),
                        fontsize=info.fontsize)
                    ax1d.set_xticks(info.ticks[info.native_index])
                    ax1d.set_xticklabels(
                        ['%.{0}g'.format(info.decimal) % s
                         for s in info.ticks[info.native_index]],
                        fontsize=info.ticksize)
                    ax1d.axis([info.x_range[info.native_index][0],
                               info.x_range[info.native_index][1],
                               0, 1.05])

                    ax1d.plot(
                        info.interp_grid, info.interp_hist,
                        lw=info.line_width, ls='-')

        # mean likelihood (optional, if comparison, it will not be printed)
        # The color cycle has to be reset, before
        if conf.plot_2d:
            ax2d.set_color_cycle(conf.cm)
        if conf.plot:
            ax1d.set_color_cycle(conf.cm)
        if conf.mean_likelihood:
            for info in information_instances:
                if not info.ignore_param:
                    try:
                        lkl_mean, _ = np.histogram(
                            info.chain[:, info.native_index+2],
                            bins=info.bin_edges,
                            normed=False,
                            weights=np.exp(
                                conf.min_minus_lkl-info.chain[:, 1])*info.chain[:, 0])
                        lkl_mean /= lkl_mean.max()
                        interp_lkl_mean, interp_grid = cubic_interpolation(
                            info, lkl_mean, info.bincenters)
                        if conf.plot_2d:
                            ax2d.plot(interp_grid, interp_lkl_mean,
                                      ls='--', lw=conf.line_width)
                        if conf.plot:
                            ax1d.plot(interp_grid, interp_lkl_mean,
                                      ls='--', lw=conf.line_width)
                    except:
                        print 'could not find likelihood contour for ',
                        print info.ref_parameters[info.native_index]

        if conf.subplot is True:
            if conf.plot_2d:
                extent2d = ax2d.get_window_extent().transformed(
                    fig2d.dpi_scale_trans.inverted())
                fig2d.savefig(os.path.join(
                    conf.folder, 'plots', file_name+'.'+conf.extension),
                    bbox_inches=extent2d.expanded(1.1, 1.4))
            if conf.plot:
                extent1d = ax1d.get_window_extent().transformed(
                    fig1d.dpi_scale_trans.inverted())
                fig1d.savefig(os.path.join(
                    conf.folder, 'plots', file_name+'.'+conf.extension),
                    bbox_inches=extent1d.expanded(1.1, 1.4))
            # Store the function in a file
            for info in information_instances:
                if not info.ignore_param:
                    hist_file_name = os.path.join(
                        info.folder, 'plots',
                        info.basename+'_%s.hist' % (
                            standard_name))
                    write_histogram(hist_file_name,
                                    info.interp_grid, info.interp_hist)

        # Now do the rest of the triangle plot
        if conf.plot_2d:
            for second_index in xrange(index):
                second_name = plotted_parameters[second_index]
                for info in information_instances:
                    if not info.ignore_param:
                        try:
                            info.native_second_index = info.ref_names.index(
                                plotted_parameters[second_index])
                            info.has_second_param = True
                            second_standard_name = info.backup_names[
                                info.native_second_index]
                        except ValueError:
                            info.has_second_param = False
                    else:
                        info.has_second_param = False
                ax2dsub = fig2d.add_subplot(
                    len(plotted_parameters),
                    len(plotted_parameters),
                    (index)*len(plotted_parameters)+second_index+1)
                for info in information_instances:
                    if info.has_second_param:
                        info.n, info.xedges, info.yedges = np.histogram2d(
                            info.chain[:, info.native_index+2],
                            info.chain[:, info.native_second_index+2],
                            weights=info.chain[:, 0],
                            bins=(info.bins, info.bins),
                            normed=False)
                        info.extent = [
                            info.x_range[info.native_second_index][0],
                            info.x_range[info.native_second_index][1],
                            info.x_range[info.native_index][0],
                            info.x_range[info.native_index][1]]
                        info.x_centers = 0.5*(info.xedges[1:]+info.xedges[:-1])
                        info.y_centers = 0.5*(info.yedges[1:]+info.yedges[:-1])

                        # plotting contours, using the ctr_level method (from Karim
                        # Benabed). Note that only the 1 and 2 sigma contours are
                        # displayed (due to the line with info.levels[:2])
                        try:
                            if info.contours_only:
                                contours = ax2dsub.contour(
                                    info.y_centers, info.x_centers, info.n,
                                    extent=info.extent, levels=ctr_level(
                                        info.n, info.levels[:2]),
                                    zorder=4, colors=info.cm[info.id],
                                    linewidths=info.line_width)
                            else:
                                contours = ax2dsub.contourf(
                                    info.y_centers, info.x_centers, info.n,
                                    extent=info.extent, levels=ctr_level(
                                        info.n, info.levels[:2]),
                                    zorder=4, cmap=info.cmaps[info.id],
                                    alpha=info.alphas[info.id])
                        except Warning:
                            warnings.warn(
                                "The routine could not find the contour of the " +
                                "'%s-%s' 2d-plot" % (
                                    info.plotted_parameters[info.native_index],
                                    info.plotted_parameters[info.native_second_index]))

                        ax2dsub.set_xticks(info.ticks[info.native_second_index])
                        if index == len(plotted_parameters)-1:
                            ax2dsub.set_xticklabels(
                                ['%.{0}g'.format(info.decimal) % s for s in
                                 info.ticks[info.native_second_index]],
                                fontsize=info.ticksize)
                            if conf.legend_style == 'sides':
                                ax2dsub.set_xlabel(
                                    info.tex_names[info.native_second_index],
                                    fontsize=info.fontsize)
                        else:
                            ax2dsub.set_xticklabels([''])

                        ax2dsub.set_yticks(info.ticks[info.native_index])
                        if second_index == 0:
                            ax2dsub.set_yticklabels(
                                ['%.{0}g'.format(info.decimal) % s for s in
                                 info.ticks[info.native_index]],
                                fontsize=info.ticksize)
                        else:
                            ax2dsub.set_yticklabels([''])

                        if conf.legend_style == 'sides':
                            if second_index == 0:
                                ax2dsub.set_ylabel(
                                    info.tex_names[info.native_index],
                                    fontsize=info.fontsize)

                if conf.subplot is True:
                    # Store the individual 2d plots.
                    if conf.plot_2d:
                        area = ax2dsub.get_window_extent().transformed(
                            fig2d.dpi_scale_trans.inverted())
                        # Pad the saved area by 10% in the x-direction and 20% in
                        # the y-direction
                        fig2d.savefig(os.path.join(
                            conf.folder, 'plots',
                            file_name+'_2d_%s-%s.%s' % (
                                standard_name, second_standard_name,
                                conf.extension)),
                            bbox_inches=area.expanded(1.4, 1.4))

                    # store the coordinates of the points for further
                    # plotting.
                    store_contour_coordinates(
                        conf, standard_name, second_standard_name, cotours)

                    for info in information_instances:
                        if not info.ignore_param and info.has_second_param:
                            info.hist_file_name = os.path.join(
                                info.folder, 'plots',
                                '{0}_2d_{1}-{2}.hist'.format(
                                    info.basename,
                                    standard_name,
                                    second_standard_name))
                    write_histogram_2d(
                        info.hist_file_name, info.x_centers, info.y_centers,
                        info.extent, info.n)

    print '-----------------------------------------------'
    if conf.plot:
        print '--> Saving figures to .{0} files'.format(info.extension)
        plot_name = '-vs-'.join([os.path.split(elem.folder)[-1]
                                for elem in information_instances])
        if conf.plot_2d:
            if len(legends) > 1:
                try:
                    fig2d.legend(legends, legend_names, 'upper right',
                                 fontsize=info.legendsize)
                except TypeError:
                    fig2d.legend(legends, legend_names, 'upper right',
                                 prop={'fontsize': info.legendsize})
            fig2d.tight_layout()
            fig2d.savefig(
                os.path.join(
                    conf.folder, 'plots', '{0}_triangle.{1}'.format(
                        plot_name, info.extension)),
                bbox_inches=0, )
        if conf.plot:
            fig1d.tight_layout()
            fig1d.savefig(
                os.path.join(
                    conf.folder, 'plots', '{0}_1d.{1}'.format(
                        plot_name, info.extension)),
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


def minimum_credible_intervals(info):
    """
    Extract minimum credible intervals (method from Jan Haman) FIXME
    """
    histogram = info.hist
    bincenters = info.bincenters
    levels = info.levels

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
                warnings.warn(
                    "please try running longer chains or reducing " +
                    "the number of bins with --bins BINS (default: 20)")
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
                failed = True
                break

        # min
        if failed:
            bounds[j][0] = np.nan
        elif indices[0] > 0:
            bounds[j][0] = bincenters[indices[0]] - delta*(histogram[indices[0]]-water_level)/(histogram[indices[0]]-histogram[indices[0]-1])
        else:
            if (left_edge > water_level):
                bounds[j][0] = bincenters[0]-0.5*delta
            else:
                bounds[j][0] = bincenters[indices[0]] - 0.5*delta*(histogram[indices[0]]-water_level)/(histogram[indices[0]]-left_edge)

        # max
        if failed:
            bounds[j][1] = np.nan
        elif indices[-1] < (len(histogram)-1):
            bounds[j][1] = bincenters[indices[-1]] + delta*(histogram[indices[-1]]-water_level)/(histogram[indices[-1]]-histogram[indices[-1]+1])
        else:
            if (right_edge > water_level):
                bounds[j][1] = bincenters[-1]+0.5*delta
            else:
                bounds[j][1] = bincenters[indices[-1]] + \
                    0.5*delta*(histogram[indices[-1]]-water_level) / \
                    (histogram[indices[-1]]-right_edge)

        j += 1

    for elem in bounds:
        for j in (0, 1):
            elem[j] -= info.mean[info.native_index]
    return bounds


def write_h(info_file, indices, name, string, quantity, modifiers=None):
    """
    Write one horizontal line of output

    """
    info_file.write('\n '+name+'\t: ')
    for i in indices:
        info_file.write(string % quantity[i]+'\t')


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


def separate_files(files):
    """
    Separate the input files in folder

    Given all input arguments to the command line files entry, separate them in
    a list of lists, grouping them by folders. The number of identified folders
    will determine the number of information instances to create
    """
    final_list = []
    temp = [files[0]]
    folder = (os.path.dirname(files[0]) if os.path.isfile(files[0])
              else files[0])
    if len(files) > 1:
        for elem in files[1:]:
            new_folder = (os.path.dirname(elem) if os.path.isfile(elem)
                          else elem)
            if new_folder == folder:
                temp.append(elem)
            else:
                folder = new_folder
                final_list.append(temp)
                temp = [elem]
    final_list.append(temp)

    return final_list


def recover_folder_and_files(files):
    """
    Distinguish the cases when analyze is called with files or folder

    Note that this takes place chronologically after the function
    `separate_files`"""
    # The following list defines the substring that a chain should contain for
    # the code to recognise it as a proper chain.
    substrings = ['.txt', '__']
    limit = 10
    # If the first element is a folder, grab all chain files inside
    if os.path.isdir(files[0]):
        folder = os.path.normpath(files[0])
        files = [os.path.join(folder, elem) for elem in os.listdir(folder)
                 if not os.path.isdir(os.path.join(folder, elem))
                 and not os.path.getsize(os.path.join(folder, elem)) < limit
                 and all([x in elem for x in substrings])]
    # Otherwise, extract the folder from the chain file-name.
    else:
        # If the name is completely wrong, say it
        if not os.path.exists(files[0]):
            raise io_mp.AnalyzeError(
                "You provided a non-existant folder/file to analyze")
        folder = os.path.relpath(
            os.path.dirname(os.path.realpath(files[0])), os.path.curdir)
        files = [os.path.join(folder, elem) for elem in os.listdir(folder)
                 if os.path.join(folder, elem) in np.copy(files)
                 and not os.path.isdir(os.path.join(folder, elem))
                 and not os.path.getsize(os.path.join(folder, elem)) < limit
                 and all([x in elem for x in substrings])]
    basename = os.path.basename(folder)
    return folder, files, basename


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
                        boundaries.append([
                            None if elem == 'None' or (isinstance(elem, int)
                                                       and elem == -1)
                            else elem for elem in array[1:3]])
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
    info.backup_names = backup_names
    info.scales = scales
    # Beware, the following two numbers are different. The first is the total
    # number of parameters stored in the chain, whereas the second is for
    # plotting purpose only.
    info.number_parameters = len(ref_names)
    info.plotted_parameters = plotted_parameters


def find_maximum_of_likelihood(info):
    """
    Finding the global maximum of likelihood

    min_minus_lkl will be appended with all the maximum likelihoods of files,
    then will be replaced by its own maximum. This way, the global
    maximum likelihood will be used as a reference, and not each chain's
    maximum.
    """
    min_minus_lkl = []
    for chain_file in info.files:
        # cheese will brutally contain everything (- log likelihood) in the
        # file chain_file being scanned.
        # This could potentially be faster with pandas, but is already quite
        # fast
        #
        # This would read the chains including comment lines:
        #cheese = (np.array([float(line.split()[1].strip())
        #                    for line in open(chain_file, 'r')]))
        #
        # This reads the chains excluding comment lines:
        with open(chain_file, 'r') as f:
            cheese = (np.array([float(line.split()[1].strip())
                                for line in ifilterfalse(iscomment,f)]))

        try:
            min_minus_lkl.append(cheese[:].min())
        except ValueError:
            pass
    # beware, it is the min because we are talking about
    # '- log likelihood'
    # Selecting only the true maximum.
    try:
        min_minus_lkl = min(min_minus_lkl)
    except ValueError:
        raise io_mp.AnalyzeError(
            "No decently sized chain was found in the desired folder. " +
            "Please wait to have more accepted point before trying " +
            "to analyze it.")

    info.min_minus_lkl = min_minus_lkl


def remove_bad_points(info):
    """
    Create an array with all the points from the chains, after removing non-markovian, burn-in and fixed fraction

    """
    # spam will brutally contain all the chains with sufficient number of
    # points, after the burn-in was removed.
    spam = list()

    # Recover the longest file name, for pleasing display
    max_name_length = max([len(e) for e in info.files])

    # Total number of steps done:
    steps = 0
    accepted_steps = 0

    # Open the log file
    log = open(info.log_path, 'w')

    for index, chain_file in enumerate(info.files):
        # To improve presentation, and print only once the full path of the
        # analyzed folder, we recover the length of the path name, and
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
        #
        # This would read the chains including comment lines:
        #cheese = (np.array([[float(elem) for elem in line.split()]
        #                    for line in open(chain_file, 'r')]))
        #
        # This read the chains excluding comment lines:
        with open(chain_file, 'r') as f:
            cheese = (np.array([[float(elem) for elem in line.split()]
                                for line in ifilterfalse(iscomment,f)]))
        # If the file contains a broken line with a different number of
        # elements, the previous array generation might fail, and will not have
        # the correct shape. Hence the following command will fail. To avoid
        # that, the error is caught.
        try:
            local_min_minus_lkl = cheese[:, 1].min()
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
        number_of_steps = cheese[:, 0].sum()
        log.write("%s\t " % os.path.basename(chain_file))
        log.write(" Number of steps:%d\t" % number_of_steps)
        log.write(" Steps accepted:%d\t" % line_count)
        log.write(" acc = %.2g\t" % (float(line_count)/number_of_steps))
        log.write("min(-loglike) = %.2f\n" % local_min_minus_lkl)
        steps += number_of_steps
        accepted_steps += line_count

        # check if analyze() is called directly by the user, or by the mcmc loop during an updating phase
        try:
            # command_line.update is defined when called by the mcmc loop
            info.update
        except:
            # in case it was not defined (i.e. when analyze() is called directly by user), set it to False
            info.update = 0

        # Removing non-markovian part, burn-in, and fraction= (1 - keep-fraction)
        start = 0
        markovian=0
        try:
            # Read all comments in chains about times when proposal was updated
            # The last of these comments gives the number of lines to be skipped in the files
            if info.markovian and not info.update:
                with open(chain_file, 'r') as f:
                    for line in ifilter(iscomment,f):
                        start = int(line.split()[2])
                markovian = start

            # Remove burn-in, defined as all points until the likelhood reaches min_minus_lkl+LOG_LKL_CUTOFF
            while cheese[start, 1] > info.min_minus_lkl+LOG_LKL_CUTOFF:
                start += 1
            burnin = start-markovian

            # Remove fixed fraction as requested by user (usually not useful if non-markovian is also removed)
            if info.keep_fraction < 1:
                start = start + (1-info.keep_fraction)*(line_count - start)

            print ": Removed",
            if info.markovian:
                print "%d non-markovian points," % markovian,
            print "%d points of burn-in," % burnin,
            if info.keep_fraction < 1:
                print "and first %.0f percent," % (100.*(1-info.keep_fraction)),
            print "keep %d steps" % (line_count-start)

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

    # Applying now new rules for scales, if the name is contained in the
    # referenced names
    for name in info.new_scales.iterkeys():
        try:
            index = info.ref_names.index(name)
            for i in xrange(len(spam)):
                spam[i][:, index+2] *= 1./info.scales[index, index]
        except ValueError:
            # there is nothing to do if the name is not contained in ref_names
            pass

    info.steps = steps
    info.accepted_steps = accepted_steps

    return spam


def compute_mean(mean, spam, total):
    """
    """
    for i in xrange(np.shape(mean)[1]):
        for j in xrange(len(spam)):
            submean = np.sum(spam[j][:, 0]*spam[j][:, i+2])
            mean[j+1, i] = submean / total[j+1]
            mean[0, i] += submean
        mean[0, i] /= total[0]


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


def compute_covariance_matrix(info):
    """
    """
    covar = np.zeros((len(info.ref_names), len(info.ref_names)))
    for i in xrange(len(info.ref_names)):
        for j in xrange(i, len(info.ref_names)):
            covar[i, j] = (
                info.chain[:, 0]*(
                    (info.chain[:, i+2]-info.mean[i]) *
                    (info.chain[:, j+2]-info.mean[j]))).sum()
            if i != j:
                covar[j, i] = covar[i, j]
    covar /= info.total

    # Removing scale factors in order to store true parameter covariance
    covar = np.dot(info.scales.T, np.dot(covar, info.scales))

    return covar


def adjust_ticks(param, information_instances):
    """
    """
    if len(information_instances) == 1:
        return
    # Recovering all x_range and ticks entries from the concerned information
    # instances
    x_ranges = []
    ticks = []
    for info in information_instances:
        if not info.ignore_param:
            x_ranges.append(info.x_range[info.native_index])
            ticks.append(info.ticks[info.native_index])

    # The new x_range and tick should min/max all the existing ones
    new_x_range = np.array(
        [min([e[0] for e in x_ranges]), max([e[1] for e in x_ranges])])
    temp_ticks = np.array(
        [min([e[0] for e in ticks]), max([e[-1] for e in ticks])])

    new_ticks = np.linspace(temp_ticks[0],
                            temp_ticks[1],
                            info.ticknumber)

    for info in information_instances:
        if not info.ignore_param:
            info.x_range[info.native_index] = new_x_range
            info.ticks[info.native_index] = new_ticks


def store_contour_coordinates(info, name1, name2, contours):
    """docstring"""
    file_name = os.path.join(
        info.folder, 'plots', '{0}_2d_{1}-{2}.dat'.format(
            info.basename, name1, name2))

    with open(file_name, 'w') as plot_file:
        plot_file.write(
            '# contour for confidence level {0}\n'.format(
                info.levels[1]))
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
                info.levels[0]))
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

def iscomment(s):
    """
    Define what we call a comment in MontePython chain files
    """
    return s.startswith('#')

class Information(object):
    """
    Hold all information for analyzing runs

    """
    # Counting the number of instances, to choose the color map
    _ids = count(0)
    # Flag checking the absence or presence of the interp1d function
    has_interpolate_module = False

    # Global colormap for the 1d plots. Colours will get chosen from this.
    # Some old versions of matplotlib do not have CMRmap, so the colours will
    # be harcoded
    # Note that, as with the other customisation options, you can specify new
    # values for this in the extra plot_file.
    cm = [
        (0.,      0.,      0.,      1.),
        (0.30235, 0.15039, 0.74804, 1.),
        (0.99843, 0.25392, 0.14765, 1.),
        (0.90000, 0.75353, 0.10941, 1.)]

    # Define colormaps for the contour plots
    cmaps = [plt.cm.gray_r, plt.cm.Purples, plt.cm.Reds_r, plt.cm.Greens]
    alphas = [1.0, 0.8, 0.6, 0.4]

    def __init__(self, command_line, other=None):
        """
        The following initialization creates the three tables that can be
        customized in an extra plot_file (see :mod:`parser_mp`).

        Parameters
        ----------
        command_line : Namespace
            it contains the initialised command line arguments
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
        # Assign a unique id to this instance
        self.id = self._ids.next()

        # Defining the sigma contours (1, 2 and 3-sigma)
        self.levels = np.array([68.26, 95.4, 99.7])/100.

        # Follows a bunch of initialisation to provide default members
        self.ref_names, self.backup_names = [], []
        self.scales, self.plotted_parameters = [], []
        self.spam = []

        # Store directly all information from the command_line object into this
        # instance, except the protected members (begin and end with __)
        for elem in dir(command_line):
            if elem.find('__') == -1:
                setattr(self, elem, getattr(command_line, elem))

        # initialize the legend size to be the same as fontsize, but can be
        # altered in the extra file
        self.legendsize = self.fontsize
        self.legendnames = []

        # Read a potential file describing changes to be done for the parameter
        # names, and number of paramaters plotted (can be let empty, all will
        # then be plotted), but also the style of the plot. Note that this
        # overrides the command line options
        if command_line.optional_plot_file:
            plot_file_vars = {'info': self}
            execfile(command_line.optional_plot_file, plot_file_vars)

        # check and store keep_fraction
        if command_line.keep_fraction<=0 or command_line.keep_fraction>1:
            raise io_mp.AnalyzeError("after --keep-fraction you should pass a float >0 and <=1")
        self.keep_fraction = command_line.keep_fraction

    def remap_parameters(self, spam):
        """
        Perform substitutions of parameters for analyzing

        .. note::

            for arbitrary combinations of parameters, the prior will not
            necessarily be flat.

        """
        if hasattr(self, 'redefine'):
            for key, value in self.redefine.iteritems():
                # Check that the key was an original name
                if key in self.backup_names:
                    print ' /|\  Transforming', key, 'into', value
                    # We recover the indices of the key
                    index_to_change = self.backup_names.index(key)+2
                    print('/_o_\ The new variable will be called ' +
                          self.ref_names[self.backup_names.index(key)])
                    # Recover all indices of all variables present in the
                    # remapping
                    variable_names = [elem for elem in self.backup_names if
                                      value.find(elem) != -1]
                    indices = [self.backup_names.index(name)+2 for name in
                               variable_names]
                    # Now loop over all files in spam
                    for i in xrange(len(spam)):
                        # Assign variables to their values
                        for index, name in zip(indices, variable_names):
                            exec("%s = spam[i][:, %i]" % (name, index))
                        # Assign to the desired index the combination
                        exec("spam[i][:, %i] = %s" % (index_to_change, value))

    def define_ticks(self):
        """
        """
        self.max_values = self.chain[:, 2:].max(axis=0)
        self.min_values = self.chain[:, 2:].min(axis=0)
        self.span = (self.max_values-self.min_values)
        # Define the place of ticks, given the number of ticks desired, stored
        # in conf.ticknumber
        self.ticks = np.array(
            [np.linspace(self.min_values[i]+self.span[i]*0.1,
                         self.max_values[i]-self.span[i]*0.1,
                         self.ticknumber) for i in range(len(self.span))])
        # Define the x range (ticks start not exactly at the range boundary to
        # avoid display issues)
        self.x_range = np.array((self.min_values, self.max_values)).T

        # In case the exploration hit a boundary (as defined in the parameter
        # file), at the level of precision defined by the number of bins, the
        # ticks and x_range should be altered in order to display this
        # meaningful number instead.
        for i in range(np.shape(self.ticks)[0]):
            x_range = self.x_range[i]
            bounds = self.boundaries[i]
            # Left boundary
            if bounds[0] is not None:
                if abs(x_range[0]-bounds[0]) < self.span[i]/self.bins:
                    self.ticks[i][0] = bounds[0]
                    self.x_range[i][0] = bounds[0]
            # Right boundary
            if bounds[-1] is not None:
                if abs(x_range[-1]-bounds[-1]) < self.span[i]/self.bins:
                    self.ticks[i][-1] = bounds[-1]
                    self.x_range[i][-1] = bounds[-1]

    def write_information_files(self):

        # Store in info_names only the tex_names that were plotted, for this
        # instance, and in indices the corresponding list of indices. It also
        # removes the $ signs, for clarity
        self.info_names = [
            name for index, name in enumerate(self.tex_names) if
            self.ref_names[index] in self.plotted_parameters]
        self.indices = [self.tex_names.index(name) for name in self.info_names]
        self.tex_names = [name for index, name in enumerate(self.tex_names) if
            self.ref_names[index] in self.plotted_parameters]
        self.info_names = [name.replace('$', '') for name in self.info_names]

        # Define the bestfit array
        self.bestfit = np.zeros(len(self.ref_names))
        for i in xrange(len(self.ref_names)):
            self.bestfit[i] = self.chain[self.sorted_indices[0], :][2+i]

        # Write down to the .h_info file all necessary information
        self.write_h_info()
        self.write_v_info()
        self.write_tex()

    def write_h_info(self):

        with open(self.h_info_path, 'w') as h_info:
            h_info.write(' param names\t:  ')
            for name in self.info_names:
                h_info.write("%-14s" % name)

            write_h(h_info, self.indices, 'R-1 values', '% .6f', self.R)
            write_h(h_info, self.indices, 'Best Fit  ', '% .6e', self.bestfit)
            write_h(h_info, self.indices, 'mean      ', '% .6e', self.mean)
            write_h(h_info, self.indices, 'sigma     ', '% .6e',
                    (self.bounds[:, 0, 1]-self.bounds[:, 0, 0])/2.)
            h_info.write('\n')
            write_h(h_info, self.indices, '1-sigma - ', '% .6e',
                    self.bounds[:, 0, 0])
            write_h(h_info, self.indices, '1-sigma + ', '% .6e',
                    self.bounds[:, 0, 1])
            write_h(h_info, self.indices, '2-sigma - ', '% .6e',
                    self.bounds[:, 1, 0])
            write_h(h_info, self.indices, '2-sigma + ', '% .6e',
                    self.bounds[:, 1, 1])
            write_h(h_info, self.indices, '3-sigma - ', '% .6e',
                    self.bounds[:, 2, 0])
            write_h(h_info, self.indices, '3-sigma + ', '% .6e',
                    self.bounds[:, 2, 1])

            # bounds
            h_info.write('\n')
            write_h(h_info, self.indices, '1-sigma > ', '% .6e',
                    self.mean+self.bounds[:, 0, 0])
            write_h(h_info, self.indices, '1-sigma < ', '% .6e',
                    self.mean+self.bounds[:, 0, 1])
            write_h(h_info, self.indices, '2-sigma > ', '% .6e',
                    self.mean+self.bounds[:, 1, 0])
            write_h(h_info, self.indices, '2-sigma < ', '% .6e',
                    self.mean+self.bounds[:, 1, 1])
            write_h(h_info, self.indices, '3-sigma > ', '% .6e',
                    self.mean+self.bounds[:, 2, 0])
            write_h(h_info, self.indices, '3-sigma < ', '% .6e',
                    self.mean+self.bounds[:, 2, 1])

    def write_v_info(self):
        """Write vertical info file"""
        with open(self.v_info_path, 'w') as v_info:
            v_info.write('%-15s\t:  %-11s' % ('param names', 'R-1'))
            v_info.write(' '.join(['%-11s' % elem for elem in [
                'Best fit', 'mean', 'sigma', '1-sigma -', '1-sigma +',
                '2-sigma -', '2-sigma +', '1-sigma >', '1-sigma <',
                '2-sigma >', '2-sigma <']]))
            for index, name in zip(self.indices, self.info_names):
                v_info.write('\n%-15s\t: % .4e' % (name, self.R[index]))
                v_info.write(' '.join(['% .4e' % elem for elem in [
                    self.bestfit[index], self.mean[index],
                    (self.bounds[index, 0, 1]-self.bounds[index, 0, 0])/2.,
                    self.bounds[index, 0, 0], self.bounds[index, 0, 1],
                    self.bounds[index, 1, 0], self.bounds[index, 1, 1],
                    self.mean[index]+self.bounds[index, 0, 0],
                    self.mean[index]+self.bounds[index, 0, 1],
                    self.mean[index]+self.bounds[index, 1, 0],
                    self.mean[index]+self.bounds[index, 1, 1]]]))

    def write_tex(self):
        """Write a tex table containing the main results """
        with open(self.tex_path, 'w') as tex:
            tex.write("\\begin{tabular}{|l|c|c|c|c|} \n \\hline \n")
            tex.write("Param & best-fit & mean$\pm\sigma$ ")
            tex.write("& 95\% lower & 95\% upper \\\\ \\hline \n")
            for index, name in zip(self.indices, self.tex_names):
                tex.write("%s &" % name)
                tex.write("$%.4g$ & $%.4g_{%.2g}^{+%.2g}$ " % (
                    self.bestfit[index], self.mean[index],
                    self.bounds[index, 0, 0], self.bounds[index, 0, 1]))
                tex.write("& $%.4g$ & $%.4g$ \\\\ \n" % (
                    self.mean[index]+self.bounds[index, 1, 0],
                    self.mean[index]+self.bounds[index, 1, 1]))

            tex.write("\\hline \n \\end{tabular} \\\\ \n")
            tex.write("$-\ln{\cal L}_\mathrm{min} =%.6g$, " % (
                self.min_minus_lkl))
            tex.write("minimum $\chi^2=%.4g$ \\\\ \n" % (
                self.min_minus_lkl*2.))
