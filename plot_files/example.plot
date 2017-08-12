# example of what you can do with the --extra files:

# use this to replace one column with a new parameter,
#defined as a function of one or more exititing parameters
info.redefine = {'omega_cdm': '(0.01*omega_b+omega_cdm)/(H0/100.)**2'}

# use this to rename a parameter (e.g. to make it look better
# in the labels). You can use plain latex here.
info.to_change = {'rs_d': 'r^s_d', 'omega_cdm': 'Omega_m'}

# use this to change the scale factor normalising the displayed paramneters
info.new_scales = {'r^s_d': 100}

# use this to plot just a selection of parameters
info.to_plot = ['omega_b', 'Omega_m', 'H0', 'Omega_Lambda', 'r^s_d']

# use this to customise the legend
# (one array entry for each plotted directory)
# The order here refers to the order in which you pass
# the directories to analyse.
info.legendnames = ['Hubble']

# use this to customise thew sequence of colors.
# Any matplotlib colormap is allowed, common ones are
# Greys,Purples,Blues,Greens,Oranges,Reds,
# but there are many others!
# (same order as for legendnames)
info.cmaps = [plt.cm.Greens, plt.cm.Oranges, plt.cm.Blues]

# adjust the transparency of the lines and filled contours
# (same order as for legendnames)
info.alphas = [1.0, 0.8, 0.6, 0.4, 0.2]

# use this to customise the ticks.
# you can write any plain python here, it will be executed as extra
# lines of codes, but only at a precise point in the initialisation of the
# "Information" class
info.ticknumber = 5
info.ticksize = 10
