# example of what you can do with the --extra files:

# use this to replace one column with a new parameter,
#defined as a function of one or more exititing parameters
info.redefine = {'omega_cdm': '(0.01*omega_b+omega_cdm)/(H0/100.)**2'}

# use this to rename a parameter (e.g. to make it look better in the labels).
# If you don't use dollars '$...$' the code will try automatically
# to convert your entry to latex (intepreting _, ^, greek letters...)
# but if you use '$...$' with plain latex in between, it will leave
# your input unchanged.
info.to_change = {'rs_d': 'r^s_d', 'omega_cdm': '$\Omega_\mathrm{m}$'}

# use this to change the scale factor normalising the parameters
info.new_scales = {'r^s_d': 100}

# use this to plot just a selection of parameters (if you have
# changed the names with 'info.to_change', you must put the new names here!)
info.to_plot = ['omega_b', 'Omega_m', 'H0', 'Omega_Lambda', 'r^s_d']

# decide whether to plot a legend or not
# (if not None, the legend will be added if
# there is more than one directoryto analyse)
info.plot_legend_1d = True
info.plot_legend_2d = False

# use this to customise the legend
# (one array entry for each plotted directory)
# The order here refers to the order in which you pass
# the directories to analyse.
info.legendnames = ['Hubble']

# If you want to change the order of colors
# (same order as for legendnames)
info.MP_color_cycle = [info.MP_color['Green'], info.MP_color['Orange'], info.MP_color['Blue']]

# You may actually even redefine the colors
# (each pair stands for [95% contour color,68% contour color]
# info.MP_color = {'Red':['#E37C80','#CE121F'],'Blue':['#7A98F6','#1157EF'],'Green':['#88B27A','#297C09'],'Orange':['#F3BE82','#ED920F'],'Grey':['#ABABAB','#737373'],'Purple':['#B87294','#88004C']}

# adjust the transparency of the lines and filled contours
# (same order as for legendnames)
info.alphas = [1.0, 0.8, 0.6, 0.4, 0.2]

# use this to control the boundaries of 1d and 2d plots
info.force_limits = {'H0':[60:70],'z_reio':[5:15]}
# use this to customise the ticks.
info.ticknumber = 5
info.ticksize = 10

# add list of python scripts for customisation of 1d or 2d plots,
# that are executed before plotting the probability lines or contours
info.custom1d = []
info.custom2d = ['add_h_contour.py','add_sigma8_Omegam_contour.py']
# any other lines of plain python can be written here without
# speacial formatting, they will be executed as extra
# lines of codes, but only at a precise point in the
# initialisation of the "Information" class
