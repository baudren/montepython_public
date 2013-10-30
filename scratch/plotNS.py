import sys, os
import pymultinest
import matplotlib.pyplot as plt
plt.clf()

print "Use as 'python plotNS.py <folder> <parameter1> ...'"

basename = str(sys.argv[1])
if basename[-1] != "/":
    basename += "/"

parameters = []
for i in range(2,len(sys.argv)) :
    parameters.append(sys.argv[i])
n_params = len(parameters)

a = pymultinest.Analyzer(n_params = n_params, outputfiles_basename = basename)
s = a.get_stats()

# Here we will plot all the marginals and whatnot, just to show off
# You may configure the format of the output here, or in matplotlibrc
# All pymultinest does is filling in the data of the plot.

# Copy and edit this file, and play with it.

p = pymultinest.PlotMarginalModes(a)
plt.figure(figsize=(5*n_params, 5*n_params))
#plt.subplots_adjust(wspace=0, hspace=0)
for i in range(n_params):
        plt.subplot(n_params, n_params, n_params * i + i + 1)
        p.plot_marginal(i, with_ellipses = True, with_points = False, grid_points=50)
        plt.ylabel("Probability")
        plt.xlabel(parameters[i])
        
        for j in range(i):
                plt.subplot(n_params, n_params, n_params * j + i + 1)
                #plt.subplots_adjust(left=0, bottom=0, right=0, top=0, wspace=0, hspace=0)
                p.plot_conditional(i, j, with_ellipses = False, with_points = True, grid_points=30)
                plt.xlabel(parameters[i])
                plt.ylabel(parameters[j])

plt.savefig("marginals_multinest.pdf") #, bbox_inches='tight')
plt.show("marginals_multinest.pdf")

for i in range(n_params):
        outfile = '%s-mode-marginal-%d.pdf' % (a.outputfiles_basename,i)
        p.plot_modes_marginal(i, with_ellipses = True, with_points = False)
        plt.ylabel("Probability")
        plt.xlabel(parameters[i])
        plt.savefig(outfile, format='pdf', bbox_inches='tight')
        plt.close()
        
        outfile = '%s-mode-marginal-cumulative-%d.pdf' % (a.outputfiles_basename,i)
        p.plot_modes_marginal(i, cumulative = True, with_ellipses = True, with_points = False)
        plt.ylabel("Cumulative probability")
        plt.xlabel(parameters[i])
        plt.savefig(outfile, format='pdf', bbox_inches='tight')
        plt.close()

print("Take a look at the pdf files in chains/") 
