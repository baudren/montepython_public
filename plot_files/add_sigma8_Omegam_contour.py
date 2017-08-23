# add 68% and 95% contours corresponding to
# a gaussian prior/likelihood on
# sigma_8 (Omega_m/a)^b

if ((second_name == 'sigma_8') and (name == 'Omega_\mathrm{m}')) or (name == 'sigma_8') and (second_name == 'Omega_\mathrm{m}'):

    center = 0.782
    sigma = 0.010
    a = 0.27
    b = 0.30
    contour_color = info.MP_color['Purple']
    contour_alpha = 0.5

    # add contours when sigma_8 is on the x axis
    if (second_name == 'sigma_8') and (name == 'Omega_\mathrm{m}'):
        x1 = info.extent[0]
        x2 = info.extent[1]
        xx = np.arange(x1,x2,(x2-x1)/20.)
        yy68_1 = a*np.power((center-sigma)/xx,1./b)
        yy68_2 = a*np.power((center+sigma)/xx,1./b)
        yy95_1 = a*np.power((center-2.*sigma)/xx,1./b)
        yy95_2 = a*np.power((center+2.*sigma)/xx,1./b)
        ax2dsub.fill_between(xx,
                             yy95_1,
                             yy95_2,
                             facecolor=contour_color[0],
                             edgecolor=contour_color[1],
                             linewidth=1,
                             alpha=contour_alpha)
        ax2dsub.fill_between(xx,
                             yy68_1,
                             yy68_2,
                             color=contour_color[1],
                             alpha=contour_alpha)

    # add contours when sigma_8 is on the y axis
    if (name == 'sigma_8') and (second_name == 'Omega_\mathrm{m}'):
        x1 = info.extent[0]
        x2 = info.extent[1]
        xx = np.arange(x1,x2,(x2-x1)/10.)
        yy68_1 = (center-sigma)/np.power(xx/a,b)
        yy68_2 = (center+sigma)/np.power(xx/a,b)
        yy95_1 = (center-2.*sigma)/np.power(xx/a,b)
        yy95_2 = (center+2.*sigma)/np.power(xx/a,b)
        ax2dsub.fill_between(xx,
                             yy95_1,
                             yy95_2,
                             facecolor=contour_color[0],
                             edgecolor=contour_color[1],
                             linewidth=1,
                             alpha=contour_alpha)
        ax2dsub.fill_between(xx,
                             yy68_1,
                             yy68_2,
                             color=contour_color[1],
                             alpha=contour_alpha)
