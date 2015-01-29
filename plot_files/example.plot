info.redefine = {'omega_cdm': '(0.01*omega_b+omega_cdm)/(H0/100.)**2'}

info.to_change = {'rs_d': 'r^s_d', 'omega_cdm': 'Omega_m'}
info.new_scales = {'r^s_d': 100}
info.to_plot = ['omega_b', 'Omega_m', 'H0', 'Omega_Lambda', 'r^s_d']

info.legendnames = ['Hubble']

info.ticknumber = 5
info.ticksize = 10
