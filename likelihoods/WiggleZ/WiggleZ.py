# likelihood for WiggleZ
# The structure differs slightly from other likelihoods, to keep a simple
# structure. Each redshift bin in WiggleZ contains a .dataset file with
# information on this redshift bin. The structure to read this has been endoded
# in the class likelihood_mpk. It will be used also with the next data release
# of SDSS.
# The whole WiggleZ is then made out of the four likelihood_mpk, WiggleZ_a, b, c
# and d
from likelihood_class import likelihood_mpk,likelihood

class WiggleZ_a(likelihood_mpk):
  pass
class WiggleZ_b(likelihood_mpk):
  pass
class WiggleZ_c(likelihood_mpk):
  pass
class WiggleZ_d(likelihood_mpk):
  pass

class WiggleZ(likelihood):
  
  def __init__(self,path,data,command_line,log_flag,default):

    likelihood.__init__(self,path,data,command_line,log_flag,default)

    # Initialize one after the other the four independant redshift bins (note:
    # the order in the array self.redshift_bins_files) must be respected !
    self.wigglez_a = WiggleZ_a(self.data_directory+self.redshift_bins_files[0],data,command_line,log_flag,default)
    self.wigglez_b = WiggleZ_b(self.data_directory+self.redshift_bins_files[1],data,command_line,log_flag,default)
    self.wigglez_c = WiggleZ_c(self.data_directory+self.redshift_bins_files[2],data,command_line,log_flag,default)
    self.wigglez_d = WiggleZ_d(self.data_directory+self.redshift_bins_files[3],data,command_line,log_flag,default)

  def loglkl(self,_cosmo,data):
    # Simply add all the sublikelihoods
    loglkl = 0
    loglkl += self.wigglez_a.loglkl(_cosmo,data)
    loglkl += self.wigglez_b.loglkl(_cosmo,data)
    loglkl += self.wigglez_c.loglkl(_cosmo,data)
    loglkl += self.wigglez_d.loglkl(_cosmo,data)
    return loglkl

