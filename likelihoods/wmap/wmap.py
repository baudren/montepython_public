from likelihood_class import likelihood

class wmap(likelihood):

  def __init__(self,path,data,command_line,log_flag,default):

    # Standard initialization, reads the .data
    likelihood.__init__(self,path,data,command_line,log_flag,default)

    # Extra definitions....
    # Extra needed Class paramters
    self.need_Class_arguments(data,{'output':'mPk'})

    # In case of a comparison, stop here
    if not default:
      return

    # More definitions that do not need to be compared between two calls
    # (independent of the experiment
    pass
  
  def loglkl(self,_cosmo,data):
    chi2 = 0.
    return -chi2/2
