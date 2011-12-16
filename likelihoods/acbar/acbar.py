# there is no specific lieklihood code for this experiment, because it falls in the category of CMB experiments described in the "newdat" format. The class below inherits the properties of a general class "likelihood_newdat", which knows how to deal with all experiments in "newdat" format.

from likelihood_class import likelihood_newdat

class acbar(likelihood_newdat):

  def __init__(self,path,data,command_line=False):
    likelihood_newdat.__init__(self,path,data,command_line)
    #self._need_Class_args(data,{'l_max_scalars':3500})
    self._need_nuisance_parameters(data,['A_SZ'])

  def _loglkl(self,_cosmo,data):
    cl = likelihood_newdat._get_cl(_cosmo)
    
    # things with nuisance

    return likelihood_newdat._compute_loglkl(cl,_cosmo,data)

