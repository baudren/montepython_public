# there is no specific lieklihood code for this experiment, because it
# falls in the category of CMB experiments described in the "newdat"
# format. The class below inherits the properties of a general class
# "Likelihood_newdat", which knows how to deal with all experiments in
# "newdat" format.

from montepython.likelihood_class import Likelihood_newdat


class boomerang(Likelihood_newdat):
    pass
