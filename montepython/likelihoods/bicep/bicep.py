# there is no specific lieklihood code for this experiment, because it
# falls in the category of CMB experiments described in the "newdat"
# format. The class below inherits the properties of a general class
# "likelihood_newdat", which knows how to deal with all experiments in
# "newdat" format.

from likelihood_class import likelihood_newdat


class bicep(likelihood_newdat):
    pass
