import os
from likelihood_class import likelihood_prior


class hst(likelihood_prior):

    # initialisation of the class is done within the parent likelihood_prior. For
    # this case, it does not differ, actually, from the __init__ method in
    # likelihood class.
    def loglkl(self, cosmo, data):

        h = cosmo.h()
        loglkl = -0.5 * (h - self.h) ** 2 / (self.sigma ** 2)
        return loglkl
