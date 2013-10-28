"""
.. module:: prior
    :synopsis: Define the prior class

.. moduleauthor:: Benjamin Audren <benjamin.audren@epfl.ch>
"""
import random as rd
import numpy as np
import io_mp 
from copy import deepcopy

class prior(object):
    """
    Store the type of prior associated to a parameter

    """

    def __init__(self, array):
        """
        It takes as an optional input argument the array of the input
        :data:`parameters` defined in the parameter file.
        
        The current implemented types are 'flat' (default), and 'gaussian',
        which expect also a mean and sigma. Possible extension would take a
        'external', needing to read an external file to read for the
        definition.

        The entry 'prior' of the dictionary :data:`mcmc_parameters` will hold
        an instance of this class. It defines one main function, called
        :func:`draw_from_prior`, that returns a number within the prior volume.

        """

        rd.seed()

        # Test the length of the array, and initialize the type.
        if len(array) == 6:
            # Default behaviour, flat prior
            self.prior_type = 'flat'
        else:
            self.prior_type = array[6].lower()
            # in case of a gaussian prior, one expects two more entries, mu and
            # sigma
            if self.prior_type == 'gaussian':
                try:
                    self.mu = array[7]
                    self.sigma = array[8]
                except IndexError:
                    io_mp.message(
                        "You asked for a gaussian prior, but provided no \
                        mean nor sigma. Please add them in the parameter file",
                        "error")

        # Store boundaries for convenient access later
        self.prior_range = deepcopy(array[1:3])
        # Put all fields that are -1 to None to avoid confusion later on.
        for index, value in enumerate(self.prior_range):
            if str(value) == str(-1):
                value = None


    def draw_from_prior(self):
        """
        Draw a random point from the prior range considering the prior type

        :Returns:
            - **value** (`float`) - a random sample inside the prior region

        """

        if self.prior_type == 'flat':
            return rd.uniform(self.prior_range[0], self.prior_range[1])

        elif self.prior_type == 'gaussian':
            within_bounds = False

            while not within_bounds:
                value = rd.gauss(self.mu, self.sigma)
                # Check for boundaries problem 
                within_bounds = calue_within_prior_range(value)

            return value
                
    def value_within_prior_range(self, value):
        """
        Check for a value being in or outside the prior range

        """
        flag = 0
        if (self.prior_range[0] is not None and value < self.prior_range[0]):
            flag += 1
        elif (self.prior_range[1] is not None and value > self.prior_range[1]):
            flag += 1

        if flag > 0:
            return False
        else:
            return True
        



