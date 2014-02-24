# there is no specific likelihood code for this experiment, because it
# falls in the category of CMB experiments described in the "mock CMB"
# format. The class below inherits the properties of a general class
# "likelihood_mock_cmb", which knows how to deal with all experiments in
# "mock CMB" format.

from likelihood_class import likelihood_mock_cmb


class fake_planck_bluebook(likelihood_mock_cmb):
    pass
