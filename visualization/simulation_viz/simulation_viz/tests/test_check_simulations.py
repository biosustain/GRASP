import os
import unittest

from simulation_viz.check_simulations import check_for_negative_concentrations


class TestImportSimulationData(unittest.TestCase):

    def setUp(self):
        this_dir, this_filename = os.path.split(__file__)

    def test_check_for_negative_concentrations(self):
        return 0
