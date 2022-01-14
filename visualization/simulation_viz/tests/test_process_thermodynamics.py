import io
import os
import unittest.mock

import pandas as pd

from simulation_viz.process_thermodynamics import get_fluxes_and_dGs, calculate_dG


class TestThermodynamicsChecks(unittest.TestCase):

    def setUp(self):
        this_dir, this_filename = os.path.split(__file__)
        self.test_folder = os.path.join(this_dir, 'test_files', 'test_thermodynamics_checks')

    def test_get_fluxes_and_dGs_modelv1(self):
        true_res_dG = pd.read_hdf(os.path.join(self.test_folder, 'true_res_get_fluxes_and_dGs_modelv1_dG.h5'), 'df')
        true_res_flux = pd.read_hdf(os.path.join(self.test_folder, 'true_res_get_fluxes_and_dGs_modelv1_flux.h5'), 'df')

        data_dict = pd.read_excel(os.path.join(self.test_folder, 'model_v1_base.xlsx'), sheet_name=None, index_col=0)
        flux_df, dG_df = get_fluxes_and_dGs(data_dict)

        self.assertTrue(dG_df.equals(true_res_dG))
        self.assertTrue(flux_df.equals(true_res_flux))

    def test_get_fluxes_and_dGs_HMP(self):
        true_res_dG = pd.read_hdf(os.path.join(self.test_folder, 'true_res_get_fluxes_and_dGs_HMP_dG.h5'), 'df')
        true_res_flux = pd.read_hdf(os.path.join(self.test_folder, 'true_res_get_fluxes_and_dGs_HMP_flux.h5'), 'df')

        data_dict = pd.read_excel(os.path.join(self.test_folder, 'toy_model_promiscuous.xlsx'), sheet_name=None, index_col=0)
        flux_df, dG_df = get_fluxes_and_dGs(data_dict)

        self.assertTrue(dG_df.equals(true_res_dG))
        self.assertTrue(flux_df.equals(true_res_flux))

    def test_calculate_dG(self):
        true_res_dG = pd.read_hdf(os.path.join(self.test_folder, 'true_res_calculate_dG.h5'), 'df')

        temperature = 298  # in K
        gas_constant = 8.314 * 10 ** -3  # in kJ K^-1 mol^-1

        data_dict = pd.read_excel(os.path.join(self.test_folder, 'model_v2_manual.xlsx'), sheet_name=None, index_col=0)
        dG_df = calculate_dG(data_dict, gas_constant, temperature)

        self.assertTrue(dG_df.equals(true_res_dG))
