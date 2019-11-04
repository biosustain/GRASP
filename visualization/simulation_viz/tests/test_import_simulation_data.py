import os
import unittest
import scipy.io
import pickle
import pandas as pd
from simulation_viz.import_simulation_data import gather_sim_data, get_met_rxn_names, get_time_series_quantiles, import_ref_conc


class TestImportSimulationData(unittest.TestCase):

    def setUp(self):
        self.this_dir, this_filename = os.path.split(__file__)
        self.raw_data_dir = os.path.join(self.this_dir, 'test_files')
        self.model_name = 'toy_model'

    def test_get_met_rxn_names(self):
        true_met_names = ['m_m5', 'm_m6', 'm_m7', 'm_m8', 'm_m9', 'm_m10', 'm_m11']
        true_rxn_names = ['r1', 'r2', 'r3', 'r4_1', 'r4_2', 'r5', 'r6', 'r7', 'r8', 'r9', 'r10', 'r11', 'r12', 'r13']

        met_names, rxn_names = get_met_rxn_names(self.raw_data_dir, self.model_name)

        self.assertListEqual(true_met_names, met_names)
        self.assertListEqual(true_rxn_names, rxn_names)

    def test_import_ref_conc(self):
        with open(os.path.join(self.this_dir, 'test_files', 'true_res_ref_conc_dic.pkl'), 'rb') as f_in:
            true_res = pickle.load(f_in)

        n_models = 5

        file_in = os.path.join(self.raw_data_dir, f'{self.model_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        n_mets = len(mat['ensemble']['mets'][0][0])
        all_met_names = [mat['ensemble']['mets'][0][0][met_i][0][0].replace('m_m_', '') for met_i in range(n_mets)]

        ref_conc_dic = import_ref_conc(mat, n_models, all_met_names)

        self.assertDictEqual(true_res, ref_conc_dic)

    def test_gather_sim_data(self):

        true_res_conc_interp = pd.read_hdf(os.path.join(self.this_dir, 'test_files', 'true_res_gather_sim_data_conc_interp.h5',), 'df')
        true_res_flux_interp = pd.read_hdf(os.path.join(self.this_dir, 'test_files', 'true_res_gather_sim_data_flux_interp.h5',), 'df')

        n_models = 5

        order_of_magnitude = [-9, -8, -7, -6, -5, -4, -3, -2, -1]
        time_points = [coefficient * 10 ** exponent for exponent in order_of_magnitude for coefficient in range(1, 10)]
        time_points.extend([1])
        time_points.insert(0, 0)

        met_names, rxn_names = get_met_rxn_names(self.raw_data_dir, self.model_name)

        file_in = os.path.join(self.raw_data_dir, f'{self.model_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        n_mets = len(mat['ensemble']['mets'][0][0])
        all_met_names = [mat['ensemble']['mets'][0][0][met_i][0][0].replace('m_m_', '') for met_i in range(n_mets)]

        ref_conc_dic = import_ref_conc(mat, n_models, all_met_names)

        simulation_name = f'{self.model_name}'
        file_in = os.path.join(self.raw_data_dir, f'simulation_{simulation_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        conc, conc_interp, flux, flux_interp = gather_sim_data(mat, met_names, rxn_names, n_models, time_points,
                                                               save_concs=False, save_fluxes=False,
                                                               ref_conc_dic=ref_conc_dic)

        self.assertTrue(pd.DataFrame.equals(true_res_conc_interp, conc_interp))
        self.assertTrue(pd.DataFrame.equals(true_res_flux_interp, flux_interp))

    def test_gather_sim_data_save_all(self):
        true_res_conc = pd.read_hdf(os.path.join(self.this_dir, 'test_files', 'true_res_gather_sim_data_save_all_conc.h5'), 'df')
        true_res_conc_interp = pd.read_hdf(os.path.join(self.this_dir, 'test_files', 'true_res_gather_sim_data_save_all_conc_interp.h5'), 'df')
        true_res_flux = pd.read_hdf(os.path.join(self.this_dir, 'test_files', 'true_res_gather_sim_data_save_all_flux.h5'), 'df')
        true_res_flux_interp = pd.read_hdf(os.path.join(self.this_dir, 'test_files', 'true_res_gather_sim_data_save_all_flux_interp.h5'), 'df')

        n_models = 5

        order_of_magnitude = [-9, -8, -7, -6, -5, -4, -3, -2, -1]
        time_points = [coefficient * 10 ** exponent for exponent in order_of_magnitude for coefficient in range(1, 10)]
        time_points.extend([1])
        time_points.insert(0, 0)

        met_names, rxn_names = get_met_rxn_names(self.raw_data_dir, self.model_name)

        file_in = os.path.join(self.raw_data_dir, f'{self.model_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        n_mets = len(mat['ensemble']['mets'][0][0])
        all_met_names = [mat['ensemble']['mets'][0][0][met_i][0][0].replace('m_m_', '') for met_i in range(n_mets)]

        ref_conc_dic = import_ref_conc(mat, n_models, all_met_names)

        simulation_name = f'{self.model_name}'
        file_in = os.path.join(self.raw_data_dir, f'simulation_{simulation_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        conc, conc_interp, flux, flux_interp = gather_sim_data(mat, met_names, rxn_names, n_models, time_points,
                                                               save_concs=True, save_fluxes=True,
                                                               ref_conc_dic=ref_conc_dic)

        self.assertTrue(pd.DataFrame.equals(true_res_conc, conc))
        self.assertTrue(pd.DataFrame.equals(true_res_conc_interp, conc_interp))
        self.assertTrue(pd.DataFrame.equals(true_res_flux, flux))
        self.assertTrue(pd.DataFrame.equals(true_res_flux_interp, flux_interp))

    def test_get_time_series_quantiles(self):
        true_res_conc = pd.read_hdf(os.path.join(self.this_dir, 'test_files', 'true_res_get_time_series_quantiles_conc.h5'), 'df')
        true_res_flux = pd.read_hdf(os.path.join(self.this_dir, 'test_files', 'true_res_get_time_series_quantiles_flux.h5'), 'df')

        n_models = 5

        order_of_magnitude = [-9, -8, -7, -6, -5, -4, -3, -2, -1]
        time_points = [coefficient * 10 ** exponent for exponent in order_of_magnitude for coefficient in range(1, 10)]
        time_points.extend([1])
        time_points.insert(0, 0)

        met_names, rxn_names = get_met_rxn_names(self.raw_data_dir, self.model_name)

        file_in = os.path.join(self.raw_data_dir, f'{self.model_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        n_mets = len(mat['ensemble']['mets'][0][0])
        all_met_names = [mat['ensemble']['mets'][0][0][met_i][0][0].replace('m_m_', '') for met_i in range(n_mets)]

        ref_conc_dic = import_ref_conc(mat, n_models, all_met_names)

        simulation_name = f'{self.model_name}'
        file_in = os.path.join(self.raw_data_dir, f'simulation_{simulation_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        conc, conc_interp, flux, flux_interp = gather_sim_data(mat, met_names, rxn_names, n_models, time_points,
                                                               save_concs=True, save_fluxes=True,
                                                               ref_conc_dic=ref_conc_dic)

        data_type = 'conc'
        conc_interp_quantiles_df = get_time_series_quantiles(conc_interp, time_points, data_type, met_names)

        data_type = 'flux'
        flux_interp_quantiles_df = get_time_series_quantiles(flux_interp, time_points, data_type, rxn_names)

        self.assertTrue(pd.DataFrame.equals(true_res_conc, conc_interp_quantiles_df))
        self.assertTrue(pd.DataFrame.equals(true_res_flux, flux_interp_quantiles_df))