import os
import pickle
import unittest

import pandas as pd
import scipy.io

from simulation_viz.simulation_viz.import_simulation_data import gather_sim_data, get_met_rxn_names, \
    get_time_series_quantiles, import_ref_conc, load_simulation, import_ref_flux


class TestImportSimulationData(unittest.TestCase):

    def setUp(self):
        self.this_dir, this_filename = os.path.split(__file__)
        self.raw_data_dir = os.path.join(self.this_dir, 'test_files')
        self.model_name = 'toy_model'

    def test_get_met_rxn_names(self):
        true_met_names = ['m5', 'm6', 'm7', 'm8', 'm9', 'm10', 'm11']
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

        ref_conc_dic = import_ref_conc(mat, n_models)

        self.assertDictEqual(true_res, ref_conc_dic)

    def test_import_ref_flux(self):
        with open(os.path.join(self.this_dir, 'test_files', 'true_res_ref_flux_dic.pkl'), 'rb') as f_in:
            true_res = pickle.load(f_in)

        file_in = os.path.join(self.raw_data_dir, f'{self.model_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        ref_flux_dic = import_ref_flux(mat)

        self.assertDictEqual(true_res, ref_flux_dic)

    def test_gather_sim_data(self):
        true_res_conc_interp = pd.read_hdf(
            os.path.join(self.this_dir, 'test_files', 'true_res_gather_sim_data_conc_interp.h5', ), 'df')
        true_res_flux_interp = pd.read_hdf(
            os.path.join(self.this_dir, 'test_files', 'true_res_gather_sim_data_flux_interp.h5', ), 'df')

        n_models = 5

        order_of_magnitude = [-9, -8, -7, -6, -5, -4, -3, -2, -1]
        time_points = [coefficient * 10 ** exponent for exponent in order_of_magnitude for coefficient in range(1, 10)]
        time_points.extend([1])
        time_points.insert(0, 0)

        met_names, rxn_names = get_met_rxn_names(self.raw_data_dir, self.model_name)

        file_in = os.path.join(self.raw_data_dir, f'{self.model_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        ref_conc_dic = import_ref_conc(mat, n_models)
        ref_flux_dic = import_ref_flux(mat)

        simulation_name = f'{self.model_name}'
        file_in = os.path.join(self.raw_data_dir, f'simulation_{simulation_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        conc, conc_interp, flux, flux_interp = gather_sim_data(mat, met_names, rxn_names, n_models, time_points,
                                                               save_concs=False, save_fluxes=False,
                                                               ref_conc_dic=ref_conc_dic, ref_flux_dic=ref_flux_dic)

        self.assertTrue(pd.DataFrame.equals(true_res_conc_interp, conc_interp))
        self.assertTrue(pd.DataFrame.equals(true_res_flux_interp, flux_interp))

    def test_gather_sim_data_save_all(self):
        true_res_conc = pd.read_hdf(
            os.path.join(self.this_dir, 'test_files', 'true_res_gather_sim_data_save_all_conc.h5'), 'df')
        true_res_conc_interp = pd.read_hdf(
            os.path.join(self.this_dir, 'test_files', 'true_res_gather_sim_data_save_all_conc_interp.h5'), 'df')
        true_res_flux = pd.read_hdf(
            os.path.join(self.this_dir, 'test_files', 'true_res_gather_sim_data_save_all_flux.h5'), 'df')
        true_res_flux_interp = pd.read_hdf(
            os.path.join(self.this_dir, 'test_files', 'true_res_gather_sim_data_save_all_flux_interp.h5'), 'df')

        n_models = 5

        order_of_magnitude = [-9, -8, -7, -6, -5, -4, -3, -2, -1]
        time_points = [coefficient * 10 ** exponent for exponent in order_of_magnitude for coefficient in range(1, 10)]
        time_points.extend([1])
        time_points.insert(0, 0)

        met_names, rxn_names = get_met_rxn_names(self.raw_data_dir, self.model_name)

        file_in = os.path.join(self.raw_data_dir, f'{self.model_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        ref_conc_dic = import_ref_conc(mat, n_models)
        ref_flux_dic = import_ref_flux(mat)

        simulation_name = f'{self.model_name}'
        file_in = os.path.join(self.raw_data_dir, f'simulation_{simulation_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        conc, conc_interp, flux, flux_interp = gather_sim_data(mat, met_names, rxn_names, n_models, time_points,
                                                               save_concs=True, save_fluxes=True,
                                                               ref_conc_dic=ref_conc_dic, ref_flux_dic=ref_flux_dic)

        self.assertTrue(pd.DataFrame.equals(true_res_conc, conc))
        self.assertTrue(pd.DataFrame.equals(true_res_conc_interp, conc_interp))
        self.assertTrue(pd.DataFrame.equals(true_res_flux, flux))
        self.assertTrue(pd.DataFrame.equals(true_res_flux_interp, flux_interp))

    def test_get_time_series_quantiles_abs(self):
        true_res_conc = pd.read_hdf(
            os.path.join(self.this_dir, 'test_files', 'true_res_get_time_series_quantiles_conc_abs.h5'), 'df')
        true_res_flux = pd.read_hdf(
            os.path.join(self.this_dir, 'test_files', 'true_res_get_time_series_quantiles_flux_abs.h5'), 'df')

        n_models = 5

        order_of_magnitude = [-9, -8, -7, -6, -5, -4, -3, -2, -1]
        time_points = [coefficient * 10 ** exponent for exponent in order_of_magnitude for coefficient in range(1, 10)]
        time_points.extend([1])
        time_points.insert(0, 0)

        met_names, rxn_names = get_met_rxn_names(self.raw_data_dir, self.model_name)

        file_in = os.path.join(self.raw_data_dir, f'{self.model_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        ref_conc_dic = import_ref_conc(mat, n_models)
        ref_flux_dic = import_ref_flux(mat)

        simulation_name = f'{self.model_name}'
        file_in = os.path.join(self.raw_data_dir, f'simulation_{simulation_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        conc, conc_interp, flux, flux_interp = gather_sim_data(mat, met_names, rxn_names, n_models, time_points,
                                                               save_concs=True, save_fluxes=True,
                                                               ref_conc_dic=ref_conc_dic, ref_flux_dic=ref_flux_dic)

        data_type = 'conc_abs'
        conc_interp_quantiles_df = get_time_series_quantiles(conc_interp, time_points, data_type, met_names)

        data_type = 'flux_abs'
        flux_interp_quantiles_df = get_time_series_quantiles(flux_interp, time_points, data_type, rxn_names)

        self.assertTrue(pd.DataFrame.equals(true_res_conc, conc_interp_quantiles_df))
        self.assertTrue(pd.DataFrame.equals(true_res_flux, flux_interp_quantiles_df))

    def test_get_time_series_quantiles_rel(self):
        true_res_conc = pd.read_hdf(
            os.path.join(self.this_dir, 'test_files', 'true_res_get_time_series_quantiles_conc_rel.h5'), 'df')
        true_res_flux = pd.read_hdf(
            os.path.join(self.this_dir, 'test_files', 'true_res_get_time_series_quantiles_flux_rel.h5'), 'df')

        n_models = 5

        order_of_magnitude = [-9, -8, -7, -6, -5, -4, -3, -2, -1]
        time_points = [coefficient * 10 ** exponent for exponent in order_of_magnitude for coefficient in range(1, 10)]
        time_points.extend([1])
        time_points.insert(0, 0)

        met_names, rxn_names = get_met_rxn_names(self.raw_data_dir, self.model_name)

        file_in = os.path.join(self.raw_data_dir, f'{self.model_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        ref_conc_dic = import_ref_conc(mat, n_models)
        ref_flux_dic = import_ref_flux(mat)

        simulation_name = f'{self.model_name}'
        file_in = os.path.join(self.raw_data_dir, f'simulation_{simulation_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        conc, conc_interp, flux, flux_interp = gather_sim_data(mat, met_names, rxn_names, n_models, time_points,
                                                               save_concs=True, save_fluxes=True,
                                                               ref_conc_dic=ref_conc_dic, ref_flux_dic=ref_flux_dic)

        data_type = 'conc_rel'
        conc_interp_quantiles_df = get_time_series_quantiles(conc_interp, time_points, data_type, met_names)

        data_type = 'flux_rel'
        flux_interp_quantiles_df = get_time_series_quantiles(flux_interp, time_points, data_type, rxn_names)

        self.assertTrue(pd.DataFrame.equals(true_res_conc, conc_interp_quantiles_df))
        self.assertTrue(pd.DataFrame.equals(true_res_flux, flux_interp_quantiles_df))

    def test_load_simulations(self):
        true_res_conc = pd.read_hdf(
            os.path.join(self.this_dir, 'test_files', 'true_res_gather_sim_data_save_all_conc.h5'), 'df')
        true_res_conc_interp = pd.read_hdf(
            os.path.join(self.this_dir, 'test_files', 'true_res_gather_sim_data_save_all_conc_interp.h5'), 'df')
        true_res_flux = pd.read_hdf(
            os.path.join(self.this_dir, 'test_files', 'true_res_gather_sim_data_save_all_flux.h5'), 'df')
        true_res_flux_interp = pd.read_hdf(
            os.path.join(self.this_dir, 'test_files', 'true_res_gather_sim_data_save_all_flux_interp.h5'), 'df')

        true_met_names = ['m5', 'm6', 'm7', 'm8', 'm9', 'm10', 'm11']
        true_rxn_names = ['r1', 'r2', 'r3', 'r4_1', 'r4_2', 'r5', 'r6', 'r7', 'r8', 'r9', 'r10', 'r11', 'r12', 'r13']

        n_models = 5

        order_of_magnitude = [-9, -8, -7, -6, -5, -4, -3, -2, -1]
        time_points = [coefficient * 10 ** exponent for exponent in order_of_magnitude for coefficient in range(1, 10)]
        time_points.extend([1])
        time_points.insert(0, 0)

        simulation_name = f'{self.model_name}'

        conc, conc_interp, flux, flux_interp, met_names, rxn_names = load_simulation(self.raw_data_dir, self.model_name,
                                                                                     simulation_name, n_models,
                                                                                     time_points, save_concs=True,
                                                                                     save_fluxes=True)

        self.assertListEqual(true_met_names, met_names)
        self.assertListEqual(true_rxn_names, rxn_names)
        self.assertTrue(pd.DataFrame.equals(true_res_conc, conc))
        self.assertTrue(pd.DataFrame.equals(true_res_conc_interp, conc_interp))
        self.assertTrue(pd.DataFrame.equals(true_res_flux, flux))
        self.assertTrue(pd.DataFrame.equals(true_res_flux_interp, flux_interp))
