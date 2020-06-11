import os
import unittest

import scipy.io

from simulation_viz.simulation_viz.import_simulation_data import gather_sim_data, get_met_rxn_names, \
    get_time_series_quantiles, import_ref_conc
from simulation_viz.simulation_viz.visualize_simulations import plot_ensemble, plot_model


class TestVisualizeSimulations(unittest.TestCase):

    def setUp(self):
        self.this_dir, this_filename = os.path.split(__file__)
        raw_data_dir = os.path.join(self.this_dir, 'test_files')
        model_name = 'toy_model'

        n_models = 5

        order_of_magnitude = [-9, -8, -7, -6, -5, -4, -3, -2, -1]
        time_points = [coefficient * 10 ** exponent for exponent in order_of_magnitude for coefficient in range(1, 10)]
        time_points.extend([1])
        time_points.insert(0, 0)

        self.met_names, self.rxn_names = get_met_rxn_names(raw_data_dir, model_name)

        file_in = os.path.join(raw_data_dir, f'{model_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        n_mets = len(mat['ensemble']['mets'][0][0])
        all_met_names = [mat['ensemble']['mets'][0][0][met_i][0][0].replace('m_m_', '') for met_i in range(n_mets)]

        ref_conc_dic = import_ref_conc(mat, n_models, all_met_names)

        simulation_name = f'{model_name}'
        file_in = os.path.join(raw_data_dir, f'simulation_{simulation_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        self.conc, self.conc_interp, self.flux, self.flux_interp = gather_sim_data(mat, self.met_names, self.rxn_names,
                                                                                   n_models,
                                                                                   time_points, save_concs=True,
                                                                                   save_fluxes=True,
                                                                                   ref_conc_dic=ref_conc_dic)

        data_type = 'conc'
        self.conc_interp_quantiles = get_time_series_quantiles(self.conc_interp, time_points, data_type, self.met_names)

        data_type = 'flux'
        self.flux_interp_quantiles = get_time_series_quantiles(self.flux_interp, time_points, data_type, self.rxn_names)

    def test_plot_ensemble_conc_1row(self):
        output_file = os.path.join(self.this_dir, 'test_files', 'plot_ensemble_conc_viz.png')
        plot_ensemble(self.conc_interp_quantiles, quant_type='conc', selected_data=[self.met_names],
                      x_scale='linear', y_scale='log', x_lim=None, y_lim=None, fig_size=None, save_plot=True,
                      output_file=output_file)

        self.assertTrue(os.path.isfile(output_file))
        os.remove(output_file)

    def test_plot_ensemble_conc_2rows(self):
        output_file = os.path.join(self.this_dir, 'test_files', 'plot_ensemble_conc_viz2.png')
        met_list1 = ['m_m5', 'm_m6']
        met_list2 = ['m_m10', 'm_m11']

        plot_ensemble(self.conc_interp_quantiles, quant_type='conc', selected_data=[met_list1, met_list2],
                      x_scale='linear', y_scale='log', x_lim=None, y_lim=None, fig_size=None, save_plot=True,
                      output_file=output_file)

        self.assertTrue(os.path.isfile(output_file))
        os.remove(output_file)

    def test_plot_ensemble_flux_1row(self):
        output_file = os.path.join(self.this_dir, 'test_files', 'plot_ensemble_flux_viz.png')
        plot_ensemble(self.flux_interp_quantiles, quant_type='flux', selected_data=[self.rxn_names],
                      x_scale='linear', y_scale='log', x_lim=None, y_lim=None, fig_size=None, save_plot=True,
                      output_file=output_file)

        self.assertTrue(os.path.isfile(output_file))
        os.remove(output_file)

    def test_plot_ensemble_flux_2rows(self):
        output_file = os.path.join(self.this_dir, 'test_files', 'plot_ensemble_flux_viz2.png')
        rxn_list1 = ['r1', 'r2']
        rxn_list2 = ['r8', 'r9']

        plot_ensemble(self.flux_interp_quantiles, quant_type='flux', selected_data=[rxn_list1, rxn_list2],
                      x_scale='linear', y_scale='log', x_lim=None, y_lim=None, fig_size=None, save_plot=True,
                      output_file=output_file)

        self.assertTrue(os.path.isfile(output_file))
        os.remove(output_file)

    def test_plot_model_conc(self):
        output_file = os.path.join(self.this_dir, 'test_files', 'plot_model_conc_viz.png')

        plot_model(self.conc, self.conc_interp, model_i=1, quant_type='conc', selected_data=self.met_names,
                   x_scale='linear', y_scale='linear', x_lim=None, y_lim=None, fig_size=None,
                   save_plot=True, output_file=output_file)
        self.assertTrue(os.path.isfile(output_file))
        os.remove(output_file)

    def test_plot_model_flux(self):
        output_file = os.path.join(self.this_dir, 'test_files', 'plot_model_flux_viz.png')

        plot_model(self.flux, self.flux_interp, model_i=1, quant_type='flux', selected_data=self.rxn_names,
                   x_scale='linear', y_scale='linear', x_lim=None, y_lim=None, fig_size=None,
                   save_plot=True, output_file=output_file)

        self.assertTrue(os.path.isfile(output_file))
        os.remove(output_file)
