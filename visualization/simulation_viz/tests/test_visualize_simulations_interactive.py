import os
import unittest
import scipy.io

from simulation_viz.visualize_simulations_interactive import plot_ensemble_interactive, plot_model_interactive
from simulation_viz.import_simulation_data import gather_sim_data, get_met_rxn_names, get_time_series_quantiles, import_ref_conc


class TestVisualizeSimulationsInteractive(unittest.TestCase):

    def setUp(self):
        self.this_dir, this_filename = os.path.split(__file__)
        raw_data_dir = os.path.join(self.this_dir, 'test_files')
        model_name = 'toy_model'

        n_models = 5

        order_of_magnitude = [-9, -8, -7, -6, -5, -4, -3, -2, -1]
        self.time_points = [coefficient * 10 ** exponent for exponent in order_of_magnitude for coefficient in range(1, 10)]
        self.time_points.extend([1])
        self.time_points.insert(0, 0)

        self.met_names, self.rxn_names = get_met_rxn_names(raw_data_dir, model_name)

        file_in = os.path.join(raw_data_dir, f'{model_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        n_mets = len(mat['ensemble']['mets'][0][0])
        all_met_names = [mat['ensemble']['mets'][0][0][met_i][0][0].replace('m_m_', '') for met_i in range(n_mets)]

        ref_conc_dic = import_ref_conc(mat, n_models, all_met_names)

        simulation_name = f'{model_name}'
        file_in = os.path.join(raw_data_dir, f'simulation_{simulation_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        self.conc, self.conc_interp, self.flux, self.flux_interp = gather_sim_data(mat, self.met_names, self.rxn_names, n_models,
                                                                                   self.time_points, save_concs=True, save_fluxes=True,
                                                                                   ref_conc_dic=ref_conc_dic)

        quant_type = 'conc'
        self.conc_interp_quantiles = get_time_series_quantiles(self.conc_interp, self.time_points, quant_type, self.met_names)

        quant_type = 'flux'
        self.flux_interp_quantiles = get_time_series_quantiles(self.flux_interp, self.time_points, quant_type, self.rxn_names)

    def test_plot_ensemble_interactive_conc(self):
        plot_ensemble_interactive(self.conc_interp_quantiles, quant_type='conc', selected_data=self.met_names,
                                  x_scale='linear', y_scale='linear', x_lim=None, y_lim=None)

    def test_plot_ensemble_interactive_flux(self):
        plot_ensemble_interactive(self.flux_interp_quantiles, quant_type='flux', selected_data=self.rxn_names,
                                  x_scale='linear', y_scale='linear', x_lim=None, y_lim=None)

    def test_plot_model_interactive_conc(self):
        plot_model_interactive(self.conc, self.conc_interp, model_i=2, quant_type='conc',
                               selected_data=self.met_names, x_scale='linear',
                               y_scale='linear', x_lim=None, y_lim=None)

    def test_plot_model_interactive_flux(self):
        plot_model_interactive(self.flux, self.flux_interp, model_i=2, quant_type='flux',
                               selected_data=self.rxn_names, x_scale='linear',
                               y_scale='linear', x_lim=None, y_lim=None)

    def test_plot_model_interactive_conc_interp_only(self):
        plot_model_interactive(None, self.conc_interp, model_i=2, quant_type='conc',
                               selected_data=self.met_names, x_scale='linear',
                               y_scale='linear', x_lim=None, y_lim=None)

    def test_plot_model_interactive_flux_interp_only(self):
        plot_model_interactive(None, self.flux_interp, model_i=2, quant_type='flux',
                               selected_data=self.rxn_names, x_scale='linear',
                               y_scale='linear', x_lim=None, y_lim=None)

