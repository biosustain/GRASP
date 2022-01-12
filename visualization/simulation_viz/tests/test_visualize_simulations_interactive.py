import os
import unittest

import scipy.io

from simulation_viz.import_simulation_data import gather_conc_data, gather_flux_data,\
    get_met_rxn_names, aggregate_time_series, import_ref_conc, import_ref_flux
from simulation_viz.visualize_simulations_interactive import plot_ensemble_interactive, \
    plot_model_interactive


class TestVisualizeSimulationsInteractive(unittest.TestCase):

    def setUp(self):
        self.this_dir, this_filename = os.path.split(__file__)
        raw_data_dir = os.path.join(self.this_dir, 'test_files')
        model_name = 'toy_model'

        n_models = 5

        order_of_magnitude = [-9, -8, -7, -6, -5, -4, -3, -2, -1]
        self.time_points = [coefficient * 10 ** exponent for exponent in order_of_magnitude for coefficient in
                            range(1, 10)]
        self.time_points.extend([1])
        self.time_points.insert(0, 0)

        self.met_names, self.rxn_names = get_met_rxn_names(raw_data_dir, model_name)

        file_in = os.path.join(raw_data_dir, f'{model_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        ref_conc_dic = import_ref_conc(mat, n_models)
        ref_flux_dic = import_ref_flux(mat)

        simulation_name = f'{model_name}'
        file_in = os.path.join(raw_data_dir, f'simulation_{simulation_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        self.conc_abs = gather_conc_data(mat, self.met_names, n_models, 'conc_abs', ref_conc_dic=ref_conc_dic)
        self.conc_rel = gather_conc_data(mat, self.met_names, n_models, 'conc_rel', ref_conc_dic=ref_conc_dic)

        self.flux_abs = gather_flux_data(mat, self.rxn_names, n_models, 'flux_abs', ref_flux_dic=ref_flux_dic)
        self.flux_rel = gather_flux_data(mat, self.rxn_names, n_models, 'flux_rel', ref_flux_dic=ref_flux_dic)

        quant_type = 'conc_rel'
        self.conc_interp_quantiles = aggregate_time_series(self.conc_rel, quant_type, self.met_names)

        quant_type = 'flux_abs'
        self.flux_interp_quantiles = aggregate_time_series(self.flux_abs, quant_type, self.rxn_names)

    def test_plot_ensemble_interactive_conc_rel(self):
        plot_ensemble_interactive(self.conc_interp_quantiles, quant_type='conc_rel', selected_data=self.met_names,
                                  x_scale='linear', y_scale='linear', x_lim=None, y_lim=None)

    def test_plot_ensemble_interactive_conc_abs(self):
        plot_ensemble_interactive(self.conc_interp_quantiles, quant_type='conc_abs', selected_data=self.met_names,
                                  x_scale='linear', y_scale='linear', x_lim=None, y_lim=None)

    def test_plot_ensemble_interactive_flux_rel(self):
        plot_ensemble_interactive(self.flux_interp_quantiles, quant_type='flux_rel', selected_data=self.rxn_names,
                                  x_scale='linear', y_scale='linear', x_lim=None, y_lim=None)

    def test_plot_ensemble_interactive_flux_abs(self):
        plot_ensemble_interactive(self.flux_interp_quantiles, quant_type='flux_abs', selected_data=self.rxn_names,
                                  x_scale='linear', y_scale='linear', x_lim=None, y_lim=None)

    def test_plot_model_interactive_conc_rel(self):
        plot_model_interactive(self.conc_rel, model_i=2, quant_type='conc_rel',
                               selected_data=self.met_names, x_scale='linear',
                               y_scale='linear', x_lim=None, y_lim=None)

    def test_plot_model_interactive_conc_abs(self):
        plot_model_interactive(self.conc_abs, model_i=2, quant_type='conc_abs',
                               selected_data=self.met_names, x_scale='linear',
                               y_scale='linear', x_lim=None, y_lim=None)

    def test_plot_model_interactive_flux_rel(self):
        plot_model_interactive(self.flux_rel, model_i=2, quant_type='flux_rel',
                               selected_data=self.rxn_names, x_scale='linear',
                               y_scale='linear', x_lim=None, y_lim=None)

    def test_plot_model_interactive_flux_abs(self):
        plot_model_interactive(self.flux_abs, model_i=2, quant_type='flux_abs',
                               selected_data=self.rxn_names, x_scale='linear',
                               y_scale='linear', x_lim=None, y_lim=None)

