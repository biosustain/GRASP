import os
import unittest

import scipy.io

from simulation_viz.import_simulation_data import gather_conc_data, gather_flux_data,\
    get_met_rxn_names, aggregate_time_series, import_ref_conc, import_ref_flux
from simulation_viz.visualize_simulations import plot_ensemble, plot_model


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

        ref_conc_dic = import_ref_conc(mat, n_models)
        ref_flux_dic = import_ref_flux(mat)

        simulation_name = f'{model_name}'
        file_in = os.path.join(raw_data_dir, f'simulation_{simulation_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        self.conc_abs = gather_conc_data(mat, self.met_names, n_models, 'conc_abs', ref_conc_dic=ref_conc_dic)
        self.conc_rel = gather_conc_data(mat, self.met_names, n_models, 'conc_rel', ref_conc_dic=ref_conc_dic)
        self.flux_abs = gather_flux_data(mat, self.rxn_names, n_models, 'flux_abs', ref_flux_dic=ref_flux_dic)
        self.flux_rel = gather_flux_data(mat, self.rxn_names, n_models, 'flux_rel', ref_flux_dic=ref_flux_dic)

        data_type = 'conc_rel'
        self.conc_interp_quantiles = aggregate_time_series(self.conc_rel, data_type, self.met_names)

        data_type = 'flux_abs'
        self.flux_interp_quantiles = aggregate_time_series(self.flux_abs, data_type, self.rxn_names)

    def test_plot_ensemble_conc_1row_rel(self):
        output_file = os.path.join(self.this_dir, 'test_files', 'plot_ensemble_conc_viz_rel.png')
        plot_ensemble(self.conc_interp_quantiles, quant_type='conc_rel', selected_data=[self.met_names],
                      x_scale='linear', y_scale='log', x_lim=None, y_lim=None, fig_size=None, save_plot=True,
                      output_file=output_file)

        self.assertTrue(os.path.isfile(output_file))
        os.remove(output_file)

    def test_plot_ensemble_conc_1row_abs(self):
        output_file = os.path.join(self.this_dir, 'test_files', 'plot_ensemble_conc_viz_abs.png')
        plot_ensemble(self.conc_interp_quantiles, quant_type='conc_abs', selected_data=[self.met_names],
                      x_scale='linear', y_scale='log', x_lim=None, y_lim=None, fig_size=None, save_plot=True,
                      output_file=output_file)

        self.assertTrue(os.path.isfile(output_file))
        os.remove(output_file)

    def test_plot_ensemble_conc_2rows_rel(self):
        output_file = os.path.join(self.this_dir, 'test_files', 'plot_ensemble_conc_viz2_rel.png')
        met_list1 = ['m_m5', 'm_m6']
        met_list2 = ['m_m10', 'm_m11']

        plot_ensemble(self.conc_interp_quantiles, quant_type='conc_rel', selected_data=[met_list1, met_list2],
                      x_scale='linear', y_scale='log', x_lim=None, y_lim=None, fig_size=None, save_plot=True,
                      output_file=output_file)

        self.assertTrue(os.path.isfile(output_file))
        os.remove(output_file)

    def test_plot_ensemble_conc_2rows_abs(self):
        output_file = os.path.join(self.this_dir, 'test_files', 'plot_ensemble_conc_viz2_abs.png')
        met_list1 = ['m_m5', 'm_m6']
        met_list2 = ['m_m10', 'm_m11']

        plot_ensemble(self.conc_interp_quantiles, quant_type='conc_abs', selected_data=[met_list1, met_list2],
                      x_scale='linear', y_scale='log', x_lim=None, y_lim=None, fig_size=None, save_plot=True,
                      output_file=output_file)

        self.assertTrue(os.path.isfile(output_file))
        os.remove(output_file)

    def test_plot_ensemble_flux_1row_rel(self):
        output_file = os.path.join(self.this_dir, 'test_files', 'plot_ensemble_flux_viz_rel.png')
        plot_ensemble(self.flux_interp_quantiles, quant_type='flux_rel', selected_data=[self.rxn_names],
                      x_scale='linear', y_scale='log', x_lim=None, y_lim=None, fig_size=None, save_plot=True,
                      output_file=output_file)

        self.assertTrue(os.path.isfile(output_file))
        os.remove(output_file)

    def test_plot_ensemble_flux_1row_abs(self):
        output_file = os.path.join(self.this_dir, 'test_files', 'plot_ensemble_flux_viz_abs.png')
        plot_ensemble(self.flux_interp_quantiles, quant_type='flux_abs', selected_data=[self.rxn_names],
                      x_scale='linear', y_scale='log', x_lim=None, y_lim=None, fig_size=None, save_plot=True,
                      output_file=output_file)

        self.assertTrue(os.path.isfile(output_file))
        os.remove(output_file)

    def test_plot_ensemble_flux_2rows_rel(self):
        output_file = os.path.join(self.this_dir, 'test_files', 'plot_ensemble_flux_viz2_rel.png')
        rxn_list1 = ['r1', 'r2']
        rxn_list2 = ['r8', 'r9']

        plot_ensemble(self.flux_interp_quantiles, quant_type='flux_rel', selected_data=[rxn_list1, rxn_list2],
                      x_scale='linear', y_scale='log', x_lim=None, y_lim=None, fig_size=None, save_plot=True,
                      output_file=output_file)

        self.assertTrue(os.path.isfile(output_file))
        os.remove(output_file)

    def test_plot_ensemble_flux_2rows_abs(self):
        output_file = os.path.join(self.this_dir, 'test_files', 'plot_ensemble_flux_viz2_abs.png')
        rxn_list1 = ['r1', 'r2']
        rxn_list2 = ['r8', 'r9']

        plot_ensemble(self.flux_interp_quantiles, quant_type='flux_abs', selected_data=[rxn_list1, rxn_list2],
                      x_scale='linear', y_scale='log', x_lim=None, y_lim=None, fig_size=None, save_plot=True,
                      output_file=output_file)

        self.assertTrue(os.path.isfile(output_file))
        os.remove(output_file)

    def test_plot_model_conc_rel(self):
        output_file = os.path.join(self.this_dir, 'test_files', 'plot_model_conc_viz_rel.png')

        plot_model(self.conc_rel, model_i=1, quant_type='conc_rel', selected_data=self.met_names,
                   x_scale='linear', y_scale='linear', x_lim=None, y_lim=None, fig_size=None,
                   save_plot=True, output_file=output_file)
        self.assertTrue(os.path.isfile(output_file))
        os.remove(output_file)

    def test_plot_model_conc_abs(self):
        output_file = os.path.join(self.this_dir, 'test_files', 'plot_model_conc_viz_abs.png')

        plot_model(self.conc_abs, model_i=1, quant_type='conc_abs', selected_data=self.met_names,
                   x_scale='linear', y_scale='linear', x_lim=None, y_lim=None, fig_size=None,
                   save_plot=True, output_file=output_file)
        self.assertTrue(os.path.isfile(output_file))
        os.remove(output_file)

    def test_plot_model_flux_rel(self):
        output_file = os.path.join(self.this_dir, 'test_files', 'plot_model_flux_viz_rel.png')

        plot_model(self.flux_rel, model_i=1, quant_type='flux_rel', selected_data=self.rxn_names,
                   x_scale='linear', y_scale='linear', x_lim=None, y_lim=None, fig_size=None,
                   save_plot=True, output_file=output_file)

        self.assertTrue(os.path.isfile(output_file))
        os.remove(output_file)

    def test_plot_model_flux_abs(self):
        output_file = os.path.join(self.this_dir, 'test_files', 'plot_model_flux_viz_abs.png')

        plot_model(self.flux_abs, model_i=1, quant_type='flux_abs', selected_data=self.rxn_names,
                   x_scale='linear', y_scale='linear', x_lim=None, y_lim=None, fig_size=None,
                   save_plot=True, output_file=output_file)

        self.assertTrue(os.path.isfile(output_file))
        os.remove(output_file)
