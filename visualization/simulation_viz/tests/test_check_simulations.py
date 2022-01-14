import io
import os
import unittest.mock

import scipy.io

from simulation_viz.check_simulations import check_for_negative_concentrations, \
    get_non_converging_models, remove_models
from simulation_viz.import_simulation_data import gather_conc_data, gather_flux_data, get_met_rxn_names,\
    import_ref_conc, import_ref_flux


class TestImportSimulationData(unittest.TestCase):

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

    def test_get_non_converging_models_abs(self):
        n_models = 5
        quant_type = 'conc_abs'
        non_converging_models = get_non_converging_models(self.conc_abs, n_models, self.met_names, quant_type,
                                                          rel_tol=5 * 10 ** -3)

        self.assertEqual(non_converging_models, [2])

    def test_get_non_converging_models_rel(self):
        n_models = 5
        quant_type = 'conc_rel'
        non_converging_models = get_non_converging_models(self.conc_rel, n_models, self.met_names, quant_type,
                                                          rel_tol=5 * 10 ** -3)

        self.assertEqual(non_converging_models, [2])

    def test_remove_models(self):
        conc_interp, flux_interp = remove_models(self.conc_abs, self.flux_abs, [2, 3])

        models_left_conc = list(conc_interp['model'].unique())
        models_left_flux = list(flux_interp['model'].unique())

        self.assertListEqual([0, 1, 4], models_left_conc)
        self.assertListEqual([0, 1, 4], models_left_flux)

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_check_for_negative_concentrations_rel(self, mock_stdout):
        true_res = ('All concentrations are above the treshold -1e-08 :)\n')

        quant_type = 'conc_rel'
        check_for_negative_concentrations(self.conc_rel, quant_type, threshold=-10 ** -8)

        self.assertEqual(true_res, mock_stdout.getvalue())

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_check_for_negative_concentrations_abs(self, mock_stdout):
        true_res = ('All concentrations are above the treshold -1e-08 :)\n')

        quant_type = 'conc_abs'
        check_for_negative_concentrations(self.conc_abs, quant_type, threshold=-10 ** -8)

        self.assertEqual(true_res, mock_stdout.getvalue())
