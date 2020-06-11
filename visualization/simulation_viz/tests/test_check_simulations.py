import io
import os
import unittest.mock

import scipy.io

from simulation_viz.simulation_viz.check_simulations import check_for_negative_concentrations, \
    get_non_converging_models, remove_models
from simulation_viz.simulation_viz.import_simulation_data import gather_sim_data, get_met_rxn_names, import_ref_conc


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

        n_mets = len(mat['ensemble']['mets'][0][0])
        all_met_names = [mat['ensemble']['mets'][0][0][met_i][0][0].replace('m_m_', '') for met_i in range(n_mets)]

        ref_conc_dic = import_ref_conc(mat, n_models)

        simulation_name = f'{model_name}'
        file_in = os.path.join(raw_data_dir, f'simulation_{simulation_name}.mat')
        mat = scipy.io.loadmat(file_in, squeeze_me=False)

        self.conc, self.conc_interp, self.flux, self.flux_interp = gather_sim_data(mat, self.met_names, self.rxn_names,
                                                                                   n_models,
                                                                                   self.time_points, save_concs=True,
                                                                                   save_fluxes=True,
                                                                                   ref_conc_dic=ref_conc_dic)

    def test_get_non_converging_models(self):
        n_models = 5
        non_converging_models = get_non_converging_models(self.conc_interp, n_models, self.met_names,
                                                          rel_tol=5 * 10 ** -3)

        self.assertEqual(non_converging_models, [])

    def test_remove_models(self):
        conc_interp, flux_interp = remove_models(self.conc_interp, self.flux, [2, 3])

        models_left_conc = list(conc_interp['model'].unique())
        models_left_flux = list(flux_interp['model'].unique())

        self.assertListEqual([0, 1, 4], models_left_conc)
        self.assertListEqual([0, 1, 4], models_left_flux)

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_check_for_negative_concentrations(self, mock_stdout):
        true_res = ('All concentrations are above the treshold -1e-08 :)\n')

        check_for_negative_concentrations(self.conc_interp, scaled=False, threshold=-10 ** -8)

        self.assertEqual(true_res, mock_stdout.getvalue())
