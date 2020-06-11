import io
import os
import unittest.mock

import pandas as pd

from simulation_viz.simulation_viz.thermodynamics_checks import check_thermodynamic_feasibility, calculate_dG


class TestThermodynamicsChecks(unittest.TestCase):

    def setUp(self):
        this_dir, this_filename = os.path.split(__file__)
        self.test_folder = os.path.join(this_dir, 'test_files', 'test_thermodynamics_checks')

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_check_thermodynamic_feasibility_putida(self, mock_stdout):
        true_res = ('\nChecking if fluxes and Gibbs energies are compatible.\n\n' +
                    'The flux and ∆G range seem to be incompatible for reaction R_G6PDH2\n' +
                    'The flux and ∆G range seem to be incompatible for reaction R_RPE\n' +
                    'The flux and ∆G range seem to be incompatible for reaction R_TALA\n' +
                    'The flux and ∆G range seem to be incompatible for reaction R_PGI\n' +
                    'The flux and ∆G range seem to be incompatible for reaction R_FBA\n' +
                    'The flux and ∆G range seem to be incompatible for reaction R_TPI\n' +
                    'The flux and ∆G range seem to be incompatible for reaction R_GAPD\n' +
                    'The flux and ∆G range seem to be incompatible for reaction R_ENO\n' +
                    'The flux and ∆G range seem to be incompatible for reaction R_EX_pyr\n' +
                    'The flux and ∆G range seem to be incompatible for reaction R_EX_pep\n')

        data_dict = pd.read_excel(os.path.join(self.test_folder, 'model_v1_base.xlsx'), sheet_name=None, index_col=0)
        flag, flux_df, dG_df = check_thermodynamic_feasibility(data_dict)
        self.assertEqual(True, flag)
        self.assertEqual(true_res, mock_stdout.getvalue())

    @unittest.mock.patch('sys.stdout', new_callable=io.StringIO)
    def test_check_thermodynamic_feasibility_HMP(self, mock_stdout):
        true_res = ('\nChecking if fluxes and Gibbs energies are compatible.\n\n' +
                    'Everything seems to be OK.\n')

        data_dict = pd.read_excel(os.path.join(self.test_folder, 'HMP2360_r0_t0.xlsx'), sheet_name=None, index_col=0)
        flag, dG_df, flux_df = check_thermodynamic_feasibility(data_dict)

        self.assertEqual(False, flag)
        self.assertEqual(true_res, mock_stdout.getvalue())

    def test_calculate_dG(self):
        true_res_ma = pd.read_pickle(os.path.join(self.test_folder, 'true_res_ma.pkl'))
        true_res_dG = pd.read_pickle(os.path.join(self.test_folder, 'true_res_dG.pkl'))
        true_res_dG_Q = pd.read_pickle(os.path.join(self.test_folder, 'true_res_dG_Q.pkl'))

        temperature = 298  # in K
        gas_constant = 8.314 * 10 ** -3  # in kJ K^-1 mol^-1

        data_dict = pd.read_excel(os.path.join(self.test_folder, 'model_v2_manual.xlsx'), sheet_name=None, index_col=0)
        ma_df, dG_Q_df, dG_df = calculate_dG(data_dict, gas_constant, temperature)

        self.assertTrue(ma_df.equals(true_res_ma))
        self.assertTrue(dG_Q_df.equals(true_res_dG_Q))
        self.assertTrue(dG_df.equals(true_res_dG))
