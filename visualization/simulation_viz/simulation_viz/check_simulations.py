import numpy as np
import pandas as pd

def get_converging_models_option1(conc_df_interp: pd.DataFrame, n_models: int) -> list:
    non_converging_models = []

    for model_i in range(n_models):
        last_conc_values = \
        conc_df_interp[(conc_df_interp['model'] == model_i) & (conc_df_interp['time_point'].between(0.75, 1.05))][
            'concentration'].values

        if len(last_conc_values) == 0 or np.any(np.abs(last_conc_values) > 1.05) or np.any(np.abs(last_conc_values) < 0.95):
            non_converging_models.append(model_i)

    return non_converging_models


def get_converging_models_option2(conc_df_interp: pd.DataFrame, n_models: int, met_names: list,
                                  rel_tol: float = 5 * 10**-3) -> list:

    non_converging_models = []

    for model_i in range(n_models):

        for met in met_names:
            last_conc_values = conc_df_interp[(conc_df_interp['model'] == model_i) &
                                              (conc_df_interp['time_point'].between(0.75, 1.05)) &
                                              (conc_df_interp['met'] == met)]['concentration'].values
            assert len(last_conc_values) == 3
            diff1 = np.abs(last_conc_values[0] - last_conc_values[1])
            diff2 = np.abs(last_conc_values[1] - last_conc_values[2])
            abs_tol = last_conc_values[0]*rel_tol

            if diff1 > abs_tol or diff2 > abs_tol:
                non_converging_models.append(model_i)

    return non_converging_models


def remove_models(conc_df: pd.DataFrame, flux_df: pd.DataFrame, model_list: list) -> tuple:
    filtered_conc_df = conc_df.drop(conc_df[conc_df['model'].isin(model_list)].index)
    filtered_flux_df = flux_df.drop(flux_df[flux_df['model'].isin(model_list)].index)

    return filtered_conc_df, filtered_flux_df


def check_for_negative_concentrations(data_df: pd.DataFrame, scaled: bool, threshold: float = -10**-8):
    """
    Checks if all concentrations are higher than the given threshold. The idea is to make sure that there are no
    negative concentrations.

    Args:
        data_df: a pandas dataframe with a concentrations or concentrations_unscaled column.
        scaled: whether or not one wants to focus on scaled concentrations.
        threshold: value to use to check if the concentrations are higher than that.

    Returns:
        None
    """

    all_pos = np.all(data_df['concentration' if scaled else 'concentration_unscaled'].values > threshold)

    if all_pos:
        print(f'All concentrations are above the treshold {threshold} :)')
    else:
        print(f'There are some concentrations below {threshold} :(')

