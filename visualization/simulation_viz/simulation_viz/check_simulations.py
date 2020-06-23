import numpy as np
import pandas as pd


def get_non_converging_models(conc_df_interp: pd.DataFrame, n_models: int, met_names: list, quant_type: str,
                              rel_tol: float = 5 * 10**-3) -> list:
    """
    Goes through every model and every metabolite time course, gets the last three time points in the simulation and
    checks if the difference between each pair of concentration values is greater than last_concentration_value*rel_tol,
    where the relative tolerance is defined by the user. If this is false, the model is considered to not converge.

    Args:
        conc_df_interp: dataframe with the interpolated concentrations.
        n_models: number of models in the dataframe.
        met_names: name of all metabolites.
        quant_type: which quantity to extract: "flux_rel", "flux_abs", "conc_rel", "conc_abs".
        rel_tol: relative tolerance to use when checking if the model is converging.

    Returns:
        A list with all model numbers that are considered to not converge.
    """

    if quant_type not in {'flux_abs', 'flux_rel', 'conc_abs', 'conc_rel'}:
        print('quant_type should be set to one of the following values; "flux_abs", "flux_rel", "conc_abs", '
              'or "conc_rel"')
        return []

    non_converging_models = []

    for model_i in range(n_models):

        for met in met_names:
            last_conc_values = conc_df_interp[(conc_df_interp['model'] == model_i) & (conc_df_interp['met'] == met)].iloc[-3:, :][quant_type].values

            assert len(last_conc_values) == 3
            diff1 = np.abs(last_conc_values[0] - last_conc_values[1])
            diff2 = np.abs(last_conc_values[1] - last_conc_values[2])
            abs_tol = last_conc_values[2]*rel_tol

            if diff1 > abs_tol or diff2 > abs_tol:
                non_converging_models.append(model_i)

    return non_converging_models


def remove_models(conc_df: pd.DataFrame, flux_df: pd.DataFrame, model_list: list) -> tuple:
    """
    Removes all models in model_list from concentration and flux dataframes.

    Args:
        conc_df: dataframe with metabolite concentrations.
        flux_df: dataframe with reactions flux values.
        model_list: list of models to be removed from dataframes.

    Returns:
        The concentration and flux dataframes with models in model_list removed.
    """

    filtered_conc_df = conc_df.drop(conc_df[conc_df['model'].isin(model_list)].index)
    filtered_flux_df = flux_df.drop(flux_df[flux_df['model'].isin(model_list)].index)

    return filtered_conc_df, filtered_flux_df


def check_for_negative_concentrations(data_df: pd.DataFrame, quant_type: str, threshold: float = -10**-8):
    """
    Checks if all concentrations are higher than the given threshold. The idea is to make sure that there are no
    negative concentrations.

    Args:
        data_df: a pandas dataframe with a concentrations or concentrations_unscaled column.
        quant_type: which quantity to extract: "flux_rel", "flux_abs", "conc_rel", "conc_abs".
        threshold: value to use to check if the concentrations are higher than that.

    Returns:
        None
    """

    if quant_type not in {'flux_abs', 'flux_rel', 'conc_abs', 'conc_rel'}:
        print('quant_type should be set to one of the following values; "flux_abs", "flux_rel", "conc_abs", '
              'or "conc_rel"')
        return []

    all_pos = np.all(data_df[quant_type].values > threshold)

    if all_pos:
        print(f'All concentrations are above the treshold {threshold} :)')
    else:
        print(f'There are some concentrations below {threshold} :(')

