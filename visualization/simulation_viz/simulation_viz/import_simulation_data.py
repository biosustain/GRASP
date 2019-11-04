import os
import time

import numpy as np
import pandas as pd
from scipy import interpolate


def import_ref_conc(mat: dict, n_models: int, met_names: list) -> dict:
    """
    Converts a .mat file from matlab into a dictionary with reference concentrations for each model.
    Metabolites are the keys, lists of concentrations are the values, where each entry is the concentration for a model.

    Args:
        mat: ensemble .mat file from matlab in form of dictionary.
        n_models: number of models to consider.
        met_names: name of the metabolites.

    Returns:
        A dictionary with the reference metabolite concentrations for each model.
    """

    n_mets = len(met_names)
    ref_conc = []
    for model_i in range(n_models):
        model_i_data = mat['ensemble']['populations'][0][0]['models'][0][0]['metConcRef'][0][model_i].transpose()[0]
        ref_conc.append(model_i_data)

    ref_conc_df = pd.DataFrame(data=ref_conc, index=range(n_models), columns=range(n_mets))
    ref_conc_df.columns = met_names
    ref_conc_dic = ref_conc_df.to_dict(orient='list')

    return ref_conc_dic


def gather_sim_data(mat: dict, met_names: list, rxn_names: list, n_models: int, time_points_spline: list,
                    save_concs: bool = False, save_fluxes: bool = False, ref_conc_dic: dict = None) -> tuple:
    """
    Imports the data from the matlab simulation results. Since the time points for which there are concentrations
    may not be the same for all models, it interpolates concentrations and fluxes, and stores the values at certain
    time points (defined in `time_points_spline`). Also, for long simulations, this saves a lot of space.

    Setting save_concs and/or save_fluxes to true may lead to a very large dataframe that occupies a lot of memory, so
    for larger models or just a large number of models, it is advised to set these variables to false.

    If a dictionary with reference concentrations is supplied the resulting dataframe will also contain a column with
    unscaled concentrations (named concentrations_unscaled), i.e. absolute concentrations.


    Args:
        mat: matlab structures with simulation results.
        met_names: names of the simulated metabolites in the model.
        rxn_names: names of the simulated reactions in the model.
        n_models: number of models to be simulated.
        time_points_spline: time points for which concentration and flux values will be returned.
        save_concs: whether or not to save all concentration values.
        save_fluxes: whether or not to save all flux values.
        ref_conc_dic: a dictionary with the metabolites reference concentrations in each model.

    Returns:
        A pandas dataframe with all concentration values (if save_concs==True), a pandas dataframe with concentration
        values at the time points in `time_points_spline`, and a pandas dataframe with flux values at the time points
        in `time_points_spline`.
    """

    data_conc = {'model': [], 'met': [], 'time_point': [], 'conc': [], 'conc_unscaled': []}
    data_conc_interp = {'model': [], 'met': [], 'time_point': [], 'conc': [], 'conc_unscaled': []}
    data_flux = {'model': [], 'rxn': [], 'time_point': [], 'flux': []}
    data_flux_interp = {'model': [], 'rxn': [], 'time_point': [], 'flux': []}

    n_missing_models = 0
    start_total = time.time()
    for model_i in range(n_models):

        try:
            time_points = mat['simulationRes'][0][model_i]['t'][0][0].transpose()[0]
            met_concs = mat['simulationRes'][0][model_i]['conc'][0][0].transpose()
            rxn_fluxes = mat['simulationRes'][0][model_i]['flux'][0][0].transpose()

        except IndexError:
            n_missing_models += 1
            print(f'{model_i} is missing.')
            continue

        # This step can be very time and space consuming
        if save_concs:
            for met_i, met in enumerate(met_names):
                data_conc['model'].extend(np.repeat(model_i, len(time_points)))
                data_conc['time_point'].extend(time_points)
                data_conc['met'].extend(np.repeat(met, len(time_points)))
                data_conc['conc'].extend(met_concs[met_i])
                if ref_conc_dic is not None:
                    data_conc['conc_unscaled'].extend(met_concs[met_i] * ref_conc_dic[met][model_i])
                else:
                    data_conc_interp['conc_unscaled'].extend(np.repeat(np.nan, len(met_concs[met_i])))

        for met_i, met in enumerate(met_names):
            data_conc_interp['model'].extend(np.repeat(model_i, len(time_points_spline)))
            data_conc_interp['time_point'].extend(time_points_spline)
            data_conc_interp['met'].extend(np.repeat(met_names[met_i], len(time_points_spline)))

            conc_interp = interpolate.CubicSpline(time_points, met_concs[met_i])
            data_conc_interp['conc'].extend(conc_interp(time_points_spline))

            if ref_conc_dic is not None:
                data_conc_interp['conc_unscaled'].extend(
                    conc_interp(time_points_spline) * ref_conc_dic[met][model_i])
            else:
                data_conc_interp['conc_unscaled'].extend(np.repeat(np.nan, len(time_points_spline)))

        if save_fluxes:
            for rxn_i in range(len(rxn_fluxes)):
                data_flux['model'].extend(np.repeat(model_i, len(time_points)))
                data_flux['time_point'].extend(time_points)
                data_flux['rxn'].extend(np.repeat(rxn_names[rxn_i], len(time_points)))
                data_flux['flux'].extend(rxn_fluxes[rxn_i])

        for rxn_i in range(len(rxn_fluxes)):
            data_flux_interp['model'].extend(np.repeat(model_i, len(time_points_spline)))
            data_flux_interp['time_point'].extend(time_points_spline)
            data_flux_interp['rxn'].extend(np.repeat(rxn_names[rxn_i], len(time_points_spline)))

            flux_interp = interpolate.CubicSpline(time_points, rxn_fluxes[rxn_i])
            data_flux_interp['flux'].extend(flux_interp(time_points_spline))

    data_df_conc = pd.DataFrame.from_dict(data_conc) if save_concs else None
    data_df_conc_interp = pd.DataFrame.from_dict(data_conc_interp)
    data_df_flux = pd.DataFrame.from_dict(data_flux) if save_fluxes else None
    data_df_flux_interp = pd.DataFrame.from_dict(data_flux_interp)

    print(f'total time: {time.time() - start_total}')
    print(f'There were a total of {n_missing_models} missing models out of {n_models}.')

    return data_df_conc, data_df_conc_interp, data_df_flux, data_df_flux_interp


def get_met_rxn_names(raw_data_dir: str, model_name: str) -> tuple:
    """
    Gets the names of metabolites and reactions in the model.

    Args:
        raw_data_dir: path to folder with the raw data.
        model_name: named of the model.

    Returns:
        A list with the metabolite names and another with the reaction names.
    """

    file_met_names = os.path.join(raw_data_dir, f'{model_name}_metsActive.dat')
    met_names = pd.read_csv(file_met_names, sep='\n').values
    met_names = list(met_names.transpose()[0])
    met_names = [met_name.replace('m_m_', '') for met_name in met_names]

    # get reaction names
    file_rxn_names = os.path.join(raw_data_dir, f'{model_name}_rxnsActive.dat')
    rxn_names = pd.read_csv(file_rxn_names, sep='\n').values
    rxn_names = list(rxn_names.transpose()[0])
    rxn_names = [rxn_name.replace('r_', '') for rxn_name in rxn_names]

    return met_names, rxn_names


def get_time_series_quantiles(data_df: pd.DataFrame, time_points: list, quant_type: str, name_list: list) -> pd.DataFrame:
    """
    Takes in a dataframe with either fluxes or concentrations and calculates the median and respective first and third
    quartiles for each flux or metabolite at a given time point over the whole model ensemble.

    Args:
        data_df: dataframe with concentrations or fluxes for each time point and model.
        time_points: time points for which there are concentrations or fluxes.
        quant_type: whether the data frame contains concentrations ('conc'), unscaled concentrations ('conc_unscaled'),
                   or fluxes ('flux').
        name_list: reaction or metabolite names, depending on what that the dataframe contains.

    Returns:
        A pandas dataframe with the concentrations/fluxes median and 1st and 3rd quartiles for each metabolite/reaction.
    """

    data_type = 'met' if quant_type.find('conc') != -1 else 'rxn'
    n_time_points = len(time_points)

    quantiles_dic = {'time_point': [], data_type: [], 'q025': np.array([]), 'median': np.array([]),
                     'q075': np.array([])}

    q025 = data_df.groupby([data_type, 'time_point']).quantile(q=0.25)
    q050 = data_df.groupby([data_type, 'time_point']).median()
    q075 = data_df.groupby([data_type, 'time_point']).quantile(q=0.75)

    for item in name_list:
        quantiles_dic[data_type].extend(np.repeat(item, n_time_points))
        quantiles_dic['time_point'].extend(time_points)

        quantiles_dic['q025'] = np.append(quantiles_dic['q025'], q025.loc[item, :][quant_type])
        quantiles_dic['median'] = np.append(quantiles_dic['median'], q050.loc[item, :][quant_type])
        quantiles_dic['q075'] = np.append(quantiles_dic['q075'], q075.loc[item, :][quant_type])

    quantiles_df = pd.DataFrame.from_dict(quantiles_dic)

    return quantiles_df
