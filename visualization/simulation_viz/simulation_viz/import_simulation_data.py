import os
import time

import numpy as np
import pandas as pd
import scipy.io
from scipy import interpolate


def import_ref_conc(mat: dict, n_models: int) -> dict:
    """
    Converts a .mat file from matlab into a dictionary with reference concentrations for each model.
    Metabolites are the keys, lists of concentrations are the values, where each entry is the concentration for a model.

    Args:
        mat: ensemble .mat file from matlab in form of dictionary.
        n_models: number of models to consider.

    Returns:
        A dictionary with the reference metabolite concentrations for each model.
    """

    # get ALL metabolite names, also the constant ones
    n_mets = len(mat['ensemble']['mets'][0][0])
    met_names = [mat['ensemble']['mets'][0][0][met_i][0][0][2:] for met_i in range(n_mets)]

    ref_conc = []
    for model_i in range(n_models):
        model_i_data = mat['ensemble']['populations'][0][0]['models'][0][0]['metConcRef'][0][model_i].transpose()[0]
        ref_conc.append(model_i_data)

    ref_conc_df = pd.DataFrame(data=ref_conc, index=range(n_models), columns=range(n_mets))
    ref_conc_df.columns = met_names
    ref_conc_dic = ref_conc_df.to_dict(orient='list')

    return ref_conc_dic


def import_ref_flux(mat: dict) -> dict:
    """
    Get the reference fluxes into a dictionary from a .mat file from matlab.
    Reaction names are the keys and fluxes are the values.

    Args:
        mat: ensemble .mat file from matlab in form of dictionary.

    Returns:
        A dictionary with the reference reaction fluxes.
    """

    n_rxns = len(mat['ensemble']['rxns'][0][0])
    rxn_names = [mat['ensemble']['rxns'][0][0][rxn_i][0][0][2:] for rxn_i in range(n_rxns)]

    ref_flux = mat['ensemble']['fluxRef'][0][0]
    ref_flux_dic = dict((rxn, flux[0]) for rxn, flux in zip(rxn_names, ref_flux))

    return ref_flux_dic


def gather_sim_data(mat: dict, met_names: list, rxn_names: list, n_models: int, time_points_spline: list,
                    save_concs: bool = False, save_fluxes: bool = False, ref_conc_dic: dict = None,
                    ref_flux_dic: dict = None) -> tuple:
    """
    Imports the data from the matlab simulation results. Since the time points for which there are concentrations
    may not be the same for all models, it interpolates concentrations and fluxes, and stores the values at certain
    time points (defined in `time_points_spline`). Also, for long simulations, this saves a lot of space.

    Setting save_concs and/or save_fluxes to true may lead to a very large dataframe that occupies a lot of memory, so
    for larger models or just a large number of models, it is advised to set these variables to false.

    If a dictionary with reference concentrations is supplied the resulting dataframe will also contain a column with
    unscaled concentrations (named concentrations_unscaled), i.e. absolute concentrations.

    The fluxes coming from matlab are absolute while the concentrations are relative to the reference state.


    Args:
        mat: matlab structures with simulation results.
        met_names: names of the simulated metabolites in the model.
        rxn_names: names of the simulated reactions in the model.
        n_models: number of models to be simulated.
        time_points_spline: time points for which concentration and flux values will be returned.
        save_concs: whether or not to save all concentration values.
        save_fluxes: whether or not to save all flux values.
        ref_conc_dic: a dictionary with the metabolites reference concentrations in each model.
        ref_flux_dic: a dictionary with the reactions reference fluxes.

    Returns:
        A pandas dataframe with all concentration values (if save_concs==True), a pandas dataframe with concentration
        values at the time points in `time_points_spline`, and a pandas dataframe with flux values at the time points
        in `time_points_spline`.
    """

    data_conc = {'model': [], 'met': [], 'time_point': [], 'conc_rel': [], 'conc_abs': []}
    data_conc_interp = {'model': [], 'met': [], 'time_point': [], 'conc_rel': [], 'conc_abs': []}
    data_flux = {'model': [], 'rxn': [], 'time_point': [], 'flux_rel': [], 'flux_abs': []}
    data_flux_interp = {'model': [], 'rxn': [], 'time_point': [], 'flux_rel': [], 'flux_abs': []}

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
                data_conc['conc_rel'].extend(met_concs[met_i])

                if ref_conc_dic is not None:
                    data_conc['conc_abs'].extend(met_concs[met_i] * ref_conc_dic[met][model_i])
                else:
                    data_conc['conc_abs'].extend(np.repeat(np.nan, len(met_concs[met_i])))

        for met_i, met in enumerate(met_names):
            data_conc_interp['model'].extend(np.repeat(model_i, len(time_points_spline)))
            data_conc_interp['time_point'].extend(time_points_spline)
            data_conc_interp['met'].extend(np.repeat(met_names[met_i], len(time_points_spline)))

            conc_interp = interpolate.CubicSpline(time_points, met_concs[met_i])
            data_conc_interp['conc_rel'].extend(conc_interp(time_points_spline))

            if ref_conc_dic is not None:
                data_conc_interp['conc_abs'].extend(
                    conc_interp(time_points_spline) * ref_conc_dic[met][model_i])
            else:
                data_conc_interp['conc_abs'].extend(np.repeat(np.nan, len(time_points_spline)))

        if save_fluxes:
            for rxn_i, rxn in enumerate(rxn_names):
                data_flux['model'].extend(np.repeat(model_i, len(time_points)))
                data_flux['time_point'].extend(time_points)
                data_flux['rxn'].extend(np.repeat(rxn_names[rxn_i], len(time_points)))
                data_flux['flux_abs'].extend(rxn_fluxes[rxn_i])

                if ref_flux_dic is not None:
                    data_flux['flux_rel'].extend(rxn_fluxes[rxn_i] / ref_flux_dic[rxn])
                else:
                    data_flux['flux_rel'].extend(np.repeat(np.nan, len(rxn_fluxes[rxn_i])))

        for rxn_i, rxn in enumerate(rxn_names):
            data_flux_interp['model'].extend(np.repeat(model_i, len(time_points_spline)))
            data_flux_interp['time_point'].extend(time_points_spline)
            data_flux_interp['rxn'].extend(np.repeat(rxn_names[rxn_i], len(time_points_spline)))

            flux_interp = interpolate.CubicSpline(time_points, rxn_fluxes[rxn_i])
            data_flux_interp['flux_abs'].extend(flux_interp(time_points_spline))

            if ref_flux_dic is not None:
                data_flux_interp['flux_rel'].extend(flux_interp(time_points_spline) / ref_flux_dic[rxn])
            else:
                data_flux_interp['flux_rel'].extend(np.repeat(np.nan, len(time_points_spline)))

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
    met_names = [met_name[2:] for met_name in met_names]

    # get reaction names
    file_rxn_names = os.path.join(raw_data_dir, f'{model_name}_rxnsActive.dat')
    rxn_names = pd.read_csv(file_rxn_names, sep='\n').values
    rxn_names = list(rxn_names.transpose()[0])
    rxn_names = [rxn_name[2:] for rxn_name in rxn_names]

    return met_names, rxn_names


def get_time_series_quantiles(data_df: pd.DataFrame, time_points: list, quant_type: str, name_list: list) -> pd.DataFrame:
    """
    Takes in a dataframe with either fluxes or concentrations and calculates the median and respective first and third
    quartiles for each flux or metabolite at a given time point over the whole model ensemble.

    Args:
        data_df: dataframe with concentrations or fluxes for each time point and model.
        time_points: time points for which there are concentrations or fluxes.
        quant_type: whether the data frame contains relative concentrations ('conc_rel'), absolute concentrations
                    ('conc_abs'), relative fluxes ('flux_rel'), or absolute_fluxes ('flux_abs'.
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


def load_simulation(raw_data_dir: str, model_name: str, simulation_name: str, n_models: int,
                    time_points: list, save_concs: bool=False, save_fluxes: bool=False) -> tuple:
    """
    Takes in the path to the folder with the model ensemble and simulation data, loads it and returns dataframes for
    concentrations and fluxes. It always returns dataframes for the concentrations and fluxes at the selected time
    points (conc/flux_interp), and if save_concs/fluxes is set to True it also returns dataframes for the
    concentrations/fluxes for all time points.


    Args:
        raw_data_dir: path to folder with simulation and model data.
        model_name: name of the model
        simulation_name: name of the simulation
        n_models: number of models in the simulation to plot
        time_points: time points to be plotted
        save_concs: whether or not to save concentrations for all time points
        save_fluxes: whether or not to save fluxes for all time points

    Returns:
        Dataframes for concentrations, interpolated concentrations, fluxes, interpolated fluxes, and lists of
          metabolite and reaction names.
    """

    met_names, rxn_names = get_met_rxn_names(raw_data_dir, model_name)

    # import model ensemble
    file_in = os.path.join(raw_data_dir, f'{model_name}.mat')
    mat = scipy.io.loadmat(file_in, squeeze_me=False)

    ref_conc_dic = import_ref_conc(mat, n_models)
    ref_flux_dic = import_ref_flux(mat)

    # load simulation data file
    file_in = os.path.join(raw_data_dir, f'simulation_{simulation_name}.mat')
    mat = scipy.io.loadmat(file_in, squeeze_me=False)

    # import simulation data
    conc, conc_interp, flux, flux_interp = gather_sim_data(mat, met_names, rxn_names, n_models, time_points,
                                                           save_concs=save_concs, save_fluxes=save_fluxes,
                                                           ref_conc_dic=ref_conc_dic, ref_flux_dic=ref_flux_dic)

    return conc, conc_interp, flux, flux_interp, met_names, rxn_names