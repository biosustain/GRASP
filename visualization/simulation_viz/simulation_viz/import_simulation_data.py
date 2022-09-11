import os
import time

import numpy as np
import pandas as pd
import scipy.io

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
    met_names = pd.read_csv(file_met_names, sep='\t').values
    met_names = list(met_names.transpose()[0])
    met_names = [met_name[2:] for met_name in met_names]

    # get reaction names
    file_rxn_names = os.path.join(raw_data_dir, f'{model_name}_rxnsActive.dat')
    rxn_names = pd.read_csv(file_rxn_names, sep='\t').values
    rxn_names = list(rxn_names.transpose()[0])
    rxn_names = [rxn_name[2:] for rxn_name in rxn_names]

    return met_names, rxn_names


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


def gather_conc_data(mat: dict, met_names: list, n_models: int, quant_type: str, ref_conc_dic: dict = None) \
        -> pd.DataFrame:
    """
    Imports the concentration data over time from the matlab simulation results.

    The concentrations are relative to the reference state.

    If a dictionary with reference concentrations is supplied the resulting dataframe will both relative and absolute
    metabolite concentrations. The columns are named conc_rel and conc_abs, respectively.

    Pre-allocated numpy arrays are used for memory performance (vs. previous list extend).

    Args:
        mat: matlab structures with simulation results.
        met_names: names of the simulated metabolites in the model.
        n_models: number of models to be simulated.
        quant_type: whether to get relative or absolute concentrations. Possible values: conc_abs, conc_rel.
        ref_conc_dic: a dictionary with the metabolites reference concentrations in each model.

    Returns:
        A pandas dataframe with all concentration values at each time point.
    """

    if quant_type not in {'conc_abs', 'conc_rel'}:
        print('quant_type should be set to either "conc_abs" for absolute concentrations or "conc_rel" for relative '
              'concentrations')
        return pd.DataFrame()

    time_points = mat['simulationRes'][0][0]['t'][0][0].transpose()[0]
    n_time_points = len(time_points)
    n_mets = len(met_names)

    conc_model = np.full(n_time_points * n_models * n_mets, np.nan)
    conc_met = list(np.full(n_time_points * n_models * n_mets, np.nan))
    conc_time_point = np.full(n_time_points * n_models * n_mets, np.nan)
    conc = np.full(n_time_points * n_models * n_mets, np.nan)

    n_missing_models = 0
    start_total = time.time()
    for model_i in range(n_models):

        try:
            time_points = mat['simulationRes'][0][model_i]['t'][0][0].transpose()[0]
            met_concs = mat['simulationRes'][0][model_i]['conc'][0][0].transpose()

        except IndexError:
            n_missing_models += 1
            print(f'{model_i} is missing.')
            continue

        for met_i, met in enumerate(met_names):
            start_i = model_i * n_time_points * n_mets + met_i * n_time_points
            end_i = model_i * n_time_points * n_mets + (met_i + 1) * n_time_points

            conc_model[start_i:end_i] = np.repeat(model_i, len(time_points))
            conc_time_point[start_i:end_i] = time_points
            conc_met[start_i:end_i] = np.repeat(met, len(time_points))
            if quant_type == 'conc_abs':
                conc[start_i:end_i] = met_concs[met_i] * ref_conc_dic[met][model_i]
            else:
                conc[start_i:end_i] = met_concs[met_i]

    data_df_conc = pd.DataFrame(data={'model': conc_model.astype(int), 'met': conc_met, 'time_point': conc_time_point,
                                      quant_type: conc})
    data_df_conc.dropna(axis=0, inplace=True)

    print(f'total time: {time.time() - start_total}')
    print(f'There were a total of {n_missing_models} missing models out of {n_models}.')

    return data_df_conc


def gather_flux_data(mat: dict, rxn_names: list, n_models: int, quant_type: str, ref_flux_dic: dict = None) \
        -> pd.DataFrame:
    """
    Imports the flux data over time from the matlab simulation results.

    The fluxes coming from matlab are absolute.

    If a dictionary with reference fluxes is supplied the resulting dataframe will both relative and absolute
    reaction fluxes. The columns are named flux_rel and flux_abs, respectively.

    Pre-allocated numpy arrays are used for memory performance (vs. previous list extend).

    Args:
        mat: matlab structures with simulation results.
        rxn_names: names of the simulated reactions in the model.
        n_models: number of models to be simulated.
        quant_type: whether to get relative or absolute fluxes. Possible values: flux_abs, flux_rel.
        ref_flux_dic: a dictionary with the reactions reference fluxes.

    Returns:
        A pandas dataframe with all concentration and flux values at each time point.
    """

    if quant_type not in {'flux_abs', 'flux_rel'}:
        print('quant_type should be set to either "flux_abs" for absolute fluxes or "flux_rel" for relative fluxes')
        return pd.DataFrame()

    time_points = mat['simulationRes'][0][0]['t'][0][0].transpose()[0]
    n_time_points = len(time_points)
    n_rxns = len(rxn_names)

    flux_model = np.full(n_time_points * n_models * n_rxns, np.nan)
    flux_rxn = list(np.full(n_time_points * n_models * n_rxns, np.nan))
    flux_time_point = np.full(n_time_points * n_models * n_rxns, np.nan)
    flux = np.full(n_time_points * n_models * n_rxns, np.nan)

    n_missing_models = 0
    start_total = time.time()
    for model_i in range(n_models):

        try:
            time_points = mat['simulationRes'][0][model_i]['t'][0][0].transpose()[0]
            rxn_fluxes = mat['simulationRes'][0][model_i]['flux'][0][0].transpose()

        except IndexError:
            n_missing_models += 1
            print(f'{model_i} is missing.')
            continue

        for rxn_i, rxn in enumerate(rxn_names):
            start_i = model_i * n_time_points * n_rxns + rxn_i * n_time_points
            end_i = model_i * n_time_points * n_rxns + (rxn_i + 1) * n_time_points

            flux_model[start_i:end_i] = np.repeat(model_i, len(time_points))
            flux_time_point[start_i:end_i] = time_points
            flux_rxn[start_i:end_i] = np.repeat(rxn_names[rxn_i], len(time_points))

            if quant_type == 'flux_abs':
                flux[start_i:end_i] = rxn_fluxes[rxn_i]
            else:
                flux[start_i:end_i] = rxn_fluxes[rxn_i] / ref_flux_dic[rxn]

    data_df_flux = pd.DataFrame(data={'model': flux_model.astype(int), 'rxn': flux_rxn, 'time_point': flux_time_point,
                                      quant_type: flux})
    data_df_flux.dropna(axis=0, inplace=True)

    print(f'total time: {time.time() - start_total}')
    print(f'There were a total of {n_missing_models} missing models out of {n_models}.')

    return data_df_flux


def aggregate_time_series(data_df: pd.DataFrame, quant_type: str, name_list: list) -> pd.DataFrame:
    """
    Takes in a dataframe with either fluxes or concentrations and calculates the median and respective first and third
    quartiles for each flux or metabolite at a given time point over the whole model ensemble.

    Args:
        data_df: dataframe with concentrations or fluxes for each time point and model.
        quant_type: whether the data frame contains relative concentrations ('conc_rel'), absolute concentrations
                    ('conc_abs'), relative fluxes ('flux_rel'), or absolute_fluxes ('flux_abs').
        name_list: reaction or metabolite names, depending on what that the dataframe contains.

    Returns:
        A pandas dataframe with the concentrations/fluxes median and 1st and 3rd quartiles for each metabolite/reaction.
    """

    data_type = 'met' if quant_type.find('conc') != -1 else 'rxn'
    time_points = data_df['time_point'].unique()
    n_time_points = len(time_points)

    quartiles_dic = {'time_point': [], data_type: [], 'q025': np.array([]), 'median': np.array([]),
                     'q075': np.array([])}

    q025 = data_df.groupby([data_type, 'time_point']).quantile(q=0.25)
    q050 = data_df.groupby([data_type, 'time_point']).median()
    q075 = data_df.groupby([data_type, 'time_point']).quantile(q=0.75)

    for item in name_list:
        quartiles_dic[data_type].extend(np.repeat(item, n_time_points))
        quartiles_dic['time_point'].extend(time_points)

        quartiles_dic['q025'] = np.append(quartiles_dic['q025'], q025.loc[item, :][quant_type])
        quartiles_dic['median'] = np.append(quartiles_dic['median'], q050.loc[item, :][quant_type])
        quartiles_dic['q075'] = np.append(quartiles_dic['q075'], q075.loc[item, :][quant_type])

    quartiles_df = pd.DataFrame.from_dict(quartiles_dic)

    return quartiles_df


def get_time_series_quartiles(mat, names_list: list, n_models: int, ref_dic: dict, quant_type: str) -> pd.DataFrame:
    """
    Given a .mat file from matlab with the simulation results, and the type of quantity the user wants to extract:
      - "flux_rel" - relative fluxes
      - "flux_abs" - absolute fluxes
      - "conc_rel" - relative concentrations
      - "conc_abs" - absolute concentrations
    it extracts the data and calculates the median, the first and third quartile for each model, metabolite/reaction,
    and time point.

    Args:
        mat: .mat file with simulations from matlab.
        names_list: list of reaction or metabolite names.
        n_models: number of models for which to extract the simulations.
        ref_dic: dictionary with either reference concentrations or fluxes.
        quant_type: which quantity to extract: "flux_rel", "flux_abs", "conc_rel", "conc_abs".

    Returns:
        A dataframe with flux/concentrations median+first+third quartiles for each each model, metabolite/reaction,
        and time point.
    """

    if quant_type not in {'flux_abs', 'flux_rel', 'conc_abs', 'conc_rel'}:
        print('quant_type should be set to one of the following values; "flux_abs", "flux_rel", "conc_abs", '
              'or "conc_rel"')
        return pd.DataFrame()

    if quant_type.find('flux') == -1:
        conc_df = gather_conc_data(mat, names_list, n_models, quant_type, ref_dic)
        quartiles_df = aggregate_time_series(conc_df, quant_type, names_list)
    else:
        flux_df = gather_flux_data(mat, names_list, n_models, quant_type, ref_dic)
        quartiles_df = aggregate_time_series(flux_df, quant_type, names_list)

    return quartiles_df


def load_simulation(raw_data_dir: str, model_name: str, simulation_name: str, n_models: int) -> tuple:
    """
    Takes in the path to the folder with the model ensemble and simulation data, loads it and returns the simulation
    Matlab file, the reference concentration and flux dictionaries, as well as the metabolite and reaction names.

    Args:
        raw_data_dir: path to folder with simulation and model data.
        model_name: name of the model
        simulation_name: name of the simulation
        n_models: number of models in the simulation to plot

    Returns:
        The matlab simulation file, reference concentration and fluxes dictionary, and lists of metabolite and
        reaction names.
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

    return mat, ref_conc_dic, ref_flux_dic, met_names, rxn_names