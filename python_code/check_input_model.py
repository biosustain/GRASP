import numpy as np
import pandas as pd
from calculate_initial_dGs_and_fluxes import calculate_dG, get_robust_fluxes


def _check_met_rxn_order(data_dict):
    flag = 0
    rxn_list = data_dict['stoic']['rxn ID'].values
    met_list = data_dict['stoic'].columns.values[1:]

    for key in list(data_dict.keys())[2:]:
        id_list = data_dict[key].iloc[:, 0].values

        if key != 'measRates':
            met_bool = id_list == met_list
            met_bool = all(met_bool) if isinstance(met_bool, np.ndarray) else met_bool

            rxn_bool = id_list == rxn_list
            rxn_bool = all(rxn_bool) if isinstance(rxn_bool, np.ndarray) else rxn_bool

            if not met_bool and not rxn_bool:
                print('Metabolite or reaction list in sheet ', key,
                      'doesn\'t match the list in the stoichiometry matrix.')
                print('Current list:\n', id_list)
                print('Metabolite list in stoichiometry matrix:\n', met_list)
                print('Reaction list in stoichiometry matrix:\n', rxn_list)
                print('\n')
                flag = 1

    return flag


def _check_kinetics_column(data_dict, col_name):
    flag = 0
    col_data = data_dict[col_name].dropna()
    for row in col_data:
        if row.find(',') != -1 or row.find(';') != -1 or row.find('.') != -1:
            print('Make sure all metabolites are separated by a single space in column "', col_name, '" row:\n', row)
            flag = 1
    
    return flag


def _check_kinetics_met_separators(data_dict):
    flag_list = []
    for key in list(data_dict.keys()):
        if key.startswith('kinetics'):
            flag = _check_kinetics_column(data_dict[key], 'order')
            flag_list.append(flag)
            flag = _check_kinetics_column(data_dict[key], 'promiscuous')
            flag_list.append(flag)
            flag = _check_kinetics_column(data_dict[key], 'inhibitors')
            flag_list.append(flag)
            flag = _check_kinetics_column(data_dict[key], 'activators')
            flag_list.append(flag)
            flag = _check_kinetics_column(data_dict[key], 'negative effector')
            flag_list.append(flag)
            flag = _check_kinetics_column(data_dict[key], 'positive effector')
            flag_list.append(flag)

    flag = 1 if 1 in flag_list else 0
    return flag


def _check_balanced_metabolites(data_dict):
    flag = 0
    stoic_df = data_dict['stoic']
    stoic_df.index = stoic_df['rxn ID']
    stoic_df = stoic_df.drop('rxn ID', axis=1)
    mets_df = data_dict['mets']

    for i, met in enumerate(stoic_df.columns):
        if stoic_df[met].gt(0).any() and stoic_df[met].lt(0).any():
            if mets_df['balanced?'][i] == 0:
                print(met, 'is marked as not balanced but it seems to be balanced.')
                flag = 1
        else:
            if mets_df['balanced?'][i] == 1:
                print(met, 'is marked as balanced but it does not seem to be balanced.')
                flag = 1
            if mets_df['fixed?'][i] == 0:
                print(met, 'is not set as constant but maybe it should, since it does not seem to be balanced.')
                flag = 1

    return flag


def _check_dG_and_flux(file_in):
    flag = 0
    temperature = 298  # in K
    gas_constant = 8.314 * 10**-3  # in kJ K^-1 mol^-1
    ma_df, dG_Q_df, dG_df = calculate_dG(file_in, gas_constant, temperature)
    flux_df = get_robust_fluxes(file_in)

    for rxn in flux_df.index:
        if flux_df.loc[rxn,'flux'] > 0 and dG_df.loc[rxn,'∆G_min'] > 0:
            print('The flux and ∆G range seem to be incompatible for reaction', rxn)
            flag = 1

        if flux_df.loc[rxn,'flux'] < 0 and dG_df.loc[rxn,'∆G_max'] < 0:
            print('The flux and ∆G range seem to be incompatible for reaction', rxn)
            flag = 1

    return flag 


def check_input_model(file_in):
    """
    Given an excel file with the model input for GRASP it does two things
     1. gets the metabolite and reaction order in the 'stoic'
    sheet and checks if it is the same in the other sheets. If it isn't, prints out a message saying in which sheet the
    order is not the same and the lists of metabolite/reactions: both the one in the current sheet and the ones in the
    'stoic' sheet;
     2. checks all lists of metabolites in the kinetics sheet to make sure there are no commas semi-colons or dots;
     3. checks if metabolites that are balanced are marked as balanced and metabolites that are not balanced are marked as
        not balanced;
     4. checks if reaction Gibbs energies are compatible with the respective fluxes, i.e. for a given reaction, if the range 
        of Gibbs energies can be negative when the flux is positive, and if the range of Gibbs energies can be positive 
        when then flux is negative.

    :param file_in: path to excel file with model input.
    :return: 0 or 1: 0 everything's fine, 1 there was an error.
    """

    try:
        data_dict = pd.read_excel(file_in, sheet_name=None)
    except IOError:
        print('File wasn\'t found: ', file_in)
        return 1

    flag_order = _check_met_rxn_order(data_dict)
    flag_kinetics_sep = _check_kinetics_met_separators(data_dict)
    flag_balanced_mets = _check_balanced_metabolites(data_dict)
    flag_dG_flux = _check_dG_and_flux(file_in)
    
    
    if flag_order == 0 and flag_kinetics_sep == 0 and flag_balanced_mets == 0 and flag_dG_flux == 0:
        print('\nYour model input seems to be all right! (take this with a grain of salt though)')
        flag = 0
    else:
        print('\nAddress the above issues before running GRASP.')
        flag = 1
    
    return flag


def main():
    file_in = '/home/mrama/GRASP_test/GRASP/input_test/HMP1489_r1_t0.xlsx'
    check_input_model(file_in)


if __name__ == '__main__':
    main()
