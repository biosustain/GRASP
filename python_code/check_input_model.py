import numpy as np
import pandas as pd


def _check_met_rxn_order(data_dict):
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


def _check_kinetics_column(data_dict, col_name):
    col_data = data_dict[col_name].dropna()
    for row in col_data:
        if row.find(',') != -1 or row.find(';') != -1 or row.find('.') != -1:
            print('Make sure all metabolites are separated by a single space in column "', col_name, '" row:\n', row)


def _check_kinetics_met_separators(data_dict):
    for key in list(data_dict.keys()):
        if key.startswith('kinetics'):
            _check_kinetics_column(data_dict[key], 'order')
            _check_kinetics_column(data_dict[key], 'promiscuous')
            _check_kinetics_column(data_dict[key], 'inhibitors')
            _check_kinetics_column(data_dict[key], 'activators')
            _check_kinetics_column(data_dict[key], 'negative effector')
            _check_kinetics_column(data_dict[key], 'positive effector')


def check_input_model(file_in):
    """
    Given an excel file with the model input for GRASP it does two things
     1. gets the metabolite and reaction order in the 'stoic'
    sheet and checks if it is the same in the other sheets. If it isn't, prints out a message saying in which sheet the
    order is not the same and the lists of metabolite/reactions: both the one in the current sheet and the ones in the
    'stoic' sheet;
     2. checks all lists of metabolites in the kinetics sheet to make sure there are no commas semi-colons or dots.

    :param file_in: path to excel file with model input.
    :return: 0 or 1: 0 everything's fine, 1 there was an error.
    """

    try:
        data_dict = pd.read_excel(file_in, sheet_name=None)
    except IOError:
        print('File wasn\'t found: ', file_in)
        return 1

    _check_met_rxn_order(data_dict)
    _check_kinetics_met_separators(data_dict)

    return 0


def main():
    file_in = '/
    check_input_model(file_in)


if __name__ == '__main__':
    main()