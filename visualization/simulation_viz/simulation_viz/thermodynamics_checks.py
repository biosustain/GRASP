import numpy as np
import pandas as pd


def _get_dG_list(rxn_names: list, stoic_matrix: np.ndarray, sub_conc: np.ndarray, prod_conc: np.ndarray,
                 dG_std: np.ndarray, gas_constant: float, temperature: float) -> tuple:

    dG_list = []
    dG_Q_list = []
    ma_ratio_list = []

    for rxn_i in range(len(rxn_names)):

        rxn_subs_conc = np.sign(stoic_matrix[rxn_i, :]) * (sub_conc ** np.abs(stoic_matrix[rxn_i, :]))
        rxn_prods_conc = np.sign(stoic_matrix[rxn_i, :]) * (prod_conc ** np.abs(stoic_matrix[rxn_i, :]))

        subs_ind = np.where(rxn_subs_conc < 0)
        subs_conc = rxn_subs_conc[subs_ind]
        prods_ind = np.where(rxn_prods_conc > 0)
        prods_conc = rxn_prods_conc[prods_ind]

        subs_prod = np.abs(np.prod(subs_conc))
        prods_prod = np.prod(prods_conc)

        ma_ratio = prods_prod / subs_prod

        dG_Q = gas_constant * temperature * np.log(ma_ratio)
        dG = dG_std[rxn_i] + dG_Q

        ma_ratio_list.append(ma_ratio)
        dG_list.append(dG)
        dG_Q_list.append(dG_Q)

    return dG_list, dG_Q_list, ma_ratio_list


def calculate_dG(data_dict: dict, gas_constant: float, temperature: float, rxn_order: list = None) -> tuple:
    """
    Given a dictionary representing a GRASP input file, calculates the minimum and maximum reaction dGs based on the
    standard dGs in thermoRxns and metabolite concentrations in thermoMets.
    It also calculates the mass-action ratio and the part of the dG based on the mass-action ratio.

    Args:
        data_dict: a dictionary that represents the excel file with the GRASP model.
        gas_constant: the gas constant to calculate the Gibbs energy.
        temperature: the temperature to calculate the Gibbs energy.
        rxn_order: a list with the reactions order (optional).

    Returns:
        Mass action ratio dataframe, dG_Q dataframe, Gibbs energies dataframe.
    """

    dG_Q_df = pd.DataFrame()
    dG_df = pd.DataFrame()
    ma_df = pd.DataFrame()

    stoic_df = data_dict['stoic']

    mets_conc_df = data_dict['thermoMets']
    mets_conc_df['mean (M)'] = (mets_conc_df['min (M)'] + mets_conc_df['max (M)']) / 2.

    dG_std_df = data_dict['thermoRxns']
    dG_std_df['∆Gr_mean'] = (dG_std_df['∆Gr\'_min (kJ/mol)'] + dG_std_df['∆Gr\'_max (kJ/mol)']) / 2.

    rxn_names = stoic_df.index.values

    stoic_matrix = stoic_df.values

    min_met_conc = mets_conc_df['min (M)'].values
    max_met_conc = mets_conc_df['max (M)'].values

    dG_list_mean, dG_Q_list_mean, ma_ratio_list_mean = _get_dG_list(rxn_names, stoic_matrix,
                                                                    mets_conc_df['mean (M)'].values,
                                                                    mets_conc_df['mean (M)'].values,
                                                                    dG_std_df['∆Gr_mean'].values,
                                                                    gas_constant, temperature)
    dG_list_min, dG_Q_list_min, ma_ratio_list_min = _get_dG_list(rxn_names, stoic_matrix, max_met_conc, min_met_conc,
                                                                 dG_std_df['∆Gr\'_min (kJ/mol)'].values,
                                                                 gas_constant, temperature)
    dG_list_max, dG_Q_list_max, ma_ratio_list_max = _get_dG_list(rxn_names, stoic_matrix, min_met_conc, max_met_conc,
                                                                 dG_std_df['∆Gr\'_max (kJ/mol)'].values,
                                                                 gas_constant, temperature)

    ma_df['ma_min'] = ma_ratio_list_min
    ma_df['ma_mean'] = ma_ratio_list_mean
    ma_df['ma_max'] = ma_ratio_list_max

    dG_Q_df['∆G_Q_min'] = dG_Q_list_min
    dG_Q_df['∆G_Q_mean'] = dG_Q_list_mean
    dG_Q_df['∆G_Q_max'] = dG_Q_list_max

    dG_df['∆G_min'] = dG_list_min
    dG_df['∆G_mean'] = dG_list_mean
    dG_df['∆G_max'] = dG_list_max

    ma_df.index = rxn_names
    dG_Q_df.index = rxn_names
    dG_df.index = rxn_names

    if rxn_order:
        ma_df = ma_df.reindex(rxn_order)
        dG_Q_df = dG_Q_df.reindex(rxn_order)
        dG_df = dG_df.reindex(rxn_order)

    return ma_df, dG_Q_df, dG_df


def _compute_robust_fluxes(stoic_matrix: np.ndarray, meas_rates: np.ndarray, meas_rates_std: np.ndarray,
                           rxn_list: list) -> tuple:

    # Determine measured fluxes and decompose stoichiometric matrix
    id_meas = np.where(meas_rates != 0)
    id_unkn = np.where(meas_rates == 0)

    stoic_meas = stoic_matrix[:, id_meas]
    stoic_meas = np.array([row[0] for row in stoic_meas])

    stoic_unkn = stoic_matrix[:, id_unkn]
    stoic_unkn = np.array([row[0] for row in stoic_unkn])

    # Initialize final fluxes
    v_mean = np.zeros(np.size(meas_rates))
    v_std = np.zeros(np.size(meas_rates))

    # Compute estimate Rred
    Dm = np.diag(meas_rates_std[id_meas] ** 2)
    Rred = np.subtract(stoic_meas, np.matmul(np.matmul(stoic_unkn, np.linalg.pinv(stoic_unkn)), stoic_meas))
    [u, singVals, vh] = np.linalg.svd(Rred)
    singVals = np.abs(singVals)
    zero_sing_vals = np.where(singVals > 10 ** -12)

    # If the system is fully determined, compute as follows
    if len(zero_sing_vals[0]) == 0:
        v_mean[id_unkn] = -np.matmul(np.matmul(np.linalg.pinv(stoic_unkn), stoic_meas), meas_rates[id_meas])
        if len(np.where(v_mean == 0)[0]) > len(id_meas[0]):
            zero_flux_rxns = list(set(np.where(v_mean == 0)[0]).difference(set(id_meas[0])))
            raise RuntimeError('According to compute robust fluxes, there are reactions with zero flux in the model.\n'+
                                'Those reactions should be removed.\n' +
                                f'The reactions are {np.array(rxn_list)[[zero_flux_rxns]]}')

        v_mean[np.where(v_mean == 0)] = meas_rates[id_meas]
        v_std[id_unkn] = np.diag(np.matmul(
            np.matmul(np.matmul(np.matmul(np.linalg.pinv(stoic_unkn), stoic_meas), Dm), np.transpose(stoic_meas)),
            np.transpose(np.linalg.pinv(stoic_unkn))))
        v_std[np.where(v_std == 0)] = np.diag(Dm)
    else:
        print('System is not fully determined and the fluxes cannot be determined.')
        exit()

    v_std = np.sqrt(v_std)  # Compute std

    return v_mean, v_std


def _get_balanced_s_matrix(data_dict: dict) -> tuple:

    stoic_df = data_dict['stoic']
    stoic_matrix = np.transpose(stoic_df.values)
    rxn_list = stoic_df.index.values

    mets_df = data_dict['mets']
    balanced_mets_ind = np.where(mets_df['balanced?'].values == 1)

    stoic_balanced = stoic_matrix[balanced_mets_ind, :][0]

    return stoic_balanced, rxn_list


def _get_meas_rates(data_dict: dict, rxn_list: list) -> tuple:

    meas_rates_df = data_dict['measRates']
    meas_rates_ids = meas_rates_df.index.values

    meas_rates_mean = np.zeros(len(rxn_list))
    meas_rates_std = np.zeros(len(rxn_list))

    meas_rates = zip(list(np.nonzero(np.in1d(rxn_list, meas_rates_ids))[0]), list(meas_rates_ids))

    for meas_rxn_ind, meas_rxn in meas_rates:
        meas_rates_mean[meas_rxn_ind] = meas_rates_df.loc[meas_rxn, 'vref_mean (mmol/L/h)']
        meas_rates_std[meas_rxn_ind] = meas_rates_df.loc[meas_rxn, 'vref_std (mmol/L/h)']

    return meas_rates_mean, meas_rates_std


def _get_inactive_reactions(data_dict: dict) -> np.ndarray:
    reactions_df = data_dict['rxns']

    inactive_rxns_ind = np.where(reactions_df['modelled?'].values == 0)

    return inactive_rxns_ind


def get_robust_fluxes(data_dict: dict, rxn_order: list = None) -> pd.DataFrame:
    """
    Given a dictionary representing a GRASP input file, it calculates the robust fluxes (almost) as in GRASP,
    unless the system is not fully determined.

    Args:
        data_dict: path to the GRASP input file
        rxn_order: a list with the reactions order (optional)

    Returns:
        fluxes_df: dataframe with flux mean and std values
    """

    fluxes_df = pd.DataFrame()
    stoic_balanced, rxn_list = _get_balanced_s_matrix(data_dict)
    # n_reactions = len(rxn_order)

    meas_rates_mean, meas_rates_std = _get_meas_rates(data_dict, rxn_list)
    inactive_rxns_ind = _get_inactive_reactions(data_dict)

    v_mean, v_std = _compute_robust_fluxes(stoic_balanced, meas_rates_mean, meas_rates_std, rxn_list)

    v_mean[inactive_rxns_ind] = 0
    v_std[inactive_rxns_ind] = 0

    fluxes_df['vref_mean (mmol/L/h)'] = v_mean
    fluxes_df['vref_std (mmol/L/h)'] = v_std

    fluxes_df.index = rxn_list
    if rxn_order:
        fluxes_df = fluxes_df.reindex(rxn_order)
        fluxes_df = fluxes_df.reindex(rxn_order)

    return fluxes_df


def check_thermodynamic_feasibility(data_dict: dict) -> tuple:
    """
    Given a dictionary representing a GRASP input file, it checks if the reaction's dG are compatible with the
    respective fluxes. It works both when all fluxes are specified in measRates and when robust fluxes are calculated
    for a fully determined system. If the fluxes are not fully specified not the system is fully determined, it
    doesn't work.

    Args:
        data_dict: a dictionary representing a GRASP input file.

    Returns:
        Whether or not the model is thermodynamically feasible plus fluxes and Gibbs energies dataframes.

    """

    print('\nChecking if fluxes and Gibbs energies are compatible.\n')

    flag = False
    temperature = 298  # in K
    gas_constant = 8.314 * 10**-3  # in kJ K^-1 mol^-1

    stoic_df = data_dict['stoic']
    flux_df = data_dict['measRates']

    ma_df, dG_Q_df, dG_df = calculate_dG(data_dict, gas_constant, temperature)

    if len(stoic_df.index) != len(flux_df.index):
        flux_df = get_robust_fluxes(data_dict)

    for rxn in flux_df.index:
        if flux_df.loc[rxn, 'vref_mean (mmol/L/h)'] > 0 and dG_df.loc[rxn, '∆G_min'] > 0:
            print(f'The flux and ∆G range seem to be incompatible for reaction {rxn}')
            flag = True

        if flux_df.loc[rxn, 'vref_mean (mmol/L/h)'] < 0 and dG_df.loc[rxn, '∆G_max'] < 0:
            print(f'The flux and ∆G range seem to be incompatible for reaction {rxn}')
            flag = True

    if flag is False:
        print('Everything seems to be OK.')

    return flag, flux_df, dG_df
