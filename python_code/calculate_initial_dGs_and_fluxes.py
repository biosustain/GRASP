import numpy as np
import pandas as pd


def _get_dG_list(rxn_names, stoic_matrix, sub_conc, prod_conc, dG_std, gas_constant, temperature):

    dG_list = []
    dG_Q_list = []
    ma_ratio_list = []

    for rxn_i in range(len(rxn_names)):
        rxn_subs_conc = stoic_matrix[rxn_i, :] * sub_conc
        rxn_prods_conc = stoic_matrix[rxn_i, :] * prod_conc
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


def calculate_dG(file_in, gas_constant, temperature, rxn_order=None):
    """
    Given a GRASP input file, calculates the minimum and maximum reaction dGs based on the standard dGs in thermoRxns
    and metabolite concentrations in thermoMets.
    It also calculates the mass-action ratio and the part of the dG based on the mass-action ratio.

    :param file_in: path to the GRASP input file.
    :param gas_constant: the gas constant to calculate the Gibbs energy.
    :param temperature: the temperature to calculate the Gibbs energy.
    :param rxn_order: a list with the reactions order (optional).
    :return: ma_df, dG_Q_df, dG_df
    """

    dG_Q_df = pd.DataFrame()
    dG_df = pd.DataFrame()
    ma_df = pd.DataFrame()

    stoic_df = pd.read_excel(file_in, sheet_name='stoic')
    stoic_df = stoic_df.iloc[:5, :]

    mets_conc_df = pd.read_excel(file_in, sheet_name='thermoMets')
    mets_conc_df['mean (M)'] = (mets_conc_df['min (M)'] + mets_conc_df['max (M)']) / 2.

    dG_std_df = pd.read_excel(file_in, sheet_name='thermoRxns')
    dG_std_df = dG_std_df.iloc[:5, :]
    dG_std_df['∆Gr_mean'] = (dG_std_df['∆Gr\'_min (kJ/mol)'] + dG_std_df['∆Gr\'_max (kJ/mol)']) / 2.

    rxn_names = stoic_df['rxn ID'].values

    stoic_matrix = stoic_df.iloc[:, 1:].values

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


def _compute_robust_fluxes(stoic_matrix, meas_rates, meas_rates_std):
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
        v_mean[np.where(v_mean == 0)] = meas_rates[id_meas]
        v_std[id_unkn] = np.diag(np.matmul(
            np.matmul(np.matmul(np.matmul(np.linalg.pinv(stoic_unkn), stoic_meas), Dm), np.transpose(stoic_meas)),
            np.transpose(np.linalg.pinv(stoic_unkn))))
        v_std[np.where(v_std == 0)] = np.diag(Dm)
    else:
        print('System is not determined')
        exit()

    v_std = np.sqrt(v_std)  # Compute std

    return v_mean, v_std


def _get_balanced_s_matrix(file_in):
    stoic_df = pd.read_excel(file_in, sheet_name='stoic')
    stoic_matrix = np.transpose(stoic_df.iloc[:, 1:].values)
    rxn_list = stoic_df['rxn ID'].values

    mets_df = pd.read_excel(file_in, sheet_name='mets')
    balanced_mets_ind = np.where(mets_df['balanced?'].values == 1)

    stoic_balanced = stoic_matrix[balanced_mets_ind, :][0]

    return stoic_balanced, rxn_list


def _get_meas_rates(file_in, rxn_list):
    meas_rates_df = pd.read_excel(file_in, sheet_name='measRates')
    meas_rates_ids = meas_rates_df['Fluxes (umol/gCDW/h)'].values

    meas_rates_mean = np.zeros(len(rxn_list))
    meas_rates_std = np.zeros(len(rxn_list))
    for rxn_i, rxn in enumerate(rxn_list):
        for meas_rxn_i in range(len(meas_rates_df.index)):
            if rxn == meas_rates_ids[meas_rxn_i]:
                meas_rates_mean[rxn_i] = meas_rates_df['vref_mean'].values[meas_rxn_i]
                meas_rates_std[rxn_i] = meas_rates_df['vref_std'].values[meas_rxn_i]

    return meas_rates_mean, meas_rates_std


def _get_inactive_reactions(file_in):
    reactions_df = pd.read_excel(file_in, sheet_name='rxns')

    inactive_rxns_ind = np.where(reactions_df['modelled?'].values == 0)

    return inactive_rxns_ind


def get_robust_fluxes(file_in, rxn_order=None):
    """
    Give a GRASP input file, it calculates the robust fluxes (almost) as in GRASP, unless the system is not determined.

    :param file_in: path to the GRASP input file.
    :param rxn_order: a list with the reactions order (optional).
    :return: fluxes_df
    """

    fluxes_df = pd.DataFrame()
    stoic_balanced, rxn_list = _get_balanced_s_matrix(file_in)
    n_reactions = len(rxn_order)

    meas_rates_mean, meas_rates_std = _get_meas_rates(file_in, rxn_list)
    inactive_rxns_ind = _get_inactive_reactions(file_in)

    v_mean, v_std = _compute_robust_fluxes(stoic_balanced, meas_rates_mean, meas_rates_std)

    v_mean[inactive_rxns_ind] = 0
    v_std[inactive_rxns_ind] = 0
    v_mean = v_mean[:n_reactions]
    v_std = v_std[:n_reactions]

    fluxes_df['flux'] = v_mean
    fluxes_df['flux_std'] = v_std

    fluxes_df.index = rxn_list[:n_reactions]
    if rxn_order:
        fluxes_df = fluxes_df.reindex(rxn_order)

    return fluxes_df
