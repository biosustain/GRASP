import pandas as pd
import matplotlib.pyplot as plt


COLOR_LIST = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6',
              '#6a3d9a', '#ffff99', '#b15928']


def plot_ensemble(data_df: pd.DataFrame, quant_type: str, selected_data: list, x_scale: str = 'linear',
                  y_scale: str = 'linear', x_lim: tuple = None, y_lim: tuple = None, fig_size: tuple = None,
                  save_plot: bool = False, output_file: str = ''):
    """
    Does the same as plot_ensemble_interactive but using matplotlib instead of altair.
    Takes in a pandas dataframe with median and 1st and 3rd quartile concentration or flux values calculated
    over all models in the ensemble for each time point, and plots them using matplotlib.
    The median is represented by a line, and uncertainty is given by the 25% quantile and the 75% one.

    Args:
        data_df: dataframe with summarized concentration or flux data for each time point.
        quant_type: the column name for the concentrations or fluxes to be plotted.
        selected_data: a list of metabolite or reaction names whose concentrations or fluxes will be plotted.
        x_scale: whether the x-scale should be log, symlog, or linear.
        y_scale: whether the y-scale should be log, symlog, or linear.
        x_lim: the limits for the x-axis.
        y_lim: the limits for the y-axis.
        fig_size: figure size in inches.
        save_plot: whether or not to save the plot.
        output_file: path to plot file to be saved, if save_plot == True.

    Returns:
        A matplotlib plot with the median, 1st quartile, and 3rd quartile metabolite concentrations or reaction fluxes.
    """

    if quant_type in {'conc_rel', 'conc_abs', 'flux_rel', 'flux_abs'}:

        data_type = 'met' if quant_type.find('conc') != -1 else 'rxn'

        n_rows = len(selected_data)
        fig, ax_list = plt.subplots(n_rows, figsize=(10, n_rows*5) if fig_size is None else fig_size)

        if selected_data is None:
            selected_data = data_df[data_type].unique()

        for row_i, name_list in enumerate(selected_data):

            for i, name in enumerate(name_list):
                if len(selected_data) == 1:
                    ax = ax_list
                else:
                    ax = ax_list[row_i]

                try:

                    ax.plot(data_df[data_df[data_type] == name]['time_point'].values, data_df[data_df[data_type] == name]['median'].values, color=COLOR_LIST[i % 12],
                            label=name)
                except ValueError:
                    raise ValueError(f'It is likely that {name} is not in the {data_type} column of the input data frame.')

                ax.fill_between(data_df[data_df[data_type] == name]['time_point'].values, data_df[data_df[data_type] == name]['q025'],
                                data_df[data_df[data_type] == name]['q075'], color=COLOR_LIST[i % 12], alpha=0.3)

            if x_lim:
                ax.set_xlim(x_lim)
            if y_lim:
                ax.set_ylim(y_lim)

            if x_scale:
                ax.set_xscale(x_scale)
            if y_scale:
                ax.set_yscale(y_scale)

            ax.set_xlabel('time')
            ax.set_ylabel(quant_type)

            ax.grid()
            ax.legend()

        plt.tight_layout()

        if save_plot:
            plt.savefig(output_file)
            plt.close()

        else:
            return plt.show()

    else:
        print('The data_type variable should be either conc_rel, conc_abs, flux_rel, or flux_abs.')
        return


def plot_model(data_df: pd.DataFrame, model_i: int, quant_type: str, selected_data: list,  x_scale: str = 'linear',
               y_scale: str = 'linear', x_lim: tuple = None, y_lim: tuple = None, fig_size: tuple = None,
               save_plot: bool = False, output_file: str = ''):
    """
    Given a pandas dataframe with metabolite concentrations or reaction fluxes, it plots the selected ones
    (specified in selected data) for the selected model.
    The plots are made using matplotlib.

    Args:
        data_df: dataframe with metabolite concentrations or reaction fluxes.
        model_i: number of the model to plot.
        quant_type: the column name for the concentrations or fluxes to be plotted.
        selected_data: a list of metabolite or reaction names whose concentrations or fluxes will be plotted.
        x_scale: whether the x-scale should be log, symlog, or linear.
        y_scale: whether the y-scale should be log, symlog, or linear.
        x_lim: the limits for the x-axis.
        y_lim: the limits for the y-axis.
        fig_size: figure size in inches.
        save_plot: whether or not to save the plot.
        output_file: path to plot file to be saved, if save_plot == True.

    Returns:
        An altair plot with the selected metabolite concentrations or reaction fluxes plotted for a given model.
    """

    if quant_type in {'conc_rel', 'conc_abs', 'flux_rel', 'flux_abs'}:

        data_type = 'met' if quant_type.find('conc') != -1 else 'rxn'

        fig, ax = plt.subplots(1, 1, figsize=(10, 5) if fig_size is None else fig_size)

        if selected_data is None:
            selected_data = data_df[data_type].unique()

        data_interp_df = data_df[
            (data_df['model'] == model_i) & (data_df[data_type].isin(selected_data))]

        for i, name in enumerate(selected_data):
            data_temp = data_interp_df[data_interp_df[data_type] == name]
            ax.plot(data_temp['time_point'].values, data_temp[quant_type].values, color=COLOR_LIST[i % 12],
                    label=name, marker='o')

        data_df = data_df[(data_df['model'] == model_i) & (data_df[data_type].isin(selected_data))]

        for i, name in enumerate(selected_data):
            data_temp = data_df[data_df[data_type] == name]
            ax.plot(data_temp['time_point'].values, data_temp[quant_type].values,
                    color=COLOR_LIST[i % 12], label=name)

        if x_lim:
            ax.set_xlim(x_lim)
        if y_lim:
            ax.set_ylim(y_lim)

        if x_scale:
            ax.set_xscale(x_scale)
        if y_scale:
            ax.set_yscale(y_scale)

        ax.set_xlabel('time')
        ax.set_ylabel(quant_type)

        ax.grid()
        ax.legend()

        plt.tight_layout()

        if save_plot:
            plt.savefig(output_file)
            plt.close()

        else:
            return plt.show()

    else:
        print('The data_type variable should be either conc_rel, conc_abs, flux_rel, or flux_abs.')
        return
