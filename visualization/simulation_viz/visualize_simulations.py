import pandas as pd
import matplotlib.pyplot as plt


COLOR_LIST = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6',
              '#6a3d9a', '#ffff99', '#b15928']


def _plot_column(data_df, time_points, data_type, names_lists, ax_list, col_i, x_scale, y_scale, x_lim, y_lim):

    for row_i, name_list in enumerate(names_lists):
        for i, name in enumerate(name_list):
            ax_list[row_i, col_i].plot(time_points, data_df[data_df[data_type] == name]['median'], color=COLOR_LIST[i % 12], label=name)
            ax_list[row_i, col_i].fill_between(time_points, data_df[data_df[data_type] == name]['q025'],
                                      data_df[data_df[data_type] == name]['q075'], color=COLOR_LIST[i % 12], alpha=0.3)

        if data_type == 'rxn':
            ax_list[row_i, col_i].set_yscale('symlog')
            ax_list[row_i, col_i].set_ylim([-10 ** 2, 10 ** 4])
        if data_type == 'met':
            ax_list[row_i, col_i].set_yscale('log')
            ax_list[row_i, col_i].set_ylim([10**-12, 10**-2])

        if col_i == 1:
            ax_list[row_i, col_i].set_xscale('log')

        if x_lim:
            ax_list[row_i, col_i].set_xlim(x_lim[col_i])
        if y_lim:
            ax_list[row_i, col_i].set_ylim(y_lim[col_i])

        if x_scale:
            ax_list[row_i, col_i].set_xscale(x_scale)
        if y_scale:
            ax_list[row_i, col_i].set_yscale(y_scale)

        ax_list[row_i, col_i].grid()
        ax_list[row_i, col_i].legend()

    return ax_list


def plot_ensemble(data_df: pd.DataFrame, time_points: list, data_type: str, selected_data: list, x_scale: str = 'linear',
                  y_scale: str = 'linear', x_lim: tuple = None, y_lim: tuple = None, save_plot: bool = False,
                  output_file: str = ''):
    """
    Does the same as plot_ensemble_interactive but using matplotlib instead of altair.
    Takes in a pandas dataframe with median and 1st and 3rd quartile concentration or flux values calculated
    over all models in the ensemble for each time point, and plots them using matplotlib.
    The median is represented by a line, and uncertainty is given by the 25% quantile and the 75% one.

    Args:
        data_df: dataframe with summarized concentration or flux data for each time point.
        time_points: time points in the simulation to be plotted.
        data_type: the column name for the concentrations or fluxes to be plotted.
        selected_data: a list of metabolite or reaction names whose concentrations or fluxes will be plotted.
        plot_title: title for the plot file
        x_scale: whether the x-scale should be log, symlog, or linear.
        y_scale: whether the y-scale should be log, symlog, or linear.
        x_lim: the limits for the x-axis.
        y_lim: the limits for the y-axis.
        save_plot: whether or not to save the plot.
        output_file: path to plot file to be saved, if save_plot == True.

    Returns:
        A matplotlib plot with the median, 1st quartile, and 3rd quartile metabolite concentrations or reaction fluxes.
    """

    n_rows = len(selected_data)
    fig, ax_list = plt.subplots(n_rows, 2, figsize=(20, n_rows*5))

    ax_list = _plot_column(data_df, time_points, data_type, selected_data, ax_list, 0, x_scale, y_scale, x_lim, y_lim)
    ax_list = _plot_column(data_df, time_points, data_type, selected_data, ax_list, 1, x_scale, y_scale, x_lim, y_lim)

    plt.tight_layout()

    if save_plot:
        plt.savefig(output_file)
        plt.close()

    else:
        return plt.show()


def plot_model():
    return 0
