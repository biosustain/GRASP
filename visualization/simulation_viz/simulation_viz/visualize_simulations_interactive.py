import altair as alt
import pandas as pd
import os


def plot_ensemble_interactive(data_df: pd.DataFrame, data_type, selected_data: list = None, x_scale: str = 'linear',
                              y_scale: str = 'linear', x_lim: tuple = None, y_lim: tuple = None):
    """
    Takes in a pandas dataframe with median and 1st and 3rd quartile concentration or flux values calculated
    over all models in the ensemble for each time point, and plots them using altair.
    The median is represented by a line, and uncertainty is given by the 25% quantile and the 75% one.

    Args:
        data_df: dataframe with summarized concentration or flux data for each time point.
        data_type: the column name for the concentrations or fluxes to be plotted.
        selected_data: a list of metabolite or reaction names whose concentrations or fluxes will be plotted.
        x_scale: whether the x-scale should be log, symlog, or linear.
        y_scale: whether the y-scale should be log, symlog, or linear.
        x_lim: the limits for the x-axis.
        y_lim: the limits for the y-axis.

    Returns:
        An altair plot with the median, 1st quartile, and 3rd quartile metabolite concentrations or reaction fluxes.
    """

    if selected_data is None:
        selected_data = data_df[data_type].unique()

    plot1 = alt.Chart(data_df[data_df[data_type].isin(selected_data)]).mark_line().encode(
        x=alt.X('time_point:Q',
                scale=alt.Scale(domain=x_lim, type=x_scale)),
        y=alt.Y('median:Q',
                scale=alt.Scale(domain=y_lim, type=y_scale)
               ),
        color=f'{data_type}:N',
        tooltip=[f'{data_type}:N', 'median:Q']
    ).properties(
            width=600,
            height=400
    )

    plot2 = alt.Chart(data_df[data_df[data_type].isin(selected_data)]).mark_area().encode(
        x=alt.X('time_point:Q',
                scale=alt.Scale(domain=x_lim, type=x_scale)),
        y=alt.Y('q025:Q',
                scale=alt.Scale(domain=y_lim, type=y_scale)
                ),
        y2='q075:Q',
        color=f'{data_type}:N',
        tooltip=[f'{data_type}:N', 'median:Q'],
        opacity=alt.OpacityValue(0.3)
    ).properties(
            width=600,
            height=400
    )

    return alt.layer(plot1, plot2).interactive()


def plot_model_interactive(data_df: pd.DataFrame, data_interp_df: pd.DataFrame, model_i: int, data_type: str,
                           selected_data: list, scaled_conc: bool = False, x_scale: str = 'linear',
                           y_scale: str = 'linear', x_lim: tuple = None, y_lim: tuple = None):
    """
    Given a pandas dataframe with metabolite concentrations or reaction fluxes, it plots the selected ones
    (specified in selected data) for the selected model.
    The plots are made using altair.
    The data plotted depends on which dataframes are given, data_df contains the pure simulation results, while the
    data_interp_df contains the interpolated simulations. The latter is mandatory and is always plotted, while the
    former is optional.

    Args:
        data_df: dataframe with metabolite concentrations or reaction fluxes.
        data_interp_df: dataframe with interpolated metabolite concentrations or reactions fluxes.
        model_i: number of the model to plot.
        data_type: the column name for the concentrations or fluxes to be plotted.
        selected_data: a list of metabolite or reaction names whose concentrations or fluxes will be plotted.
        scaled_conc: whether concentrations are scaled or absolute.
        x_scale: whether the x-scale should be log, symlog, or linear.
        y_scale: whether the y-scale should be log, symlog, or linear.
        x_lim: the limits for the x-axis.
        y_lim: the limits for the y-axis.

    Returns:
        An altair plot with the selected metabolite concentrations or reaction fluxes plotted for a given model.
    """

    plot1 = alt.Chart(data_interp_df[data_interp_df[data_interp_df['model'] == model_i][data_type].isin(selected_data)]).mark_point().encode(
        alt.X('time_point:Q',
              scale=alt.Scale(domain=x_lim, type=x_scale)
              ),
        alt.Y('concentration:Q' if scaled_conc else 'concentration_unscaled:Q',
              scale=alt.Scale(domain=y_lim, type=y_scale)
              ),
        color='met:N',
        tooltip='met:N'
    ).properties(
        width=500,
        height=400
    )

    if data_df is not None:

        plot2 = alt.Chart(data_df[data_df['model'] == model_i]).mark_line().encode(
            alt.X('time_point:Q',
                  scale=alt.Scale(domain=x_lim, type=x_scale)
                  ),
            alt.Y('concentration:Q' if scaled_conc else 'concentration_unscaled:Q',
                  scale=alt.Scale(domain=y_lim, type=y_scale)
                  ),
            color='met:N',
            tooltip='met:N'
        ).properties(
            width=500,
            height=400
        )

        return alt.layer(plot1, plot2).interactive()

    return plot1.interactive()
