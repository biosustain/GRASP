import altair as alt
import pandas as pd


def plot_dG_and_flux(dG_df, flux_df, height=400, width=600):
    threshold_df = pd.DataFrame([{"ThresholdValue": 0, "Threshold": "hazardous"}])

    rect_dG = alt.Chart(dG_df).mark_rect().encode(
        y='rxn:N',
        x='∆G_min',
        x2='∆G_max'
    ).properties(
        height=height,
        width=width
    )

    point_dG = alt.Chart(dG_df).mark_point(size=100, color='red', filled=True).encode(
        y='rxn:N',
        x='∆G_mean',
        tooltip=['rxn', '∆G_min', '∆G_mean', '∆G_max']
    ).interactive()

    rule_dG = alt.Chart(threshold_df).mark_rule().encode(
        x='ThresholdValue:Q'
    )

    rect_flux = alt.Chart(flux_df).mark_rect().encode(
        y='rxn:N',
        x='flux_min:Q',
        x2='flux_max:Q'
    ).properties(
        height=height,
        width=width
    )

    point_flux = alt.Chart(flux_df).mark_point(size=100, color='red', filled=True).encode(
        y='rxn:N',
        x='flux',
        tooltip=['rxn', 'flux_min', 'flux', 'flux_max']
    ).interactive()

    rule_flux = alt.Chart(threshold_df).mark_rule().encode(
        x='ThresholdValue:Q'
    )

    return alt.hconcat(rect_dG + point_dG + rule_dG, rect_flux + point_flux + rule_flux)
