import pymc as pm
import xarray as xr
import pandas as pd
import numpy as np


def load_realization(ds: xr.Dataset, realization: int):
    data = pd.DataFrame(
        ds["data"][:, :-1, realization], columns=ds.attrs["analytes"][:-1]
    )
    for col in ("control", "round"):
        data.insert(0, col, ds[col])
    data["control"] = data["control"].replace({1: "control", 0: "treatment"})
    data = data.set_index(["control", "round"]).sort_index()
    return data


def init_model(df):

    with pm.Model() as model:

        # ----------------------------------
        # deployment information

        wet_feedstock_mass = pm.Normal("wet feedstock mass", mu=12800, sigma=100)
        treatment_area = pm.Normal("treatment area", mu=3200, sigma=32)

        sample_depth = pm.Gamma("sample depth", mu=0.1, sigma=0.025)
        feedstock_moisture_fraction = pm.Normal(
            "feedstock moisture fraction", mu=0.125, sigma=0.025
        )
        dry_feedstock_mass = pm.Deterministic(
            "dry feedstock mass", (1 - feedstock_moisture_fraction) * wet_feedstock_mass
        )
        application_rate = pm.Deterministic(
            "application rate", dry_feedstock_mass / treatment_area
        )
        feedstock_concentration = pm.Normal(
            "feedstock concentration", mu=[0.07, 0.05], sigma=[0.0035, 0.0025]
        )

        # ----------------------------------
        # mixing

        soil_density = pm.Normal("soil density", mu=1e3, sigma=100)
        soil_mass = pm.Deterministic("soil mass", soil_density * sample_depth)

        mixing_fraction = pm.Deterministic(
            "mixing fraction", application_rate / (application_rate + soil_mass)
        )

        mixed_concentration = pm.Deterministic(
            "mixed concentration",
            mixing_fraction * feedstock_concentration
            + (1 - mixing_fraction) * df.loc["treatment", 1].values,
        )

        # ----------------------------------
        # observed post-spreading enrichment

        enrichment = pm.Deterministic(
            "enrichment", mixed_concentration - df.loc["treatment", 1].values
        )

        enrichment_sigma = pm.Exponential("enrichment sigma", scale=1e-3, shape=(1, 2))

        pm.Normal(
            "obs enriched",
            mu=enrichment,
            sigma=enrichment_sigma,
            observed=df.loc["treatment", 2].values - df.loc["treatment", 1].values,
        )

        # ----------------------------------
        # change in control concentrations

        control_change = pm.Normal("control change", mu=0, sigma=1e-3, shape=(1, 2))
        control_change_sigma = pm.HalfNormal(
            "control change sigma", sigma=1e-3, shape=(1, 2)
        )

        pm.Normal(
            "obs control",
            mu=control_change,
            sigma=control_change_sigma,
            observed=df.loc["control", 3]
            - (df.loc["control", 1].values / 2 + df.loc["control", 2].values / 2),
        )

        # ----------------------------------
        # observed concentration loss

        norm_loss_mu = pm.Uniform("norm loss mu")
        norm_loss_sigma = pm.Beta("norm loss sigma", alpha=1, beta=6)

        norm_loss = pm.Normal("norm loss", mu=norm_loss_mu, sigma=norm_loss_sigma, shape=(1, 2))
        concentration_loss = pm.Deterministic(
            "concentration loss", norm_loss * enrichment
        )
        weathered_concentration = pm.Deterministic(
            "weathered concentration",
            df.loc["treatment", 2].values - (concentration_loss - control_change),
        )

        weathered_sigma = pm.Exponential("weathered sigma", scale=1e-3, shape=(1, 2))

        pm.Normal(
            "obs weathered",
            mu=weathered_concentration,
            sigma=weathered_sigma,
            observed=df.loc["treatment", 3].values,
        )

        # ----------------------------------
        # CDR

        cdrpot = pm.Deterministic(
            "CDR potential (per mass)",
            pm.math.sum(feedstock_concentration * np.array([2.196, 3.621])),
        )
        cdrpot_area = pm.Deterministic(
            "CDR potential (per area)", application_rate * cdrpot
        )
        cdr = pm.Deterministic(
            "CDR (per area)",
            pm.math.sum(norm_loss * feedstock_concentration * np.array([2.196, 3.621]))
            * application_rate,
        )
        pm.Deterministic("CDR [metric ton CO2]", cdr * treatment_area / 1e3)
        pm.Deterministic("CDR completion [-]", cdr / cdrpot_area)

    return model
