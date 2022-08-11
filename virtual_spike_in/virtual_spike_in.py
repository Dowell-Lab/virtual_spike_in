#!/bin/env python3
# virtual_spike_in.py --- 3' Virtual Spike-In Method
#
# Filename: virtual_spike_in.py
# Description: Implementation of the 3' Virtual Spike-In Method
# Author: Zachary Maas <zama8258@colorado.edu>
# Maintainer: Zachary Maas <zama8258@colorado.edu>
# Created: Thu Feb 27 09:46:34 2020 (-0700)
# Version: 0.1.0
#

# Commentary:
#
# This file contains the basic implementation of a new 'virtual
# spike-in' method that uses 3' ends to determine a normalization
# factor between samples in an experiment. The method is based on
# Bayesian inference, and easily scales to a large number of samples.
#
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street, Fifth
# Floor, Boston, MA 02110-1301, USA.
#
#

# Code:

import setuptools
import itertools as it
import warnings
import re
import pickle
from os import path
import time
import logging
import scipy.stats as st
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pymc as pm
import arviz as az
from sklearn.linear_model import LinearRegression
from tqdm import tqdm


def import_count_data(data_csv_filename):
    """Generic wrapper for importing counts tables as a dataframe"""
    data = pd.read_csv(data_csv_filename, sep="\t")
    sample_names = list(data.columns)[7:]
    data = data.rename(columns=lambda x: re.sub(".sorted.bam", "", x))
    data = data.drop(data.columns[list(range(0, 6))], axis=1)
    data = data.astype("float32")
    data = data[sorted(data.columns)]
    # TODO Check to see if counts are all positive (should be...)
    # Filtering for only columns that have some reads
    data["total"] = np.sum(data, axis=1)
    data = data[data["total"] > 0]
    data = data.drop(labels=["total"], axis=1)
    return data


def log2(x):
    """Wrapper to support log2 natively in pymc3"""
    return pm.math.log(x) / pm.math.log(2)


def laplace_inflated_log2_ratio(a, b):
    """Take the log2 ratio of two things with laplace smoothing"""
    return log2((a + 1) / (b + 1))


def lfc_model(
    obs_dat,
    a,
    b,
    diagnostic_prefix="DEBUG",
    plot_model=False,
    num_samples=25_000,
    burn=2_500,
):
    """The definition for our linear log-fold change model"""
    with pm.Model() as model:
        # Calculate data for priors
        # TODO Move NegBinom priors outside here for variance calculation
        a_mu = np.mean(a)
        b_mu = np.mean(b)
        data_var = np.std(np.abs(a - b))
        # Each of our data sets are negative binomial modeled as read counts
        # We use Laplace smoothing to avoid log problems
        a_dist = pm.NegativeBinomial("a_dat", mu=a_mu, alpha=0.5)
        b_dist = pm.NegativeBinomial("b_dat", mu=b_mu, alpha=0.5)
        # We expect variance to be positive, so use a Exponential. Half
        # cauchy also makes sure we won't hit invalid negative values
        # std = pm.Exponential("std_line", lam=1 / data_var)
        std = pm.InverseGamma("std_line", mu=1, sigma=1)
        # Normal or Cauchy works here, but the steepness of Cauchy
        # seems to make the model converge more effectively
        # laplace_ratio =
        mean_std = pm.InverseGamma("std", mu=1, sigma=1)
        mean = pm.Normal(
            "mean",
            mu=a_dist / b_dist,
            sigma=mean_std,
            initval=1,
        )
        intercept = pm.Normal("intercept", 0, sigma=1)
        # Fold changes are log-normally distributed. By using the log
        # of the ratio, we transform it to be normally distributed,
        # making parameter inference simpler.
        # obs = pm.Normal("obs", mu=mean, sigma=std, observed=obs, initval=0)
        obs = pm.Normal(
            "obs",
            mu=mean * a + intercept,
            sigma=std,
            observed=b,
            initval=1,
        )

        # Plot the model if we explicitly request it
        if plot_model:
            pm.model_to_graphviz(model).render("model_graph.gv")

        # Run the actual sampling process in MCMC. The discrete
        # elements of the model preclude the use of variational
        # inference.
        inferencedata = pm.sample(
            draws=num_samples,
            tune=burn,
            init="auto",
            target_accept=0.95,
            # progressbar=False,
        )
        # Generate diagnostic arviz plots
        az.plot_posterior(inferencedata)
        plt.savefig(
            f"{diagnostic_prefix}_posterior.pdf",
            format="pdf",
            orientation="landscape",
        )
        az.plot_autocorr(inferencedata)
        plt.savefig(
            f"{diagnostic_prefix}_autocorr.pdf",
            format="pdf",
            orientation="landscape",
        )
        az.plot_trace(inferencedata)
        plt.savefig(
            f"{diagnostic_prefix}_trace.pdf",
            format="pdf",
            orientation="landscape",
        )
        # az.plot_pair(inferencedata, divergences=True)
        # plt.savefig(
        #     f"{diagnostic_prefix}_divergences.pdf",
        #     format="pdf",
        #     orientation="landscape",
        # )
        plt.close("all")
        return inferencedata


def run_normalization_model(i, j, tag, data, logger):
    """Run the linear log-fold change model on our data"""

    # We want to time each run for diagnostics
    start = time.time()

    # Each run gets a unique design string
    design = f"{i}_vs_{j}"
    logger.debug(f"Running {design}")

    # We express fold-changes in the log transformed form, since
    # log-transformed fold-change values are normally distributed.
    logger.debug(f"Using Laplace Smoothing on {design}")
    data_a = data[i] + 1
    data_b = data[j] + 1
    data_ratio = np.divide(data_a, data_b)
    # np.seterr(divide="ignore")
    data_ratio = np.log2(
        data_ratio,
    )
    logger.debug(f"Number of samples: {len(data_ratio)}")

    # Execute the model
    diagnostic_prefix = f"{tag}_{design}"
    inferencedata = lfc_model(data_ratio, data_a, data_b, diagnostic_prefix)
    post = inferencedata.posterior
    burn = 2500  # Make a global var
    posterior = post.sel(draw=slice(burn, None))

    # Gather results
    # mu = 2 ** np.mean(posterior["mean"])  ## np.mean(trace["mean"]) ** 2
    ## FIXME, is this exponentiated?
    # sigma = 2 ** np.mean(posterior["std"])  ## p.mean(trace["std"]) ** 2
    mean_dat = np.mean(posterior["mean"], axis=0)
    std_dat = np.mean(posterior["std"], axis=0)
    intercept_dat = np.mean(posterior["intercept"], axis=0)
    mu = float(np.mean(mean_dat))
    sigma = float(np.mean(std_dat))
    intercept = float(np.mean(intercept_dat))

    # Generate diagnostic plot compared to data
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, xlabel=i, ylabel=j, title="Data and Model")
    ax.plot(data_a, data_b, ".", label="Observed Data")
    x = np.linspace(min(data_a), max(data_a), 200)
    # Plot the regression line
    reg_line = mu * x + intercept
    ax.plot(x, reg_line, label="Model Fit", lw=1.0, color="red")
    # hdi = az.hdi(inferencedata)
    # lower_line_hdi = float(hdi["mean"][0]) * x + intercept
    # upper_line_hdi = float(hdi["mean"][1]) * x + intercept
    # print(lower_line_hdi, upper_line_hdi)
    # plt.fill_between(x, lower_line_hdi, upper_line_hdi, label="89% CI")
    plt.legend(loc=0)
    plt.savefig(
        f"{diagnostic_prefix}_comparison.pdf",
        format="pdf",
        orientation="landscape",
    )

    # Finish timing
    end = time.time()
    logger.debug(f"{design} took: {end-start:.2f}s")

    # Output results
    return dict(
        design=design,
        mu=mu,
        sigma=sigma,
        intercept=intercept,
        ratio="deprecated",
    )


def run_linear_model(i, j, tag, data, logger):
    """Run a standard linear regression on the data"""

    # Each run gets a unique design string
    design = f"{i}_vs_{j}"
    diagnostic_prefix = f"{tag}_{design}"

    data_a = np.array(data[i]).reshape(-1, 1)
    data_b = np.array(data[j]).reshape(-1, 1)

    reg = LinearRegression().fit(data_a, data_b)
    mu = reg.coef_[0][0]

    # Output results
    return dict(design=design, mu=mu, sigma="NA")


def plot_model_output(ax, result, n=250):
    # Grab Data
    plot_title = result["design"]
    mu = result["mu"]
    sigma = result["sigma"]

    # Generate points
    x = np.linspace(0, 2, n)
    y = st.norm.pdf(x, loc=mu, scale=sigma)

    # Perform plotting
    ax.plot(x, y, color="blue", label="Estimated")
    ax.axvline(0, color="grey", linewidth=0.25)
    ax.set_title(plot_title)
    ax.set_xlabel("Log_2(Normalization Factor)")
    ax.set_ylabel("Relative Proportion")
    ax.legend()


def run_vsi(
    count_data: str,
    label: str,
    outdir: str,
    generate_plots: bool = True,
):
    # Set up our module logger
    logger = logging.getLogger("VSILogger")
    logger.setLevel(logging.DEBUG)
    fh = logging.StreamHandler()
    fh.setLevel(level=logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s %(levelname)s:%(message)s")
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    # Disable most messages for pymc3 to not clutter logs
    debug = False
    if not debug:
        pm_logger = logging.getLogger("pymc3")
        pm_logger.setLevel(logging.ERROR)
        pm_logger.propagate = False
        # Put warnings in our logs too
        logging.captureWarnings(True)
        py_warnings_logger = logging.getLogger("py.warnings")
        py_warnings_logger.addHandler(fh)
        # Set our main logger to show debug messages too
        fh.setLevel(logging.INFO)

    # Perform data import for each organism and merge the data
    seqdata = import_count_data(count_data)
    logger.info(f"Loaded {count_data}")

    # We want every pairwise combination
    sample_names = list(seqdata.columns)

    # We only do combinatorial pairwise combinations if we're working
    # on debugging the algorithm, otherwise a simple one-to-rest
    # comparison is much faster and still gives us size factors.
    do_combinatorial = False
    if do_combinatorial:
        pairwise_combinations = list(it.combinations(sample_names, 2))
    else:
        # TODO Fix
        pairwise_combinations = [(sample_names[0], x) for x in sample_names[1:]]

    # Run each model
    # We want to allow any number of datasets to be used together
    # This just requires a small amount of reparameterization throughout
    diagnostic_prefix = f"{outdir}/{label}"
    for model, model_name in (
        (run_linear_model, "linear"),
        (run_normalization_model, "vsi"),
    ):
        iter_dat = []
        logger.info(f"Performing Inference on {model_name} Models")
        for (i, j) in tqdm(
            pairwise_combinations, desc="Pairwise Comparisons", unit="models"
        ):
            ## We run j vs i to match the standard normalization multiplier
            result = model(j, i, diagnostic_prefix, seqdata, logger)
            iter_dat.append(result)
        # Write size factors out
        size_factor_path = f"{outdir}/{label}_{model_name}_size_factors.txt"
        with open(size_factor_path, "w+") as out_file_handle:
            for result in iter_dat:
                # Gt our parameters for normalization, adjusted back to normal space
                design = result["design"]
                mu = result["mu"]
                sigma = result["sigma"]
                # Write each to file
                out_file_handle.write(f"{design}\t{mu}\t{sigma}\n")

        last_run_path = f"{outdir}/{model_name}_last_run.pickle"
        with open(last_run_path, "wb") as last_run:
            pickle.dump(iter_dat, last_run)
            logger.info(f"Pickled results to {last_run_path}")

    if generate_plots:
        # Do output printing and quantitation
        logger.info("Plotting Model Output")

        # Plot result naively
        ncol = 4
        nrow = int(np.ceil(len(iter_dat) / ncol))
        fig, axs = plt.subplots(
            nrow, ncol, constrained_layout=True, figsize=(18.0, 9.0)
        )
        for idx, result in enumerate(iter_dat):
            plot_model_output(ax=axs.flat[idx], result=result)
        # lst = axs.flat[(ncol * nrow) - 1].axis("off")
        plt.savefig(
            f"{outdir}/{label}_vsi_results.pdf",
            format="pdf",
            orientation="landscape",
        )


#
# virtual_spike_in.py ends here
