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

import itertools as it
import re
import scipy.stats as st
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pymc3 as pm
import pickle
from os import path
import time

# Constants required for now, eventually these will be parameterized
count_data_human = (
    "/home/zach/dowell_lab/virtual_spike_in/dat/counts_long_ends.txt"
)
count_data_drosophila = (
    "/home/zach/dowell_lab/virtual_spike_in/dat/drosophila_counts_fix.txt"
)
# Samples can be specified in terms of a flat list or tuples of replicates
# We ought to be able to infer these from counts data at some point
samples = [
    ("wt_U2OS_b1", "wt_U2OS_b2"),
    ("dG3BP_U2OS_b1", "dG3BP_U2OS_b2"),
    ("wt_1hAs_U2OS_b1", "wt_1hAs_U2OS_b2"),
    ("dG3BP_1hAs_U2OS_b1", "dG3BP_1hAs_U2OS_b2"),
]
# Static spike-in values, eventually these will be obsolete
spike_in_dict = dict(
    wt_U2OS_b1=85780194,
    wt_U2OS_b2=93689147,
    dG3BP_U2OS_b1=58282107,
    dG3BP_U2OS_b2=71052598,
    wt_1hAs_U2OS_b1=91587100,
    wt_1hAs_U2OS_b2=66604306,
    dG3BP_1hAs_U2OS_b1=56958302,
    dG3BP_1hAs_U2OS_b2=69950273,
)
sample_cols = list(it.chain.from_iterable(samples))


def import_count_data(data_csv_filename):
    """Generic wrapper for importing counts tables as a dataframe"""
    return (
        pd.read_csv(data_csv_filename, sep="\t")
        .rename(columns=lambda x: re.sub(".sorted.bam", "", x))[sample_cols]
        .astype("float32")
    )


# Perform data import for each organism and merge the data
seqdata_human = import_count_data(count_data_human)
seqdata_drosophila = import_count_data(count_data_drosophila)
seqdata_merged = pd.concat([seqdata_human, seqdata_drosophila])


# We want every pairwise combination
pairwise_combinations = list(it.combinations(sample_cols, 2))


def lfc_model(obs, plot_model=False, num_samples=2000, burn=500):
    """The definition for our linear log-fold change model"""
    with pm.Model() as model:
        # Calculate data for priors
        data_var = np.var(obs)
        data_mean = np.mean(obs)
        # We expect variance to be positive, so use a Exponential. Half
        # cauchy also makes sure we won't hit invalid negative values
        std = pm.Exponential("std", lam=1 / data_var)
        # Normal or Cauchy works here, but the steepness of Cauchy
        # seems to make the model converge more effectively
        mean = pm.Laplace("mean", mu=data_mean, b=data_var)
        # Fold changes are log-normally distributed. By using the log
        # of the ratio, we transform it to be normally distributed,
        # making parameter inference simpler.
        obs = pm.Normal("obs", mu=mean, sigma=std, observed=obs)

        # Plot the model if we explicitly request it
        if plot_model:
            pm.model_to_graphviz(model).render("model_graph.gv")

        # Run the actual sampling process
        trace = pm.sample(draws=num_samples, tune=burn, init="auto",)
        return trace


def run_normalization_model(i, j, tag, data):
    """Run the linear log-fold change model on our data"""

    # We want to time each run for diagnostics
    start = time.time()

    # Each run gets a unique design string
    design = f"{i}_vs_{j}"

    # We express fold-changes in the log transformed form, since
    # log-transformed fold-change values are normally distributed.
    data_ratio = np.log2(np.divide(np.add(data[i], 1), np.add(data[j], 1),))
    ratio = np.log2(spike_in_dict[i] / spike_in_dict[j])

    # Execute the model
    trace = lfc_model(data_ratio)

    # Gather results
    mu = np.mean(trace["mean"])
    sigma = np.mean(trace["std"])

    # Finish timing
    end = time.time()
    print(design, "took:", end - start)

    # Output results
    return dict(design=design, mu=mu, sigma=sigma, ratio=ratio)


# Hacky way to run the model right now. Eventually this will need to
# be formalized and cleaned up
run_model = True
last_run_path = "/home/zach/dowell_lab/virtual_spike_in/dat/last_run.pickle"
iter_dat = None
if run_model or not path.exists(last_run_path):
    print("Performing Inference on Human Data")
    results_human = [
        run_normalization_model(i, j, "human", seqdata_human)
        for (i, j) in pairwise_combinations
    ]
    print("Performing Inference on Drosophila Data")
    results_drosophila = [
        run_normalization_model(i, j, "drosophila", seqdata_drosophila)
        for (i, j) in pairwise_combinations
    ]
    print("Performing Inference on Merged Data")
    results_merged = [
        run_normalization_model(i, j, "merged", seqdata_merged)
        for (i, j) in pairwise_combinations
    ]
    # Merge data
    iter_dat = list(
        map(list, zip(*[results_human, results_drosophila, results_merged]))
    )
    # Write result
    with open(last_run_path, "wb") as last_run:
        pickle.dump(iter_dat, last_run)
else:
    with open(last_run_path, "wb") as last_run:
        iter_dat = pickle.read(last_run)


# Do output printing and quantitation
print("Plotting Model Output")


def plot_model_output(ax, result_human, result_drosophila, result_merge, n=100):
    # Grab Data
    plot_title = result_human["design"]
    mu_human = result_human["mu"]
    sigma_human = result_human["sigma"]
    mu_drosophila = result_drosophila["mu"]
    sigma_drosophila = result_drosophila["sigma"]
    mu_merge = result_merge["mu"]
    sigma_merge = result_merge["sigma"]

    # Generate points
    x = np.linspace(-4, 4, n)
    y_human = st.norm.pdf(x, loc=mu_human, scale=sigma_human)
    y_drosophila = st.norm.pdf(x, loc=mu_drosophila, scale=sigma_drosophila)
    y_merge = st.norm.pdf(x, loc=mu_merge, scale=sigma_merge)

    # Compute statistics
    norm_human = np.exp(mu_human)
    norm_drosophila = np.exp(mu_drosophila)
    mean_norm = (norm_human + norm_drosophila) / 2
    dnorm = abs(norm_human - norm_drosophila) / mean_norm
    pdiff = mean_norm * dnorm

    # Perform plotting
    ax.plot(x, y_human, color="blue")
    ax.plot(x, y_drosophila, color="green")
    ax.plot(x, y_merge, color="orange")
    ax.axvline(0, color="grey", linewidth=0.25)
    ax.annotate(f"fc = {pdiff:.3f}", xy=(0.1, 0.8), xycoords="axes fraction")
    ax.set_title(plot_title)
    ax.set_xlabel("Log(Normalization Factor)")
    ax.set_ylabel("Relative Proportion")


# Plot result naively
ncol = 4
nrow = int(np.ceil(len(results_human) / ncol))
fig, axs = plt.subplots(
    nrow, ncol, constrained_layout=True, figsize=(18.0, 36.0)
)
for idx, (result_human, result_drosophila, result_merged) in enumerate(
    iter_dat
):
    plot_model_output(
        axs.flat[idx], result_human, result_drosophila, result_merged
    )
# lst = axs.flat[(ncol * nrow) - 1].axis("off")
plt.savefig(
    "/home/zach/dowell_lab/virtual_spike_in/dat/bayes_spike_in.pdf",
    format="pdf",
    orientation="landscape",
)

# Generate a n x n matrix with more detailed plots
# test_mat = np.zeros((8, 8))
# mat_size = 8
# for col in range(0, mat_size - 1):
#     for row in range(col + 1, mat_size):
#         test_mat[row, col] = -1
#         test_mat[col, row] = 1
#
# x_lin = [i[1]["mu"] for i in iter_dat]
# y_lin = [j[2]["mu"] for j in iter_dat]
# fig, ax = plt.subplots(ncols=1)
# ax.scatter(x_lin, y_lin)
# lin = np.linspace(*ax.get_xlim())
# ax.plot(lin, lin)


#
# virtual_spike_in.py ends here
