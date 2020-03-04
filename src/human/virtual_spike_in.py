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

count_data = "/home/zach/dowell_lab/virtual_spike_in/dat/counts_long_ends.txt"
# Samples can be specified in terms of a flat list or tuples of replicates
samples = [
    ("wt_U2OS_b1", "wt_U2OS_b2"),
    ("dG3BP_U2OS_b1", "dG3BP_U2OS_b2"),
    ("wt_1hAs_U2OS_b1", "wt_1hAs_U2OS_b2"),
    ("dG3BP_1hAs_U2OS_b1", "dG3BP_1hAs_U2OS_b2"),
]
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

# Flatten rows to subset our data
sample_cols = list(it.chain.from_iterable(samples))
seqdata = pd.read_csv(count_data, sep="\t").rename(
    columns=lambda x: re.sub(".sorted.bam", "", x)
)


def tpm_normalize_col(col, data):
    "TPM normalize a column of our counts data"
    rpk = np.divide(col, data["Length"])
    scale_factor = np.divide(np.sum(rpk), 1e6)
    tpm = np.divide(rpk, scale_factor)
    return tpm


tpm_seqdata = seqdata.apply(
    lambda x: tpm_normalize_col(x, seqdata) if x.name in sample_cols else x
)

# We want every pairwise combination
pairwise_combinations = [(sample_cols[0], j) for j in sample_cols[1:]]


def run_glm(i, j, data):
    "Run a GLM on our data"
    design = i + " ~  " + j
    # We express fold-changes in the log transformed form, since
    # log-transformed fold-change values are normally distributed.
    data_ratio = (
        np.log(np.divide(data[i], data[j]))
        .replace([np.inf, -np.inf], np.nan)
        .dropna()
    )
    data_var = np.var(data_ratio)
    data_mean = np.mean(data_ratio)
    ratio = np.log(spike_in_dict[i] / spike_in_dict[j])

    # Execute the model
    with pm.Model() as model:
        num_samples = 5000
        burn = 100
        # We expect variance to be positive, so use a Exponential. Half
        # cauchy also makes sure we won't hit invalid negative values
        std = pm.Exponential("std", lam=data_var ** 2)
        # Normal or Cauchy works here, but the steepness of Cauchy
        # seems to make the model converge more effectively
        mean = pm.Normal("mean", mu=data_mean, sigma=data_var ** 2)
        # Fold changes are log-normally distributed. By using the log
        # of the ratio, we transform it to be normally distributed,
        # making parameter inference simpler.
        obs = pm.Normal("obs", mu=mean, sigma=std, observed=data_ratio)
        # Plot our model first
        pm.model_to_graphviz(model).render("model_graph.gv")
        # Then, run sampling
        trace = pm.sample(num_samples, init="advi+adapt_diag")
        trace_post_burn = trace[burn:]
        mu = np.mean(trace_post_burn["mean"])
        sigma = np.mean(trace_post_burn["std"])
        return dict(design=design, mu=mu, sigma=sigma, ratio=ratio)


results = [run_glm(i, j, seqdata) for (i, j) in pairwise_combinations]


def plot_model_output(ax, result_dict, n=100):
    plot_title = result_dict["design"]
    mu = result_dict["mu"]
    sigma = result_dict["sigma"]
    ratio = result_dict["ratio"]
    # x = np.linspace(mu - 3 * sigma, mu + 3 * sigma, n)
    x = np.linspace(-2, 2, n)
    y = st.norm.pdf(x, mu, sigma)
    z = st.norm.ppf(st.norm.cdf(ratio, loc=mu, scale=sigma))
    z_color = "black" if np.abs(z) > 1 else "red"
    ax.plot(x, y)
    ax.axvline(0, color="grey", linewidth=0.25)
    ax.axvline(ratio, color=z_color)
    ax.annotate(f"z = {z:.3f}", xy=(0.1, 0.8), xycoords="axes fraction")
    ax.set_title(plot_title)
    ax.set_xlabel("Log(Normalization Factor)")
    ax.set_ylabel("Relative Proportion")


ncol = 4
nrow = int(np.ceil(len(results) / ncol))
fig, axs = plt.subplots(
    nrow, ncol, constrained_layout=True, figsize=(16.0, 9.0)
)
for idx, result in enumerate(results):
    plot_model_output(axs.flat[idx], result)
lst = axs.flat[(ncol * nrow) - 1].axis("off")
plt.savefig(
    "/home/zach/dowell_lab/virtual_spike_in/dat/bayes_spike_in.pdf",
    format="pdf",
    orientation="landscape",
)

#
# virtual_spike_in.py ends here
