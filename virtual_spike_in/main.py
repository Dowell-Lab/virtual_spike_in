# main.py --- Main wrapper for the virtual spike-in code
#
# Filename: main.py
# Author: Zach Maas
# Created: Mon Mar 14 10:18:16 2022 (-0600)
#
#

# Commentary:
#
#
# This file contains a main module to run the virtual-spike-in model
# using a well documented CLI. We use the typer library which
# generates our CLI from the type hints provided in these functions.
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Emacs.  If not, see <https://www.gnu.org/licenses/>.
#
#

# Code:

import click
from virtual_spike_in import run_vsi


@click.command()
@click.option(
    "--label", required=True, help="A prefix label for generated plots"
)
@click.option(
    "--outdir", required=True, help="The directory to output files to"
)
@click.option(
    "--count-data",
    required=True,
    help="Counts table using long gene ends",
)
@click.option("--generate-plots", is_flag=True, help="Whether to plot")
def vsi(count_data: str, label: str, outdir: str, generate_plots: bool):
    """Run the virtual spike-in"""
    print("Running spike-in!")
    run_vsi(
        count_data=count_data,
        label=label,
        outdir=outdir,
        generate_plots=generate_plots,
    )


if __name__ == "__main__":
    vsi()

#
# main.py ends here
