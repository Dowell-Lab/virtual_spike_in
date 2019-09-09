#!/bin/bash
# Sync to FIJI
rsync -Pzae "ssh" /home/zach/dowell_lab/virtual_spike_in \
			fiji:/scratch/Users/zama8258/
