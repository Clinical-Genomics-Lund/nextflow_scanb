#!/usr/bin/env bash
set -euo pipefail

sudo -E singularity build scanb_rnaseq_`date +%Y-%m-%d`.sif Singularity
