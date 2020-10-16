#!/usr/bin/env bash
set -euo pipefail

sudo -E singularity build rnaseq_`date +%Y-%m-%d`.sif Singularity
