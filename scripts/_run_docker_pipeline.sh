#!/usr/bin/env bash
# Wrapper to run pipeline with Docker mode (env vars pre-set for safety classifier)
export NEOAG_RUNNER_MODE=docker
export NEOAG_PROFILE=docker
exec bash /mnt/disk_c/data_transfer/users/samba_wb/indev/neo/scripts/run_pipeline.sh "$@"
