#!/bin/bash

GIT_DIR=$(git rev-parse --git-dir)
HOOKS_DIR="${GIT_DIR}/hooks"

TOPLEVEL_DIR=$(git rev-parse --show-toplevel)
SCRIPT_DIR="${TOPLEVEL_DIR}/scripts"

ln -sf "${SCRIPT_DIR}/pre-commit" "${HOOKS_DIR}/pre-commit"
