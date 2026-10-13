#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
set -eu

MSC_DIR=${1} # Should be set to the directory containing the Marc20XX.X and Marc20XX.X_documentation folders
MARC_VERSION=${2:-2026.1}
BUILDCMD=${3:-docker} # or buildah

if [ -z "${MSC_DIR}" ]; then echo Location of MSC install files \(\$MSC_DIR\) not set; exit 1;fi

source ./set_build_vars.sh ${BUILDCMD} # set $VERSION, $REMOTE $PULL and $CREDS

${BUILDCMD} build -f Containerfile.marc_base \
                --build-arg MARC_VERSION=${MARC_VERSION} \
                --build-arg MSC_DIR=${MSC_DIR} \
                --build-arg REMOTE=${REMOTE} \
                -t ${REMOTE}/marc_base:${MARC_VERSION} \
                ${PULL} .
${BUILDCMD} push ${CREDS} ${REMOTE}/marc_base:${MARC_VERSION}

${BUILDCMD} build -f Containerfile.marc_patched \
                --build-arg MARC_VERSION=${MARC_VERSION} \
                -t ${REMOTE}/marc_patched:${MARC_VERSION} \
                ${PULL} .
${BUILDCMD} push ${CREDS} ${REMOTE}/marc_patched:${MARC_VERSION}
