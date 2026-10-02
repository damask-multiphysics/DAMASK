#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
set -eu

BUILDCMD=${1:-docker} # or buildah

source ./set_build_vars.sh ${BUILDCMD} # set $VERSION, $REMOTE $PULL and $CREDS

${BUILDCMD} build -f Containerfile.python \
                -t ${REMOTE}/python:${VERSION} \
                ${PULL} .

${BUILDCMD} push ${CREDS} ${REMOTE}/python:${VERSION}
