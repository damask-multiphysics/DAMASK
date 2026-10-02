#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
set -eu

BUILDCMD=${1:-docker} # or buildah
source ./set_build_vars.sh ${BUILDCMD} # set $VERSION, $REMOTE $PULL and $CREDS
COMPILER=${2}  # gcc/oneapi/llvm
VERSION=${3}  # overwrite to user-defined patch version

if [[ ! "$VERSION" =~ [[:lower:]]$ ]]; then
    echo "The last letter is NOT lowercase."
fi

for VARIANT in default debug i64; do
  IMAGE=${REMOTE}/${COMPILER}_petsc-${VARIANT}:${VERSION}

  ${BUILDCMD} build -f Containerfile.${COMPILER}_petsc-${VERSION} -t ${IMAGE} \
                  --build-arg VARIANT=${VARIANT} \
                  --build-arg REMOTE=${REMOTE} \
                  ${PULL} .

  ${BUILDCMD} push ${CREDS} ${IMAGE}
done
