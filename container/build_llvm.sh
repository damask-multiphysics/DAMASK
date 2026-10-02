#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
set -eu

BUILDCMD=${1:-docker} # or buildah

source ./set_build_vars.sh ${BUILDCMD} # set $VERSION, $REMOTE $PULL and $CREDS

IMAGE=${REMOTE}/llvm_base:${VERSION}
${BUILDCMD} build -f Containerfile.llvm_base -t ${IMAGE} ${PULL} .
${BUILDCMD} push ${CREDS} ${IMAGE}

for VARIANT in default; do
  case ${VARIANT} in
    default)
      export FLAGS="-FOPTFLAGS='-O0 -g' -COPTFLAGS='-O0 -g' -CXXOPTFLAGS='-O0 -g'"
      ;;
  esac

  IMAGE=${REMOTE}/llvm_petsc-${VARIANT}:${VERSION}

  ${BUILDCMD} build -f Containerfile.llvm_petsc -t ${IMAGE} \
                  --build-arg VERSION=${VERSION} \
                  --build-arg VARIANT=${VARIANT} \
                  --build-arg REMOTE=${REMOTE} \
                  --build-arg FLAGS="${FLAGS}" \
                  ${PULL} .

  ${BUILDCMD} push ${CREDS} ${IMAGE}
done
