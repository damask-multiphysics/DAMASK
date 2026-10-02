#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
set -eu

BUILDCMD=${1:-docker} # or buildah

source ./set_build_vars.sh ${BUILDCMD} # set $VERSION, $REMOTE $PULL and $CREDS

IMAGE=${REMOTE}/oneapi_base:${VERSION}
${BUILDCMD} build -f Containerfile.oneapi_base -t ${IMAGE} ${PULL} .
${BUILDCMD} push ${CREDS} ${IMAGE}

for VARIANT in default debug i64; do
  case ${VARIANT} in
    default)
      export FLAGS="-FOPTFLAGS='-O2' -COPTFLAGS='-O2' -CXXOPTFLAGS='-O2'"
      ;;
    debug)
      export FLAGS="-FOPTFLAGS='-O0 -g' -COPTFLAGS='-O0 -g' -CXXOPTFLAGS='-O0 -g'"
      ;;
    i64)
      export FLAGS="-FOPTFLAGS='-O2' -COPTFLAGS='-O2' -CXXOPTFLAGS='-O2' --with-64-bit-indices=1"
      ;;
  esac

  IMAGE=${REMOTE}/oneapi_petsc-${VARIANT}:${VERSION}

  ${BUILDCMD} build -f Containerfile.oneapi_petsc -t ${IMAGE} \
                  --build-arg VERSION=${VERSION} \
                  --build-arg VARIANT=${VARIANT} \
                  --build-arg REMOTE=${REMOTE} \
                  --build-arg FLAGS="${FLAGS}" \
                  ${PULL} .

  ${BUILDCMD} push ${CREDS} ${IMAGE}
done
