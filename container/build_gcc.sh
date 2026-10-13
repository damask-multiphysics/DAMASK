#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
set -eu

BUILDCMD=${1:-docker} # or buildah

source ./set_build_vars.sh ${BUILDCMD} # set $VERSION, $REMOTE $PULL and $CREDS

IMAGE=${REMOTE}/gcc_base:${VERSION}
${BUILDCMD} build -f Containerfile.gcc_base -t ${IMAGE} ${PULL} .
${BUILDCMD} push ${CREDS} ${IMAGE}

for VARIANT in default debug i64; do
  case ${VARIANT} in
    default)
      export FLAGS="-FOPTFLAGS='-O2 -mtune=native' -COPTFLAGS='-O2 -mtune=native' -CXXOPTFLAGS='-O2 -mtune=native'"
      ;;
    debug)
      export FLAGS="-FOPTFLAGS='-Og -g' -COPTFLAGS='-Og -g' -CXXOPTFLAGS='-Og -g'"
      ;;
    i64)
      export FLAGS="-FOPTFLAGS='-O2 -mtune=native' -COPTFLAGS='-O2 -mtune=native' -CXXOPTFLAGS='-O2 -mtune=native'\
                   --with-64-bit-indices=1"
      ;;
  esac

  IMAGE=${REMOTE}/gcc_petsc-${VARIANT}:${VERSION}

  ${BUILDCMD} build -f Containerfile.gcc_petsc -t ${IMAGE} \
                  --build-arg VERSION=${VERSION} \
                  --build-arg VARIANT=${VARIANT} \
                  --build-arg REMOTE=${REMOTE} \
                  --build-arg FLAGS="${FLAGS}" \
                  ${PULL} .

  ${BUILDCMD} push ${CREDS} ${IMAGE}
done
