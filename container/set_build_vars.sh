# SPDX-License-Identifier: AGPL-3.0-or-later
BUILDCMD="${1:-docker}"

printf -v VERSION '%(%Y.%m.%d)T' -1
REMOTE='git.damask-multiphysics.org:5050/damask/damask'

CREDS=""
if [ "$BUILDCMD" = "buildah" ]; then
  PULL="--pull=newer"
  if [ -n "${REGISTRY_USER:-}" ] && [ -z "${REGISTRY_PW:-}" ]; then
    echo User name \(\$REGISTRY_PW\) for registry not set; exit 1
  fi
  if [ -z "${REGISTRY_USER:-}" ] && [ -n "${REGISTRY_PW:-}" ]; then
    echo Password \(\$REGISTRY_USER\) for registry not set; exit 1
  fi
  if [ -n "${REGISTRY_USER:-}" ] && [ -n "${REGISTRY_PW:-}" ]; then
    CREDS="--creds ${REGISTRY_USER}:${REGISTRY_PW}"
  fi
else
  PULL="--pull"
fi
