#!/bin/bash

set -o errexit
set -o nounset

BWA_REPO="https://github.com/lh3/bwa.git"
BWA_RELEASE="0.7.10"

if [ $# -eq 0 ]; then
    BWA_DIR="$(mktemp --dry-run bwa-build-XXXXXXXXXXXXXXXXX)"
else
    BWA_DIR="${1}"
fi
echo "Using BWA build dir $BWA_DIR"

git clone "${BWA_REPO}" "${BWA_DIR}"

cd "${BWA_DIR}"

git checkout "${BWA_RELEASE}"

# check whether we're running on Linux
if type -P uname >/dev/null && [ $(uname) == 'Linux' ]; then
    # if so, add '-fPIC' to CFLAGS
    sed -i -e 's/CFLAGS\s*=.*/& -fPIC/' Makefile
fi

# build libbwa.a
make clean
make libbwa.a

exit 0
