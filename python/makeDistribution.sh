#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $DIR

VERSION=$(head -n1 damask/VERSION)
read -p "updated version [${VERSION}]: " UP
echo ${UP:-$VERSION} > damask/VERSION

rm -rf dist/
python3 setup.py sdist bdist_wheel
python3 -m twine upload dist/*
