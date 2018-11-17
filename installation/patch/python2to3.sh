#! /usr/bin/env bash
if [ $1x != 3to2x ]; then
  echo 'python2.7 to python3'
  find . -name '*.py' -type f | xargs sed -i 's/usr\/bin\/env python2.7/usr\/bin\/env python3/g'
else
  echo 'python3 to python2.7'
  find . -name '*.py' -type f | xargs sed -i 's/usr\/bin\/env python3/usr\/bin\/env python2.7/g'
fi
