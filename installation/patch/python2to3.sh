#! /usr/bin/env bash
if [ $1x != 3to2x ]; then
  echo 'python2.7 to python'
  find . -name '*.py' -type f | xargs sed -i 's/usr\/bin\/env python2.7/usr\/bin\/env python/g'
else
  echo 'python to python2.7'
  find . -name '*.py' -type f | xargs sed -i 's/usr\/bin\/env python/usr\/bin\/env python2.7/g'
fi
