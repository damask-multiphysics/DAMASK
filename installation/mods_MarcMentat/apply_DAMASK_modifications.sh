#!/usr/bin/env bash

SCRIPTLOCATION="$( cd "$( dirname "$0" )" && pwd )"
DAMASK_ROOT=$SCRIPTLOCATION/../../
# defining set() allows to source the same file for tcsh and bash, with and without space around =
set() {
    export $1$2$3
 }
source $DAMASK_ROOT/CONFIG

if [ "x$MSC_ROOT" != "x" ]; then
 DEFAULT_DIR=$MSC_ROOT
fi
if [ "x$MARC_VERSION" != "x" ]; then
 DEFAULT_VERSION=$MARC_VERSION
fi
if [ "x$DAMASK_BIN" != "x" ]; then
 BIN_DIR=$DAMASK_BIN
fi

while [ ! -d "$SCRIPTLOCATION/$VERSION" ] || [ -z "$VERSION" ]
do
  echo "Input version of MARC/MENTAT installation: [${DEFAULT_VERSION}]"
  read VERSION
  if [ -z "$VERSION" ]; then
    VERSION=${DEFAULT_VERSION}
  fi
  [[ -d "$SCRIPTLOCATION/$VERSION" ]] || echo -e "$VERSION not supported..!\n"
done
echo "MSC version: $VERSION"

while [ ! -d "$INSTALLDIR" ] || [ -z "$INSTALLDIR" ]
do
  echo "Input path of MARC/MENTAT installation: [${DEFAULT_DIR}]"
  read INSTALLDIR
  if [ -z "$INSTALLDIR" ]; then
    INSTALLDIR=${DEFAULT_DIR}
  fi
  [[ -d "$INSTALLDIR" ]] || echo -e "$INSTALLDIR not found..!\n"
done

INSTALLDIR=${INSTALLDIR%/}               # remove trailing slash
echo "MSC installation path: $INSTALLDIR"

DEFAULT_EDITOR='vi'
EDITOR=''
while [ -z "$EDITOR" ]
do
  echo "Input command to invoke your preferred editor: [${DEFAULT_EDITOR}]"
  read EDITOR
  if [ -z "$EDITOR" ]; then
    EDITOR=${DEFAULT_EDITOR}
  fi
done
echo "Editor: $EDITOR"

# tools
echo ''
echo 'adapting Marc tools...'
theDIR=$INSTALLDIR/marc$VERSION/tools
for filename in 'comp_damask' \
                'comp_damask_l' \
                'comp_damask_h' \
                'comp_damask_mp' \
                'comp_damask_lmp' \
                'comp_damask_hmp' \
                'run_damask' \
                'run_damask_l' \
                'run_damask_h' \
                'run_damask_mp' \
                'run_damask_lmp' \
                'run_damask_hmp' \
                'include_linux64'; do
  cp $SCRIPTLOCATION/$VERSION/Marc_tools/$filename $theDIR
  echo $theDIR/$filename | xargs perl -pi -e "s:%INSTALLDIR%:${INSTALLDIR}:g"
  echo $theDIR/$filename | xargs perl -pi -e "s:%VERSION%:${VERSION}:g"
  echo $filename
done

# Mentat scripts
echo ''
echo 'adapting Mentat scripts...'
theDIR=$INSTALLDIR/mentat$VERSION/bin
for filename in 'edit_window' \
                'submit4' \
                'submit5' \
                'submit6' \
                'submit7' \
                'submit8' \
                'submit9' \
                'kill4' \
                'kill5' \
                'kill6' \
                'kill7' \
                'kill8' \
                'kill9'; do
  cp $SCRIPTLOCATION/$VERSION/Mentat_bin/$filename $theDIR
  echo $theDIR/$filename | xargs perl -pi -e "s:%INSTALLDIR%:${INSTALLDIR}:g"
  echo $theDIR/$filename | xargs perl -pi -e "s:%VERSION%:${VERSION}:g"
  echo $theDIR/$filename | xargs perl -pi -e "s:%EDITOR%:${EDITOR}:g"
  echo $filename
done

# Mentat scripts
echo -e '\nadapting Mentat menus...'
theDIR=$INSTALLDIR/mentat$VERSION/menus
for filename in 'job_run.ms'; do
  cp $SCRIPTLOCATION/$VERSION/Mentat_menus/$filename $theDIR
  echo $theDIR/$filename | xargs perl -pi -e "s:%INSTALLDIR%:${INSTALLDIR}:g"
  echo $theDIR/$filename | xargs perl -pi -e "s:%VERSION%:${VERSION}:g"
  echo $filename
done

# compile menus
echo ''
echo 'compiling Mentat menu binaries...'
$(which xvfb-run 2>/dev/null) $INSTALLDIR/mentat$VERSION/bin/mentat -compile $INSTALLDIR/mentat$VERSION/menus/linux64/main.msb
[[ $? != 0 ]] && echo '...failed. Try installing xvfb-run on your system.'

# setting access rights
echo ''
echo 'setting file access rights...'
for filename in marc$VERSION/tools/run_damask* \
                marc$VERSION/tools/comp_damask* \
                mentat$VERSION/bin/submit{4..9} \
                mentat$VERSION/bin/kill{4..9} \
                
  chmod 755 $INSTALLDIR/${filename}
done

#creating symlinks for run_damask_scripts

if [ -d "$BIN_DIR" ]; then
  echo ''
  echo "Do you want to create symlinks for run_damask scripts in ${BIN_DIR} [YES/no] ?"
  read YESNO
    if [ -z "$YESNO" ]; then
      YESNO=yes
    fi
  case $YESNO in
   y* | Y* )
    echo''
    echo 'creating symlinks ...'
    echo''
    theDIR=$INSTALLDIR/marc$VERSION/tools
    for filename in 'run_damask' \
                    'run_damask_l' \
                    'run_damask_h' \
                    'run_damask_mp' \
                    'run_damask_lmp' \
                    'run_damask_hmp'; do
      echo ${filename:4}$VERSION
      [ -f $BIN_DIR/${filename:4}$VERSION ] && rm $BIN_DIR/${filename:4}$VERSION
      ln -s $theDIR/$filename $BIN_DIR/${filename:4}$VERSION 
    done
    ;;
  esac
fi

# precompiling user subroutine
echo ''
echo 'precompiling $VERSION HYPELA2 user subroutine...'
echo 'not yet implemented..!'

echo -e '\ndone.'
