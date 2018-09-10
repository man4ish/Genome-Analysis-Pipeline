#!/bin/bash

# resolve links
PRG="$0"
while [ -h "$PRG" ] ; do
  ls=`ls -ld "$PRG"`
  link=`expr "$ls" : '.*-> \(.*\)$'`
  if expr "$link" : '.*/.*' > /dev/null; then
    PRG="$link"
  else
    PRG=`dirname "$PRG"`/"$link"
  fi
done

PRGDIR=`dirname "$PRG"`


# Build
./build.sh
if [ "$?" = "0" ]; then
    echo "Build completed."
else
    echo "Build failed!" 1>&2
    exit 1
fi


# Get version number, os and platform
MAJOR_VERSION=$(grep MAJOR_VERSION VERSION | cut -d= -f2)
MINOR_VERSION=$(grep MINOR_VERSION VERSION | cut -d= -f2)
PATCH_VERSION=$(grep PATCH_VERSION VERSION | cut -d= -f2)
IS_DEV_VERSION=$(grep IS_DEV_VERSION VERSION | cut -d= -f2)
VER_NUM=$MAJOR_VERSION.$MINOR_VERSION.$PATCH_VERSION
if [ $IS_DEV_VERSION -eq 1 ]; then
  VER_NUM=$VER_NUM-SNAPSHOT
fi

PLATFORM=`uname -i`
OS='unknown'
unamestr=`uname`
if [[ "$unamestr" == 'Linux' ]]; then
    OS='linux'
fi


# create package
PKG_NAME=`basename $PWD`-$VER_NUM
rm -Rf dist
mkdir -p dist/$PKG_NAME/bin
mkdir -p dist/$PKG_NAME/conf
cp bin/* dist/$PKG_NAME/bin
cp Scripts/ilmn_script/target/* dist/$PKG_NAME/bin
cp conf/* dist/$PKG_NAME/conf
cd dist
tar zcvf $PKG_NAME-$OS-$PLATFORM.tar.gz $PKG_NAME
rm -Rf $PKG_NAME