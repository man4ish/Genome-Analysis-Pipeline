#!/bin/sh

if [ -z "$SYNABASEROOT" ] ; then
  if [ -d "/opt/synamatix/synabase/synabase" ]; then
    export SYNABASEROOT=/opt/synamatix/synabase/synabase
    export LD_LIBRARY_PATH=$SYNABASEROOT/lib:$LD_LIBRARY_PATH
    export PATH=$SYNABASEROOT/../:$PATH
  else
    echo "The SYNABASEROOT environment variable is not defined."
    return 1
  fi
fi
if [ -z "$SXDBConfig" ] ; then
  if [ -f "/data/dbref/synabase/sxdbconfig.xml" ]; then
    export SXDBConfig=/data/dbref/synabase/sxdbconfig.xml
  else
    echo "The SXDBConfig environment variable is not defined."
    return 1
  fi
fi
if [ -z "$SYNASEARCH" ] ; then
  if [ -f "/opt/synamatix/synasearch/synasearch/bin/synasearch.sh" ]; then
    export SYNASEARCH=/opt/synamatix/synasearch/synasearch/bin/synasearch.sh
  else
    echo "The SYNASEARCH environment variable is not defined."
    return 1
  fi
fi
if [ -z "$SXDB" ] ; then
  echo "The SXDB environment variable is not defined."
  return 1
fi