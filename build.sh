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

# Build project
cd Scripts/ilmn_script
rm -Rf obj
rm -Rf target
make

### shouldn't store class files in svn, change is required
cd ../../Java_SXAppendZygo
jar cfe CGZygosity_SingleSample.jar CGZygosity_SingleSample *.class 
mv CGZygosity_SingleSample.jar ../Scripts/ilmn_script/target
