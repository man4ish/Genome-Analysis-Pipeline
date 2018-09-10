#!/bin/sh

CALL_DIR=$(pwd)
cd $1
jar cfe CGZygosity_SingleSample.jar CGZygosity_SingleSample *.class 
eval mv CGZygosity_SingleSample.jar $CALL_DIR"/target/CGZygosity_SingleSample.jar"
