#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/conf.cfg

mvn package
export MGRID_PATH=$MGRID_FILE_PATH
export APARAPI_CL_BIN_FOLDER=$DIR/../cl
export LES_LIB_PATH=$DIR/../src/main/resources/libles_ocl_gpu.so
java \
    -classpath $APARAPI_JAR_PATH:$DIR/../src/main/resources/jna-4.2.1.jar:$DIR/../target/les-map-reduce-0.1.jar \
	-Djava.library.path=$APARAPI_JNI_PATH \
	-Dcom.amd.aparapi.flowType=binary \
	-Djna.nosys=true \
	-Dcom.amd.aparapi.executionMode=GPU \
	-Dcom.amd.aparapi.enableProfiling=true \
	-Dcom.amd.aparapi.enableProfilingCSV=true \
	-Dcom.amd.aparapi.profilingFileNameFormatStr=profiling/halos \
	hoos.project.LES.Drivers.HalosDriver
