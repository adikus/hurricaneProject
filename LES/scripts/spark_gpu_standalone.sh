#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/conf.cfg

#mvn package
export MGRID_PATH=$MGRID_FILE_PATH
export APARAPI_CL_BIN_FOLDER=$DIR/../cl
export SPARK_WORKER_DIR=$WORK_DIR
export LES_LIB_PATH=$DIR/../src/main/resources/libles_ocl.so
export SPARK_JAVA_OPTS="-Xss5m -Djava.library.path=$APARAPI_JNI_PATH -Dcom.amd.aparapi.flowType=binary -Djna.nosys=true -Dcom.amd.aparapi.executionMode=GPU -Dcom.amd.aparapi.enableProfiling=true -Dcom.amd.aparapi.enableProfilingCSV=true -Dcom.amd.aparapi.profilingFileNameFormatStr=profiling/halos" 
$SPARK_PATH/bin/spark-submit --verbose \
  --class hoos.project.LES.spark.SparkHalosDriver \
  --master local[1] --conf spark.executor.memory=4G --conf spark.driver.memory=4G \
  --conf spark.local.dir=$WORK_DIR/temp \
  --jars $APARAPI_JAR_PATH,$DIR/../src/main/resources/jna-4.2.1.jar \
  $DIR/../target/les-map-reduce-0.1.jar $@
