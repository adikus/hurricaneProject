#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/conf.cfg

export MGRID_PATH=$MGRID_FILE_PATH
export APARAPI_CL_BIN_FOLDER=$DIR/../cl
export SPARK_WORKER_DIR=$WORK_DIR
export LES_LIB_PATH=$DIR/../src/main/resources/libles_ocl_gpu.so
export SPARK_JAVA_OPTS="-Djava.library.path=$APARAPI_JNI_PATH -Dcom.amd.aparapi.flowType=binary -Djna.nosys=true -Dcom.amd.aparapi.executionMode=GPU" 
$SPARK_PATH/bin/spark-class org.apache.spark.deploy.worker.Worker -m 2G --cores 1 spark://$IP:7077
