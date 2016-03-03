#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source $DIR/conf.cfg

#mvn package
export SPARK_JAVA_OPTS="-Xss5m"
$SPARK_PATH/bin/spark-submit --verbose \
  --class hoos.project.LES.spark.SparkHalosDriver \
  --master spark://$IP:7077 --conf spark.executor.memory=4G \
  --driver-memory 2g --conf spark.shuffle.spill=false \
  --jars $APARAPI_JAR_PATH,$DIR/../src/main/resources/jna-4.2.1.jar \
  $DIR/../target/les-map-reduce-0.1.jar $@
