# Level 4 CS project - Accelerating hurricane simulator using Apache Spark

This repository contains modified version of Aparapi-ucores `/aparapi-ucores`, 
the Game of Life proof of concept (in `/dummy`),
and the parallelised simulation code in `/LES`.

## Building Aparapi-ucores (`/aparapi-ucores`)

 *  Install OpenCL drivers for your system. (
[AMD APP SDK](http://developer.amd.com/tools-and-sdks/opencl-zone/amd-accelerated-parallel-processing-app-sdk/) 
used in project)
 *  Make sure environmental variable `LD_LIBRARY_PATH` contains path to OpenCL libraries
 *  Build Aparapi using `ant build` in ` src/aparapi`. Configuration in
`src/aparapi/com.amd.aparapi.jni/build.xml` might be neeeded

## Running Game of life proof of concept (`/dummy`)

 *  Download [Apache Spark](http://spark.apache.org/) (project tested with version 1.4.1)
 *  Build project using `mvn package`
 *  To run the cogroup version:
    ```
    path/to/spark/bin/spark-submit \
      --class hoos.project.dummy.DummyPartitionedJoin \
      --master local[8] \
      target/dummy-map-reduce-0.1.jar data/gol.txt data/out.txt N
    ```
    where first two arguments are the input and ouptut files and N is the number of iterations.
`local[8]` means it will run on 8 CPU threads

## Running LES (`/LES`)

 *  Build Aparapi-ucores (above)
 *  Download Apache Spark
 *  Configure `scripts/conf.cfg`
    ```
    APARAPI_PATH=/path/to/project/aparapi-ucores/src/aparapi
    SPARK_PATH=/path/to/spark/
    IP=192.168.1.64 # IP address of Spark master
    WORK_DIR=/temp/or/scratch
    MGRID_FILE_PATH=/path/to/project/GIS/Tokyo_20mgrid.txt
    
    APARAPI_JNI_PATH=$APARAPI_PATH/com.amd.aparapi.jni/dist
    APARAPI_JAR_PATH=$APARAPI_PATH/com.amd.aparapi/dist/aparapi.jar
    ```

 *  Build project using `mvn package`

 *  Combile OpenCL binaries
    ```
    cl/cl-compile cl/kernel.cl cl/hoos.project.LES.Kernels.Single.abcl
    cl/cl-compile cl/kernel_halos.cl cl/hoos.project.LES.Kernels.Halos.abcl
    ```
If the `cl-compile` binary does not work on your system, you might
need to compile it from [source](https://github.com/adikus/cl-compile).
File type of the binaries should be `.abcl` if using AMD APP, ` .ibcl` if
using Intel SDK and `.nbcl` if using NVidia.

 *  To run version without Spark you can run either of:
    ```
      path-to-project/LES$ bash scripts/single.sh N
      path-to-project/LES$ bash scripts/halos.sh N
    ```
    where halos also runs single node halo exchanges and N is the number of iterations.

 *  For local mode Spark versions, run either of:
    ```
      path-to-project/LES$ bash scripts/spark_cpu_standalone.sh N X Y
      path-to-project/LES$ bash scripts/spark_gpu_standalone.sh N X Y
    ```
    to run wither on CPU or GPU and X and Y are sizes of the node grid. 
Note: if you specify more than 1x1 node grid, these will run sequentially.

 *  To run on cluster:
     * Start Spark master
          ```
          export SPARK_MASTER_IP=192.168.1.64
          sh path/to/spark/sbin/start-master.sh
          ```
        where 192.168.1.64 is replaced with the IP address configured in `conf.cfg`
     *  Start one or more workers (run on the target node)
        ```
        path-to-project/LES$ bash scripts/start_cpu_worker.sh
        path-to-project/LES$ bash scripts/start_gpu_worker.sh
        ```
     *  Submit spark application
          ```
          path-to-project/LES$ bash scripts/spark-submit N X Y
          ```
     *  You can monitor the cluster in Spark WebUI at [http://localhost:8080/](http://localhost:8080/)
  
