java \
 -Djava.library.path=../../com.amd.aparapi.jni/dist \
 -Dcom.amd.aparapi.enableVerboseJNI=false \
 -Dcom.amd.aparapi.enableShowGeneratedOpenCL=true \
 -Dcom.amd.aparapi.executionMode=%1 \
 -classpath ../../com.amd.aparapi/dist/aparapi.jar:add.jar \
 com.amd.aparapi.sample.add.Main
