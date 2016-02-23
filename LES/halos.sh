APARAPI_JNI_PATH=../aparapi-ucores/src/aparapi/com.amd.aparapi.jni/dist

java \
    -classpath ../aparapi-ucores/src/aparapi/com.amd.aparapi/dist/aparapi.jar:src/main/resources/jna-4.2.1.jar:bin \
	-Djava.library.path=$APARAPI_JNI_PATH \
	-Dcom.amd.aparapi.flowType=binary \
	-Djna.nosys=true \
	-Dcom.amd.aparapi.executionMode=GPU \
	-Dcom.amd.aparapi.enableProfiling=true \
	-Dcom.amd.aparapi.enableProfilingCSV=true \
	-Dcom.amd.aparapi.profilingFileNameFormatStr=profiling/halos \
	hoos.project.LES.Drivers.HalosDriver
