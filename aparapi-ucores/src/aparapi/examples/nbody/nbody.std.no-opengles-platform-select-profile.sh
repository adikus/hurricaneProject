
java \
  -Djava.library.path=../../com.amd.aparapi.jni/dist.std:../third-party/jogamp \
  -Dcom.amd.aparapi.executionMode=$1 \
  -Dbodies=$2 \
  -Dcom.amd.aparapi.platformHint=$3 \
  -Dcom.amd.aparapi.enableProfiling=true \
  -Dheight=600 \
  -Dwidth=600 \
  -Djogl.disable.opengles \
  -classpath ../third-party/jogamp/jogl-all.jar:../third-party/jogamp/gluegen-rt.jar:../../com.amd.aparapi/dist/aparapi.jar:nbody.jar \
  com.amd.aparapi.examples.nbody.Main 

