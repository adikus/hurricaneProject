package hoos.project.LES;

import com.sun.jna.Native;

public interface InitLib extends com.sun.jna.Library {
    InitLib lib = (InitLib) Native.loadLibrary("/home/andrej/fortran/LES_hoos/RefactoredSources/libles_ocl_gpu.so", InitLib.class);
    
    void init_les_params(
    		float[] p2, float[] uvw, float[] uvwsum, float[] fgh, float[] fgh_old,
    		float[] rhs, float[] mask1, float[] diu, float[] sm,
    		float[] dxs, float[] dys, float[] dzs, float[] dx1, float[] dy1, float[] dzn, float[] z2, 
    		float[] cn1, float[] cn2l, float[] cn2s, float[] cn3l, float[] cn3s, float[] cn4l, float[] cn4s
	);
}
