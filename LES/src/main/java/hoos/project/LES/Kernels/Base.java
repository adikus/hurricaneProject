package hoos.project.LES.Kernels;

import hoos.project.LES.InitLib;

import com.amd.aparapi.Kernel;
import com.amd.aparapi.Range;

public abstract class Base extends Kernel {
	
	protected int ip, jp, kp;
	
	protected float[] p2;
	protected float[] uvw;
	protected float[] uvwsum;
	protected float[] fgh;
	protected float[] fgh_old;
	protected float[] rhs;
	protected float[] mask1;
	protected float[] diu;
	protected float[] sm;

	protected float[] dxs;
	protected float[] dys;
	protected float[] dzs;
	protected float[] dx1;
	protected float[] dy1;
	protected float[] dzn;
	protected float[] z2;
	
	protected float[] cn1;
	protected float[] cn2l;
	protected float[] cn2s;
	protected float[] cn3l;
	protected float[] cn3s;
	protected float[] cn4l;
	protected float[] cn4s;
	
	protected float[] val_ptr;
	protected float[] chunks_num;
	protected float[] chunks_denom;
	protected int[] n_ptr;
	protected int[] state_ptr;
	protected float[] pdt = {0.2f};
	protected int[] pim = new int[1];
	protected int[] pjm = new int[1];
	protected int[] pkm = new int[1];
	
	public void init(int ip, int jp, int kp) {		
		this.ip = ip;
		this.jp = jp;
		this.kp = kp;
		
		p2 = new float[2 * (ip+3) * (jp+3) * (kp+2)];    // f2, 3D, (0,0,0)
		uvw = new float[4 * (ip+3) * (jp+3) * (kp+3)];   // f4, 3D, (0,-1,-1)
		uvwsum = new float[4 * (ip+1) * (jp+1) * (kp+1)];// f4, 3D, (0,0,0)
		fgh = new float[4 * (ip+1) * (jp+1) * (kp+1)];   // f4, 3D, (0,0,0)
		fgh_old = new float[4 * ip * jp * kp];           // f4, 3D, (1,1,1)
		rhs = new float[(ip+2) * (jp+2) * (kp+2)];       // f1, 3D, (0,0,0)
		mask1 = new float[4 * (ip+3) * (jp+3) * (kp+2)]; // f4, 3D, (-1,0,0)
		diu = new float[16 * (ip+4) * (jp+3) * (kp+3)];  // f16, 3D, (-1,0,0)
		sm = new float[(ip+3) * (jp+3) * (kp+2)];        // f1, 3D, (-1,-1,0)

		dxs = new float[ip+1]; // 1D, 0
		dys = new float[jp+1]; // 1D, 0
		dzs = new float[kp+4]; // 1D, -1
		dx1 = new float[ip+3]; // 1D, -1
		dy1 = new float[jp+2]; // 1D, 0
		dzn = new float[kp+4]; // 1D, -1
		z2 = new float[kp+2];  // 1D, 1
		
		cn1 = new float[ip*jp*kp]; // 3D, (1,1,1)
		cn2l = new float[ip];      // 1D, 1
		cn2s = new float[ip];      // 1D, 1
		cn3l = new float[jp];      // 1D, 1
		cn3s = new float[jp];      // 1D, 1
		cn4l = new float[kp];      // 1D, 1
		cn4s = new float[kp];      // 1D, 1
		
		val_ptr = new float[256];     // 1D, 1
		chunks_num = new float[512];  // 1D, 1
		chunks_denom = new float[512];// 1D, 1
		n_ptr = new int[256];           // 1D, 1
		state_ptr = new int[256];       // 1D, 1
		
		pim[0] = ip;
		pjm[0] = jp;
		pkm[0] = kp;
		
		InitLib.lib.init_les_params(p2,uvw,uvwsum,fgh,fgh_old,rhs,mask1,diu,sm,dxs,dys,dzs,dx1,dy1,dzn,z2,cn1,cn2l,cn2s,cn3l,cn3s,cn4l,cn4s);
		
		this.setExplicit(true);
		
		this.put(p2);
		this.put(uvw);
		this.put(uvwsum);
		this.put(fgh);
		this.put(fgh_old);
		this.put(rhs);
		this.put(sm);
		this.put(diu);
		this.put(mask1);
		
		this.put(dxs);
		this.put(dys);
		this.put(dzs);
		this.put(dx1);
		this.put(dy1);
		this.put(dzn);
		this.put(z2);
		
		this.put(cn1);
		this.put(cn2l);
		this.put(cn2s);
		this.put(cn3l);
		this.put(cn3s);
		this.put(cn4l);
		this.put(cn4s);
		
		this.put(pdt);
		this.put(pim);
		this.put(pjm);
		this.put(pkm);
		
		this.executeState(States.INIT, Range.create(1, 1));
	}
	
	// A bit cheaty way to deal with any numerical arrays
	public void put(Object array){
		switch(array.getClass().getCanonicalName()){
		case "int[]":
			this.put((int[])array);
			break;
		case "long[]":
			this.put((long[])array);
			break;
		case "float[]":
			this.put((float[])array);
			break;
		case "double[]":
			this.put((double[])array);
			break;
		default:
			System.err.println("Unable to put " + array.getClass().getCanonicalName() + " into kernel");
		}
	}
	
	public void get(Object array){
		switch(array.getClass().getCanonicalName()){
		case "int[]":
			this.get((int[])array);
			break;
		case "long[]":
			this.get((long[])array);
			break;
		case "float[]":
			this.get((float[])array);
			break;
		case "double[]":
			this.get((double[])array);
			break;
		default:
			System.err.println("Unable to get " + array.getClass().getCanonicalName() + " from kernel");
		}
	}
	
	public void executeState(int state, Range range) {
		state_ptr[0] = state;
		this.put(state_ptr);
		this.execute(range);
	}
}
