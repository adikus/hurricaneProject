package hoos.project.LES.Kernels;

import hoos.project.LES.InitLib;

import com.amd.aparapi.Kernel;
import com.amd.aparapi.Range;

public class Single extends Kernel {
	
	private int ip, jp, kp;
	
	private float[] p2;
	private float[] uvw;
	private float[] uvwsum;
	private float[] fgh;
	private float[] fgh_old;
	private float[] rhs;
	private float[] mask1;
	private float[] diu;
	private float[] sm;

	private float[] dxs;
	private float[] dys;
	private float[] dzs;
	private float[] dx1;
	private float[] dy1;
	private float[] dzn;
	private float[] z2;
	
	private float[] cn1;
	private float[] cn2l;
	private float[] cn2s;
	private float[] cn3l;
	private float[] cn3s;
	private float[] cn4l;
	private float[] cn4s;
	
	private float[] val_ptr;
	private float[] chunks_num;
	private float[] chunks_denom;
	private int[] n_ptr;
	private int[] state_ptr;
	private float[] pdt = {0.2f};
	private int[] pim = new int[1];
	private int[] pjm = new int[1];
	private int[] pkm = new int[1];
	
	@Override 
	public void run() {
		float sum = p2[0] + uvw[0] + uvwsum[0] + fgh[0] + fgh_old[0] + rhs[0] + mask1[0] + diu[0] + sm[0];
		sum += dxs[0] + dys[0] + dzs[0] + dx1[0] + dy1[0] + dzn[0] + z2[0] + cn1[0] + cn2l[0] + cn2s[0] + cn3l[0] + cn3s[0] + cn4l[0] + cn4s[0];
		sum += val_ptr[0] + chunks_num[0] + chunks_denom[0] + n_ptr[0] + state_ptr[0] + pdt[0] + pim[0] + pjm[0] + pkm[0];
		state_ptr[0] = n_ptr[1];
		chunks_num[0] = sum;
		chunks_denom[0] = sum;
		uvw[0] = sum;
		uvwsum[0] = sum;
		p2[0] = sum;
	}
	
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
		
		state_ptr[0] = 0;
		this.put(state_ptr);
		this.execute(Range.create(1, 1));
	}
	
	public void run(int state) {		
		switch(state) {
		case 1:
			state_ptr[0] = 1;
			n_ptr[0] = 1;
			this.put(state_ptr);
			this.put(n_ptr);
			this.execute(Range.create((ip+1)*jp*kp));
			break;
		case 2:
			state_ptr[0] = 2;
			this.put(state_ptr);
			this.execute(Range.create(jp, jp));
			break;
		case 3:
			state_ptr[0] = 3;
			this.put(state_ptr);
			this.execute(Range.create((kp*jp)+(kp+2)*(ip+2)+(ip+3)*(jp+3)));
			break;
		case 4:
			state_ptr[0] = 4;
			this.put(state_ptr);
			this.execute(Range.create(ip*jp*kp));
			break;
		case 5:	
			state_ptr[0] = 5;
			this.put(state_ptr);
			int max_range = Math.max(Math.max(ip+3,jp+3),kp+2);
			this.execute(Range.create(max_range*max_range, max_range));
			break;
		case 6:
			state_ptr[0] = 6;
			this.put(state_ptr);
			this.execute(Range.create(ip*jp*kp));
			break;
		case 7:
			state_ptr[0] = 7;
			this.put(state_ptr);
			rhsavState();
			break;
		case 8:
			state_ptr[0] = 8;
			this.put(state_ptr);
			sorState();			
			break;
		case 9:
			state_ptr[0] = 9;
			this.put(state_ptr);
			pavState();
			break;
		case 10:
			state_ptr[0] = 10;
			this.put(state_ptr);
			this.execute(Range.create(ip*jp*kp));
			break;
		case 11:
			state_ptr[0] = 11;
			this.put(state_ptr);
			max_range = Math.max(Math.max(ip+2,jp+2),kp+2);
			this.execute(Range.create(max_range*max_range, max_range));
			break;
		}	
	}
	
	public void dispose() {
		this.dispose();
	}
	
	private void rhsavState() {		
		this.execute(Range.create(ip*kp, ip));
		
		this.get(chunks_num);
		this.get(chunks_denom);
		
		float rhsav = 0f;
        float area = 0f;
		for(int i = 0; i < ip; i++){
			rhsav += chunks_num[i];
            area += chunks_denom[i];
		}
		val_ptr[0] = rhsav / area;
		this.put(val_ptr);
		System.out.println("State 7 value: " + val_ptr[0]);
	}
	
	private void sorState() {		
		float pjuge = 100.0f;//0.0001f;
		int nmaxp = 50;
		float sor = pjuge * 1.1f;
		int iter = 0;
		
		int oclGlobalRange;
		int oclLocalRange;
		int ngroups;
		
		System.out.println("Starting SOR");
		
		while (sor > pjuge && iter < nmaxp){
			iter++;
			
			for(int i = 0; i <= 2; i++){
				if (i < 2){
					oclGlobalRange = kp*jp;
			        oclLocalRange = jp;
			        ngroups = kp;
				}else{
					oclGlobalRange = (ip+2)*(jp+2);
			        oclLocalRange = jp+2;
			        ngroups = ip+2;
				}
				
				n_ptr[0] = i;
				this.put(n_ptr);
				this.execute(Range.create(oclGlobalRange, oclLocalRange));
				
				if(i == 1){
					this.get(chunks_num);
					sor = 0f;
					for(int j = 0; j < ngroups; j++){
						sor = sor + chunks_num[j];
					}
					sor = (float) Math.sqrt(sor);
				}
			}
			
			//System.out.println(sor);
		}
	}
	
	private void pavState() {		
		this.execute(Range.create(ip*kp, ip));
		
		float pav = 0f;
        float pco = 0f;
		for(int i = 0; i < ip; i++){
			pav += chunks_num[i];
			pco += chunks_denom[i];
		}
		val_ptr[0] = pav / pco;
		this.put(val_ptr);
		System.out.println("State 9 value: " + val_ptr[0]);
	}
}
