package hoos.project.LES.Kernels;

import com.amd.aparapi.Range;

public class Single extends Base {
	
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
	
	public void run(int state) {
		System.out.println("Kernel running state: " + state);	
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
			
			System.out.println(iter + " - " + sor);
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
