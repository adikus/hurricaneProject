package hoos.project.LES.Kernels;

import java.util.Arrays;

import com.amd.aparapi.Range;

public class Halos extends Base {
	private float[] p_halo;
	private float[] uvw_halo;
	private float[] uvwsum_halo;
	private float[] fgh_halo;
	private float[] fgh_old_halo;
	private float[] diu_halo;
	private float[] rhs_halo;
	private float[] sm_halo;
	
	public float[] get_p_halo() {
		this.get(p_halo);
		return p_halo;
	}

	public void set_p_halo(float[] p_halo) {
		this.p_halo = p_halo;
	}

	public float[] get_uvw_halo() {
		return uvw_halo;
	}

	public void set_uvw_halo(float[] uvw_halo) {
		this.uvw_halo = uvw_halo;
	}

	public float[] get_uvwsum_halo() {
		return uvwsum_halo;
	}

	public void set_uvwsum_halo(float[] uvwsum_halo) {
		this.uvwsum_halo = uvwsum_halo;
	}

	public float[] get_fgh_halo() {
		return fgh_halo;
	}

	public void set_fgh_halo(float[] fgh_halo) {
		this.fgh_halo = fgh_halo;
	}

	public float[] get_fgh_old_halo() {
		return fgh_old_halo;
	}

	public void set_fgh_old_halo(float[] fgh_old_halo) {
		this.fgh_old_halo = fgh_old_halo;
	}

	public float[] get_diu_halo() {
		return diu_halo;
	}

	public void set_diu_halo(float[] diu_halo) {
		this.diu_halo = diu_halo;
	}

	public float[] get_rhs_halo() {
		return rhs_halo;
	}

	public void set_rhs_halo(float[] rhs_halo) {
		this.rhs_halo = rhs_halo;
	}

	public float[] get_sm_halo() {
		return sm_halo;
	}

	public void set_sm_halo(float[] sm_halo) {
		this.sm_halo = sm_halo;
	}
	
	@Override 
	public void run() {
		float sum = p2[0] + uvw[0] + uvwsum[0] + fgh[0] + fgh_old[0] + rhs[0] + mask1[0] + diu[0] + sm[0];
		sum += dxs[0] + dys[0] + dzs[0] + dx1[0] + dy1[0] + dzn[0] + z2[0] + cn1[0] + cn2l[0] + cn2s[0] + cn3l[0] + cn3s[0] + cn4l[0] + cn4s[0];
		sum += val_ptr[0] + chunks_num[0] + chunks_denom[0] + n_ptr[0] + state_ptr[0] + pdt[0] + pim[0] + pjm[0] + pkm[0];
		float halo = p_halo[0] + uvw_halo[0] + uvwsum_halo[0] + fgh_halo[0] + fgh_old_halo[0] + diu_halo[0] + rhs_halo[0] + sm_halo[0];
		state_ptr[0] = n_ptr[1];
		chunks_num[0] = sum;
		chunks_denom[0] = sum;
		uvw[0] = sum;
		uvwsum[0] = sum;
		p2[0] = sum;
		p_halo[0] = halo;
		uvw_halo[0] = halo;
		uvwsum_halo[0] = halo;
		fgh_halo[0] = halo;
		fgh_old_halo[0] = halo;
		diu_halo[0] = halo;
		rhs_halo[0] = halo;
		sm_halo[0] = halo;
	}
	
	@Override
	public void init(int ip, int jp, int kp) {
		// Outer Halo size: 2 x v_dim x h_w x (ip+jp + 2h_w) x kp
		// Inner Halo size: 2 x v_dim x h_w x (ip+jp - 2h_w) x kp
		p_halo = new float[2 * 2 * (ip+jp+6 + 2) * (kp+2)];
		uvw_halo = new float[2 * 4 * (ip+jp+6 + 2) * (kp+3)];
		uvwsum_halo = new float[2 * 4 * (ip+jp+2 + 2) * (kp+1)];
		fgh_halo = new float[2 * 4 * (ip+jp+2 + 2) * (kp+1)];
		fgh_old_halo = new float[2 * 4 * (ip+jp + 2) * kp];
		diu_halo = new float[2 * 16 * (ip+jp+7 + 2) * (kp+3)];
		rhs_halo = new float[2 * (ip+jp+4 + 2) * (kp+2)];
		sm_halo = new float[2 * (ip+jp+6 + 2) * (kp+2)];
		
		super.init(ip, jp, kp);
		
		this.readHalos(States.HALO_READ_ALL);
		this.get(p_halo);
		this.get(uvw_halo);
		this.get(uvwsum_halo);
		this.get(fgh_halo);
		this.get(fgh_old_halo);
		this.get(diu_halo);
		this.get(rhs_halo);
		this.get(sm_halo);
	}
	
	public void writeHalos(int state) {
		Range haloRange = Range.create((kp+3) * Math.max(ip+4, jp+3));
		this.executeState(state, haloRange);
	}
	
	public void readHalos(int state) {
		Range haloRange = Range.create((kp+3) * Math.max(ip+4, jp+3));
		this.executeState(state, haloRange);
	}
	
	public float getReductionValue() {
		return val_ptr[0];
	}
	
	public void setReductionValue(float redVal) {
		val_ptr[0] = redVal;
		this.put(val_ptr);
	}
	
	public void run(int state) {		
		switch(state) {
		case States.VELNW__BONDV1_INIT_UVW:
			n_ptr[0] = 1;
			this.put(n_ptr);

			this.writeHalos(States.HALO_WRITE_VELNW__BONDV1_INIT_UVW);
			this.executeState(States.VELNW__BONDV1_INIT_UVW, Range.create((ip+1)*jp*kp));
			this.readHalos(States.HALO_READ_VELNW__BONDV1_INIT_UVW);
			break;
		case States.BONDV1_CALC_UOUT:
			bondv1Reduction();
			break;
		case States.BONDV1_CALC_UVW:
			this.writeHalos(States.HALO_WRITE_BONDV1_CALC_UVW);
			this.executeState(States.BONDV1_CALC_UVW, Range.create((kp*jp)+(kp+2)*(ip+2)+(ip+3)*(jp+3)));
			this.readHalos(States.HALO_READ_BONDV1_CALC_UVW);
			break;
		case States.VELFG__FEEDBF__LES_CALC_SM:
			this.writeHalos(States.HALO_WRITE_VELFG__FEEDBF__LES_CALC_SM);
			this.executeState(States.VELFG__FEEDBF__LES_CALC_SM, Range.create(ip*jp*kp));
			this.readHalos(States.HALO_READ_VELFG__FEEDBF__LES_CALC_SM);
			break;
		case States.LES_BOUND_SM:
			int max_range = Math.max(Math.max(ip+3,jp+3),kp+2);
			this.writeHalos(States.HALO_WRITE_LES_BOUND_SM);
			this.executeState(States.LES_BOUND_SM, Range.create(max_range*max_range, max_range));
			this.readHalos(States.HALO_READ_LES_BOUND_SM);
			break;
		case States.LES_CALC_VISC__ADAM:
			this.writeHalos(States.HALO_WRITE_LES_CALC_VISC__ADAM);
			this.executeState(States.LES_CALC_VISC__ADAM, Range.create(ip*jp*kp));
			this.readHalos(States.HALO_READ_LES_CALC_VISC__ADAM);
			break;
		case States.PRESS_RHSAV:
			rhsavState();
			break;
		case States.PRESS_SOR:
			// rest needs to be called using pressSORIteration and getPressSORValue		
			break;
		case States.PRESS_PAV:
			pavState();
			break;
		case States.PRESS_ADJ:
			this.writeHalos(States.HALO_WRITE_PRESS_ADJ);
			this.executeState(States.PRESS_ADJ, Range.create(ip*jp*kp));
			this.readHalos(States.HALO_READ_PRESS_ADJ);
			break;
		case States.PRESS_BOUNDP:
			max_range = Math.max(Math.max(ip+2,jp+2),kp+2);
			this.writeHalos(States.HALO_WRITE_PRESS_BOUNDP);
			this.executeState(States.PRESS_BOUNDP, Range.create(max_range*max_range, max_range));
			this.readHalos(States.HALO_READ_PRESS_BOUNDP);
			//System.out.println(Arrays.toString(Arrays.copyOfRange(p_halo, 0, 50)));
			break;
		}	
	}
	
	public void dispose() {
		this.dispose();
	}
	
	private void bondv1Reduction() {
		this.writeHalos(States.HALO_WRITE_BONDV1_CALC_UOUT);
		this.executeState(States.BONDV1_CALC_UOUT, Range.create(jp, jp));
		this.readHalos(States.HALO_READ_BONDV1_CALC_UOUT);
		
		this.get(chunks_num);
		this.get(chunks_denom);
		
		float aaa = 0f;
        float bbb = 0f;
		for(int i = 0; i < jp; i++){
			aaa = Math.max(aaa, chunks_num[i]);
			bbb = Math.max(bbb, chunks_denom[i]);
		}
		val_ptr[0] = (aaa + bbb)*0.5f;
		System.out.println("State 2 value: " + val_ptr[0]);
	}
	
	private void rhsavState() {		
		this.writeHalos(States.HALO_WRITE_PRESS_RHSAV);
		this.executeState(States.PRESS_RHSAV, Range.create(ip*kp, ip));
		this.readHalos(States.HALO_READ_PRESS_RHSAV);
		
		this.get(chunks_num);
		this.get(chunks_denom);
		
		float rhsav = 0f;
        float area = 0f;
		for(int i = 0; i < ip; i++){
			rhsav += chunks_num[i];
            area += chunks_denom[i];
		}
		val_ptr[0] = rhsav / area;
		System.out.println("State 7 value: " + val_ptr[0]);
	}
	
	public void pressSORIteration(int i) {	
		int oclGlobalRange;
		int oclLocalRange;
		
		if (i < 2){
			oclGlobalRange = kp*jp;
	        oclLocalRange = jp;
		}else{
			oclGlobalRange = (ip+2)*(jp+2);
	        oclLocalRange = jp+2;
		}
		
		n_ptr[0] = i;
		this.put(n_ptr);
		this.writeHalos(States.HALO_WRITE_PRESS_SOR);
		this.executeState(States.PRESS_SOR, Range.create(oclGlobalRange, oclLocalRange));
		this.readHalos(States.HALO_READ_PRESS_SOR);
	}
	
	public float getPressSORValue() {
		float sor;
		
		this.get(chunks_num);
		sor = 0f;
		for(int j = 0; j < kp; j++){
			sor = sor + chunks_num[j];
		}
		return (float) Math.sqrt(sor);
	}
	
	private void pavState() {		
		this.writeHalos(States.HALO_WRITE_PRESS_PAV);
		this.executeState(States.PRESS_PAV, Range.create(ip*kp, ip));
		this.readHalos(States.HALO_READ_PRESS_PAV);
		
		float pav = 0f;
        float pco = 0f;
		for(int i = 0; i < ip; i++){
			pav += chunks_num[i];
			pco += chunks_denom[i];
		}
		val_ptr[0] = pav / pco;
		System.out.println("State 9 value: " + val_ptr[0]);
	}
}
