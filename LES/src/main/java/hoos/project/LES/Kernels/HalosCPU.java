package hoos.project.LES.Kernels;

import java.io.Serializable;

import com.amd.aparapi.Range;
import com.amd.aparapi.device.OpenCLDevice;

/**
 * @author      Andrej Hoos
 * Kernel host class for kernel with halo read and write capability and designed to work on the CPU
 */
public class HalosCPU extends Base implements Serializable {
	private static final long serialVersionUID = 1L;
	
	private static final int NTH = 128;
	private int computeUnits;
	
	private float[] p_halo;
	private float[] uvw_halo;
	private float[] uvwsum_halo;
	private float[] fgh_halo;
	private float[] fgh_old_halo;
	private float[] diu_halo;
	private float[] rhs_halo;
	private float[] sm_halo;
	
	private float nominator, denominator;
	
	transient private Range haloRange;
	
	/*
	 * Halo getter and setter methods
	 */
	
	
	public float[] get_p_halo() {
		this.get(p_halo);
		return p_halo;
	}

	public void set_p_halo(float[] p_halo) {
		this.p_halo = p_halo;
		this.put(p_halo);
	}

	public float[] get_uvw_halo() {
		this.get(uvw_halo);
		return uvw_halo;
	}

	public void set_uvw_halo(float[] uvw_halo) {
		this.uvw_halo = uvw_halo;
		this.put(uvw_halo);
	}

	public float[] get_uvwsum_halo() {
		this.get(uvwsum_halo);
		return uvwsum_halo;
	}

	public void set_uvwsum_halo(float[] uvwsum_halo) {
		this.uvwsum_halo = uvwsum_halo;
		this.put(uvwsum_halo);
	}

	public float[] get_fgh_halo() {
		this.get(fgh_halo);
		return fgh_halo;
	}

	public void set_fgh_halo(float[] fgh_halo) {
		this.fgh_halo = fgh_halo;
		this.put(fgh_halo);
	}

	public float[] get_fgh_old_halo() {
		this.get(fgh_old_halo);
		return fgh_old_halo;
	}

	public void set_fgh_old_halo(float[] fgh_old_halo) {
		this.fgh_old_halo = fgh_old_halo;
		this.put(fgh_old_halo);
	}

	public float[] get_diu_halo() {
		this.get(diu_halo);
		return diu_halo;
	}

	public void set_diu_halo(float[] diu_halo) {
		this.diu_halo = diu_halo;
		this.put(diu_halo);
	}

	public float[] get_rhs_halo() {
		this.get(rhs_halo);
		return rhs_halo;
	}

	public void set_rhs_halo(float[] rhs_halo) {
		this.rhs_halo = rhs_halo;
		this.put(rhs_halo);
	}

	public float[] get_sm_halo() {
		this.get(sm_halo);
		return sm_halo;
	}

	public void set_sm_halo(float[] sm_halo) {
		this.sm_halo = sm_halo;
		this.put(sm_halo);
	}
	
	/**
	 * Dummy Aparapi kernel function to create correct OpenCL method signature
	 */
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
	
	/**
	 * Create all parameter arrays and initialise their values and also create halo buffers
	 * @param ip x size of the domain on a single node 
	 * @param jp y size of the domain on a single node 
	 * @param kp z size of the domain on a single node 
	 */
	@Override
	public void init(int ip, int jp, int kp) {
		System.out.println("Initialising kernel..");
		
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
		
		this.haloRange = Range.create((kp+3) * Math.max(ip+4, jp+3));
		
		this.computeUnits = ((OpenCLDevice) (getExecutionMode().equals(EXECUTION_MODE.GPU) ? OpenCLDevice.best() : (OpenCLDevice)OpenCLDevice.firstCPU())).getMaxComputeUnits();
		
		super.init(ip, jp, kp);
		
		this.executeState(States.HALO_READ_ALL, haloRange);
		this.get(p_halo);
		this.get(uvw_halo);
		this.get(uvwsum_halo);
		this.get(fgh_halo);
		this.get(fgh_old_halo);
		this.get(diu_halo);
		this.get(rhs_halo);
		this.get(sm_halo);

		System.out.println("Finished initialising kernel..");
	}
	
	public float getReductionValue() {
		return val_ptr[0];
	}
	
	public float getReductionNominator() {
		return this.nominator;
	}
	
	public float getReductionDenominator() {
		return this.denominator;
	}
	
	public void setReductionValue(float redVal) {
		val_ptr[0] = redVal;
		System.out.println("Reduction value: " + redVal);
		this.put(val_ptr);
	}
	
	/**
	 * Runs the given step of the simulation
	 * @param state id of the step to be executed
	 */
	public void run(int state) {		
		System.out.println("Kernel running state: " + state);
		long time = System.nanoTime();
		switch(state) {
		case States.VELNW__BONDV1_INIT_UVW:
			n_ptr[0] = 1;
			this.put(n_ptr);

			this.executeState(States.HALO_WRITE_VELNW__BONDV1_INIT_UVW, haloRange);
			this.executeState(States.VELNW__BONDV1_INIT_UVW, Range.create((ip+1)*jp*kp));
			this.executeState(States.HALO_READ_VELNW__BONDV1_INIT_UVW, haloRange);
			break;
		case States.BONDV1_CALC_UOUT:
			bondv1Reduction();
			break;
		case States.BONDV1_CALC_UVW:
			this.executeState(States.HALO_WRITE_BONDV1_CALC_UVW, haloRange);
			this.executeState(States.BONDV1_CALC_UVW, Range.create((kp*jp)+(kp+2)*(ip+2)+(ip+3)*(jp+3)));
			this.executeState(States.HALO_READ_BONDV1_CALC_UVW, haloRange);
			break;
		case States.VELFG__FEEDBF__LES_CALC_SM:
			this.executeState(States.HALO_WRITE_VELFG__FEEDBF__LES_CALC_SM, haloRange);
			this.executeState(States.VELFG__FEEDBF__LES_CALC_SM, Range.create(ip*jp*kp));
			this.executeState(States.HALO_READ_VELFG__FEEDBF__LES_CALC_SM, haloRange);
			break;
		case States.LES_BOUND_SM:
			this.executeState(States.HALO_WRITE_LES_BOUND_SM, haloRange);
			this.executeState(States.LES_BOUND_SM, Range.create((jp+3)*(kp+2) + (kp+2)*(ip+2) + (jp+3)*(ip+2)));
			this.executeState(States.HALO_READ_LES_BOUND_SM, haloRange);
			break;
		case States.LES_CALC_VISC__ADAM:
			this.executeState(States.HALO_WRITE_LES_CALC_VISC__ADAM, haloRange);
			this.executeState(States.LES_CALC_VISC__ADAM, Range.create(ip*jp*kp));
			this.executeState(States.HALO_READ_LES_CALC_VISC__ADAM, haloRange);
			break;
		case States.PRESS_RHSAV:
			rhsavState();
			break;
		case States.PRESS_SOR:
			// needs to be called using pressSORIteration and getPressSORValue		
			break;
		case States.PRESS_PAV:
			pavState();
			break;
		case States.PRESS_ADJ:
			this.executeState(States.HALO_WRITE_PRESS_ADJ, haloRange);
			this.executeState(States.PRESS_ADJ, Range.create(ip*jp*kp));
			this.executeState(States.HALO_READ_PRESS_ADJ, haloRange);
			break;
		case States.PRESS_BOUNDP:
			this.executeState(States.HALO_WRITE_PRESS_BOUNDP, haloRange);
			this.executeState(States.PRESS_BOUNDP, Range.create((jp+2)*(kp+2) + (kp+2)*(ip+2) + (jp+2)*(ip+2)));
			this.executeState(States.HALO_READ_PRESS_BOUNDP, haloRange);
			//System.out.println(Arrays.toString(Arrays.copyOfRange(p_halo, 0, 50)));
			break;
		}	
		System.out.println("Kernel running state: " + state + " took: " + (System.nanoTime() - time)/1000000 + " ms");
	}

	private void bondv1Reduction() {
		this.executeState(States.HALO_WRITE_BONDV1_CALC_UOUT, haloRange);
		this.executeState(States.BONDV1_CALC_UOUT, Range.create(NTH*computeUnits, NTH));
		this.executeState(States.HALO_READ_BONDV1_CALC_UOUT, haloRange);
		
		this.get(chunks_num);
		this.get(chunks_denom);
		
		nominator = 0f;
        denominator = 0f;
		for(int i = 0; i < computeUnits; i++){
			nominator = Math.max(nominator, chunks_num[i]);
			denominator = Math.max(denominator, chunks_denom[i]);
		}
		//val_ptr[0] = (nominator + denominator)*0.5f;
		//System.out.println("State 2 value: " + val_ptr[0]);
	}
	
	private void rhsavState() {		
		this.executeState(States.HALO_WRITE_PRESS_RHSAV, haloRange);
		this.executeState(States.PRESS_RHSAV, Range.create(NTH*computeUnits, NTH));
		this.executeState(States.HALO_READ_PRESS_RHSAV, haloRange);
		
		this.get(chunks_num);
		this.get(chunks_denom);
		
		nominator = 0f;
		denominator = 0f;
		for(int i = 0; i < computeUnits; i++){
			nominator += chunks_num[i];
			denominator += chunks_denom[i];
		}
		//val_ptr[0] = nominator / denominator;
		//System.out.println("State 7 value: " + val_ptr[0]);
	}
	
	/**
	 * Rans a single iteration of the 8th step
	 * @param i the current subiteration
	 */
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
		this.executeState(States.HALO_WRITE_PRESS_SOR, haloRange);
		this.executeState(States.PRESS_SOR, Range.create(oclGlobalRange, oclLocalRange));
		this.executeState(States.HALO_READ_PRESS_SOR, haloRange);
	}
	
	public float getPressSORValue() {		
		this.get(chunks_num);
		float sor = 0f;
		for(int j = 0; j < kp; j++){
			sor += chunks_num[j];
		}
		return sor;
		//return (float) Math.sqrt(sor);
	}
	
	private void pavState() {
		this.executeState(States.HALO_WRITE_PRESS_PAV, haloRange);
		this.executeState(States.PRESS_PAV, Range.create(NTH*computeUnits, NTH));
		this.executeState(States.HALO_READ_PRESS_PAV, haloRange);
		
		nominator = 0f;
		denominator = 0f;
		for(int i = 0; i < computeUnits; i++){
			nominator += chunks_num[i];
			denominator += chunks_denom[i];
		}
		//val_ptr[0] = nominator / denominator;
		//System.out.println("State 9 value: " + val_ptr[0]);
	}
}
