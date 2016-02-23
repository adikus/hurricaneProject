package hoos.project.LES.Kernels;

public class States {
	public final static int INIT = 0;
	public final static int HALO_READ_ALL = 38;
	public final static int HALO_WRITE_ALL = 39;
	
	public final static int VELNW__BONDV1_INIT_UVW = 1;
	public final static int HALO_WRITE_VELNW__BONDV1_INIT_UVW = 13;
	public final static int HALO_READ_VELNW__BONDV1_INIT_UVW = 14;
	
	public final static int BONDV1_CALC_UOUT = 2;
	public final static int HALO_WRITE_BONDV1_CALC_UOUT = 15;
	public final static int HALO_READ_BONDV1_CALC_UOUT = 16;
	
	public final static int BONDV1_CALC_UVW = 3;
	public final static int HALO_WRITE_BONDV1_CALC_UVW = 17;
	public final static int HALO_READ_BONDV1_CALC_UVW = 18;
	
	public final static int VELFG__FEEDBF__LES_CALC_SM = 4;
	public final static int HALO_WRITE_VELFG__FEEDBF__LES_CALC_SM = 19;
	public final static int HALO_READ_VELFG__FEEDBF__LES_CALC_SM = 20;
	
	public final static int LES_BOUND_SM = 5;
	public final static int HALO_WRITE_LES_BOUND_SM = 21;
	public final static int HALO_READ_LES_BOUND_SM = 22;
	
	public final static int LES_CALC_VISC__ADAM = 6;
	public final static int HALO_WRITE_LES_CALC_VISC__ADAM = 23;
	public final static int HALO_READ_LES_CALC_VISC__ADAM = 24;
	
	public final static int PRESS_RHSAV = 7;
	public final static int HALO_WRITE_PRESS_RHSAV = 25;
	public final static int HALO_READ_PRESS_RHSAV = 26;
	
	public final static int PRESS_SOR = 8;
	public final static int HALO_WRITE_PRESS_SOR = 27;
	public final static int HALO_READ_PRESS_SOR = 28;
	
	public final static int PRESS_PAV = 9;
	public final static int HALO_WRITE_PRESS_PAV = 29;
	public final static int HALO_READ_PRESS_PAV = 33;
	
	public final static int PRESS_ADJ = 10;
	public final static int HALO_WRITE_PRESS_ADJ = 34;
	public final static int HALO_READ_PRESS_ADJ = 35;
	
	public final static int PRESS_BOUNDP = 11;
	public final static int HALO_WRITE_PRESS_BOUNDP = 36;
	public final static int HALO_READ_PRESS_BOUNDP = 37;
	
	public final static int BONDV1_CALC_UVW__VELFG__FEEDBF__LES_CALC_SM = 30;
	public final static int VELFG = 31;
	public final static int FEEDBF__LES_CALC_SM = 32;
	
	public final static int DONE = 12;
}
