package hoos.project.LES.Drivers;

import hoos.project.LES.HaloExchange.HaloConstructor;
import hoos.project.LES.Kernels.HalosCPU;
import hoos.project.LES.Kernels.States;

public class HalosDriver{

	public static void main(String[] _args) {
		HalosCPU kernel = new HalosCPU();
		
		int ip = 150;
		int jp = 150;
		int kp = 90;
		
		kernel.init(ip, jp, kp);
		
		long previousTime = kernel.getAccumulatedExecutionTime();
		System.out.println("First execution time: " + previousTime);
		System.out.println("Conversion execution time: " + kernel.getConversionTime());

		System.out.println("Running kernel..");
		
		int iter = 0;
		while(iter < 5){			
			// 1
			kernel.run(States.VELNW__BONDV1_INIT_UVW);
			float [] p_halo = kernel.get_p_halo();
			float [] uvw_halo = kernel.get_uvw_halo();
			float [] uvwsum_halo = kernel.get_uvwsum_halo();
			
			p_halo = exchange(ip+3, jp+3, kp+2, 2, p_halo);
			uvw_halo = exchange(ip+3, jp+3, kp+3, 4, uvw_halo);
			uvwsum_halo = exchange(ip+1, jp+1, kp+1, 4, uvwsum_halo);
			
			kernel.set_p_halo(p_halo);
			kernel.set_uvw_halo(uvw_halo);
			kernel.set_uvwsum_halo(uvwsum_halo);
			
			// 2
			kernel.run(States.BONDV1_CALC_UOUT);
			uvw_halo = kernel.get_uvw_halo();
			float redVal = (kernel.getReductionNominator() + kernel.getReductionDenominator())*0.5f;
			
			// reduction
			uvw_halo = exchange(ip+3, jp+3, kp+3, 4, uvw_halo);
			
			kernel.set_uvw_halo(uvw_halo);
			kernel.setReductionValue(redVal);
			
			// 3
			kernel.run(States.BONDV1_CALC_UVW);
			uvw_halo = kernel.get_uvw_halo();

			uvw_halo = exchange(ip+3, jp+3, kp+3, 4, uvw_halo);
			
			kernel.set_uvw_halo(uvw_halo);
			
			// 4
			kernel.run(States.VELFG__FEEDBF__LES_CALC_SM);
			uvw_halo = kernel.get_uvw_halo();
			uvwsum_halo = kernel.get_uvwsum_halo();
			float [] fgh_halo = kernel.get_fgh_halo();
			float [] diu_halo = kernel.get_diu_halo();
			float [] sm_halo = kernel.get_sm_halo();

			uvw_halo = exchange(ip+3, jp+3, kp+3, 4, uvw_halo);
			uvwsum_halo = exchange(ip+1, jp+1, kp+1, 4, uvwsum_halo);
			fgh_halo = exchange(ip+1, jp+1, kp+1, 4, fgh_halo);
			diu_halo = exchange(ip+4, jp+3, kp+3, 16, diu_halo);
			sm_halo = exchange(ip+3, jp+3, kp+2, 1, sm_halo);

			kernel.set_uvw_halo(uvw_halo);
			kernel.set_uvwsum_halo(uvwsum_halo);
			kernel.set_fgh_halo(fgh_halo);
			kernel.set_diu_halo(diu_halo);
			kernel.set_sm_halo(sm_halo);
			
			// 5
			kernel.run(States.LES_BOUND_SM);
			sm_halo = kernel.get_sm_halo();

			sm_halo = exchange(ip+3, jp+3, kp+2, 1, sm_halo);
			
			kernel.set_sm_halo(sm_halo);
			
			// 6
			kernel.run(States.LES_CALC_VISC__ADAM);
			fgh_halo = kernel.get_fgh_halo();
			float [] fgh_old_halo = kernel.get_fgh_old_halo();

			fgh_halo = exchange(ip+1, jp+1, kp+1, 4, fgh_halo);
			fgh_old_halo = exchange(ip, jp, kp, 4, fgh_old_halo);
			
			kernel.set_fgh_halo(fgh_halo);
			kernel.set_fgh_old_halo(fgh_old_halo);
			
			// 7
			kernel.run(States.PRESS_RHSAV);
			fgh_halo = kernel.get_fgh_halo();
			float [] rhs_halo = kernel.get_rhs_halo();
			redVal = kernel.getReductionNominator() / kernel.getReductionDenominator();
			
			// reduction
			fgh_halo = exchange(ip+1, jp+1, kp+1, 4, fgh_halo);
			rhs_halo = exchange(ip+2, jp+2, kp+2, 1, rhs_halo);
			
			kernel.set_fgh_halo(fgh_halo);
			kernel.set_rhs_halo(rhs_halo);
			kernel.setReductionValue(redVal);
			
			// 8
			pressSORState(kernel, ip, jp, kp);
			
			// 9
			kernel.run(States.PRESS_PAV);
			p_halo = kernel.get_p_halo();
			redVal = kernel.getReductionNominator() / kernel.getReductionDenominator();
			
			// reduction
			p_halo = exchange(ip+3, jp+3, kp+2, 2, p_halo);
			
			kernel.set_p_halo(p_halo);
			kernel.setReductionValue(redVal);

			// 10
			kernel.run(States.PRESS_ADJ);
			p_halo = kernel.get_p_halo();
			
			p_halo = exchange(ip+3, jp+3, kp+2, 2, p_halo);
			
			kernel.set_p_halo(p_halo);

			// 11
			kernel.run(States.PRESS_BOUNDP);
			p_halo = kernel.get_p_halo();
			
			p_halo = exchange(ip+3, jp+3, kp+2, 2, p_halo);
			
			kernel.set_p_halo(p_halo);
			
			iter++;
		}
	}
	
	private static void pressSORState(HalosCPU kernel, int ip, int jp, int kp) {
		float pjuge = 0.0001f;
		int nmaxp = 50;
		float sor = pjuge * 1.1f;
		int iter = 0;
		
		System.out.println("Starting SOR");

		while (sor > pjuge && iter < nmaxp){
			iter++;
			
			for(int i = 0; i <= 2; i++){
				kernel.pressSORIteration(i);
				float [] p_halo = kernel.get_p_halo();
				
				p_halo = exchange(ip+3, jp+3, kp+2, 2, p_halo);
				
				kernel.set_p_halo(p_halo);
				
				if(i == 1){
					float redVal = kernel.getPressSORValue();
					// Reduction					
					sor = (float) Math.sqrt(redVal);
				}
			}
			
			System.out.println(iter + " " + sor);
		}
	}
	
	private static float [] exchange(int im, int jm, int km, int v_dim, float [] innerHalo) {
		float [] outerHalo = new float[2 * v_dim * (im+jm + 2) * km];
		float [] northHalo = new float[v_dim * km * im];
		float [] southHalo = new float[v_dim * km * im];
		float [] westHalo = new float[v_dim * km * jm];
		float [] eastHalo = new float[v_dim * km * jm];
		float [] northWestHalo = new float[v_dim * km];
		float [] northEastHalo = new float[v_dim * km];
		float [] southWestHalo = new float[v_dim * km];
		float [] southEastHalo = new float[v_dim * km];
		
		HaloConstructor hc = new HaloConstructor();
		
		hc.deconstruct(innerHalo, northHalo, southHalo, westHalo, eastHalo, northWestHalo, northEastHalo, southWestHalo, southEastHalo, v_dim, im, jm, km);
		hc.construct(outerHalo, southHalo, northHalo, eastHalo, westHalo, southEastHalo, southWestHalo, northEastHalo, northWestHalo, v_dim, im, jm, km);
		
		return outerHalo;
	}
}