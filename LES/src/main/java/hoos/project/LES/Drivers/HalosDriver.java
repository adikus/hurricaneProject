package hoos.project.LES.Drivers;

import hoos.project.LES.HaloExchange.HaloExchanger;
import hoos.project.LES.Kernels.Halos;
import hoos.project.LES.Kernels.States;

public class HalosDriver{

	public static void main(String[] _args) {
		Halos kernel = new Halos();

		System.out.printf("Initialising kernel..");
		
		int ip = 150;
		int jp = 150;
		int kp = 90;
		
		kernel.init(ip, jp, kp);
		
		long previousTime = kernel.getAccumulatedExecutionTime();
		System.out.println("First execution time: " + previousTime);
		System.out.println("Conversion execution time: " + kernel.getConversionTime());

		System.out.println("Running kernel..");
		
		int iter = 0;
		while(iter < 20){			
			// 1
			kernel.run(States.VELNW__BONDV1_INIT_UVW);
			float [] p_halo = kernel.get_p_halo();
			float [] uvw_halo = kernel.get_uvw_halo();
			float [] uvwsum_halo = kernel.get_uvwsum_halo();
			
			p_halo = new HaloExchanger(ip+3, jp+3, kp+2, 2).exchange(p_halo);
			uvw_halo = new HaloExchanger(ip+3, jp+3, kp+3, 4).exchange(uvw_halo);
			uvwsum_halo = new HaloExchanger(ip+1, jp+1, kp+1, 4).exchange(uvwsum_halo);
			
			kernel.set_p_halo(p_halo);
			kernel.set_uvw_halo(uvw_halo);
			kernel.set_uvwsum_halo(uvwsum_halo);
			
			// 2
			kernel.run(States.BONDV1_CALC_UOUT);
			uvw_halo = kernel.get_uvw_halo();
			float redVal = kernel.getReductionValue();
			
			// reduction
			uvw_halo = new HaloExchanger(ip+3, jp+3, kp+3, 4).exchange(uvw_halo);
			
			kernel.set_uvw_halo(uvw_halo);
			kernel.setReductionValue(redVal);
			
			// 3
			kernel.run(States.BONDV1_CALC_UVW);
			uvw_halo = kernel.get_uvw_halo();

			uvw_halo = new HaloExchanger(ip+3, jp+3, kp+3, 4).exchange(uvw_halo);
			
			kernel.set_uvw_halo(uvw_halo);
			
			// 4
			kernel.run(States.VELFG__FEEDBF__LES_CALC_SM);
			uvw_halo = kernel.get_uvw_halo();
			uvwsum_halo = kernel.get_uvwsum_halo();
			float [] fgh_halo = kernel.get_fgh_halo();
			float [] diu_halo = kernel.get_diu_halo();
			float [] sm_halo = kernel.get_sm_halo();

			uvw_halo = new HaloExchanger(ip+3, jp+3, kp+3, 4).exchange(uvw_halo);
			uvwsum_halo = new HaloExchanger(ip+1, jp+1, kp+1, 4).exchange(uvwsum_halo);
			fgh_halo = new HaloExchanger(ip+1, jp+1, kp+1, 4).exchange(fgh_halo);
			diu_halo = new HaloExchanger(ip+4, jp+3, kp+3, 16).exchange(diu_halo);
			sm_halo = new HaloExchanger(ip+3, jp+3, kp+2, 1).exchange(sm_halo);

			kernel.set_uvw_halo(uvw_halo);
			kernel.set_uvwsum_halo(uvwsum_halo);
			kernel.set_fgh_halo(fgh_halo);
			kernel.set_diu_halo(diu_halo);
			kernel.set_sm_halo(sm_halo);
			
			// 5
			kernel.run(States.LES_BOUND_SM);
			sm_halo = kernel.get_sm_halo();

			sm_halo = new HaloExchanger(ip+3, jp+3, kp+2, 1).exchange(sm_halo);
			
			kernel.set_sm_halo(sm_halo);
			
			// 6
			kernel.run(States.LES_CALC_VISC__ADAM);
			fgh_halo = kernel.get_fgh_halo();
			float [] fgh_old_halo = kernel.get_fgh_old_halo();

			fgh_halo = new HaloExchanger(ip+1, jp+1, kp+1, 4).exchange(fgh_halo);
			fgh_old_halo = new HaloExchanger(ip, jp, kp, 4).exchange(fgh_old_halo);
			
			kernel.set_fgh_halo(fgh_halo);
			kernel.set_fgh_old_halo(fgh_old_halo);
			
			// 7
			kernel.run(States.PRESS_RHSAV);
			fgh_halo = kernel.get_fgh_halo();
			float [] rhs_halo = kernel.get_rhs_halo();
			redVal = kernel.getReductionValue();
			
			// reduction
			fgh_halo = new HaloExchanger(ip+1, jp+1, kp+1, 4).exchange(fgh_halo);
			rhs_halo = new HaloExchanger(ip+2, jp+2, kp+2, 1).exchange(rhs_halo);
			
			kernel.set_fgh_halo(fgh_halo);
			kernel.set_rhs_halo(rhs_halo);
			kernel.setReductionValue(redVal);
			
			// 8
			pressSORState(kernel, ip, jp, kp);
			
			// 9
			kernel.run(States.PRESS_PAV);
			p_halo = kernel.get_p_halo();
			redVal = kernel.getReductionValue();
			
			// reduction
			p_halo = new HaloExchanger(ip+3, jp+3, kp+2, 2).exchange(p_halo);
			
			kernel.set_p_halo(p_halo);
			kernel.setReductionValue(redVal);

			// 10
			kernel.run(States.PRESS_ADJ);
			p_halo = kernel.get_p_halo();
			
			p_halo = new HaloExchanger(ip+3, jp+3, kp+2, 2).exchange(p_halo);
			
			kernel.set_p_halo(p_halo);

			// 11
			kernel.run(States.PRESS_BOUNDP);
			p_halo = kernel.get_p_halo();
			
			p_halo = new HaloExchanger(ip+3, jp+3, kp+2, 2).exchange(p_halo);
			
			kernel.set_p_halo(p_halo);
			
			/*long time = kernel.getAccumulatedExecutionTime();
			System.out.println("Current state: " + i + " duration: " + (time - previousTime));
			System.out.println("Total kernel execution time: " + time);
			previousTime = time;*/
			
			iter++;
		}
	}
	
	private static void pressSORState(Halos kernel, int ip, int jp, int kp) {
		float pjuge = 0.0001f;
		int nmaxp = 50;
		float sor = pjuge * 1.1f;
		int iter = 0;
		
		System.out.println("Starting SOR");
		kernel.run(States.PRESS_SOR);

		while (sor > pjuge && iter < nmaxp){
			iter++;
			
			for(int i = 0; i <= 2; i++){
				kernel.pressSORIteration(i);
				float [] p_halo = kernel.get_p_halo();
				
				p_halo = new HaloExchanger(ip+3, jp+3, kp+2, 2).exchange(p_halo);
				
				kernel.set_p_halo(p_halo);
				
				if(i == 1){
					float redVal = kernel.getPressSORValue();
					// Reduction					
					sor = redVal;
				}
			}
			
			//System.out.println(iter + " " + sor);
		}
	}
}