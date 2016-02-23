package hoos.project.LES;
import java.util.Arrays;

import hoos.project.LES.HaloConstructor;

public class HaloTest {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int ip = 1;
		int jp = 1;
		int kp = 0;
		
		float [] p_halo = new float[2 * 2 * (ip+jp+6 + 2) * (kp+2)];
		float [] p_north = new float[2 * (kp+2) * (ip+3)];
		float [] p_east = new float[2 * (kp+2) * (jp+3)];
		float [] p_south = new float[2 * (kp+2) * (ip+3)];
		float [] p_west = new float[2 * (kp+2) * (jp+3)];
		float [] p_northwest = new float[2 * (kp+2)];
		float [] p_northeast = new float[2 * (kp+2)];
		float [] p_southwest = new float[2 * (kp+2)];
		float [] p_southeast = new float[2 * (kp+2)];
		
		for(int i = 0; i < p_halo.length; i++) {
			p_halo[i] = i;
		}
		System.out.println(Arrays.toString(p_halo));
		
		HaloConstructor hc = new HaloConstructor();
		hc.deconstruct(p_halo, p_north, p_south, p_west, p_east, p_northwest, p_northeast, p_southwest, p_southeast, 2, ip+3, jp+3, kp+2);
		
		System.out.println(Arrays.toString(p_north));
		System.out.println(Arrays.toString(p_south));
		System.out.println(Arrays.toString(p_northwest));
		System.out.println(Arrays.toString(p_northeast));
		System.out.println(Arrays.toString(p_southwest));
		System.out.println(Arrays.toString(p_southeast));
		System.out.println(Arrays.toString(p_west));
		System.out.println(Arrays.toString(p_east));
		
		hc.construct(p_halo, p_north, p_south, p_west, p_east, p_northwest, p_northeast, p_southwest, p_southeast, 2, ip+3, jp+3, kp+2);
		
		System.out.println(Arrays.toString(p_halo));
	}

}
