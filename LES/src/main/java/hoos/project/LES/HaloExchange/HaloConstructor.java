package hoos.project.LES.HaloExchange;

/**
 * @author      Andrej Hoos
 * Class for constructing and deconstructing halos
 */
public class HaloConstructor {
	
	// Combines north, south, west, east, northwest, northeast, southwest, southeast into outer halo
	public void construct(
		float [] halo, 
		float [] north, float [] south, float [] west, float [] east,
		float [] northwest, float [] northeast, float [] southwest, float [] southeast,
		int v_dim, int im, int jm, int km
	) {
		int n = 0;
		
		for(int v = 0; v < v_dim; v++){
			for(int k = 0; k < km; k++){
				//NW
				halo[n++] = northwest[v*km + k];		
				//N
				System.arraycopy(north, v*km*im + k*im, halo, n, im);
				n += im;
				//NE
				halo[n++] = northeast[v*km + k];
			}
			for(int k = 0; k < km; k++){				
				//SW
				halo[n++] = southwest[v*km + k];			
				//S
				System.arraycopy(south, v*km*im + k*im, halo, n, im);
				n += im;
				//SE
				halo[n++] = southeast[v*km + k];
			}
			System.arraycopy(west, v*km*jm, halo, n, km*jm);
			n += km*jm;
			System.arraycopy(east, v*km*jm, halo, n, km*jm);
			n += km*jm;
		}
	}
	
	// Deconstruct inner halo into north, south, west, east, northwest, northeast, southwest, southeast
	public void deconstruct(
		float [] halo, 
		float [] north, float [] south, float [] west, float [] east,
		float [] northwest, float [] northeast, float [] southwest, float [] southeast,
		int v_dim, int im, int jm, int km
	) {
		for(int v = 0; v < v_dim; v++){
			for(int k = 0; k < km; k++){
				System.arraycopy(halo, v*2*(im+jm-2)*km + k*im, north, v*km*im + k*im, im);
				System.arraycopy(halo, v*2*(im+jm-2)*km + km*im + k*im, south, v*km*im + k*im, im);
				
				northwest[v*km + k] = halo[v*2*(im+jm-2)*km + k*im];
				northeast[v*km + k] = halo[v*2*(im+jm-2)*km + k*im + im - 1];
				
				southwest[v*km + k] = halo[v*2*(im+jm-2)*km + km*im + k*im];
				southeast[v*km + k] = halo[v*2*(im+jm-2)*km + km*im + k*im + im - 1];
				
				System.arraycopy(halo, v*2*(im+jm-2)*km + km*im*2 + k*(jm-2), west, 1+v*km*jm + k*jm, jm-2);
				System.arraycopy(halo, v*2*(im+jm-2)*km + km*im*2 + km*(jm-2) + k*(jm-2), east, 1+v*km*jm + k*jm, jm-2);

				west[v*km*jm + k*jm] = northwest[v*km + k];
				west[v*km*jm + k*jm + jm - 1] = southwest[v*km + k];
				east[v*km*jm + k*jm] = northeast[v*km + k];
				east[v*km*jm + k*jm + jm - 1] = southeast[v*km + k];
			}
		}
	}
}
