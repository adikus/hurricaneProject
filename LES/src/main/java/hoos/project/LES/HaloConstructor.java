package hoos.project.LES;

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
				for(int i = 0; i < im; i++){
					halo[n++] = north[v*km*im + k*im + i];
				}	
				//NE
				halo[n++] = northeast[v*km + k];
			}
			for(int k = 0; k < km; k++){				
				//SW
				halo[n++] = southwest[v*km + k];			
				//S
				for(int i = 0; i < im; i++){
					halo[n++] = south[v*km*im + k*im + i];
				}			
				//SE
				halo[n++] = southeast[v*km + k];
			}
			for(int k = 0; k < km; k++){				
				//W
				for(int j = 0; j < im; j++){
					halo[n++] = west[v*km*jm + k*jm + j];
				}
			}
			for(int k = 0; k < km; k++){				
				//E
				for(int j = 0; j < im; j++){
					halo[n++] = east[v*km*jm + k*jm + j];
				}
			}
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
				for(int i = 0; i < im; i++){
					north[v*km*im + k*im + i] = halo[v*2*(im+jm-2)*km + k*im + i];
					south[v*km*im + k*im + i] = halo[v*2*(im+jm-2)*km + km*im + k*im + i];
				}
				
				northwest[v*km + k] = halo[v*2*(im+jm-2)*km + k*im];
				northeast[v*km + k] = halo[v*2*(im+jm-2)*km + k*im + im - 1];
				
				southwest[v*km + k] = halo[v*2*(im+jm-2)*km + km*im + k*im];
				southeast[v*km + k] = halo[v*2*(im+jm-2)*km + km*im + k*im + im - 1];
				
				for(int j = 0; j < jm-2; j++){
					west[1+v*km*jm + k*jm + j] = halo[v*2*(im+jm-2)*km + km*im*2 + k*(jm-2) + j];
					east[1+v*km*jm + k*jm + j] = halo[v*2*(im+jm-2)*km + km*im*2 + km*(jm-2) + k*(jm-2) + j];
				}
				west[v*km*jm + k*jm] = northwest[v*km + k];
				west[v*km*jm + k*jm + jm - 1] = southwest[v*km + k];
				east[v*km*jm + k*jm] = northeast[v*km + k];
				east[v*km*jm + k*jm + jm - 1] = southeast[v*km + k];
			}
		}
	}
}
