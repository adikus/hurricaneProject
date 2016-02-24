package hoos.project.LES.HaloExchange;

public class HaloExchanger {
	int im, jm, km, v_dim;
	
	public HaloExchanger(int im, int jm, int km, int v_dim) {
		this.im = im;
		this.jm = jm;
		this.km = km;
		this.v_dim = v_dim;
	}
	
	public float [] exchange(float [] innerHalo) {
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
		
		//long time = System.nanoTime();
		hc.deconstruct(innerHalo, northHalo, southHalo, westHalo, eastHalo, northWestHalo, northEastHalo, southWestHalo, southEastHalo, v_dim, im, jm, km);
		hc.construct(outerHalo, southHalo, northHalo, eastHalo, westHalo, southEastHalo, southWestHalo, northEastHalo, northWestHalo, v_dim, im, jm, km);
		//System.out.println("Halo exchange (" + im + ", " + jm + ", " + km + ", " + v_dim + ") took: " + (System.nanoTime() - time));
		
		return outerHalo;
	}
}
