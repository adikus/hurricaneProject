package hoos.project.LES.HaloExchange;

import java.util.ArrayList;

import scala.Tuple2;

public class HaloExchanger {
	int im, jm, km, v_dim, X, Y;
	
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

	public Iterable<Tuple2<Integer, Neighbour>> deconstructIntoPairs(float[] innerHalo, Integer i, int X, int Y) {
		this.X = X;
		this.Y = Y;
		
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
		
		ArrayList<Tuple2<Integer, Neighbour>> neighbourPairs = new ArrayList<Tuple2<Integer, Neighbour>>();
		
		int x = i % X;
		int y = i / Y;
		
		neighbourPairs.add(new Tuple2<Integer, Neighbour>(key(x, y-1), new Neighbour(northHalo, "N")));
		neighbourPairs.add(new Tuple2<Integer, Neighbour>(key(x, y+1), new Neighbour(southHalo, "S")));
		neighbourPairs.add(new Tuple2<Integer, Neighbour>(key(x-1, y), new Neighbour(westHalo, "W")));
		neighbourPairs.add(new Tuple2<Integer, Neighbour>(key(x+1, y), new Neighbour(eastHalo, "E")));
		
		neighbourPairs.add(new Tuple2<Integer, Neighbour>(key(x-1, y-1), new Neighbour(northWestHalo, "NW")));
		neighbourPairs.add(new Tuple2<Integer, Neighbour>(key(x+1, y-1), new Neighbour(northEastHalo, "NE")));
		neighbourPairs.add(new Tuple2<Integer, Neighbour>(key(x-1, y+1), new Neighbour(southWestHalo, "SW")));
		neighbourPairs.add(new Tuple2<Integer, Neighbour>(key(x+1, y+1), new Neighbour(southEastHalo, "SE")));
		
		return neighbourPairs;
	}

	public float[] constructFromPairs(Iterable<Neighbour> neighbours) {
		float [] outerHalo = new float[2 * v_dim * (im+jm + 2) * km];
		float [] northHalo = new float[v_dim * km * im];
		float [] southHalo = new float[v_dim * km * im];
		float [] westHalo = new float[v_dim * km * jm];
		float [] eastHalo = new float[v_dim * km * jm];
		float [] northWestHalo = new float[v_dim * km];
		float [] northEastHalo = new float[v_dim * km];
		float [] southWestHalo = new float[v_dim * km];
		float [] southEastHalo = new float[v_dim * km];
		
		for(Neighbour neighbour : neighbours){
			switch(neighbour.getName()){
			case "N":
				southHalo = neighbour.getArray();
				System.out.println("Got southHalo");
				break;
			case "S":
				northHalo = neighbour.getArray();
				System.out.println("Got northHalo");
				break;
			case "W":
				eastHalo = neighbour.getArray();
				System.out.println("Got eastHalo");
				break;
			case "E":
				westHalo = neighbour.getArray();
				System.out.println("Got westHalo");
				break;
			case "NW":
				southEastHalo = neighbour.getArray();
				System.out.println("Got southEastHalo");
				break;
			case "NE":
				southWestHalo = neighbour.getArray();
				System.out.println("Got southWestHalo");
				break;
			case "SW":
				northEastHalo = neighbour.getArray();
				System.out.println("Got northEastHalo");
				break;
			case "SE":
				northWestHalo = neighbour.getArray();
				System.out.println("Got northWestHalo");
				break;
			}
		}
		
		HaloConstructor hc = new HaloConstructor();
		hc.construct(outerHalo, southHalo, northHalo, eastHalo, westHalo, southEastHalo, southWestHalo, northEastHalo, northWestHalo, v_dim, im, jm, km);
		
		return outerHalo;
	}
	
	private Integer key(int x, int y) {
		int dx = (x % X);
		int dy = (y % Y);
		if(dx < 0)dx += X;
		if(dy < 0)dy += Y;
		return (dx) + (dy)*X;
	}
}
