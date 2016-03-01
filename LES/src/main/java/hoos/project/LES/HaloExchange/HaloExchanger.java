package hoos.project.LES.HaloExchange;

import hoos.project.LES.Kernels.Halos;

import java.util.ArrayList;

import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.PairFlatMapFunction;

import scala.Tuple2;

public class HaloExchanger {
	int im, jm, km, v_dim, X, Y;
	
	public HaloExchanger(int im, int jm, int km, int v_dim) {
		this.im = im;
		this.jm = jm;
		this.km = km;
		this.v_dim = v_dim;
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
				break;
			case "S":
				northHalo = neighbour.getArray();
				break;
			case "W":
				eastHalo = neighbour.getArray();
				break;
			case "E":
				westHalo = neighbour.getArray();
				break;
			case "NW":
				southEastHalo = neighbour.getArray();
				break;
			case "NE":
				southWestHalo = neighbour.getArray();
				break;
			case "SW":
				northEastHalo = neighbour.getArray();
				break;
			case "SE":
				northWestHalo = neighbour.getArray();
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
	
	public static JavaPairRDD<Integer, Halos> exchangeHalos(JavaPairRDD<Integer, Halos> kernels, String halos, final int ip, final int jp, final int kp, final int X, final int Y) {
		return exchangeHalos(kernels, halos, ip, jp, kp, X, Y, false);
	}
	
	public static JavaPairRDD<Integer, Halos> exchangeHalos(JavaPairRDD<Integer, Halos> kernels, String halos, final int ip, final int jp, final int kp, final int X, final int Y, Boolean checkpoint) {
		final String[] haloNames = halos.split(",");
		
		for(final String haloName : haloNames){
			JavaPairRDD<Integer, Neighbour> neighbours = kernels.flatMapToPair(new PairFlatMapFunction<Tuple2<Integer,Halos>, Integer, Neighbour>() {
				public Iterable<Tuple2<Integer, Neighbour>> call(Tuple2<Integer, Halos> pair) throws Exception {
					System.out.println("Deconstructing and sending halos for " + haloName);
					switch(haloName){
					case "p":
						return new HaloExchanger(ip+3, jp+3, kp+2, 2).deconstructIntoPairs(pair._2().get_p_halo(), pair._1(), X, Y);
					case "uvw":
						return new HaloExchanger(ip+3, jp+3, kp+3, 4).deconstructIntoPairs(pair._2().get_uvw_halo(), pair._1(), X, Y);
					case "uvwsum":
						return new HaloExchanger(ip+1, jp+1, kp+1, 4).deconstructIntoPairs(pair._2().get_uvwsum_halo(), pair._1(), X, Y);	
					case "fgh":
						return new HaloExchanger(ip+1, jp+1, kp+1, 4).deconstructIntoPairs(pair._2().get_fgh_halo(), pair._1(), X, Y);
					case "fgh_old":
						return new HaloExchanger(ip, jp, kp, 4).deconstructIntoPairs(pair._2().get_fgh_old_halo(), pair._1(), X, Y);
					case "diu":
						return new HaloExchanger(ip+4, jp+3, kp+3, 16).deconstructIntoPairs(pair._2().get_diu_halo(), pair._1(), X, Y);
					case "rhs":
						return new HaloExchanger(ip+2, jp+2, kp+2, 1).deconstructIntoPairs(pair._2().get_rhs_halo(), pair._1(), X, Y);
					case "sm":
						return new HaloExchanger(ip+3, jp+3, kp+2, 1).deconstructIntoPairs(pair._2().get_sm_halo(), pair._1(), X, Y);
					default:
						return null;							
					}
				}
			});
			
			JavaPairRDD<Integer, Halos> newKernels = kernels.cogroup(neighbours).mapValues(new Function<Tuple2<Iterable<Halos>,Iterable<Neighbour>>, Halos>() {
				public Halos call(Tuple2<Iterable<Halos>, Iterable<Neighbour>> pair){
					Halos kernel = pair._1().iterator().next();
					System.out.println("Receiveing and constructing halos for " + haloName);
					switch(haloName){
					case "p":
						kernel.set_p_halo(new HaloExchanger(ip+3, jp+3, kp+2, 2).constructFromPairs(pair._2()));
						break;
					case "uvw":
						kernel.set_uvw_halo(new HaloExchanger(ip+3, jp+3, kp+3, 4).constructFromPairs(pair._2()));
						break;
					case "uvwsum":
						kernel.set_uvwsum_halo(new HaloExchanger(ip+1, jp+1, kp+1, 4).constructFromPairs(pair._2()));
						break;	
					case "fgh":
						kernel.set_fgh_halo(new HaloExchanger(ip+1, jp+1, kp+1, 4).constructFromPairs(pair._2()));
						break;
					case "fgh_old":
						kernel.set_fgh_old_halo(new HaloExchanger(ip, jp, kp, 4).constructFromPairs(pair._2()));
						break;
					case "diu":
						kernel.set_diu_halo(new HaloExchanger(ip+4, jp+3, kp+3, 16).constructFromPairs(pair._2()));
						break;
					case "rhs":
						kernel.set_rhs_halo(new HaloExchanger(ip+2, jp+2, kp+2, 1).constructFromPairs(pair._2()));
						break;
					case "sm":
						kernel.set_sm_halo(new HaloExchanger(ip+3, jp+3, kp+2, 1).constructFromPairs(pair._2()));
						break;				
					}
					
					
					return kernel;
				}
				
			});
			newKernels.cache();
			if(checkpoint)newKernels.checkpoint();
			newKernels.count();
			kernels.unpersist();
			kernels = newKernels;
		}
				
		return kernels;
	}
}
