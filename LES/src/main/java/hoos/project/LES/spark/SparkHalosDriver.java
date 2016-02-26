package hoos.project.LES.spark;

import hoos.project.LES.HaloExchange.HaloExchanger;
import hoos.project.LES.HaloExchange.Neighbour;
import hoos.project.LES.Kernels.Halos;
import hoos.project.LES.Kernels.States;

import java.util.Arrays;

import org.apache.spark.api.java.*;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.PairFlatMapFunction;
import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.SparkConf;

import scala.Tuple2;

public class SparkHalosDriver {
	
	public static void main(String[] _args) {
		SparkConf conf = new SparkConf().setAppName("LES Spark Simulation");
		JavaSparkContext sc = new JavaSparkContext(conf);
		sc.setCheckpointDir("temp");
		
		// Number of nodes
		final int X = 2;
		final int Y = 1;
		Integer[] a = new Integer[X*Y];
		for(int i = 0; i < X*Y; i++){ a[i]=i; }
		JavaRDD<Integer> indexes = sc.parallelize(Arrays.asList(a));
		JavaPairRDD<Integer, Halos> kernels = indexes.mapToPair(new PairFunction<Integer, Integer, Halos>() {
			public Tuple2<Integer, Halos> call(Integer index) { 
				return new Tuple2<Integer, Halos>(index, new Halos());
			}
		}).partitionBy(new KeyPartitioner(X*Y)).cache();
		
		final int ip = 150;
		final int jp = 150;
		final int kp = 90;

		System.out.println("Running first map");
		
		kernels = kernels.mapValues(new Function<Halos, Halos>() {
			public Halos call(Halos kernel){
				kernel.init(ip, jp, kp);
				
				kernel.run(States.VELNW__BONDV1_INIT_UVW);
				
				return kernel;
			}
		});
		
		kernels.cache();
		
		JavaPairRDD<Integer, Neighbour> p_neighbourhs = kernels.flatMapToPair(new PairFlatMapFunction<Tuple2<Integer,Halos>, Integer, Neighbour>() {
			public Iterable<Tuple2<Integer, Neighbour>> call(Tuple2<Integer, Halos> pair) throws Exception {
				float [] p_halo = pair._2().get_p_halo();
				return new HaloExchanger(ip+3, jp+3, kp+2, 2).deconstructIntoPairs(p_halo, pair._1(), X, Y);
			}
		});
		
		kernels = kernels.cogroup(p_neighbourhs).mapValues(new Function<Tuple2<Iterable<Halos>,Iterable<Neighbour>>, Halos>() {
			public Halos call(Tuple2<Iterable<Halos>, Iterable<Neighbour>> pair){
				Halos kernel = pair._1().iterator().next();
				
				float[] outerHalo = new HaloExchanger(ip+3, jp+3, kp+2, 2).constructFromPairs(pair._2());
				kernel.set_p_halo(outerHalo);
				return kernel;
			}
			
		});

		System.out.println(kernels.count());
		
		
		sc.close();
	}
	
}
