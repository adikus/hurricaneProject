package hoos.project.LES.spark;

import hoos.project.LES.Kernels.Halos;

import java.util.Arrays;

import org.apache.spark.api.java.*;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.SparkConf;

import scala.Tuple2;

public class SparkHalosDriver {
	
	public static void main(String[] _args) {
		SparkConf conf = new SparkConf().setAppName("LES Spark Simulation");
		JavaSparkContext sc = new JavaSparkContext(conf);
		sc.setCheckpointDir("temp");
		
		// Number of nodes
		int X = 1;
		int Y = 1;
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
				System.out.println("Working Directory = " + System.getProperty("user.dir"));
				kernel.init(ip, jp, kp);
				return kernel;
			}
		});

		System.out.println(kernels.count());
		
		
		sc.close();
	}
	
}
