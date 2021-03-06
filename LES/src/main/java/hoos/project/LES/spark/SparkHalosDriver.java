package hoos.project.LES.spark;

import hoos.project.LES.HaloExchange.HaloExchanger;
import hoos.project.LES.Kernels.Halos;
import hoos.project.LES.Kernels.States;

import java.util.Arrays;

import org.apache.spark.api.java.*;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.SparkConf;

import scala.Tuple2;

/**
 * @author      Andrej Hoos
 * Class that runs LES inside Spark
 */
public class SparkHalosDriver {
	
	public static void main(String[] args) {
		if(args.length < 3){
			System.err.println("Please provide the number of iterations and the domain size");
			return;
		}
		final int iterationsNum = new Integer(args[0]);
		
		SparkConf conf = new SparkConf().setAppName("LES Spark Simulation (" + iterationsNum + ")");
		JavaSparkContext sc = new JavaSparkContext(conf);
		sc.setCheckpointDir("temp");
		
		// Number of nodes
		final int X = new Integer(args[1]);
		final int Y = new Integer(args[2]);
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
		
		// Init
		kernels = kernels.mapValues(new Function<Halos, Halos>() {
			public Halos call(Halos kernel){
				kernel.init(ip, jp, kp);
				
				return kernel;
			}
		});
		
		for(int i = 1; i <= iterationsNum; i++){
			System.out.println("Start time step " + i);
			// 1
			kernels = executeKernelStep(kernels, States.VELNW__BONDV1_INIT_UVW);
			kernels = HaloExchanger.exchangeHalos(kernels, "p,uvw,uvwsum", ip, jp, kp, X, Y);
			
			// 2
			kernels = executeKernelStep(kernels, States.BONDV1_CALC_UOUT);
			kernels = HaloExchanger.exchangeHalos(kernels, "uvw", ip, jp, kp, X, Y);

			// 3
			kernels = executeKernelStep(kernels, States.BONDV1_CALC_UVW);
			kernels = HaloExchanger.exchangeHalos(kernels, "uvw", ip, jp, kp, X, Y);
			
			// 4
			kernels = executeKernelStep(kernels, States.VELFG__FEEDBF__LES_CALC_SM);
			kernels = HaloExchanger.exchangeHalos(kernels, "uvw,uvwsum,fgh,diu,sm", ip, jp, kp, X, Y);
					
			// 5
			kernels = executeKernelStep(kernels, States.LES_BOUND_SM);
			kernels = HaloExchanger.exchangeHalos(kernels, "sm", ip, jp, kp, X, Y);
			
			// 6
			kernels = executeKernelStep(kernels, States.LES_CALC_VISC__ADAM);
			kernels = HaloExchanger.exchangeHalos(kernels, "fgh,fgh_old", ip, jp, kp, X, Y);
			
			// 7
			kernels = executeKernelStep(kernels, States.PRESS_RHSAV);
			kernels = HaloExchanger.exchangeHalos(kernels, "rhs,fgh", ip, jp, kp, X, Y);
			kernels = divisionReduction(kernels);
			
			// 8
			kernels = pressSORState(kernels, ip, jp, kp, X, Y);
			
			// 9
			kernels = executeKernelStep(kernels, States.PRESS_PAV);
			kernels = HaloExchanger.exchangeHalos(kernels, "p", ip, jp, kp, X, Y);
			kernels = divisionReduction(kernels);
			
			// 10
			kernels = executeKernelStep(kernels, States.PRESS_ADJ);
			kernels = HaloExchanger.exchangeHalos(kernels, "p", ip, jp, kp, X, Y);
					
			// 11
			kernels = executeKernelStep(kernels, States.PRESS_BOUNDP);
			kernels = HaloExchanger.exchangeHalos(kernels, "p", ip, jp, kp, X, Y, true);
		}
		
		// Force final Spark computation
		kernels.count();
		
		sc.close();
	}
	
	/**
	 * Helper method for running the 8th step of the simulation
	 * @param kernels the kernel collection
	 * @param ip x size of the domain on a single node 
	 * @param jp y size of the domain on a single node 
	 * @param kp z size of the domain on a single node 
	 * @param X node grid x size
	 * @param Y node grid y size
	 * @return the kernel collection
	 */
	private static JavaPairRDD<Integer, Halos> pressSORState(JavaPairRDD<Integer, Halos> kernels, int ip, int jp, int kp, int X, int Y) {
		float pjuge = 0.0001f;
		int nmaxp = 50;
		float sor = pjuge * 1.1f;
		int iter = 0;

		while (sor > pjuge && iter < nmaxp){
			iter++;
			
			for(int i = 0; i <= 2; i++){
				final int iteration = iter;
				final int ii = i;
				JavaPairRDD<Integer, Halos> newKernels = kernels.mapValues(new Function<Halos, Halos>() {
					public Halos call(Halos kernel){
						kernel.pressSORIteration(ii);
						return kernel;
					}
				});
				newKernels.cache();
				newKernels.count();
				kernels.unpersist();
				if(iter % 20 == 0 && i == 0){
					kernels = HaloExchanger.exchangeHalos(newKernels, "p", ip, jp, kp, X, Y, true);
				}else{
					kernels = HaloExchanger.exchangeHalos(newKernels, "p", ip, jp, kp, X, Y);
				}
				
				if(i == 1){
					sor = kernels.values().map(new Function<Halos, Float>() {
						public Float call(Halos kernel){
							float sor = kernel.getPressSORValue(); 
							System.out.println("SOR: " + iteration + " " + Math.sqrt(sor));
							return sor;
						}
					}).reduce(new Function2<Float, Float, Float>() {
						public Float call(Float value1, Float value2){ return value1 + value2; }
					});
					sor = (float) Math.sqrt(sor);
				}
			}
		}
		
		return kernels;
	}
	
	/**
	 * Helper method to execute a single simulation step
	 * @param kernels the kernel collection
	 * @param state id of the step to be executed
	 * @return the kernel collection
	 */
	private static JavaPairRDD<Integer, Halos> executeKernelStep(JavaPairRDD<Integer, Halos> kernels, final int state) {
		JavaPairRDD<Integer, Halos> newKernels = kernels.mapValues(new Function<Halos, Halos>() {
			public Halos call(Halos kernel){
				kernel.run(state);
				return kernel;
			}
		});
		newKernels.cache();
		newKernels.count();
		kernels.unpersist();
		return newKernels;
	}
	
	/**
	 * Helper method to run a reduction over the global domain
	 * @param kernels the kernel collection
	 * @return the kernel collection
	 */
	private static JavaPairRDD<Integer, Halos> divisionReduction(JavaPairRDD<Integer, Halos> kernels) {
		Tuple2<Float, Float> reductionPair = kernels.values().map(new Function<Halos, Tuple2<Float, Float>>() {
			public Tuple2<Float, Float> call(Halos kernel) {
				return new Tuple2<Float, Float>(kernel.getReductionNominator(), kernel.getReductionDenominator());
			}			
		}).reduce(new Function2<Tuple2<Float, Float>, Tuple2<Float, Float>, Tuple2<Float, Float>>() {
			public Tuple2<Float, Float> call(Tuple2<Float, Float> pair1, Tuple2<Float, Float> pair2) {
				return new Tuple2<Float, Float>(pair1._1() + pair2._1(), pair1._2() + pair2._2());
			}
		});
		final float reductionValue = reductionPair._1() / reductionPair._2();
		JavaPairRDD<Integer, Halos> newKernels = kernels.mapValues(new Function<Halos, Halos>(){
			public Halos call(Halos kernel){
				kernel.setReductionValue(reductionValue);
				return kernel;
			}
		});
		newKernels.cache();
		newKernels.count();
		kernels.unpersist();
		return newKernels;
	}
	
}
