package hoos.project.dummy;

import org.apache.spark.api.java.*;
import org.apache.spark.HashPartitioner;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.api.java.function.PairFlatMapFunction;
import org.apache.spark.api.java.function.PairFunction;
import org.apache.spark.storage.StorageLevel;

import scala.Tuple2;
import hoos.project.dummy.datastore.ArrayDataStore;
import hoos.project.dummy.gol.GameOfLife;

import java.lang.Iterable;

public class DummyPartitionedJoin {
	public static void main(String[] args) {
		DataManipulationAlgorithm algorithm = new GameOfLife(args[0]);
		Chunk startingChunk = new Chunk(algorithm, new ArrayDataStore());
		
		final int N = 8;
		
		startingChunk.multiply(100);
		final int width = startingChunk.getWidth();
		final int height = startingChunk.getHeight();

		SparkConf conf = new SparkConf().setAppName("Dummy Map Reduce " + args[2]);
		JavaSparkContext sc = new JavaSparkContext(conf);
		sc.setCheckpointDir("temp");
		
		int widthCount = new Double(Math.ceil(Math.sqrt(N))).intValue();
		while(N % widthCount > 0)widthCount++;
		int heightCount = N / widthCount;
		
		if(height > width){
			int tempCount = widthCount;
			widthCount = heightCount;
			heightCount = tempCount;
		}
		
		JavaRDD<Chunk> startingChunks = sc.parallelize(startingChunk.split(widthCount, heightCount));
		JavaPairRDD<String, Chunk> chunks = startingChunks.mapToPair(new PairFunction<Chunk, String, Chunk>() {
			public Tuple2<String, Chunk> call(Chunk chunk) { return new Tuple2<String, Chunk>(chunk.key(), chunk); }
		}).partitionBy(new SpatialPartitioner(widthCount, heightCount)).cache();
		
		for(int i = 0; i < new Integer(args[2]); i++){
			JavaPairRDD<String, Chunk> neighbourChunks = chunks.flatMapToPair(new PairFlatMapFunction<Tuple2<String, Chunk>, String, Chunk>() {
			 	public Iterable<Tuple2<String, Chunk>> call(Tuple2<String, Chunk> pair) { return pair._2().getNeighbourPairs(width, height); }
			});
			
			chunks = chunks.cogroup(neighbourChunks).mapValues(new Function<Tuple2<Iterable<Chunk>,Iterable<Chunk>>, Chunk>() {
				public Chunk call(Tuple2<Iterable<Chunk>,Iterable<Chunk>> pair) {
					Chunk result = null;
					for(Chunk chunk : pair._1()){
						result = result == null ? chunk : chunk.combine(result);
					}
					for(Chunk chunk : pair._2()){
						result = result == null ? chunk : chunk.combine(result);
					}
					
					result.merge();
					
					DataManipulationAlgorithm algorithm = new GameOfLife("");
					result.setAlgorithm(algorithm);
					return result.step();				
				}
			});
			
			/*newChunks = newChunks.persist(StorageLevel.MEMORY_ONLY());
			chunks.unpersist();
			chunks = newChunks;*/
			
			//chunks = chunks.persist(StorageLevel.MEMORY_ONLY());
			chunks = chunks.cache();
			
			/*if(i % 10 == 0 && i > 0){
				chunks.checkpoint();
				System.out.println(chunks.count());
			}*/
		}
		
		Chunk result = chunks.values().reduce(new Function2<Chunk, Chunk, Chunk>(){
			public Chunk call(Chunk chunk1, Chunk chunk2) { return chunk1.combine(chunk2); }
		});
		
		result.setAlgorithm(algorithm);
		result.saveToFile(args[1]);
		
		sc.close();
	}
}
