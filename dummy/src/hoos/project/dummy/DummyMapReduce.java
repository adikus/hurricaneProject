package hoos.project.dummy;

import org.apache.spark.api.java.*;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.api.java.function.PairFlatMapFunction;
import scala.Tuple2;
import hoos.project.dummy.datastore.ArrayDataStore;
import hoos.project.dummy.gol.GameOfLife;

import java.lang.Iterable;

public class DummyMapReduce {
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

		JavaRDD<Chunk> chunks = sc.parallelize(startingChunk.split(widthCount, heightCount));
		
		for(int i = 0; i < new Integer(args[2]); i++){
			JavaPairRDD<String, Chunk> mappedChunks = chunks.flatMapToPair(new PairFlatMapFunction<Chunk, String, Chunk>() {
			 	public Iterable<Tuple2<String, Chunk>> call(Chunk chunk) { return chunk.splitIntoTuples(width, height); }
			});
			
			JavaRDD<Chunk> reducedChunks = mappedChunks.reduceByKey(new Function2<Chunk, Chunk, Chunk>(){
				public Chunk call(Chunk chunk1, Chunk chunk2) { return chunk1.combine(chunk2); }
			}).values();
			
			chunks = reducedChunks.map(new Function<Chunk, Chunk>(){
				public Chunk call(Chunk chunk) {
					DataManipulationAlgorithm algorithm = new GameOfLife("");
					chunk.setAlgorithm(algorithm);
					return chunk.step();
				}
			});
			
			if(i % 20 == 0 && i > 0){
				chunks.checkpoint();
				System.out.println(chunks.count());
			}
		}
		
		Chunk result = chunks.reduce(new Function2<Chunk, Chunk, Chunk>(){
			public Chunk call(Chunk chunk1, Chunk chunk2) { return chunk1.combine(chunk2); }
		});
		
		result.merge();
		result.setAlgorithm(algorithm);
		result.saveToFile(args[1]);
		
		sc.close();
	}
}
