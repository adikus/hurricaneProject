package hoos.project.dummy;

import org.apache.spark.api.java.*;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.function.Function;
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.api.java.function.PairFlatMapFunction;
import scala.Tuple2;
import hoos.project.dummy.gol.GoLDataSet;

import java.lang.Iterable;

public class DummyMapReduce {
	public static void main(String[] args) {
		DataSet ds = new GoLDataSet(args[0]);
		ds.initialize();

		SparkConf conf = new SparkConf().setAppName("Dummy Map Reduce");
		JavaSparkContext sc = new JavaSparkContext(conf);

		JavaRDD<DataSetChunk> chunks = sc.parallelize(ds.getChunkArray());
		
		for(int i = 0; i < new Integer(args[2]); i++){
			JavaPairRDD<String, DataSetChunk> mappedChunks = chunks.flatMapToPair(new PairFlatMapFunction<DataSetChunk, String, DataSetChunk>() {
			 	public Iterable<Tuple2<String, DataSetChunk>> call(DataSetChunk chunk) { return chunk.splitIntoTuples(); }
			});
			
			JavaRDD<DataSetChunk> reducedChunks = mappedChunks.reduceByKey(new Function2<DataSetChunk, DataSetChunk, DataSetChunk>(){
				public DataSetChunk call(DataSetChunk chunk1, DataSetChunk chunk2) { return chunk1.combine(chunk2); }
			}).values();
			
			chunks = reducedChunks.map(new Function<DataSetChunk, DataSetChunk>(){
				public DataSetChunk call(DataSetChunk chunk)  { return chunk.step(); }
			});
		}
		
		DataSetChunk result = chunks.reduce(new Function2<DataSetChunk, DataSetChunk, DataSetChunk>(){
			public DataSetChunk call(DataSetChunk chunk1, DataSetChunk chunk2) { return chunk1.combine(chunk2); }
		});
		
		result.outputToFile(args[1]);
		
		sc.close();
	}
}
