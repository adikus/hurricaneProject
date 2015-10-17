import org.apache.spark.api.java.*;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.function.PairFlatMapFunction;
import scala.Tuple2;
import java.lang.Iterable;

public class DummyMapReduce {
	public static void main(String[] args) {
		DataSet ds = new DataSet(args[0]);
		ds.initialize();

		SparkConf conf = new SparkConf().setAppName("Dummy Map Reduce");
		JavaSparkContext sc = new JavaSparkContext(conf);

		JavaRDD<DataSetChunk> chunks = sc.parallelize(ds.getChunkArray());
        
		JavaPairRDD<String, DataSetChunk> mappedChunks = chunks.flatMapToPair(new PairFlatMapFunction<DataSetChunk, String, DataSetChunk>() {
		 	public Iterable<Tuple2<String, DataSetChunk>> call(DataSetChunk chunk) { return chunk.splitIntoTuples(); }
		});
	}
}
