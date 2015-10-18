package hoos.project.dummy;
import scala.Tuple2;
import java.lang.Iterable;

public interface DataSetChunk {
	public abstract Iterable<Tuple2<String, DataSetChunk>> splitIntoTuples();

	public abstract DataSetChunk combineAndCompute(DataSetChunk result, DataSetChunk chunk);
}
