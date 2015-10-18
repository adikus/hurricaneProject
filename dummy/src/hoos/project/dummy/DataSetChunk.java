package hoos.project.dummy;

import scala.Tuple2;
import java.lang.Iterable;
import java.util.HashMap;

public interface DataSetChunk extends java.io.Serializable {
	public abstract Iterable<Tuple2<String, DataSetChunk>> splitIntoTuples();

	public abstract DataSetChunk combine(DataSetChunk result);
	
	public abstract String signature();

	public abstract int getX();
	public abstract int getY();
	public abstract int getWidth();
	public abstract int getHeight();
	
	public abstract int getOriginalX();
	public abstract int getOriginalY();
	public abstract int getOriginalW();
	public abstract int getOriginalH();

	public abstract HashMap<Tuple2<Integer, Integer>, Integer> getMapData();

	public abstract DataSetChunk step();
	
	public abstract void outputToFile(String fileName);
}
