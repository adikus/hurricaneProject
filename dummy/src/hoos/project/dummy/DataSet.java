package hoos.project.dummy;

import java.util.List;

public interface DataSet {
	public abstract void initialize();

	public abstract List<DataSetChunk> getChunkArray();
}
