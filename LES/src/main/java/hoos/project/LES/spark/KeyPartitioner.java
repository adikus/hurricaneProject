package hoos.project.LES.spark;

import org.apache.spark.Partitioner;

public class KeyPartitioner extends Partitioner {
	private static final long serialVersionUID = 1L;
	
	private int size;
	
	public KeyPartitioner(int size) {
		this.size = size;
	}

	@Override
	public int getPartition(Object key) {
		if(key instanceof Integer){
			return ((Integer)key) % size;
		}
		return 0;
	}

	@Override
	public int numPartitions() {
		return size;
	}

}
