package hoos.project.dummy;

import org.apache.spark.Partitioner;

public class SpatialPartitioner extends Partitioner {
	private static final long serialVersionUID = 3517776862294784152L;
	
	private int width;
	private int height;
	
	public SpatialPartitioner(int width, int height) {
		this.width = width;
		this.height = height;
	}

	@Override
	public int getPartition(Object key) {
		if(key instanceof String){
			String[] coordinates = ((String) key).split("_");
			int x = new Integer(coordinates[0]);
			int y = new Integer(coordinates[1]);
			
			return (x + y * width) % (width * height);
		}
		return 0;
	}

	@Override
	public int numPartitions() {
		return width * height;
	}

}
