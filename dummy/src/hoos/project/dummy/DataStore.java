package hoos.project.dummy;

import java.io.Serializable;

public interface DataStore extends Serializable {	
	public void put(int i, int j, int value);
	public Integer get(int i, int j);
	
	public int getWidth();
	public int getHeight();
}
