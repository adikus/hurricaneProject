package hoos.project.LES.HaloExchange;

import java.io.Serializable;

/**
 * @author      Andrej Hoos
 * Class that encapsulates neighbour data array and the direction it came from
 */
public class Neighbour implements Serializable{
	private static final long serialVersionUID = 1L;
	private float[] array;
	private String name;
	
	public Neighbour(float[] array, String name){
		this.array = array;
		this.name = name;
	}
	
	public String getName() {
		return this.name;
	}
	
	public float[] getArray() {
		return this.array;
	}
}
