package hoos.project.dummy.datastore;

import hoos.project.dummy.DataStore;

import java.util.HashMap;

import scala.Tuple2;

public class HashDataStore implements DataStore{

	private static final long serialVersionUID = 8557096785566236961L;
	
	private HashMap<Tuple2<Integer, Integer>, Integer> mapData = new HashMap<Tuple2<Integer, Integer>, Integer>();
	private int w = 0;
	private int h = 0;
	
	public HashDataStore(){}

	public void put(int x, int y, int value) {
		mapData.put(coordinate(x, y), value);
		if(x >= w)w = x + 1;
		if(y >= h)h = y + 1;
	}
	
	public Integer get(int x, int y) {
		return mapData.get(coordinate(x, y));
	}
	
	public static Tuple2<Integer, Integer> coordinate(int i, int j){
		return new Tuple2<Integer, Integer>(i, j);
	}

	public int getWidth() {
		return w;
	}

	public int getHeight() {
		return h;
	}
}
