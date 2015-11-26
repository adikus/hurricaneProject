package hoos.project.dummy.datastore;

import hoos.project.dummy.DataStore;

import java.util.ArrayList;

public class ArrayDataStore implements DataStore{

	private static final long serialVersionUID = 8557096785566236961L;
	
	private ArrayList<ArrayList<Integer>> mapData = new ArrayList<ArrayList<Integer>>();
	
	public ArrayDataStore() {}

	public void put(int x, int y, int value) {
		if(mapData.size() <= x){
			while(mapData.size() <= x) mapData.add(new ArrayList<Integer>(y));
		}
		ArrayList<Integer> column = mapData.get(x);
		
		if(column.size() <= y){
			while(column.size() <= y) column.add(null);
		}
		column.set(y, value);
	}
	
	public Integer get(int x, int y) {
		if(x < 0 || y < 0)return null;
		if(mapData.size() <= x)return null;
		ArrayList<Integer> row = mapData.get(x);
		
		if(row.size() <= y)return null;
		return row.get(y);
	}

	public int getWidth() {
		return mapData.size();
	}

	public int getHeight() {
		if(mapData.size() == 0)return 0;
		return mapData.get(0).size();
	}
}
