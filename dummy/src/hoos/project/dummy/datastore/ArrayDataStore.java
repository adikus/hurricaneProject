package hoos.project.dummy.datastore;

import hoos.project.dummy.DataStore;

import java.util.ArrayList;

public class ArrayDataStore implements DataStore{

	private static final long serialVersionUID = 8557096785566236961L;
	
	private int width = 0;
	private int height = 0;
	
	private ArrayList<ArrayList<Integer>> mapData = new ArrayList<ArrayList<Integer>>();
	
	public ArrayDataStore() {}

	public void put(int x, int y, int value) {
		if(mapData.size() <= x){
			this.width = x + 1;
			while(mapData.size() < width) mapData.add(new ArrayList<Integer>(Math.max(height, y + 1)));
		}
		ArrayList<Integer> column = mapData.get(x);
		
		if(column.size() <= y){
			this.height = y + 1;
			while(column.size() < height) column.add(null);
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
		return width;
	}

	public int getHeight() {
		return height;
	}
}
