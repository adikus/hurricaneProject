package hoos.project.dummy;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import scala.Tuple2;

public class Chunk implements Serializable, DataStore{
	private static final long serialVersionUID = 5686116991330689263L;
	
	transient private DataManipulationAlgorithm algorithm;
	private HashMap<Tuple2<Integer, Integer>, DataStore> storeMap = new HashMap<Tuple2<Integer, Integer>, DataStore>();
	
	private int x, y, w, h, ix, iy;
	private int sizeFactor = 1;
	private int originalX, originalY, originalW, originalH;
	
	public Chunk(DataManipulationAlgorithm algorithm, DataStore dataStore) {
		this.algorithm = algorithm;
		
		this.algorithm.initialise(dataStore);
		
		this.x = 0;
		this.y = 0;
		this.ix = 0;
		this.iy = 0;
		this.w = dataStore.getWidth();
		this.h = dataStore.getHeight();
		
		this.storeMap.put(coordinate(x, y), dataStore);
		setOriginalValues(0, 0, 0, 0);
	}
	
	public Chunk(int x, int y, int ix, int iy, DataManipulationAlgorithm algorithm, DataStore dataStore) {
		this.x = x;
		this.y = y;
		this.ix = ix;
		this.iy = iy;
		
		this.algorithm = algorithm;
		
		this.w = dataStore.getWidth();
		this.h = dataStore.getHeight();

		this.storeMap.put(coordinate(x, y), dataStore);
		setOriginalValues(x, y, w, h);
	}
	
	private void setOriginalValues(int x, int y, int w, int h){
		this.originalX = x;
		this.originalY = y;
		this.originalW = w;
		this.originalH = h;
	}
	
	private int getBaseWidth() {
		return w;
	}
	
	private int getBaseHeight() {
		return h;
	}
	
	public int getWidth() {
		return w * sizeFactor;
	}
	
	public int getHeight() {
		return h * sizeFactor;
	}
	
	private int getOriginalWidth() {
		return originalW;
	}
	
	private int getOriginalHeight() {
		return originalH;
	}
	
	public int getX() {
		return x;
	}
	
	public int getY() {
		return y;
	}
	
	public int getOriginalX() {
		return originalX;
	}
	
	public int getOriginalY() {
		return originalY;
	}
	
	public HashMap<Tuple2<Integer, Integer>, DataStore> getStoreMap() {
		return this.storeMap;
	}
	
	public void setAlgorithm(DataManipulationAlgorithm algorithm) {
		this.algorithm = algorithm;
	}
	
	public void put(int x, int y, int value) {
		int finalX = x % getBaseWidth();
		int finalY = y % getBaseHeight();
		
		// Put into the first one
		storeMap.values().iterator().next().put(finalX, finalY, value);
	}
	
	public Integer get(int x, int y) {
		int finalX = x % getBaseWidth();
		int finalY = y % getBaseHeight();
		
		if(storeMap.size() == 1) {
			storeMap.values().iterator().next().get(finalX, finalY);
		}
		
		for(Map.Entry<Tuple2<Integer, Integer>, DataStore> entry : storeMap.entrySet()){
			Tuple2<Integer, Integer> key = entry.getKey();
			DataStore dataStore = entry.getValue();
			int sx = key._1() - this.x;
			int sy = key._2() - this.y;
			int width = dataStore.getWidth();
			int height = dataStore.getHeight();
			
			if(finalX >= sx && finalY >= sy && finalX - sx < width && finalY - sy < height){
				if(dataStore.get(finalX - sx, finalY - sy) == null)System.out.println(finalX + "," + finalY + " - " + sx + ", " + sy) ;
				return dataStore.get(finalX - sx, finalY - sy);
			}
		}
		return null;
	}

	public void multiply(int i) {
		this.sizeFactor = i;
	}
	
	public DataStore getSubset(int x1, int y1, int x2, int y2){
		DataStore subset;
		try {
			subset = storeMap.values().iterator().next().getClass().newInstance();
			
			for(int k = x1; k < x2; k++){
				for(int l = y1; l < y2; l++){
					subset.put(k - x1, l - y1, get(k, l));
				}
			}
			
			return subset;
		} catch (InstantiationException | IllegalAccessException e) {
			e.printStackTrace();
		}
		return null;
	}

	public List<Chunk> split(int n) {
		int width = getWidth();
		int height = getHeight();
		
		int widthCount = new Double(Math.ceil(Math.sqrt(n))).intValue();
		while(n % widthCount > 0)widthCount++;
		int heightCount = n / widthCount;
		
		if(height > width){
			int tempCount = widthCount;
			widthCount = heightCount;
			heightCount = tempCount;
		}
		
		double widthPart = new Double(width) / widthCount;
		double heightPart = new Double(height) / heightCount;
		System.out.println("Ideal chunk size: " + widthPart + " x " + heightPart);
		
		ArrayList<Chunk> chunks = new ArrayList<Chunk>();
		
		for(double i = 0; Math.ceil(i) < width; i += widthPart){
			for(double j = 0; Math.ceil(j) < height; j += heightPart){
				int x = new Double(Math.ceil(i)).intValue();
				int y = new Double(Math.ceil(j)).intValue();
				int x2 = Math.min(width, new Double(Math.ceil(i + widthPart)).intValue());
				int y2 = Math.min(height, new Double(Math.ceil(j + heightPart)).intValue());
				
				DataStore chunkData = getSubset(x, y, x2, y2);
				
				Chunk chunk = new Chunk(x, y, (int)Math.round(i/widthPart), (int)Math.round(j/heightPart), algorithm, chunkData);
				System.out.println("New " + x + ", " + y + " -- " + (x2 - x) + "x" + (y2 - y));
				chunks.add(chunk);
			}
		}
		
		return chunks;
	}
	
	public Iterable<Tuple2<String, Chunk>> splitIntoTuples(int width, int height) {
		ArrayList<Tuple2<String, Chunk>> chunkTuples = new ArrayList<Tuple2<String, Chunk>>();
		chunkTuples.add(tuple(key(ix, iy), this));
		
		// Add borders
		if(y - 1     > 0     )chunkTuples.add(createBorderChunk( 0, -1)); // Top 
		if(y + h + 1 < height)chunkTuples.add(createBorderChunk( 0,  1)); // Bottom
		if(x - 1     > 0     )chunkTuples.add(createBorderChunk(-1,  0)); // Left
		if(x + w + 1 < width )chunkTuples.add(createBorderChunk( 1,  0)); // Right
		
		// Add Corners
		if(x - 1     > 0     && y - 1     > 0     )chunkTuples.add(createBorderChunk(-1, -1)); // TopLeft
		if(x + w + 1 < width && y - 1     > 0     )chunkTuples.add(createBorderChunk( 1, -1)); // TopRight
		if(x - 1     > 0     && y + h + 1 < height)chunkTuples.add(createBorderChunk(-1,  1)); // BottomLeft
		if(x + w + 1 < width && y + h + 1 < height)chunkTuples.add(createBorderChunk( 1,  1)); // BottomRight
		
		System.out.println("Split " + key(ix, iy));
		
		return chunkTuples;
	}
	
	private Tuple2<String, Chunk> createBorderChunk(int dx, int dy) {
		boolean corner = dy != 0 && dx != 0;
		int dw = corner || dy == 0 ? 1 : w;
		int dh = corner || dx == 0 ? 1 : h;
		
		int x1 = dx > 0 ? x + w - 1: x;
		int x2 = dx > 0 ? x + w    : x + dw;
		int y1 = dy > 0 ? y + h - 1: y;
		int y2 = dy > 0 ? y + h    : y + dh;
		
		//System.out.println("Border " + " " + (x2 - x1) + "x" + (y2 - y1) + " " + dx + ", " + dy);
		
		DataStore chunkData = getSubset(x1, y1, x2, y2);
		Chunk chunk = new Chunk(x1, y1, ix + dx, iy + dy, algorithm, chunkData);
		return tuple(key(ix + dx, iy + dy), chunk);
	}
	
	public Chunk combine(Chunk chunk) {
		//String initialSigntatures = "Combined " + this.signature() + " with " + chunk.signature();
		System.out.println("Combining " + " " + w + "x" + h + " " + x + ", " + y + " + " + chunk.getWidth() + "x" + chunk.getHeight() + " " + chunk.getX() + ", " + chunk.getY());
		
		if(chunk.getOriginalWidth() * chunk.getOriginalHeight() > originalW * originalH){
			setOriginalValues(chunk.getOriginalX(), chunk.getOriginalY(), chunk.getOriginalWidth(), chunk.getOriginalHeight());
		}
		
		int oldX = x;
		int oldY = y;
		this.x = Math.min(x, chunk.getX());
		this.y = Math.min(y, chunk.getY());
		this.w = Math.max(oldX + w, chunk.getX() + chunk.getWidth()) - x;
		this.h = Math.max(oldY + h, chunk.getY() + chunk.getHeight()) - y;
	
		for(Map.Entry<Tuple2<Integer, Integer>, DataStore> entry : chunk.getStoreMap().entrySet()){
			Tuple2<Integer, Integer> key = entry.getKey();
			DataStore dataStore = entry.getValue();
			this.storeMap.put(key, dataStore);
		}
		
		System.out.println("Result " + key() + " "+ w + "x" + h + " " + x + ", " + y);
		
		//System.out.println(initialSigntatures + " => " + signature() + " ( " + originalSignature() + " )");
		
		return this;
	}
	
	public void merge() {
		DataStore newStore;
		try {
			System.out.println("Merge " + storeMap.size() + " stores");
			newStore = storeMap.values().iterator().next().getClass().newInstance();
			
			for(Map.Entry<Tuple2<Integer, Integer>, DataStore> entry : storeMap.entrySet()){
				Tuple2<Integer, Integer> key = entry.getKey();
				DataStore dataStore = entry.getValue();
				int sx = key._1() - this.x;
				int sy = key._2() - this.y;
				int width = dataStore.getWidth();
				int height = dataStore.getHeight();
				
				for(int i = sx; i < sx + width; i++){
					for(int j = sy; j < sy + height; j++){
						newStore.put(i, j, dataStore.get(i - sx, j - sy));
					}
				}
			}
			
			storeMap.clear();
			System.out.println("New store: " + newStore.getWidth() + "x" + newStore.getHeight());
			this.storeMap.put(coordinate(x, y), newStore);
		} catch (InstantiationException | IllegalAccessException e) {
			e.printStackTrace();
		}
	}
	
	public Chunk step() {
		DataStore newStore;
		try {
			newStore = storeMap.values().iterator().next().getClass().newInstance();
			
			algorithm.step(this, newStore, originalX - x, originalY - y, originalX - x + originalW, originalY - y + originalH);
			
			System.out.println("Step " + key());
			
			return new Chunk(originalX, originalY, ix, iy, algorithm, newStore);
		} catch (InstantiationException | IllegalAccessException e) {
			e.printStackTrace();
		}
		
		return null;
	}
	
	public void saveToFile(String fileName) {
		algorithm.outputToFile(fileName, this);
	}
	
	private Tuple2<Integer, Integer> coordinate(int x, int y) {
		return new Tuple2<Integer, Integer>(x, y);
	}
	
	private Tuple2<String, Chunk> tuple(String key, Chunk chunk) {
		return new Tuple2<String, Chunk>(key, chunk);
	}
	
	public String key(int x, int y) {
		return x + "_" + y;		
	}
	
	public String key() {
		return ix + "_" + iy;		
	}
}
