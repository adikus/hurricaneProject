package hoos.project.dummy.gol;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import scala.Tuple2;
import hoos.project.dummy.DataSetChunk;

public class GoLDataSetChunk implements DataSetChunk {

	private static final long serialVersionUID = 596892733849066651L;
	private int x;
	private int y;
	private int height;
	private int width;
	private int idealSize;
	
	private int originalX;
	private int originalY;
	private int originalW;
	private int originalH;
	
	private HashMap<Tuple2<Integer, Integer>, Integer> mapData;

	public GoLDataSetChunk(int x, int y, int w, int h, int s, HashMap<Tuple2<Integer, Integer>, Integer> chunkData) {
		this.x = x;
		this.y = y;
		this.width = w;
		this.height = h;
		this.idealSize = s;
		this.mapData = chunkData;
		
		setOriginalValues(x, y, w, h);
	}

	@Override
	public Iterable<Tuple2<String, DataSetChunk>> splitIntoTuples() {
		ArrayList<Tuple2<String, DataSetChunk>> chunkTuples = new ArrayList<Tuple2<String, DataSetChunk>>();
		chunkTuples.add(tuple(key(x, y), this));
		System.out.println("Main " + signature() + " key: " + key(x, y));
		
		// Add borders
		if(y - 1 > 0)chunkTuples.add(createBorderChunk(x, x + width, y, y + 1, key(x, y - 1))); // Top 
		chunkTuples.add(createBorderChunk(x, x + width, y + height - 1	, y + height, key(x, y + height + 1))); // Bottom
		if(x - 1 > 0)chunkTuples.add(createBorderChunk(x, x + 1, y, y + height, key(x - 1, y))); // Left
		chunkTuples.add(createBorderChunk(x + width - 1, x + width, y, y + height, key(x + width + 1, y))); // Right
		
		// Add Corners
		if(x - 1 > 0 && y - 1 > 0)chunkTuples.add(createBorderChunk(x, x + 1, y, y + 1, key(x - 1, y - 1))); // TopLeft
		if(y - 1 > 0)chunkTuples.add(createBorderChunk(x + width - 1, x + width, y, y + 1, key(x + width + 1, y - 1))); // TopRight
		if(x - 1 > 0)chunkTuples.add(createBorderChunk(x, x + 1, y + height - 1, y + height, key(x - 1, y + height + 1))); // BottomLeft
		chunkTuples.add(createBorderChunk(x + width - 1, x + width, y + height - 1, y + height, key(x + width + 1, y + height + 1))); // BottomRight
		
		return chunkTuples;
	}

	@Override
	public DataSetChunk combine(DataSetChunk chunk) {
		String initialSigntatures = "Combined " + this.signature() + " with " + chunk.signature();
		
		if(chunk.getOriginalW() * chunk.getOriginalH() > originalW * originalH){
			setOriginalValues(chunk.getOriginalX(), chunk.getOriginalY(), chunk.getOriginalW(), chunk.getOriginalH());
		}
		
		this.x = Math.min(x, chunk.getX());
		this.y = Math.min(y, chunk.getY());
		this.width = Math.max(x + width, chunk.getX() + chunk.getWidth()) - x;
		this.height = Math.max(y + height, chunk.getY() + chunk.getHeight()) - y;
		
		Iterator<Tuple2<Integer, Integer>> keySetIterator = chunk.getMapData().keySet().iterator(); 
		while(keySetIterator.hasNext()){ 
			Tuple2<Integer, Integer> key = keySetIterator.next();
			mapData.put(key, chunk.getMapData().get(key)); 
		}
		
		System.out.println(initialSigntatures + " => " + signature() + " ( " + originalSignature() + " )");
		
		return this;
	}

	@Override
	public DataSetChunk step() {
		Iterator<Tuple2<Integer, Integer>> keySetIterator = mapData.keySet().iterator();
		
		HashMap<Tuple2<Integer, Integer>, Integer> newMapData = new HashMap<Tuple2<Integer, Integer>, Integer>();
		
		while(keySetIterator.hasNext()){ 
			Tuple2<Integer, Integer> key = keySetIterator.next();
			
			int y = key._1();
			int x = key._2();
			
			int cellsAlive = 0;
			boolean isAlive = mapData.get(key) > 0;
			
			for(int i = y - 1; i <= y + 1; i++){
				for(int j = x - 1; j <= x + 1; j++){
					Integer alive = mapData.get(GoLDataSet.coordinate(i, j));
					if(alive != null){
						cellsAlive += alive;
					}
				}
			}
			
			if(x >= originalX && y >= originalY && x < originalX + originalW && y < originalY + originalH){				
				if(isAlive && (cellsAlive >= 3 && cellsAlive <= 4)){ newMapData.put(key, 1); }
				else if(!isAlive && cellsAlive == 3){ newMapData.put(key, 1); }
				else newMapData.put(key, 0);
			}
		}
		
		DataSetChunk chunk = new GoLDataSetChunk(originalX, originalY, originalW, originalH, idealSize, newMapData);
		
		System.out.println("Computed " + signature() + " => " + chunk.signature());
		
		return chunk;
	}

	@Override
	public void outputToFile(String fileName) {
		PrintWriter writer;
		try {
			writer = new PrintWriter(fileName, "UTF-8");
			
			for(int i = 0; i < height; i ++){
				String line = "";
				
				for(int j = 0; j < width; j ++){
					line += mapData.get(GoLDataSet.coordinate(i, j)) > 0 ? 'X' : '.';
				}
				
				writer.println(line);
			}
			
			writer.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (UnsupportedEncodingException e) {
			e.printStackTrace();
		}
	}

	@Override
	public int getX() {
		return x;
	}

	@Override
	public int getY() {
		return y;
	}

	@Override
	public int getWidth() {
		return width;
	}

	@Override
	public int getHeight() {
		return height;
	}
	
	@Override
	public int getOriginalX() {
		return originalX;
	}
	
	@Override
	public int getOriginalY() {
		return originalY;
	}
	
	@Override
	public int getOriginalW() {
		return originalW;
	}
	
	@Override
	public int getOriginalH() {
		return originalH;
	}
	
	@Override
	public HashMap<Tuple2<Integer, Integer>, Integer> getMapData() {
		return mapData;
	}
	
	public String key(int x, int y) {
		return Math.round(Math.floor(x/idealSize)) + "_" + Math.round(Math.floor(y/idealSize));		
	}
	
	public String signature() {
		return width + " x " + height + " chunk created at (" + x + "," + y + ")";
	}
	
	public String originalSignature() {
		return originalW + " x " + originalH + " chunk created at (" + originalX + "," + originalY + ")";
	}
	
	private Tuple2<String, DataSetChunk> createBorderChunk(int x1, int x2, int y1, int y2, String key) {
		HashMap<Tuple2<Integer, Integer>, Integer> chunkData = new HashMap<Tuple2<Integer, Integer>, Integer>();
		
		for(int i = y1; i < y2; i++){
			for(int j = x1; j < x2; j++){
				chunkData.put(GoLDataSet.coordinate(i, j), mapData.get(GoLDataSet.coordinate(i, j)));
			}
		}
		
		GoLDataSetChunk chunk = new GoLDataSetChunk(x1, y1, x2 - x1, y2 - y1, idealSize, chunkData);
		System.out.println("New border " + chunk.signature() + " key: " + key);
		return tuple(key, chunk);
	}
	
	private Tuple2<String, DataSetChunk> tuple(String key, DataSetChunk chunk) {
		return new Tuple2<String, DataSetChunk>(key, chunk);
	}
	
	private void setOriginalValues(int x, int y, int w, int h){
		this.originalX = x;
		this.originalY = y;
		this.originalW = w;
		this.originalH = h;
	}
}
