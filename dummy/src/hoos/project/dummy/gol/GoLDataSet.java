package hoos.project.dummy.gol;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import scala.Tuple2;

import hoos.project.dummy.DataSet;
import hoos.project.dummy.DataSetChunk;

public class GoLDataSet implements DataSet {
	
	private String fileName;
	private HashMap<Tuple2<Integer, Integer>, Integer> mapData = new HashMap<Tuple2<Integer, Integer>, Integer>();
	private int height;
	private int width;
	
	final private int N = 4;
	
	public GoLDataSet(String fileName) {
		this.fileName = fileName;
	}

	@Override
	public void initialize() {
        try {
			BufferedReader in = new BufferedReader(new FileReader(fileName));
	        String line = "";
	        this.height = 0;
	        this.width = 0;
	        while ((line = in.readLine()) != null) {
	        	for(int i = 0; i < line.length(); i++){
	        		char c = line.charAt(i);
	        		mapData.put(coordinate(this.height, i), c == 'X' ? 1 : 0);
	        		this.width = i;
	        	}
	        	this.height++;
	        }
	        this.width++;
			in.close();
			
			System.out.println("Dataset size: " + width + " x " + height);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	@Override
	public List<DataSetChunk> getChunkArray() {
		int idealSize = (int) Math.round(Math.sqrt(height*width/N));
		System.out.println("Ideal chunk size: " + idealSize + " x " + idealSize);
		
		ArrayList<DataSetChunk> chunks = new ArrayList<DataSetChunk>();
		
		for(int i = 0; i < height; i += idealSize){
			for(int j = 0; j < width; j += idealSize){
				HashMap<Tuple2<Integer, Integer>, Integer> chunkData = new HashMap<Tuple2<Integer, Integer>, Integer>();
				int w = 0;
				int h = 0;
				for(int k = i; k < height && k < i + idealSize; k++){
					for(int l = j; l < width && l < j + idealSize; l++){
						chunkData.put(coordinate(k, l), mapData.get(coordinate(k, l)));
						w = l - j;
						h = k - i;
					}
				}
				DataSetChunk chunk = new GoLDataSetChunk(j, i, ++w, ++h, idealSize, chunkData);
				System.out.println("New " + chunk.signature());
				chunks.add(chunk);
			}
		}
		
		return chunks;
	}
	
	public static Tuple2<Integer, Integer> coordinate(int i, int j){
		return new Tuple2<Integer, Integer>(i, j);
	}
}
