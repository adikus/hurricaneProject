package hoos.project.dummy;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import com.google.common.collect.Lists;

import scala.Tuple2;

import hoos.project.dummy.datastore.ArrayDataStore;
import hoos.project.dummy.gol.GameOfLife;

public class Test {
	
	private static double serializedSize(Chunk chunk) {
		try {
			ByteArrayOutputStream fileOut = new ByteArrayOutputStream();
			ObjectOutputStream out = new ObjectOutputStream(fileOut);
			out.writeObject(chunk);
			out.close();
			fileOut.close();
			String res = fileOut.toString();
			
			return res.length() / 1024.0 / 1024;
		}catch(IOException i) {
			i.printStackTrace();
		}
		return 0;
	}

	public static void main(String[] args) {
		DataManipulationAlgorithm algorithm = new GameOfLife("data/gol.txt");
		Chunk startingChunk = new Chunk(algorithm, new ArrayDataStore());
		
		startingChunk.multiply(50);
		
		int width = startingChunk.getWidth();
		int height = startingChunk.getHeight();
		
		List<Chunk> chunks = startingChunk.split(8);
		
		for(int k = 0; k < 1; k++){
			ArrayList<Tuple2<String, Chunk>> chunkTuples = new ArrayList<Tuple2<String, Chunk>>();
			HashMap<String, ArrayList<Chunk>> groupedChunks = new HashMap<String, ArrayList<Chunk>>();
			
			for(Chunk chunk : chunks){
				System.out.println();
				System.out.println("Size: " + serializedSize(chunk));
				chunkTuples.addAll(Lists.newArrayList(chunk.splitIntoTuples(width, height)));
				System.out.println();
			}
			
			for(Tuple2<String, Chunk> tuple : chunkTuples){
				String key = tuple._1();
				Chunk chunk = tuple._2();
				
				if(groupedChunks.get(key) == null)groupedChunks.put(key, new ArrayList<Chunk>());
				groupedChunks.get(key).add(chunk);
			}
			
			chunks.clear();
			
		    for (ArrayList<Chunk> chunksToCombine : groupedChunks.values()) {
		        Chunk result = null;

				System.out.println();
		    	for(Chunk chunk : chunksToCombine){
		        	result = result == null ? chunk : chunk.combine(result);
		        }

		    	result.merge();
		    	result = result.step();
		    	
		    	chunks.add(result);
		    	
				System.out.println();
		    }
		}		
	    
	    Chunk finalChunk = null;
	    for(Chunk chunk : chunks){
	    	finalChunk = finalChunk == null ? chunk : chunk.combine(finalChunk);
	    }
	    finalChunk.merge();
	    finalChunk.saveToFile("temp/out.txt");
	}
}
