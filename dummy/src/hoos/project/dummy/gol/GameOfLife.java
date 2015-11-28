package hoos.project.dummy.gol;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Scanner;

import hoos.project.dummy.DataManipulationAlgorithm;
import hoos.project.dummy.DataStore;

public class GameOfLife implements DataManipulationAlgorithm {
	
	private String fileName;
	
	public GameOfLife(String fileName) {
		this.fileName = fileName;
	}
	
	@Override
	public void initialise(DataStore dataStore) {
		try {
        	Scanner sc = new Scanner(new File(fileName));
	        String line = "";
	        int j = 0;
	        
	        while (sc.hasNextLine()) {
	        	line = sc.nextLine();
	        	for(int i = 0; i < line.length(); i++){
	        		char c = line.charAt(i);
	        		dataStore.put(i, j, c == 'X' ? 1 : 0);
	        	}
	        	j++;
	        }
			sc.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	@Override
	public void step(DataStore dataStore, DataStore newStore, int x1, int y1, int x2, int y2) {
		for(int i = x1; i < x2; i++){
			for(int j = y1; j < y2; j++) {
				int cellsAlive = 0;
				boolean isAlive = dataStore.get(i, j) > 0;
				
				for(int x = i - 1; x <= i + 1; x++){
					for(int y = j - 1; y <= j + 1; y++){
						Integer alive = dataStore.get(x, y);
						if(alive != null){
							cellsAlive += alive;
						}
					}
				}
				
				if(isAlive && (cellsAlive >= 3 && cellsAlive <= 4)){ newStore.put(i - x1, j - y1, 1); }
				else if(!isAlive && cellsAlive == 3){ newStore.put(i - x1, j - y1, 1); }
				else newStore.put(i - x1, j - y1, 0);
			}
		}
	}
	
	@Override
	public void outputToFile(String fileName, DataStore dataStore) {
		PrintWriter writer;
		
		System.out.println("Write to " + fileName);
		
		try {
			writer = new PrintWriter(fileName, "UTF-8");
			
			for(int j = 0; j < dataStore.getHeight(); j ++){
				String line = "";
				
				for(int i = 0; i < dataStore.getWidth(); i ++){
					Integer val = dataStore.get(i, j);
					line += val == null ? '-' : (val > 0 ? 'X' : '.');
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
	
}
