package hoos.project.dummy;

public interface DataManipulationAlgorithm {

	void initialise(DataStore dataStore);

	void step(DataStore dataStore, DataStore newStore, int x1, int y1, int x2, int y2);

	void outputToFile(String fileName, DataStore dataStore);

}
