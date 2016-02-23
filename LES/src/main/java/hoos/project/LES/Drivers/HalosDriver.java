package hoos.project.LES.Drivers;
import hoos.project.LES.Kernels.Halos;
import hoos.project.LES.Kernels.States;

public class HalosDriver{

	public static void main(String[] _args) {
		Halos kernel = new Halos();

		System.out.printf("Initialising kernel..");
		
		kernel.init(150, 150, 90);
		
		long previousTime = kernel.getAccumulatedExecutionTime();
		System.out.println("First execution time: " + previousTime);
		System.out.println("Conversion execution time: " + kernel.getConversionTime());

		System.out.println("Running kernel..");
		
		int iter = 0;
		while(iter < 20){
			for(int i = States.VELNW__BONDV1_INIT_UVW; i <= States.PRESS_BOUNDP; i++){
				kernel.run(i);
				long time = kernel.getAccumulatedExecutionTime();
				System.out.println("Current state: " + i + " duration: " + (time - previousTime));
				System.out.println("Total kernel execution time: " + time);
				previousTime = time;
			}
			
			iter++;
		}
	}
}