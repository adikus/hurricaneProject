package hoos.project.LES.Drivers;
import hoos.project.LES.Kernels.Single;

public class SingleDriver{

	public static void main(String[] _args) {
		Single kernel = new Single();

		System.out.printf("Initialising kernel..");
		
		kernel.init(150, 150, 90);
		
		System.out.println("First execution time: " + kernel.getAccumulatedExecutionTime());
		System.out.println("Conversion execution time: " + kernel.getConversionTime());

		System.out.printf("Running kernel..");
		
		for(int i = 1; i <= 11; i++){
			kernel.run(i);
			System.out.println("Total kernel execution time: " + kernel.getAccumulatedExecutionTime());
		}
		
	}
}