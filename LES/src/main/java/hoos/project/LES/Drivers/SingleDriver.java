package hoos.project.LES.Drivers;
import hoos.project.LES.Kernels.SingleCPU;

public class SingleDriver{

	public static void main(String[] args) {
		if(args.length < 1){
			System.err.println("Please provide the number of iterations");
			return;
		}
		final int iterationsNum = new Integer(args[0]);
		
		 SingleCPU kernel = new SingleCPU();

		System.out.printf("Initialising kernel..");
		
		kernel.init(150, 150, 90);
		
		System.out.println("First execution time: " + kernel.getAccumulatedExecutionTime());
		System.out.println("Conversion execution time: " + kernel.getConversionTime());

		System.out.printf("Running kernel..");
		
		int iter = 0;
		while(iter < iterationsNum){
			for(int i = 1; i <= 11; i++){
				kernel.run(i);
				System.out.println("Total kernel execution time: " + kernel.getAccumulatedExecutionTime());
			}
			
			iter++;
		}
		
	}
}