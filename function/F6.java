package jp.ac.anan_nct.pso.function;

public class F6 implements Function{ //Rastrigin Function

    private final static double[] RANGE = {-5.12, 5.12, -5.12, 2.0};
    
    public double[] get_range(){
	return RANGE;
    }

    public double criterion(double[] positions){
	double sum = 0;

	for(int i = 0; i < positions.length; i++){
	    double x = positions[i];
	    sum+= x*x -10*Math.cos(2*Math.PI * x)+10; 
	}
	
	return -1*sum;
    }
}
