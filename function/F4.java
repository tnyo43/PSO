package jp.ac.anan_nct.pso.function;

public class F4 implements Function{ //Noise
    
    private final static double[] RANGE = {-1.28, 1.28, -1.28, 0.64};
    //RANGE[0] - RANGE[1] is search range
    //RANGE[1] - RANGE[2] is initial range
    
    public double[] getRange(){
	return RANGE;
    }

    public double criterion(double[] positions){
	double sum = 0;

	for(double x : positions){
	    sum += x*x*x*x + Math.random();
	}
	
	return -1*sum;
    }
}
