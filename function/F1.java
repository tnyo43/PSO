package jp.ac.anan_nct.pso.function;

public class F1 implements Function{ //Sphere

    private final static double[] RANGE = {-100, 100, -100, 50};
    //RANGE[0] - RANGE[1] is search range
    //RANGE[1] - RANGE[2] is initial range
    
    public double[] getRange(){
	return RANGE;
    }

    public double criterion(double[] positions){
	double sum = 0;

	for(double x : positions){
	    sum += x*x;
	}
	
	return -1*sum;
    }
}
