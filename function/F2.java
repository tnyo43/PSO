package jp.ac.anan_nct.pso.function;

public class F2 implements Function{ //Sphere

    private final static double[] RANGE = {-10, 10, -10, 5};
    //RANGE[0] - RANGE[1] is search range
    //RANGE[1] - RANGE[2] is initial range
    
    public double[] getRange(){
	return RANGE;
    }

    public double criterion(double[] positions){
	double sum = 0;
	double prod = 1;

	for(double x : positions){
	    sum += Math.abs(x);
	    prod *= x;
	}
	prod = Math.abs(prod);
	
	return -1*(sum + prod);
    }
}
