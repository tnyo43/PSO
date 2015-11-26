package jp.ac.anan_nct.pso.function;

public class F9 implements Function{ //Penalized1
    
    private final static double[] RANGE = {-50, 50, -50, 25};
    //RANGE[0] - RANGE[1] is search range
    //RANGE[1] - RANGE[2] is initial range
    
    public double[] getRange(){
	return RANGE;
    }

    public double criterion(double[] positions){
	double sum = 0;
	double prod = 1;

	for(int i = 0; i < positions.length; i++){
	    double x = positions[i];

	    sum += x*x;
	    prod *= Math.cos(x/Math.sqrt(i));
	}
	
	return -( 1/4000*sum -prod +1);
    }
}
