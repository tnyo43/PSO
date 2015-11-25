package jp.ac.anan_nct.pso.function;

public class F5 implements Function{ //Schwefel2.26

    private final static double[] RANGE = {-500, 500, -500, 500};
    //RANGE[0] - RANGE[1] is search range
    //RANGE[1] - RANGE[2] is initial range
    
    public double[] getRange(){
	return RANGE;
    }

    public double criterion(double[] positions){
	double sum = 0;

	for(int i = 0; i < positions.length; i++){
	    double x = positions[i];
	    sum+= 418.98288727243369 - x * Math.sin(Math.sqrt(Math.abs(x)));
	}
	
	return -1*sum;
    }
}
