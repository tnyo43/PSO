package jp.ac.anan_nct.pso.function;

public class F7 implements Function{ //Ackley

    private final static double[] RANGE = {-32, 32, -32, 16};
    //RANGE[0] - RANGE[1] is search range
    //RANGE[1] - RANGE[2] is initial range
    
    public double[] getRange(){
	return RANGE;
    }

    public double criterion(double[] positions){
	double sum1 = 0;
	double sum2 = 0;
	
	for(int i = 0; i < positions.length; i++){
	    double x = positions[i];
	    sum1 += x*x;
	    sum2 += Math.cos(2*Math.PI*x);
	}
	double ave1 = sum1 /positions.length;
	double ave2 = sum2 /positions.length;

	return 20*Math.exp(-0.2*Math.sqrt(ave1)) +Math.exp(ave2) -20 +Math.E;   
    }
}
