package jp.ac.anan_nct.pso.function;

public class F3 implements Function{ //Rosenbrock

    private final static double[] RANGE = {-10, 10, -10, 10};
    //RANGE[0] - RANGE[1] is search range
    //RANGE[1] - RANGE[2] is initial range
    
    public double[] getRange(){
	return RANGE;
    }

    public double criterion(double[] positions){
	double sum = 0;

	double xi = positions[0]; 
	for(int i = 0; i < positions.length-1; i++){
	    double x = xi;
	    xi = positions[i+1];
	    
	    sum += 100*Math.pow((xi - x*x),2) + (x-1)*(x-1);
	}

	double x = xi;
	xi = positions[0];

	sum += 100*Math.pow((xi - x*x),2) + (x-1)*(x-1);
	
	return -1*sum;
    }
}
