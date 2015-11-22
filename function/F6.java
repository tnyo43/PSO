package jp.ac.anan_nct.pso.function;

class F6 implements Function{ //Rastrigin Function

    public double[] get_range(){
	double[] RANGE =  {-5.12, 5.12, 5.12, 2};
	return RANGE;
    }

    public double criterion(double[] positions){
	double sum = 0;
	
	for(double x : positions){
	    sum+= x*x -10*Math.cos(2*Math.PI * x)+10; 
	}
	return -1*sum;
    }
}
