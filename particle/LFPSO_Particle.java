package jp.ac.anan_nct.pso.particle;

import jp.ac.anan_nct.pso.particle.Particle;
import jp.ac.anan_nct.pso.function.*;

import java.util.Random;
import org.apache.commons.math3.special.*;

public class LFPSO_Particle extends Particle{
    
    static int LIMIT = 10;
    private int trial;
    
    private double beta;
    private double sigma_u;
    private double sigma_v;

    public LFPSO_Particle(Function function, int dimension){
	super(function, dimension);
	
	trial = 0;
	beta = -2*Math.random()+2; // beta <- (0,2]
	sigma_u = Math.pow(    (Math.exp(Gamma.logGamma(1+beta))*Math.sin(Math.PI*beta/2))
			       /
			       (Math.exp(Gamma.logGamma((1+beta)/2))*beta*Math.pow(2, ((beta-1)  /   2)  ))
			       ,1/beta);
	sigma_v = 1;
	//sigma_u, _v are standard deviations
    }
    
    protected LFPSO_Particle(Particle original){
	super(original);
    }

    @Override
    public LFPSO_Particle clone(){
	return new LFPSO_Particle(this);
    }

    @Override
    public LFPSO_Particle creatNew(){
	return new LFPSO_Particle(function, DIMENSION);
    }

    @Override
    public int getAmount(){
	return 40;
    }

    @Override
    public void updatePBest(){
	double score = function.criterion(positions);
	if(score > this.score){
	    this.score = score;
	    pBestPositions = positions.clone();
	    trial = 0;
	}else{
	    trial++;
	}
    }

    protected void leavyFlight(Particle gBest){
	
	double u = rand.nextGaussian() * sigma_u;
	double v = rand.nextGaussian() * sigma_v;

	double[] gBestPosition = gBest.getPBestPositions().clone();
	double step = 0.01 * u / (Math.pow(Math.abs(v), 1/beta));
	for(int i = 0; i < DIMENSION; i++){
	    double s = (step * (positions[i] - gBestPosition[i]) )%width;
	    double p = positions[i] + s;

	    if     (p > RANGE[1]) positions[i] = p - width;
	    else if(p < RANGE[0]) positions[i] = p + width;
	    else                  positions[i] = p;

	    if(positions[i] > RANGE[1] || positions[i] < RANGE[0]){
		System.out.println("positions[" + i + "] is " + positions[i] + "; s = "+ s);
		System.exit(-1);
	    }
	}
    }

    @Override
    public void update(Particle gbest, int iter, Particle[] particles, int index){
	if(trial < LIMIT){
	    updateVelocity(gbest, iter);
	    updatePosition();
	}else{
	    leavyFlight(gbest);
	    trial = 0;
	}
	updatePBest();
    }

    protected void updateVelocity(Particle gbest, int iter){
	double ro_max = 2;
	
	double rand1 = rand.nextDouble() * ro_max;
	double rand2 = rand.nextDouble() * ro_max;

	double w = 1-iter/200000.0;

	for(int i = 0; i < DIMENSION; i++){
	    velocities[i] = w*velocities[i] + rand1*(pBestPositions[i] - positions[i]) + rand2*(gbest.getPBestPosition(i) - positions[i]);
	    double v = velocities[i];

	    if     (v < -1*width*0.2) velocities[i] = -1*width*0.2;
	    else if(v >    width*0.2) velocities[i] =    width*0.2;
	}
    }

    @Override
    public void printType(){
	System.out.println("Levy Flight PSO");
    }
}
