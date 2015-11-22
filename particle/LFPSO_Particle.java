package jp.ac.anan_nct.pso.particle;

import java.util.Random;
import jp.ac.anan_nct.pso.particle.Particle;
import org.apache.commons.math3.special.*;

public class LFPSO_Particle extends Particle{

    final static int LIMIT = 10;
    private int trial;

    private double beta;
    private double sigma_u;
    private double sigma_v;

    
    private double         width;
    private double initial_width;
    Random rand;
    
    private void init(){
	DIMENSION = 30;
	RANGE[0] = -500;
	RANGE[1] =  500;

	INITIAL_RANGE[0] = -500;
	INITIAL_RANGE[1] =  500;
    }
    
    public LFPSO_Particle(){
	rand = new Random();
	
	RANGE = new double[2];
	INITIAL_RANGE = new double[2];
	init();
	trial = 0;
	
	beta = -2*Math.random()+2;
	sigma_u = Math.pow(    (Math.exp(Gamma.logGamma(1+beta))*Math.sin(Math.PI*beta/2))
			       /
			       (Math.exp(Gamma.logGamma((1+beta)/2))*beta*Math.pow(2, ((beta-1)  /   2)  ))
			       ,1/beta);
	sigma_v = 1;
	
	width         =         RANGE[1] -         RANGE[0];
	initial_width = INITIAL_RANGE[1] - INITIAL_RANGE[0];

	positions = new double[DIMENSION];
	for(int i = 0; i < DIMENSION; i++){
	    positions[i] = rand.nextDouble()*(initial_width)+INITIAL_RANGE[0];
	}
	velocities = new double[DIMENSION];

	score = criterion();

	pbest_positions = positions.clone();
    }

    private LFPSO_Particle(LFPSO_Particle original){
	this.positions = original.positions.clone();
	this.velocities = original.velocities.clone();
	this.score = original.get_score();
	this.pbest_positions = original.get_pbest_positions();
    }

    public LFPSO_Particle clone(){
	return new LFPSO_Particle(this);
    }
    
    protected double criterion(){
	double sum = 0;
	
	for(int i = 0; i < DIMENSION; i++){
	    double x = positions[i];

	    
	    //f5
	    sum += 418.98288727243369 - x * Math.sin(Math.sqrt(Math.abs(x)));
	    /*
	    //f6
	    sum+= x*x -10*Math.cos(2*Math.PI * x)+10; 
	    */
	}
	return -1*sum;
    }

    @Override
    public void update_pbest(){
	double score = criterion();
	if(score > this.score){
	    this.score = score;
	    pbest_positions = positions.clone();
	    trial = 0;
	}else{
	    trial++;
	}
    }
    

    private void leavy_flight(Particle gbest){
	
	double u = rand.nextGaussian() * sigma_u;
	double v = rand.nextGaussian() * sigma_v;

	double[] gbest_position = gbest.get_pbest_positions().clone();
	double step = 0.01 * u / (Math.pow(Math.abs(v), 1/beta));
	for(int i = 0; i < DIMENSION; i++){
	    double s = (step * (positions[i] - gbest_position[i]) )%width;
	    double p = positions[i] + s;

	    if     (p > RANGE[1]) positions[i] = p - width;
	    else if(p < RANGE[0]) positions[i] = p + width;
	    else                  positions[i] = p;

	    
	    if(positions[i] > RANGE[1] || positions[i] < RANGE[0]){
		System.out.println("positions[" + i + "] is " + positions[i] + "; s = "+ s);
		System.exit(-1);
	    }
	    velocities[i] = 0;
	}
    }

    @Override
    public void update(Particle gbest, int iter){
	if(trial < LIMIT){
	    update_velocity(gbest, iter);
	    update_position();
	}else{
	    leavy_flight(gbest);
	    trial = 0;
	}
	update_pbest();
    }

    @Override
    protected void update_velocity(Particle gbest, int iter){
	double ro_max = 2;
	
	double rand1 = rand.nextDouble() * ro_max;
	double rand2 = rand.nextDouble() * ro_max;

	double w = 1-iter/200000.0;
	
	for(int i = 0; i < DIMENSION; i++){
	    velocities[i] = w*velocities[i] + rand1*(pbest_positions[i] - positions[i]) + rand2*(gbest.get_position(i) - positions[i]);
	    
	    //if     (velocities[i] > RANGE[1]*0.2) velocities[i] = RANGE[1]*0.2;
	    //else if(velocities[i] < RANGE[0]*0.2) velocities[i] = RANGE[0]*0.2;
	}
    }

    @Override
    protected void update_position(){
	for(int i = 0; i < DIMENSION; i++){
	    double v = velocities[i]%(width*0.2);

	    double p = positions[i] + v;
	    if     (p > RANGE[1]) positions[i] = p - width;
	    else if(p < RANGE[0]) positions[i] = p + width;
	    else                  positions[i] = p;
	    
	    if(positions[i] > RANGE[1] || positions[i] < RANGE[0]){
		System.out.println("positions[" + i + "] is " + positions[i] + "; v = "+ v);
		System.exit(-1);
	    }
	}
    }
}
