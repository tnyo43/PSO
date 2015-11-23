package jp.ac.anan_nct.pso.particle;

import jp.ac.anan_nct.pso.particle.Particle;
import jp.ac.anan_nct.pso.function.*;

import java.util.Random;
import org.apache.commons.math3.special.*;

public class MyLFPSO_Particle extends Particle{

    final static int LIMIT = 10;
    private int trial;

    private double beta;
    private double sigma_u;
    private double sigma_v;

    Random rand;
    
    private void init(Function function){
	double[] range = function.get_range();
	
	DIMENSION = 30;
	RANGE[0] = range[0];
	RANGE[1] = range[1];

	INITIAL_RANGE[0] = range[2];
	INITIAL_RANGE[1] = range[3];
	
	width         =         RANGE[1] -         RANGE[0];
	initial_width = INITIAL_RANGE[1] - INITIAL_RANGE[0];
    }
    
    public MyLFPSO_Particle(Function function){
	this.function = function;
	
	rand = new Random();
	
	RANGE = new double[2];
	INITIAL_RANGE = new double[2];
	init(function);
	
	trial = 0;
	beta = -2*Math.random()+2;
	sigma_u = Math.pow(    (Math.exp(Gamma.logGamma(1+beta))*Math.sin(Math.PI*beta/2))
			       /
			       (Math.exp(Gamma.logGamma((1+beta)/2))*beta*Math.pow(2, ((beta-1)  /   2)  ))
			       ,1/beta);
	sigma_v = 1;
	

	positions = new double[DIMENSION];
	for(int i = 0; i < DIMENSION; i++){
	    positions[i] = rand.nextDouble()*(initial_width)+INITIAL_RANGE[0];
	}
	velocities = new double[DIMENSION];

	score = criterion();

	pbest_positions = positions.clone();
    }

    private MyLFPSO_Particle(MyLFPSO_Particle original){
	this.positions = original.positions.clone();
	this.velocities = original.velocities.clone();
	this.score = original.get_score();
	this.pbest_positions = original.get_pbest_positions();
    }

    @Override
    public MyLFPSO_Particle clone(){
	return new MyLFPSO_Particle(this);
    }

    @Override
    public void update_pbest(){
	double score = function.criterion(positions);
	if(score > this.score){
	    this.score = score;
	    pbest_positions = positions.clone();
	    trial = 0;
	}else{
	    trial++;
	}
    }

    private int look_around(Particle[] particles, int index){
	int trial = 0;
	
	for(int i = index+1; i < 30; i++){
	    
	    double d = 0;
	    double[] p = particles[i].get_positions();
	    for(int j = 0; j < DIMENSION; j++){
		d += Math.pow((positions[j]-p[j]),2);
	    }
	    d = Math.sqrt(d);

	    trial += (int)(DIMENSION/d/2);
	}

	return trial;
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
	}
    }

    public void update(Particle gbest, int iter, Particle[] particles, int index){
	if(look_around(particles, index) < LIMIT){
	    update_velocity(gbest, iter);
	    update_position();
	}else{
	    leavy_flight(gbest);
	    trial = 0;
	}
	update_pbest();
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
	    velocities[i] = w*velocities[i] + rand1*(pbest_positions[i] - positions[i]) + rand2*(gbest.get_pbest_position(i) - positions[i]);
	    double v = velocities[i];

	    if     (v < -1*width*0.2) velocities[i] = -1*width*0.2;
	    else if(v >    width*0.2) velocities[i] =    width*0.2;
	}
    }
}