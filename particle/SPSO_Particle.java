package jp.ac.anan_nct.pso.particle;

import java.util.Random;
import jp.ac.anan_nct.pso.particle.Particle;

public class SPSO_Particle extends Particle{

    Random rand;
    
    private void init(){
	DIMENSION = 30;
	RANGE[0] = -5.12;
	RANGE[1] =  5.12;

	INITIAL_RANGE[0] = -5.12;
	INITIAL_RANGE[1] =  5.12;
    }
    
    public SPSO_Particle(){
	rand = new Random();
	
	RANGE = new double[2];
	INITIAL_RANGE = new double[2];
	init();
	
	double         width =         RANGE[1] -         RANGE[0];
	double initial_width = INITIAL_RANGE[1] - INITIAL_RANGE[0];

	positions = new double[DIMENSION];
	for(int i = 0; i < DIMENSION; i++){
	    positions[i] = rand.nextDouble()*(initial_width)+INITIAL_RANGE[0];
	}
	velocities = new double[DIMENSION];

	score = criterion();

	pbest_positions = positions.clone();
    }

    private SPSO_Particle(SPSO_Particle original){
	this.positions = original.positions.clone();
	this.velocities = original.velocities.clone();
	this.score = original.get_score();
	this.pbest_positions = original.get_pbest_positions();
    }

    public SPSO_Particle clone(){
	return new SPSO_Particle(this);
    }

    protected double criterion(){
	double sum = 0;
	
	for(int i = 0; i < DIMENSION; i++){
	    double x = positions[i];

	    /*
	    //f5
	    sum -= x * Math.sin(Math.sqrt(Math.abs(x)));
	    */
	    //f6
	    sum+= x*x -10*Math.cos(2*Math.PI * x)+10; 
	    
	}
	return -1*sum;
    }

    public void update(Particle gbest, int iter){
	update_velocity(gbest, iter);
	update_position();
	update_pbest();
    }

    protected void update_velocity(Particle gbest, int iter){
	double ro_max = 1.1931;
	
	double rand1 = rand.nextDouble() * ro_max;
	double rand2 = rand.nextDouble() * ro_max;

	double w = 0.7213;

	for(int i = 0; i < DIMENSION; i++){
	    velocities[i] = w*velocities[i] + rand1*(pbest_positions[i] - positions[i]) + rand2*(gbest.get_position(i) - positions[i]);
	}
    }

    protected void update_position(){
	for(int i = 0; i < DIMENSION; i++){
	    positions[i] = positions[i] + velocities[i];
	    if(positions[i] < RANGE[0]) positions[i] = RANGE[0];
	    if(positions[i] > RANGE[1]) positions[i] = RANGE[1];
	}
    }
}
