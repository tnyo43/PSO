package jp.ac.anan_nct.pso.particle;

import java.util.Random;
import jp.ac.anan_nct.pso.particle.Particle;
import jp.ac.anan_nct.pso.function.*;

public class SPSO_Particle extends Particle{
    
    Random rand;
    
    private void init(Function function){
	double[] range = function.get_range();
	
	DIMENSION = 30;
	RANGE[0] = range[0];
	RANGE[1] = range[1];

	INITIAL_RANGE[0] = range[2];
	INITIAL_RANGE[1] = range[3];
	
	        width =         RANGE[1] -         RANGE[0];
	initial_width = INITIAL_RANGE[1] - INITIAL_RANGE[0];
    }
    
    public SPSO_Particle(Function function){
	this.function = function;
	
	rand = new Random();
	
	RANGE = new double[2];
	INITIAL_RANGE = new double[2];
	init(function);

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

    @Override
    public SPSO_Particle clone(){
	return new SPSO_Particle(this);
    }

    @Override
    public void update(Particle gbest, int iter){
	update_velocity(gbest, iter);
	update_position();
	update_pbest();
    }

    @Override
    protected void update_velocity(Particle gbest, int iter){
	double ro_max = 1.1931;
	
	double rand1 = rand.nextDouble() * ro_max;
	double rand2 = rand.nextDouble() * ro_max;

	double w = 0.7213;
        
        double max = 0;
	for(int i = 0; i < DIMENSION; i++){
	    velocities[i] = w*velocities[i] + rand1*(pbest_positions[i] - positions[i]) + rand2*(gbest.get_position(i) - positions[i]);
	    double v = Math.abs(velocities[i]);
	    
            if(v > max) max = v;
	    /*
            if     (v < -1*width*0.2) velocities[i] = -1*width*0.2;
	    else if(v >    width*0.2) velocities[i] =    width*0.2;
            */
	}

        if(max > RANGE[1]*0.2){
            for(int i = 0; i < DIMENSION; i++){
                velocities[i] *= RANGE[1]*0.2/max;
            }
        }
    }
}
