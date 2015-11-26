package jp.ac.anan_nct.pso.particle;

import jp.ac.anan_nct.pso.particle.Particle;
import jp.ac.anan_nct.pso.function.*;

import java.util.Random;
import org.apache.commons.math3.special.*;

public class MyLFPSO_Particle extends LFPSO_Particle{

    static int LIMIT = 10;
    static int THRESHOLD = 10;
    private int trial;

    private double beta;
    private double sigma_u;
    private double sigma_v;
    
    public MyLFPSO_Particle(Function function, int dimension){
	super(function, dimension);
    }

    private MyLFPSO_Particle(MyLFPSO_Particle original){
	super(original);
    }

    @Override
    public MyLFPSO_Particle clone(){
	return new MyLFPSO_Particle(this);
    }
    
    @Override
    public MyLFPSO_Particle creatNew(){
	return new MyLFPSO_Particle(function, DIMENSION);
    }

    public void setThreshold(int t){
	THRESHOLD = t;
    }

    public void setLimit(int l){
	LIMIT = l;
    }

    @Override
    public void updatePBest(){
	double score = function.criterion(positions);
	if(score > this.score){
	    this.score = score;
	    pBestPositions = positions.clone();
	}else{
	    
	}
    }

    @Override
    protected void addTrial(){
	trial++;
    }

    private void lookAround(Particle[] particles, int index){
	double[] p0 = particles[index].getCurrentPositions();

	for(int i = index+1; i < particles.length; i++){
	    double[] p = particles[i].getCurrentPositions();

	    double sum = 0;
	    for(int j = 0; j < DIMENSION; j++){
		sum += Math.pow((p0[j]-p[j]),2);
	    }

	    if(!(sum < THRESHOLD)){
		trial++;
		particles[i].addTrial();
	    }
	}
    }

    @Override
    public void update(Particle gbest, double w, Particle[] particles, int index){
       	lookAround(particles, index);

	if(trial < LIMIT){
	    updateVelocity(gbest, w, 2);
	    updatePosition();
	}else{
	    leavyFlight(gbest);
	}
	trial = 0;
	updatePBest();
    }
}
