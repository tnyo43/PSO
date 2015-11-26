package jp.ac.anan_nct.pso.particle;

import jp.ac.anan_nct.pso.particle.Particle;
import jp.ac.anan_nct.pso.function.*;

import java.util.Random;
import org.apache.commons.math3.special.*;

public class MyLFPSO_Particle extends LFPSO_Particle{

    final static int LIMIT = 10;
    final static int DISTANCE = 30;
    final static int THRESHOLD = 20
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

    @Override
    public void updatePBest(){
	double score = function.criterion(positions);
	if(score > this.score){
	    this.score = score;
	    pBestPositions = positions.clone();
	}else{
	    
	}
    }

    protected void addTrial(){
	trial++;
    }

    private void lookAround(Particle[] particles, int index){
	double[] p0 = particles[index].getCurrentPositions();

	for(int i = index+1; i < particles.length; i++){
	    double[] p = particles[i].getCurrentPositions();

	    double sum = 0;
	    for(int j = 0; j < DIMENSION; i++){
		sum += Math.pow((p0[j]-p[j]),2);
	    }

	    if(!(sum < THRESHOLD)){
		trial++;
		particles[i].addTrial();
	    }
	}
    }

    public void update(Particle gbest, int iter, Particle[] particles, int index){
	if(trial < LIMIT){
	    updateVelocity(gbest, iter);
	    updatePosition();
	}else{
	    leavyFlight(gbest);
	}
	trial = 0;
	updatePBest();
    }
}
