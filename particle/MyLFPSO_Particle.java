package jp.ac.anan_nct.pso.particle;

import jp.ac.anan_nct.pso.particle.Particle;
import jp.ac.anan_nct.pso.function.*;

import java.util.Random;
import org.apache.commons.math3.special.*;

public class MyLFPSO_Particle extends LFPSO_Particle{

    final static int LIMIT = 10;
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
	    trial = 0;
	}else{
	    trial++;
	}
    }

    private int lookAround(Particle[] particles, int index){
	int trial = 0;
	
	for(int i = 0; i < 30; i++){
	    
	    double d = 0;
	    double[] p = particles[i].getCurrentPositions();
	    for(int j = 0; j < DIMENSION; j++){
		d += Math.pow((positions[j]-p[j]),2);
	    }
	    d = Math.sqrt(d);

	    trial += (int)(DIMENSION/d/2); // 01
	    //trial += (int)(DIMENSION/d); // 02
	}
	return trial;
    }

    public void update(Particle gbest, int iter, Particle[] particles, int index){
	if(lookAround(particles, index) < LIMIT){
	    updateVelocity(gbest, iter);
	    updatePosition();
	}else{
	    leavyFlight(gbest);
	    trial = 0;
	}
	updatePBest();
    }
}
