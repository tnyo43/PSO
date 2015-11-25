package jp.ac.anan_nct.pso.particle;

import jp.ac.anan_nct.pso.function.*;

import java.util.Random;

public class Particle implements ParticleInterface{
    protected static Function function;

    protected static int DIMENSION;
    
    protected static double[] RANGE;
    protected static double[] INITIAL_RANGE;
    protected static double width;
    protected static double initial_width;
    
    protected double[] positions;
    protected double[] velocities;
    protected double score;

    protected double[] pBestPositions;

    protected Random rand;

    public Particle(Function function, int dimension){
	this.function = function;
	
	rand = new Random();
	
	RANGE = new double[2];
	INITIAL_RANGE = new double[2];
	double[] range = function.getRange();
	
	DIMENSION = dimension;
	RANGE[0] = range[0];
	RANGE[1] = range[1];
	
	INITIAL_RANGE[0] = range[2];
	INITIAL_RANGE[1] = range[3];
	
	width =         RANGE[1] -         RANGE[0];
	initial_width = INITIAL_RANGE[1] - INITIAL_RANGE[0];
	
	positions = new double[DIMENSION];
	for(int i = 0; i < DIMENSION; i++){
	    positions[i] = rand.nextDouble()*(initial_width)+INITIAL_RANGE[0];
	}
	velocities = new double[DIMENSION];
	
	score = criterion();
	
	pBestPositions = positions.clone();
    }
    
    protected Particle(Particle original){
	this.positions = original.positions.clone();
	this.velocities = original.velocities.clone();
	this.score = original.getScore();
	this.pBestPositions = original.getPBestPositions();
    }

    @Override
    public Particle clone(){
	return new Particle(this);
    }

    public Particle creatNew(){
	return new Particle(function, DIMENSION);
    }

    public int getAmount(){
	return 10+(int)(2*Math.sqrt(DIMENSION));
    }

    @Override
    public double criterion(){
	return function.criterion(positions);
    }

    @Override
    public void update(Particle gbest, int iter, Particle[] particles, int index){
	updateVelocity(gbest, iter);
	updatePosition();
	updatePBest();
    }

    private void updateVelocity(Particle gbest, int iter){
	double ro_max = 1.1931;
	
	double rand1 = rand.nextDouble() * ro_max;
	double rand2 = rand.nextDouble() * ro_max;

	double w = 0.7213;
        
	for(int i = 0; i < DIMENSION; i++){
	    velocities[i] = w*velocities[i] + rand1*(pBestPositions[i] - positions[i]) + rand2*(gbest.getPBestPosition(i) - positions[i]);
	    double v = velocities[i];
	    
            if     (v < -1*width*0.2) velocities[i] = -1*width*0.2;
	    else if(v >    width*0.2) velocities[i] =    width*0.2;
            
	}	
    }

    protected void updatePosition(){
	for(int i = 0; i < DIMENSION; i++){
	    double p = positions[i] + velocities[i];
	    
	    if     (p > RANGE[1]) positions[i] = p - width;
	    else if(p < RANGE[0]) positions[i] = p + width;
	    else                  positions[i] = p;
	    
	    if(positions[i] > RANGE[1] || positions[i] < RANGE[0]){
		System.out.println("positions[" + i + "] is " + positions[i] + "; v = "+ velocities[i] + "; w = " + width);
		System.exit(-1);
	    }
	}
    }

    @Override
    public void updatePBest(){
	double score = criterion();
	if(score > this.score){
	    this.score = score;
	    pBestPositions = positions.clone();
	}
    }

    @Override
    public void printCurrentPosition(){
	System.out.println("current position:");
	for(int i = 0; i < DIMENSION; i++){
	    System.out.println("X" + i + ":" + positions[i]);
	}
    }

    @Override
    public String toString(){
	String pos = "";
	for(int i = 0; i < DIMENSION; i++){
	    pos += ((i+1) + " " + positions[i] + "\n");
	}
	return pos;
    }

    @Override
    public void print(){
	printPBestPosition();
	printScore();
    }

    @Override
    public void printPBestPosition(){
	System.out.println("best position:");
	for(int i = 0; i < DIMENSION; i++){
	    System.out.println("X" + i + ":" + pBestPositions[i]);
	}
    }

    @Override
    public void printScore(){
	System.out.println("best score:" + score + "\n");
    }

    public void printType(){
	System.out.println("state-of-the-art PSO");
    }

    @Override
    public int getDimension(){
	return DIMENSION;
    }
    @Override
    public double getCurrentPosition(int i){
	return positions[i];
    }
    @Override
    public double[] getCurrentPositions(){
	return positions;
    }
    @Override
    public double getScore(){
	return score;
    }
    @Override
    public double getPBestPosition(int i){
	return pBestPositions[i];
    }
    @Override
    public double[] getPBestPositions(){
	return pBestPositions;
    }
}
