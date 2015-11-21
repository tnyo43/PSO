package jp.ac.anan_nct.pso.particle;

import java.util.Random;

public abstract class Particle{
    protected static int DIMENSION;
    protected static double[] RANGE;
    protected static double[] INITIAL_RANGE;
    protected double[] positions;
    protected double[] velocities;
    protected double score;

    protected double[] pbest_positions;

    protected Random rand;
    

    public abstract Particle clone();
    
    protected abstract double criterion();
    
    public abstract void update(Particle gbest, int iter);
    
    protected abstract void update_position();
    protected abstract void update_velocity(Particle gbest, int iter);

    public void update_pbest(){
	double score = criterion();
	if(score > this.score){
	    this.score = score;
	    pbest_positions = positions.clone();
	}
    }
    
    public void print_position(){
	System.out.println("current position:");
	for(int i = 0; i < DIMENSION; i++){
	    System.out.println("X" + i + ":" + positions[i]);
	}
    }

    public String position_toString(){
	String pos = "";
	for(int i = 0; i < DIMENSION; i++){
	    pos += ((i+1) + " " + positions[i] + "\n");
	}
	return pos;
    }

    public void print(){
	print_pbest_position();
	print_pbest_score();
    }
    
    public void print_pbest_position(){
	System.out.println("best position:");
	for(int i = 0; i < DIMENSION; i++){
	    System.out.println("X" + i + ":" + pbest_positions[i]);
	}
    }

    public void print_pbest_score(){
	System.out.println("best score:" + score + "\n");
    }
    
    public int get_dimension(){
	return DIMENSION;
    }
    
    public double get_position(int i){
	return positions[i];
    }
    public double[] get_positions(){
	return positions;
    }
    public double get_score(){
	return score;
    }

    public double get_pbest_position(int i){
	return pbest_positions[i];
    }
    public double[] get_pbest_positions(){
	return pbest_positions;
    }
}
