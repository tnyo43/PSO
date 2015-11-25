package jp.ac.anan_nct.pso.pso;

import jp.ac.anan_nct.pso.particle.*;
import jp.ac.anan_nct.pso.function.*;

import java.util.Random;
import java.io.*;

class PSO{
    Function function;
    
    final static int N = 20; //amount of particles
    final static int FEs = 200000; //number of function evaluations
    final static int DIMENSION = 30; //dimension of the search area
    final static int RUN = 30; //number of running the algorithms
    final static Function[] functions = {new F5(), new F6()};
    
    Particle[] particles;
    Particle gBest; //global best

    double[] results; //sum of the score of gBest for each evaluation
    double[] finals; //final score of gBest in each runs
    
    private PSO(){
	function = functions[1];
	
	Random r = new Random();

	init();
	results = new double[FEs];
	finals = new double[RUN];
    }
    
    void init(){
	particles = new Particle[N];
	gBest = null;
	
	for(int i = 0; i < N; i++){
	    particles[i] = new LFPSO_Particle(function, DIMENSION);
	    if(gBest == null){
		gBest = particles[i].clone();
	    }
	    else if(particles[i].getScore() > gBest.getScore()){
		gBest = particles[i].clone();
	    }
	}
    }

    protected void updateGBest(){
	double bestScore = -Double.MAX_VALUE;
	for(Particle p : particles){
	    if(p.getScore() > bestScore){
		bestScore = p.getScore();
		gBest = p.clone();
	    }
	}
    }

    public static double std(double[] items){
	int n = items.length;
	double ave = 0;
	for(double item : items){
	    ave += item;
	}
	ave /= n;
	
	double sum = 0;
	for(double item : items){
	    sum += Math.pow(  item-ave , 2  );
	}

	return Math.sqrt( sum/n );
    }

    void fileWrite(String filename){
	File file = new File(filename);
	try{
	    PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(file)));

	    for(int t = 0; t < FEs; t++){
		double r = results[t]/RUN;
		pw.println((t+1) + " " + r);	
	    }
	    pw.println("\nstd = " + std(results));
	    
	    pw.close();
	}catch(IOException e){
	    System.err.println(e);
	}
    }
    
    protected void execute(){
	for(int r = 0; r < RUN; r++){
	    for(int t = 0; t < FEs; t++){
		for(int i = 0; i < particles.length; i++){
		    particles[i].update(gBest, t, particles, i);
		}
		updateGBest();
		results[t] += gBest.getScore();
		System.out.println((t/1000) + "\t" + gBest.getScore());
	    }
	    finals[r] = gBest.getScore();

	    init();
	}
    }
    
    public static void main(String[] args){
	PSO pso = new PSO();
	pso.execute();
	pso.fileWrite("PSOtest.dat");
    }
}

