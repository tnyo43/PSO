package jp.ac.anan_nct.pso.pso;

import jp.ac.anan_nct.pso.particle.*;
import jp.ac.anan_nct.pso.function.*;

import java.util.Random;
import java.io.*;

class PSOMain{
    Function function;
    
    static int N; //amount of particles
    final static int FEs = 200000; //number of function evaluations
    final static int DIMENSION = 30; //dimension of the search area
    final static int RUN = 30; //number of running the algorithms
    final static Function[] functions = {new F1(), new F2(), new F3(), new F4(), new F5(), new F6(), new F7()}; //array of functions, and you can select one designating its index
 
    Particle[] particles;
    Particle gBest; //global best

    double[] results; //sum of the score of gBest for each evaluation
    double[] finals; //final score of gBest in each runs
      
    void init(Particle p){
	particles = new Particle[N];
	gBest = null;
	
	for(int i = 0; i < N; i++){
	    particles[i] = p.creatNew();
	    if(gBest == null){
		gBest = particles[i].clone();
	    }
	    else if(particles[i].getScore() > gBest.getScore()){
		gBest = particles[i].clone();
	    }
	}
    }
    
    private PSOMain(int funcNum, int partiIdx){
	function = functions[funcNum-1];

	Particle p;
	switch(partiIdx){
	case 1:
	    p = new Particle(function, DIMENSION);
	    break;
	case 2:
	    p = new LFPSO_Particle(function, DIMENSION);
	    break;
	case 3:
	    p = new MyLFPSO_Particle(function, DIMENSION);
	    break;
	default:
	    p = new Particle(function, DIMENSION);
	}
	N = p.getAmount();
	
	Random r = new Random();

	init(p);
	results = new double[FEs];
	finals = new double[RUN];
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
	    pw.println("\nstd = " + std(finals));
	    
	    pw.close();
	}catch(IOException e){
	    System.err.println(e);
	}
    }
    
    protected void execute(){
	for(int r = 0; r < RUN; r++){
	    for(int t = 0; t < FEs; t++){
		for(int i = 0; i < particles.length; i++){
		    particles[i].update(gBest, 1.0-(i/FEs), particles, i);
		}
		updateGBest();
		results[t] += gBest.getScore();
		System.out.println((r+1) + "\t" + (t/10000) + "\t" + gBest.getScore());
	    }
	    finals[r] = gBest.getScore();

	    init(gBest);
	}
    }
    
    public static void main(String[] args){
	PSOMain pso = new PSOMain(6, 3);
	//1st argument -> index of the function(1~7)
	//2nd argument -> index of the particle(1~3)
	
	pso.execute();
	pso.fileWrite("PSOtest.dat");
    }
}

