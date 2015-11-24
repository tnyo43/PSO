package jp.ac.anan_nct.pso.pso;

import jp.ac.anan_nct.pso.function.*;
import jp.ac.anan_nct.pso.particle.*;

import java.util.Random;
import java.text.*;
import java.io.*;

class MyLFPSO{
    static Function function;
    
    final static int N = 40; //number of particles
    final static int T = 200000; //number of roops
    final static int ROOP = 30;

    double[] result;
    double[] last;

    MyLFPSO_Particle[] particles;
    MyLFPSO_Particle gbest; //global best

    Random r;
    File file;
    PrintWriter pw;
    
    MyLFPSO(){
	function = new F6();
	
	r = new Random();
	
	result = new double[T];
	last = new double[ROOP];
        init();
    }
    
    void init(){
	particles = new MyLFPSO_Particle[N];
	gbest = null;
	for(int i = 0; i < N; i++){
	    particles[i] = new MyLFPSO_Particle(function);
	    if(gbest == null){
		gbest = particles[i].clone();
	    }
	    else if(particles[i].get_score() > gbest.get_score()){
		gbest = particles[i].clone();
	    }
	}
    }

    void execute(String filename){
	for(int r = 0; r < ROOP; r++){
	    java.util.Calendar calendar = java.util.Calendar.getInstance();
	    System.out.println(calendar.getTime().toString());
	    
	    System.out.println(r);
	    
	    gbest.update_pbest();
	    
	    for(int t = 0; t < T; t++){
		double gb_score = gbest.get_score();
		for(int i = 0; i < N; i++){
		    MyLFPSO_Particle p = particles[i];
		    p.update(gbest, t+1, particles, i);
		    
		    double p_score = p.get_score();
		    if(p_score > gb_score){
			gbest = p.clone();
			gb_score = p_score;
		    }
		}

		result[t] += gbest.get_score();
		System.out.println(gbest.get_score());
	    }

	    init();
	}
	
	file = new File(file_name);
	
	try{
	    pw = new PrintWriter(new BufferedWriter(new FileWriter(file)));
	    for(int i = 0; i < T; i++){
		double r = result[i]/ROOP;
		pw.println((i+1) + " " + r);
	    }
	}catch(Exception e){
	    System.err.println(e);
	    System.exit(0);
	}
	pw.close();	
    }

    private 
    public static void main(String[] args){
	MyLFPSO lfpso = new MyLFPSO();
	
        java.util.Calendar calendar = java.util.Calendar.getInstance();
        System.out.println(calendar.getTime().toString());
	lfpso.execute("MyLFPSO_test1.dat");

	lfpso = new MyLFPSO();
	
	lfpso.execute("MyLFPSO_test2.dat");
	calendar = java.util.Calendar.getInstance();
        System.out.println(calendar.getTime().toString());
    }
}

