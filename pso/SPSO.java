package jp.ac.anan_nct.pso.pso;

import java.util.Random;
import java.text.*;
import java.io.*;
import jp.ac.anan_nct.pso.particle.*;
import jp.ac.anan_nct.pso.function.*;

class SPSO{
    final static int N = 20; //number of particles
    final static int T = 200000; //number of roops
    final static int ROOP = 30;

    Function function;

    double[] result;
    double[] last;

    SPSO_Particle[] particles;
    SPSO_Particle gbest; //global best

    
    Random r;
    File file;
    PrintWriter pw;

    
    SPSO(){
	function = new F6();
	
	r = new Random();

	result = new double[T];
	last = new double[ROOP];
        init();
    }

    void init(){
	particles = new SPSO_Particle[N];
	gbest = null;
	for(int i = 0; i < N; i++){
	    particles[i] = new SPSO_Particle(function);
	    if(gbest == null){
		gbest = particles[i].clone();
	    }
	    else if(particles[i].get_score() > gbest.get_score()){
		gbest = particles[i].clone();
	    }
	}
    }
    
    void execute(String x){
	for(int i = 0; i < ROOP; i++){
	    java.util.Calendar calendar = java.util.Calendar.getInstance();
	    System.out.println(calendar.getTime().toString());
	    
	    System.out.println(i);
	    
	    gbest.update_pbest();
	    
	    for(int t = 0; t < T; t++){
		double gb_score = gbest.get_score();
		for(SPSO_Particle p : particles){
		    p.update(gbest, t);
		    
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
	
	
	String file_name = x;
	file = new File(file_name);
	
	try{
	    pw = new PrintWriter(new BufferedWriter(new FileWriter(file)));
	    for(int i = 0; i < T; i++){
		pw.println((i+1) + " " + result[i]/ROOP);
	    }
	}catch(Exception e){
	    System.err.println(e);
	    System.exit(0);
	}
	pw.close();
	
    }
    
    public static void main(String[] args){
	SPSO spso = new SPSO();

	
        java.util.Calendar calendar = java.util.Calendar.getInstance();
        System.out.println(calendar.getTime().toString());
	spso.execute("SPSO_test1.dat");

	spso = new SPSO();
	
	spso.execute("SPSO_test2.dat");
	calendar = java.util.Calendar.getInstance();
        System.out.println(calendar.getTime().toString());
    }
}

