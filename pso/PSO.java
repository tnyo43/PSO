package jp.ac.anan_nct.pso.pso;

import java.util.Random;
import java.text.*;
import java.io.*;
import jp.ac.anan_nct.pso.particle.*;

class PSO{
    final static int N = 100; //number of particles
    final static int T = 10000; //number of roops
    final static double x_min = -5, x_max = 5;
    final static double y_min = -5, y_max = 5;

    Particle[] particles;
    Particle gbest; //global best

    
    Random r;
    File file;
    PrintWriter pw;
    
    PSO(){
	r = new Random();
	
	particles = new Particle[N];	
	for(int i = 0; i < N; i++){
	    particles[i] = new Particle(r.nextDouble()*(x_max-x_min)+x_min, r.nextDouble()*(y_max-y_min)+y_min);
	    if(gbest == null){
		gbest = particles[i].clone();
	    }
	    else if(particles[i].get_score() > gbest.get_score()){
		gbest = particles[i].clone();
	    }//initialize gbest

	    //  particles[i].set_area(x_min, x_max, y_min, y_max);

	}
    }
    
    void execute(){
	for(int t = 0; t < T; t++){
	    try{
		NumberFormat nf = NumberFormat.getInstance();
		nf.setMinimumIntegerDigits(4);
		nf.setGroupingUsed(false);
		
		String file_name = "test" + nf.format(t) + ".dat";
		file = new File(file_name);
		pw = new PrintWriter(new BufferedWriter(new FileWriter(file)));
		
		double w = 0.5;
		
		for(Particle p : particles){
		    p.update(gbest, t);
		    
		    pw.println(p.print_position());
		}
		pw.close();
	    }catch(Exception e){
		System.err.println(e);
	    }
	    
	}
	
	gbest.print();
    }
    
    public static void main(String[] args){
	PSO pso = new PSO();
	pso.execute();
    }
}

