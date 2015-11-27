package jp.ac.anan_nct.pso.pso;

import jp.ac.anan_nct.pso.particle.MyLFPSO_Particle;
import jp.ac.anan_nct.pso.function.*;

import java.util.Random;
import java.io.*;

class MyLFPSOTest extends PSOMain{
    private static final int partiIdx = 3;
    private static final int limitNum = 10;
    private static final int thresholdNum = 10;
    private static final int LOOP = 10;
    
    private MyLFPSO_Particle[] particles;
    private MyLFPSO_Particle gBest; //global best

    private double[][] evaluationValues;

    MyLFPSOTest(int funcNum, int partiIdx, int dimension){
	super(funcNum, partiIdx, dimension);
	RUN = 10;
	evaluationValues = new double[limitNum][thresholdNum];
    }

    private static int getThreshold(int t){
	return t*5 + 10;
    }

    private static int getLimit(int l){
	return l*2 + 5;
    }

    private void init(MyLFPSO_Particle p, int limit, int threshold){
	particles = new MyLFPSO_Particle[P];
	gBest = null;
	
	for(int i = 0; i < P; i++){
	    particles[i] = p.creatNew();
	    particles[i].setLimit(limit);
	    particles[i].setThreshold(threshold);
	    
	    if(gBest == null){
		gBest = particles[i].clone();
	    }
	    else if(particles[i].getScore() > gBest.getScore()){
		gBest = particles[i].clone();
	    }
	}
    }

    private double runEvaluation(int limit, int threshold){
	
	init(new MyLFPSO_Particle(function, DIMENSION), limit, threshold);
	
	for(int t = 0; t < FEs; t++){
	    for(int i = 0; i < P; i++){
		particles[i].update(gBest, 1.0-((i+1)/FEs), particles, i);
	    }
	    updateGBest();
	}
	return gBest.getScore();
    }

    private void fileWrite(){
	File file = new File("MyLFPSOTest.dat");
	try{
	    PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(file)));
	    for(int j = 0; j < limitNum; j++){
		String str = Integer.toString(j+1);
		for(int i = 0; i < thresholdNum; i++){
		    str += " " + evaluationValues[j][i];
		}
		pw.println(str);
	    }
	    pw.close();
	}catch(Exception e){
	    System.err.println(e);
	}
    }

    private void execute(){
	for(int l = 0; l < limitNum; l++){
	    int limit = getLimit(l);
	    
	    for(int t = 0; t < thresholdNum; t++){
		int threshold = getThreshold(t);

		double s = 0;
		for(int lo = 1; lo <= LOOP; lo++){
		    s += runEvaluation(limit, threshold);
		}
		evaluationValues[l][t] = s/LOOP;
		System.out.println("limit = " + limit + ": threshold = " + threshold);
		System.out.println("  " + evaluationValues[l][t]);
	    }
	}
	fileWrite();
    }
    
    public static void main(String[] args){
	MyLFPSOTest ml = new MyLFPSOTest(6, partiIdx, 30);
	ml.execute();
    }
}
