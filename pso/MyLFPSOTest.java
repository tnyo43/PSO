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
	return t*2 + 10;
    }

    private static int getLimit(int l){
	return l*2 + 10;
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
	init(gBest, limit, threshold);
	
	for(int t = 0; t < FEs; t++){
	    for(int i = 0; i < P; i++){
		particles[i].update(gBest, 1.0-((i+1)/FEs), particles, i);
	    }
	    updateGBest();
	}
	return gBest.getScore();
    }

    private void execute(){
	for(int l = 0; l < limitNum; l++){
	    double limit = getLimit(l);
	    
	    for(int t = 0; t < thresholdNum; t++){
		double threshold = getThreshold(t);

		for(int lo = 1; lo <= LOOP; lo++){
		    evaluationValues[l][t] += runEvaluation(l, t);
		}
		evaluationValues[l][t] /= LOOP;
	    }
	}
    }
    
    public static void main(String[] args){
	MyLFPSOTest ml = new MyLFPSOTest(6, partiIdx, 30);
	ml.execute();
    }
}
