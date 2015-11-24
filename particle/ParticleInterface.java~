package jp.ac.anan_nct.pso.particle;

import jp.ac.anan_nct.pso.function.*;

public interface ParticleInterface{
    public double criterion();
    public double getScore();
    
    public void update(Particle gbest, int iter);
    public void updatePBest();
    
    public void printCurrentPosition();
    public void printPBestPosition();
    public void printScore();
    public void print();
    public String toString();

    public int getDimension();
    
    public double getCurrentPosition(int i);
    public double[] getCurrentPositions();
    public double getPBestPosition(int i);
    public double[] getPBestPositions();
}
