package jp.ac.anan_nct.pso.function;

public interface Function{
    
    public double[] getRange();
    public double criterion(double[] positions);
    
}
