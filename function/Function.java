package jp.ac.anan_nct.pso.function;

interface Function{
    
    public double[] get_range();
    public double criterion(double[] positions);
    
}
