import org.apache.commons.math3.special.*;
import java.util.Random;

class test{

    static Random rand = new Random();
    
    public static void main(String[] args){
	for(int i = 0; i < 10000; i++){
	    double v = rand.nextGaussian() * 0.6;
	    System.out.print(v+",");
	}
	System.out.println();
    }
    
    
}
