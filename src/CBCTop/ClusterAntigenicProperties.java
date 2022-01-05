/*
 * fitur merupakan hasil PCA dari fitur yang digunakan pada paper sebelumnya
 * reduksi dimensi dilakukan 
 * and open the template in the editor.
 */
package CBCTop;

import ComplexStructure.Complex;
import SpatialAnalysis.Cluster;
import java.util.ArrayList;
import jsat.classifiers.DataPoint;

/**
 *
 * @author dell
 */
public class ClusterAntigenicProperties {
    private ArrayList<String> vertexs;
    private int vertexNum;
    private DataPoint p;
    private double[] fiturSum;
    private double[] fiturSquare;
    private double[] fiturAve;
    private double[] fiturMin;
    private double[] fiturMax;

    public DataPoint getP() {
        return p;
    }

    public void setP(DataPoint p) {
        this.p = p;
    }
    
    
    public double[] getfiturResume(){
        double[] resumeFitur = new double[6];
        for(int i=0;i<6;i++){
            resumeFitur[i]=0;
        }
        resumeFitur[5]= vertexNum;
        for (int i=0; i< fiturSum.length;i++){
            
            resumeFitur[0]= resumeFitur[0]+ Math.log(fiturSum[i]);
            resumeFitur[1]= resumeFitur[1]+ Math.log(fiturSquare[i]);
            resumeFitur[2]= resumeFitur[2]+ fiturAve[i];
            resumeFitur[3]= resumeFitur[3]+ fiturMin[i];
            resumeFitur[4]= resumeFitur[4]+ fiturMax[i];          
        }
        System.out.println("resumeFitur: "+resumeFitur[0]+","+resumeFitur[1]+", "+resumeFitur[2]+", "+resumeFitur[3]+", "+resumeFitur[4]+","+ resumeFitur[5]);
        return resumeFitur;
    }
    public ArrayList<String> getVertexs() {
        return vertexs;
    }

    public void setVertexs(ArrayList<String> vertexs) {
        this.vertexs = vertexs;
    }

    public int getVertexNum() {
        return vertexNum;
    }

    public void setVertexNum(int vertexNum) {
        this.vertexNum = vertexNum;
    }

    public double[] getFiturSum() {
        return fiturSum;
    }

    public void setFiturSum(double[] fiturSum) {
        this.fiturSum = fiturSum;
    }

    public double[] getFiturSquare() {
        return fiturSquare;
    }

    public void setFiturSquare(double[] fiturSquare) {
        this.fiturSquare = fiturSquare;
    }

    public double[] getFiturAve() {
        return fiturAve;
    }

    public void setFiturAve(double[] fiturAve) {
        this.fiturAve = fiturAve;
    }

    public double[] getFiturMin() {
        return fiturMin;
    }

    public void setFiturMin(double[] fiturMin) {
        this.fiturMin = fiturMin;
    }

    public double[] getFiturMax() {
        return fiturMax;
    }

    public void setFiturMax(double[] fiturMax) {
        this.fiturMax = fiturMax;
    }
    
    
}
