/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ComplexFeature;

import java.util.ArrayList;

/**
 *
 * @author OpenGress
 */
public class AAIndex {
    private int size;
    private ArrayList<double[]> aaindex;
    private double[][]aaiAA ;
     //index aaindex yang digunakan pada Qi et all (seppa 2)
    private int[]aaiSun ={27,77,87,103,117,138,147,175,264,309,316,318,353,376,385,429,433,458,485,491,528};
    private int[]aaiWithNanData={473,474,475,476,477,478,479,480,481,482,521,522,524,525};
    private int []aaiRen={7,37,40,111,142,143,144,321,425,426,427,529};

    public double[] getAaiRen(int aa) {
        double[] aaitemp  = this.getAAIndexAA(aa);
        double[]aaiInRen = new double[aaiRen.length];
        for(int i=0;i<aaiRen.length;i++){
            aaiInRen[i]=aaitemp[aaiRen[i]];
        }
        return aaiInRen;
        
        
    }

    
    public double[] getAaiWithoutNanData(int aa) {
        double[]aaiwithoutNan = this.getAAIndexAA(aa);
        for(int i=0; i<aaiWithNanData.length;i++){
            aaiwithoutNan[aaiWithNanData[i]]=0;
        }
        return aaiwithoutNan;
    }
    
    public double[] getAaiSun(int aa) {
        double[] aaitemp  = this.getAAIndexAA(aa);
        double[]aaibysun = new double[aaiSun.length];
        for(int i=0;i<aaiSun.length;i++){
            aaibysun[i]=aaitemp[aaiSun[i]];
        }
        return aaibysun;
    }
    public String printArray1D(double []array1D){
        
        String ss ="";
        for(int i=0;i<array1D.length;i++){
            ss += array1D[i]+",";
        }
        return ss;
    }
    public AAIndex(){
        aaindex = new ArrayList();
    }
    public void setAindex(ArrayList<double[]> aai){
        this.aaindex=aai;
    }
    public void addAAIndexComponent(double[]idx){
        aaindex.add(idx);
    }

    public ArrayList<double[]> getAaindex() {
        return aaindex;
    }
    public double[] getAAIndexAA(int aa){
        if(aa!=-1){
            double[]aaindex= this.aaindex.get(aa);
            return aaindex;
        }else
            return null;
        
            
    }
    public void composeAAI(){
        double[][] aai= new double[544][20];
        double[]aaik;
        for(int i=0;i<544;i++){
            aaik = aaindex.get(i);
            for(int j=0;j<20;j++){
               aai[i][j]=aaik[j];
            }
                
        }
        
        ArrayList<double[]> newAAI = new ArrayList();
        for(int k=0;k<20;k++){
            double[] aaim= new double[544];
            for(int m=0;m<544;m++){
                aaim[m]=aai[m][k];
            }
            newAAI.add(k, aaim);
            
        }
        this.aaindex = newAAI;
    }
    
}
