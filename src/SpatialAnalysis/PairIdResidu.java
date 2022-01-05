/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package SpatialAnalysis;

/**
 *
 * @author dell
 */
public class PairIdResidu{
        private String idResiduF;
        private String idResiduS;
        private double weight;

    public void setWeight(double weight) {
        this.weight = weight;
    }
    public String[]getPairIdResidu(){
        String[]idr = new String[2];
        idr[0]=this.idResiduF;
        idr[1]=this.idResiduS;
        return idr;
    }
    
    public String getIdResiduF() {
        return idResiduF;
    }

    public String getIdResiduS() {
        return idResiduS;
    }

    public double getWeight() {
        return weight;
    }
        
        public PairIdResidu(String s1, String s2, double w) {
            this.idResiduF =s1;
            this.idResiduS=s2;
            this.weight = w;
        }
        public String getPairResidu(){
            String p = this.idResiduF + "," + this.idResiduS + ":"+ this.weight;
            return p;
        }
        public String getPairResiduInTab(){
            String p = this.idResiduF + "\t" + this.idResiduS + "\t"+ this.weight;
            return p;
        }
        
    }

