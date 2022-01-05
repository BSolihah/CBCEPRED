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
public class PSAI {
    private String psaiId;
    private double[]psaiParam;
    public PSAI(){
        psaiParam = new double[23];
    }

    public String getPsaiId() {
        return psaiId;
    }

    public void setPsaiId(String psaiId) {
        this.psaiId = psaiId;
    }
    public ArrayList<Double> getListPSAI(){
        ArrayList<Double> ps = new ArrayList();
        for(int i=0;i<psaiParam.length;i++){
            ps.add(psaiParam[i]);
        }
        return ps;
    }
    public double[] getPsaiParam() {
        return psaiParam;
    }
    public double getASAFromPSAI(){
        return psaiParam[0];
    }
    public void setPsaiParam(double[] psaiParam) {
        this.psaiParam = psaiParam;
    }
    public String printPSAIParam(){
        String psai="";
        for(int i=0; i<psaiParam.length;i++){
            psai += psaiParam[i]+",";
        }
        return psai;
    }
}
