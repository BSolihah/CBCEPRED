/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package StructureAnalysis;

import java.util.ArrayList;

/**
 *
 * @author OpenGress
 */
public class ParameterASA {
    private int[]efse;
    private int[]sfse;
    private int[] fse;
    private int[] qse;
    private int[] hse;
    private int cn;
    private int cnfromHSE;

    public String printCNHSE(){
        String shse = String.valueOf(cn)+","+String.valueOf(hse[0])+","+ String.valueOf(hse[1])+",";
        return shse;
    }
    public String printHSE(){
        String shse = String.valueOf(hse[0])+","+ String.valueOf(hse[1])+",";
        return shse;
    }
    public String printQSE(){
        String sqse="";
        for(int i=0;i<qse.length;i++){
            sqse = sqse+String.valueOf(qse[i])+",";
        }
        return sqse;
    }
    public String printFSE(){
        String sfse="";
        for(int j=0;j<fse.length;j++){
            sfse = sfse+ String.valueOf(fse[j])+",";
        }
        return sfse;
    }
    public String printSFSE(){
        String ssfse="";
        for(int j=0;j<sfse.length;j++){
            ssfse = ssfse+ String.valueOf(sfse[j])+",";
        }
        return ssfse;
    }
    public String printEFSE(){
        String sefse="";
        for(int j=0;j<efse.length;j++){
            sefse = sefse+ String.valueOf(efse[j])+",";
        }
        return sefse;
    }
    
    public ArrayList<Double> getParamASAD(){
        String[] asa = this.printParamASA().split(",");
        ArrayList<Double> d= new ArrayList();
        for(int i=0;i<asa.length;i++){
            d.add(Double.valueOf(asa[i]));
        }
        return d;
    }
    public String  printParamASA(){
        String pasa = String.valueOf(cn)+"," + String.valueOf(this.cnfromHSE)+ ",";
        String shse = String.valueOf(hse[0])+","+ String.valueOf(hse[1])+",";
        String stringqse=""; 
        for(int i=0;i<qse.length;i++){
            stringqse = stringqse+String.valueOf(qse[i])+",";
        }
        String stringfse="";
        for(int j=0;j<fse.length;j++){
            stringfse = stringfse+ String.valueOf(fse[j])+",";
        }
        String stringsfse="";
        for(int j=0;j<sfse.length;j++){
            stringsfse = stringsfse+ String.valueOf(sfse[j])+",";
        }
        String stringefse="";
        for(int j=0;j<efse.length;j++){
            stringefse = stringefse+ String.valueOf(efse[j])+",";
        }
        return pasa+shse+stringqse+stringfse+stringsfse+stringefse;
    }

    public int[] getEfse() {
        return efse;
    }

    public void setEfse(int[] efse) {
        this.efse = efse;
    }

    public int[] getSfse() {
        return sfse;
    }

    public void setSfse(int[] sfse) {
        this.sfse = sfse;
    }
    
    public int getCn() {
        return cn;
    }

    public void setCn(int cn) {
        this.cn = cn;
    }

    public int[] getFse() {
        return fse;
    }

    public void setFse(int[] fse) {
        this.fse = fse;
    }

    public int[] getHse() {
        return hse;
    }

    public void setHse(int[] hse) {
        this.hse = hse;
    }

    public int[] getQse() {
        return qse;
    }

    public void setQse(int[] qse) {
        this.qse = qse;
    }

    public int getCnfromHSE() {
        return cnfromHSE;
    }

    public void setCnfromHSE(int cnfromHSE) {
        this.cnfromHSE = cnfromHSE;
    }
    
    
}
