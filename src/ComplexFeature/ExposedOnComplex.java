/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ComplexFeature;

import ComplexStructure.Aminoacid;
import ComplexStructure.ExposedResidu;
import java.util.ArrayList;

/**
 *
 * @author dell
 */
public class ExposedOnComplex {
    private String pdbid;
    private String chain;
    private ArrayList<ExposedResidu> listExposedR;
    private ArrayList<Aminoacid> listAA;

    public void printEOC(){
        System.out.println(pdbid+"_"+chain+ listExposedR.size());
        for(ExposedResidu e:listExposedR){
            System.out.println(e.getIdResidue()+": "+e.isIsEpitope()+"; ");
        }
        
    }

    public ArrayList<Aminoacid> getListAA() {
        return listAA;
    }

    public void setListAA(ArrayList<Aminoacid> listAA) {
        this.listAA = listAA;
    }
    
    public String getPdbid() {
        return pdbid;
    }

    public void setPdbid(String pdbid) {
        this.pdbid = pdbid;
    }

    public String getChain() {
        return chain;
    }

    public void setChain(String chain) {
        this.chain = chain;
    }

    public ArrayList<ExposedResidu> getListExposedR() {
        return listExposedR;
    }

    public void setListExposedR(ArrayList<ExposedResidu> listExposedR) {
        this.listExposedR = listExposedR;
    }
    
    
}
