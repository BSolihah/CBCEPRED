/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ComplexFeature;

import ComplexStructure.Aminoacid;
import ComplexStructure.Atom;
import ComplexStructure.Complex;
import StructureAnalysis.StructureBasedAnalysis;

import java.util.ArrayList;
import java.util.HashSet;

/**
 *
 * @author dell
 */
public class AminoAcidPairByCN {
    private Complex complex;
    public AminoAcidPairByCN(Complex c){
        this.complex = c;
    }
    public int calcPairByCN(String resid_a, String resid_b){
        
        Aminoacid a = this.complex.getAminoById(resid_a);
        //System.out.println("a: "+a.getresidueID());
        Aminoacid b = this.complex.getAminoById(resid_b);
        //System.out.println("b: "+b.getresidueID());
        StructureBasedAnalysis sba = new StructureBasedAnalysis();
        ArrayList<Integer> ca_a= sba.ContactIdOfNeighbor(a, this.complex);
        ArrayList<Integer> ca_b= sba.ContactIdOfNeighbor(b, this.complex);
        HashSet<Integer> pair_ca = new HashSet();
        for(int i=0;i<ca_a.size();i++){
            int ca = ca_a.get(i);
            for(int j=0;j<ca_b.size();j++){
                int cb = ca_b.get(j);
                if (ca==cb)
                    pair_ca.add(cb);
            }
        }
        ArrayList<Integer> list_capair = new ArrayList(pair_ca);
        int sum = list_capair.size();
        return sum;
    }
    
    
    
    
}
