/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ComplexStructure;

import java.util.ArrayList;

/**
 *
 * @author OpenGress
 */
public class Epitope {
    private String pdbid;
    private String chain;
    private String no_chain;
    private String sequence;
    private String epipos;
    private ArrayList<String> sequenceAsEpitope;
    private ArrayList<Integer> idxsequenceAsEpitope;
    private ArrayList<String> idxEpi;

    public ArrayList<String> getIdxEpi() {
        return idxEpi;
    }

    public void setIdxEpi(ArrayList<String> idxEpi) {
        this.idxEpi = idxEpi;
    }
    
    
    public void setPdbid(String pdbid) {
        this.pdbid = pdbid;
    }

    public void setChain(String chain) {
        this.chain = chain;
    }

    
    public String getEpipos() {
        return epipos;
    }

    public void setEpipos(String epipos) {
        this.epipos = epipos;
    }

    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    
    
    public String getSegmentOfEpi(){
        String segmen ="";
        int prev = -1;
        for(int i=0;i<this.idxsequenceAsEpitope.size();i++){
            int j = this.idxsequenceAsEpitope.get(i);
            String temp = this.sequenceAsEpitope.get(i);
            if(prev == -1){
                segmen= temp;
                
            }else if (prev == j-1){
                segmen = segmen + temp;
            }else{
                segmen = segmen + "\n" + temp;
            }
            prev = j;
        }
        return segmen;
    }
    public String getChain() {
        return chain;
    }

    public ArrayList<Integer> getIdxsequenceAsEpitope() {
        return idxsequenceAsEpitope;
    }
    public void printEpitope(){
        System.out.print(this.pdbid+", "+ this.chain+": ");
        System.out.print(this.getSequenceAsEpitope().size()+"\t");
        System.out.print(this.getIdxsequenceAsEpitope().size()+"\n");
        /*
        for(int i=0;i<this.getIdxsequenceAsEpitope().size();i++){
            System.out.print(this.getIdxsequenceAsEpitope().get(i));
            System.out.print("\t");
            System.out.print(this.getSequenceAsEpitope().get(i));
            System.out.print("\t");
        }*/
    }
    public String getPdbid() {
        return pdbid;
    }
    public void setIdxsequenceAsEpitope(ArrayList<Integer> idxepi){
        this.idxsequenceAsEpitope = idxepi;
    }
    public void setSequenceAsEpitope(ArrayList<String> sepi){
        this.sequenceAsEpitope = sepi;
    }
    public ArrayList<String> getSequenceAsEpitope() {
        return sequenceAsEpitope;
    }
    public Epitope(){
    
    }
    public Epitope(String[] epitopefromfile){
        if (epitopefromfile.length==5){
        this.pdbid = epitopefromfile[0];
        this.chain = epitopefromfile[1];
        this.no_chain = epitopefromfile[2];
        this.sequence = epitopefromfile[3];
        this.epipos = epitopefromfile[4];
        extractEpitope(epitopefromfile[3],epitopefromfile[4]);
        }
        else {
        this.pdbid = epitopefromfile[0];
        this.chain = epitopefromfile[1];
        this.sequence = epitopefromfile[2];
        this.epipos = epitopefromfile[3];
        extractEpitope(epitopefromfile[2],epitopefromfile[3]);
        }
    }

    public String getNo_chain() {
        return no_chain;
    }

    public void setNo_chain(String no_chain) {
        this.no_chain = no_chain;
    }
    
    private void extractEpitope(String seq, String asepi){
        ArrayList<String> sequence = new ArrayList();
        ArrayList<Integer> idx = new ArrayList();
        if(seq.length()==asepi.length()){
            for(int i=0;i<asepi.length();i++){
              //  System.out.print(asepi.charAt(i)+"\t");
                if(asepi.charAt(i) == ("1").charAt(0)|asepi.charAt(i)==("E").charAt(0)){
                    
                    sequence.add(String.valueOf(seq.charAt(i)));
                    idx.add(i);
                }
            }
        }
        this.sequenceAsEpitope = sequence;
        this.idxsequenceAsEpitope = idx;
    }
}
