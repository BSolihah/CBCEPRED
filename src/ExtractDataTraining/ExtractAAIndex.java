/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ExtractDataTraining;

import ComplexCalculate.UnboundDSSP;
import ComplexFeature.AAIndex;
import ComplexStructure.Aminoacid;
import ComplexStructure.Chain;
import ComplexStructure.Complex;
import ComplexStructure.Epitope;
import ExtractFromFile.IO;
import bpred.AnalysisOnSurfaceResidue;
import java.io.*;
import java.util.ArrayList;

/**
 *
 * @author OpenGress
 */
public class ExtractAAIndex {
   
    public ExtractAAIndex() throws IOException{
         //ekstrak aaindex
    
    //
     ExtractDataBenchmark edb = new ExtractDataBenchmark();
        ArrayList<Epitope> epilist = edb.loadBenchmarkData("epitopelist.txt");
    //buat obyek UnboundDSSP untuk setiap pdbid
       // System.out.println("amino acid vector length");
        for(Epitope e:epilist){
            String pdbid = e.getPdbid();
            String chain = e.getChain();
            String chain_num = e.getNo_chain();
            Complex complex = new Complex(pdbid, chain);
            //complex.getChainByChianID(chain.charAt(0)).setEpitopeStatus(e, chain);
            complex.getChainByChianID(chain.charAt(0)).setAAAsEpitope(e);
            int numaa = complex.getChainByChianID(chain.charAt(0)).getaminoacidNum();
            //System.out.println (pdbid+"\t"+complex.getChainByChianID(chain.charAt(0)).getAminoacidVector().size());
           // complex.getChainByChianID(chainid.charAt(0)).getAminoacidVector().size();
          //  System.out.println("panjang sequence"+ e.getSequence().length());
           // complex.setEpitopeStatus(e, chainid.charAt(0));
            UnboundDSSP unbound_dssp = new UnboundDSSP(complex, chain);
            
            this.writeAAIndexToFile(complex, chain.charAt(0), pdbid);
        }
    }
     public void  writeAAIndexToFile(Complex c, char chainId,String filename) throws IOException {
      IO io = new IO();
    AAIndex aai= io.extractAAIndexFromFile();
    aai.composeAAI();
      Chain chain = c.getChainByChianID(chainId);
      
      File fout = new File("src/datatraining/AAI/"+"aai_"+filename+".txt");
      FileOutputStream fos = new FileOutputStream(fout);
      BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fos));
      for (int i=0; i<chain.getaminoacidNum();i++){
          Aminoacid aa= chain.getAminoacidbyorder(i);
          
          
          if(aa.getRsaRose1985()>0.05 || aa.getRsaMiller1987()>0.05 ||aa.getRsaRostSander2005()>0.05 ||aa.getRsaTien2013()>0.05){
          String dataline= c.getComplexName()+"_"+chainId+"_"+aa.getresidueID()+"_"+aa.getresidueName1()+"\t"+aa.isAsEpitope()+"\t"+io.DAtoString(aai.getAAIndexAA(aa.letterToNum()))+"\n";
          bw.write(dataline);
          }
          
          
      }
     
      bw.close();
  }
}
