/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package bpred;

import ComplexCalculate.UnboundDSSP;
import ComplexFeature.AAIndex;
import ComplexStructure.*;
import ExtractFromFile.IO;
import StructureAnalysis.ParameterASA;
import StructureAnalysis.StructureBasedAnalysis;
import java.util.ArrayList;

/**
 *
 * @author OpenGress
 */
public class AnalysisOnSurfaceResidue {
    private Boolean epiandexposedAAGenerated = false;
    private ArrayList<Chain> complexchainlist;
    private ArrayList<Aminoacid> exposedAA;
    private ArrayList<Complex> complexList;
    //private ArrayList<Complex> listComplex;
    //analisis nilai cn, hse, qse dan fse residu permukaan 
    //1.dari ekstraksi epitop ambil daftar pdbid dan antigen chain
        //baca file pdb untuk membentuk struktur complex
    //  get dssp dari pdbid dalam daftar
    //2. dari setiap complex, ambil residu permukaan dengan nilai RSA > 0.05
    //3. untuk setiap residu permukaan  epitope hitung fse, qse, hse dan cn
    //4. untuk setiap residu permukaan nonepitope hitung fse, qse, hse dan cn
    //simpan hasilnya dalam file :

    public ArrayList<Aminoacid> getExposedAA() {
        return exposedAA;
    }

    public ArrayList<Chain> getComplexchainlist() {
        return complexchainlist;
    }

   
    
    
    public void generateAAIndexFromEpi(ArrayList<Aminoacid>AAepi){
       double[]aaiAA;
        IO io = new IO();
       AAIndex aai= io.extractAAIndexFromFile();
       aai.composeAAI();
       for(int i=0;i<AAepi.size();i++){
           aaiAA = aai.getAAIndexAA(AAepi.get(i).letterToNum());
           //write to file
       }
       
    }
    public int[] getAANumListOfNonEpiExpRes(){
        
        int[]aanumlist = new int[this.exposedAA.size()];
        //dari daftar exposed AA, 
        //ambil nomor exposedAA
        int numofepi =0;
         for (int i=0;i<this.exposedAA.size();i++){
             if(this.exposedAA.get(i).isAsEpitope()){
             numofepi = numofepi +1;
             aanumlist[i]= -1;
             }
             aanumlist[i] = this.exposedAA.get(i).letterToNum();
         }
         return aanumlist;
         
    }
    
    public ArrayList<Aminoacid> getAAFromEpilist(){
         IO io = new IO();
        ArrayList<Epitope> epilist = io.getEpitopeListFromFile();
        ArrayList<Aminoacid>allAaepi= new ArrayList();
        for(int i=0;i<epilist.size();i++){
            //ambil pdbid dan chain
            String pdbid = epilist.get(i).getPdbid();
            String chain = epilist.get(i).getChain();
            //buat struktur complexnya
            Complex complex = new Complex(pdbid, chain);
            complex.getChainByChianID(chain.charAt(0)).setAAAsEpitope(epilist.get(i).getIdxsequenceAsEpitope());
            ArrayList<Aminoacid> epi = complex.getChainByChianID(chain.charAt(0)).getAAAsEpitope();
            allAaepi.addAll(epi);
        }
        return allAaepi;
    }
    
    public int[] getAANum(){
        
        ArrayList<Aminoacid> aa= this.getAAFromEpilist();
        int[]aanum=new int[aa.size()];
        for(int i=0;i<aa.size();i++){
            aanum[i]=aa.get(i).letterToNum();
            /*
            if(aanum[i]==-1){
                System.out.print("aanum -1: "+aa.get(i).getresidueName());
            }
            * 
            */
        }
        return aanum;
    }
    public ArrayList<NearestNeighboor> distanceAnalysisOnSurfaceResidu(){
        ArrayList<NearestNeighboor> nna= new ArrayList();
        IO io = new IO();
        ArrayList<Epitope> epilist = io.getEpitopeListFromFile();
        for(int i=0;i<epilist.size();i++){
            String pdbid = epilist.get(i).getPdbid();
            String chain = epilist.get(i).getChain();
            //buat struktur complexnya
            Complex complex = new Complex(pdbid, chain);
            complex.getChainByChianID(chain.charAt(0)).setAAAsEpitope(epilist.get(i).getIdxsequenceAsEpitope()); 
            UnboundDSSP unbound_dssp = new UnboundDSSP(complex, chain);
            ArrayList<Aminoacid> exposedaminoacid = unbound_dssp.getExposedAA(complex, chain.charAt(0),1);
            if(exposedaminoacid.size()>0){
            for(int j=0;j<exposedaminoacid.size();j++){
                StructureBasedAnalysis sba = new StructureBasedAnalysis();
                ArrayList<Atom> cnearneighbor= sba.listNearNeighboorAtom(exposedaminoacid.get(j), complex);
                String nid="";
                for(int k=0;k<cnearneighbor.size();k++){
                    nid=nid+";"+cnearneighbor.get(k).getResidueID();
                }
                NearestNeighboor nn=new NearestNeighboor(exposedaminoacid.get(j).isAsEpitope(),cnearneighbor.size(),exposedaminoacid.get(j).getresidueID(),nid);
                nna.add(nn);
            }}
        }
        return nna;
    }
    public void doAnalysisOnSurfaceResidu(){
        //pdbid dan antigenchain dapat diambil dari epilist
        exposedAA= new ArrayList();
        complexList = new ArrayList();
        IO io = new IO();
        ArrayList<Epitope> epilist = io.getEpitopeListFromFile();
        //System.out.print("contoh hasil ekstraksi"+ epilist.get(epilist.size()-1).getPdbid()+epilist.get(epilist.size()-1).getChain());
        complexchainlist = new ArrayList();
        for(int i=0;i<epilist.size();i++){
            //ambil pdbid dan chain
            String pdbid = epilist.get(i).getPdbid();
            String chain = epilist.get(i).getChain();
            //buat struktur complexnya
            Complex complex = new Complex(pdbid, chain);
            //set asam amino sebagai epitope berdasarkan epilist
            //System.out.print("complex: "+complex.getcomplexName()+"\n");
            //cek jumlah epitope dari complex
            //System.out.print("jumlah epitope: "+pdbid+": "+epilist.get(i).getIdxsequenceAsEpitope().size()+"\n");
           /*
            for(int j=0;j<epilist.get(i).getIdxsequenceAsEpitope().size();j++)
                System.out.println(epilist.get(i).getIdxsequenceAsEpitope().get(j));
            */
            complex.getChainByChianID(chain.charAt(0)).setAAAsEpitope(epilist.get(i).getIdxsequenceAsEpitope());
            //cek apakah proses seting berhasil
            
            ArrayList<Aminoacid> epi = complex.getChainByChianID(chain.charAt(0)).getAAAsEpitope();
           //System.out.print("jumlah epitope: "+ epi.size());
            exposedAA.addAll(epi);
            
            //ambil residu permukaan dengan RSA >0.05
         
          //  System.out.println(complex.getAntigenChainIDs());
            //System.out.println(complex.getcomplexName()+"_"+chain);
            UnboundDSSP unbound_dssp = new UnboundDSSP(complex, chain);
           
           ArrayList<Aminoacid> nonepi = unbound_dssp.getExposedAA(complex, chain.charAt(0),1);
           complexchainlist.add(complex.getChainByChianID(chain.charAt(0)));
           exposedAA.addAll(nonepi);
           this.complexList.add(complex);
            
            //hitung nilai fse
            //hitung nilai qse
            //hitung nilai hse
            //hitung nilai cn
           
            if(exposedAA.size()>0){
            for(int j=0;j<exposedAA.size();j++){
                StructureBasedAnalysis sba = new StructureBasedAnalysis();
                //cek disini untuk print koordinat
                int[]efse = sba.calculateEightSubAreaFSE(exposedAA.get(j), complex);
                
                int[]sfse = sba.calculateSixSubAreaFSE(exposedAA.get(j), complex);
                int[] fse= sba.calculateFSE(exposedAA.get(j), complex);
                int[] qse = sba.calculateQSE(exposedAA.get(j), complex);
                int[] hse = sba.calculateHSE(qse);
                int cn = sba.calculateCNFromHSE(hse);
                
                ParameterASA pAsa = new ParameterASA();
                pAsa.setCn(cn);
                pAsa.setFse(fse);
                pAsa.setHse(hse);
                pAsa.setQse(qse);
                pAsa.setEfse(efse);
                pAsa.setSfse(sfse);
                exposedAA.get(j).setpASA(pAsa);
                
                
            }
            //simpan hasilnya dalam file
            String filename = "paramasa.txt";
          // io.writeAnalysisParamASAOfExpossedResidue(exposedAA,filename);
           //epiandexposedAAGenerated = true;
            
           
            }
            
         
        }
        
        
            
    }
    //tulis hasil analisis residu ke layar atau ke file
    
}
