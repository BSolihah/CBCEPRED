/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CBCTop;

import ComplexFeature.ExposedOnComplex;
import ComplexStructure.ExposedResidu;
import ComplexStructure.PredictedEpitope;
import SpatialAnalysis.GraphPredictedEpitope;
import SpatialAnalysis.PairOfPredictedEpitope;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author dell
 */
public class TestSpatialClusterAnalysis {
     public static void main(String[] args) throws IOException {
         String list_file ="list_antigen_test_seppa3.txt";
         SpatialClusterAnalysis sca = new SpatialClusterAnalysis (list_file);
         System.out.println("mulai klasifikasi");
         sca.classifyAllEOC();
         ArrayList<ExposedOnComplex>eocList = sca.getEocList();
         //ArrayList<ExposedResidu> erList = sca.getEocList().get(0).getListExposedR();
         System.out.println("eocList"+ eocList.size());
          BufferedWriter bw = null;
         FileWriter fw;
         try{
              for(ExposedOnComplex eoc: eocList){
                String writeToFile = eoc.getPdbid()+"_"+eoc.getChain()+"\n";
                System.out.println(writeToFile);
             String filename=eoc.getPdbid()+"_"+eoc.getChain();
            // fw = new FileWriter(new File("src/dataset/analysis/cluster/"+filename+"_"+"cluster_anal"+".txt"),true);
             fw = new FileWriter(new File("src/dataset/seppa3/test/cluster/"+"cluster_anal_all"+".txt"),true);
              bw = new BufferedWriter(fw);
            ArrayList<PredictedEpitope> peList= sca.extractSpatialFeatureOnResiduTpOrFp(eoc.getListExposedR());
            writeToFile +=  "jumlah: "+ peList.size()+"\n";
           ArrayList<PairOfPredictedEpitope> listPe= sca.distanceAnalysis(peList);
           writeToFile+= "jumlah PPE"+ listPe.size()+"\n";
           bw.write(writeToFile);
           System.out.println(writeToFile);
           writeToFile="";
           for(PairOfPredictedEpitope ppe: listPe){
               writeToFile += ppe.getPe1().getIdResidue()+", "+ppe.getPe2().getIdResidue()+","+
                       ppe.getPe1().getPredictedState()+","+ppe.getPe2().getPredictedState()+":"+String.valueOf(ppe.getWeight())+"\n";
           }
           bw.write(writeToFile);
           System.out.println(writeToFile);
           ArrayList<GraphPredictedEpitope> gpeList=sca.doSpatialClusterOnPairOfPE(listPe);
           writeToFile= sca.cekSubClusterCondition(gpeList);
           System.out.println("first development phase:");
           System.out.println(writeToFile);
          ArrayList<GraphPredictedEpitope> newGpeList = sca.removeSubClusterWithRareEdge(gpeList);
          
          // System.out.println("jumlah sub custer ="+ gpeListV2.size());
           /*
           for(GraphPredictedEpitope gpe: gpeList){
              System.out.println("aromatic:"+ gpe.isAromaticCluster());
           }*/
            //sca.clusterOnPredictedEpitope(peList);
            writeToFile= sca.cekSubClusterCondition(newGpeList);
            
             bw.write(writeToFile);
             System.out.println("after filtering:");
             System.out.println(writeToFile);
            }
              
             
         }catch(Exception except){
                System.out.println(except.getMessage());
            }
         finally{
             bw.close();
         }
        
             
//System.out.println(er.getIdResidue()+": "+er.isIsEpitope()+" vs "+er.getPredictedAs().mostLikely());
     }
     public static void writeClusterAnalysisResult(String filename, String toWrite) throws IOException{
         BufferedWriter bw;
         FileWriter fw;
         try{
            //   fw = new FileWriter(new File("src/datatraining/seppa3/tes/"+filename+"_"+threshold+".arff"),true);
              fw = new FileWriter(new File("src/dataset/analysis/cluster/"+filename+"_"+"cluster_anal"+".txt"),true);
              bw = new BufferedWriter(fw);
              bw.write(toWrite);
            
             bw.close();
         }     
         catch(Exception except){
                System.out.println(except.getMessage());
            }
     }
     
}
