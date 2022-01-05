/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package SpatialAnalysis;

import ComplexStructure.PredictedEpitope;
import java.util.ArrayList;

/**
 *
 * @author dell
 */
public class GraphPredictedEpitope {
        private ArrayList<PredictedEpitope> peList;
        private ArrayList<PairIdResidu> edgeList;
       
        public GraphPredictedEpitope(PairOfPredictedEpitope pope) {
            this.peList = new ArrayList();
            peList.add(pope.getPe1());
            peList.add(pope.getPe2());
            this.edgeList = new ArrayList();
            PairIdResidu pir = new PairIdResidu(pope.getPe1().getIdResidue(),pope.getPe2().getIdResidue(),pope.getWeight());
            edgeList.add(pir);            
        }
        public boolean isAromaticCluster(){
            boolean aromatic=false;
            for(PredictedEpitope pe: peList){
                String residuName= pe.getAmino().getresidueName();
                 //cek aromatik pada cluster
                //F or W or Y or H
                if(residuName.equals("F")||residuName.equals("PHE")||residuName.equals("W")||residuName.equals("TRP")||residuName.equals("Y")||residuName.equals("TYR")||residuName.equals("H")||residuName.equals("HIS")){
                    aromatic=true;
                }
    
                
            }
            return aromatic;
        }
        public ArrayList<PredictedEpitope> getPeList() {
            return peList;
        }       
        public ArrayList<PairIdResidu> getEdgeList() {
            return edgeList;
        }
        public void addEdge(PairIdResidu pir){
            this.edgeList.add(pir);
        }
        public void addPe(PredictedEpitope pe){
            this.peList.add(pe);
        }
        public void addSubGraph(GraphPredictedEpitope gpe){
            
                this.peList.addAll(gpe.getPeList());
                this.edgeList.addAll(gpe.getEdgeList());
            
            
        }
        public void printGraphPE(){
            System.out.println("daftar pe:");
            for(int i=0;i<this.peList.size();i++){
                System.out.print(this.getPeList().get(i).getIdResidue()+",");
            }
            System.out.println("daftar edge");
            for(int i=0;i<this.getEdgeList().size();i++){
                System.out.println(this.getEdgeList().get(i).getPairResidu());
            }
            
        }
        public boolean insertNewPair(PairOfPredictedEpitope pope){
            boolean inserted = false;
            //jika salah satu node ada pada daftar node maka tambahkan node yang belum ada dan tambahkan edge
            //jika keduanya ada maka tambahkan edge saja
            if(peList.contains(pope.getPe1())&&peList.contains(pope.getPe2())) {
                //tambahkan edge
                inserted = true;
            }else if (peList.contains(pope.getPe1())){
                inserted = true;
            }
            else if(peList.contains(pope.getPe2())){
                inserted = true;
            }
            return inserted;
        }
    }
