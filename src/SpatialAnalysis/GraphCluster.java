/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package SpatialAnalysis;

import java.util.ArrayList;

/**
 *
 * @author dell
 */
public class GraphCluster {

        String pdbid;
        String id="";
        Integer num_epi =0;
        Integer num_ne =0;
        double persen_epi=0;
        ArrayList<String> id_epi;
        ArrayList<String> id_nonepi;
        boolean epi_cluster;
        double[] features;
        double[] antigen_features;
        
        public GraphCluster (String pdbid){
            this.pdbid=pdbid;
            this.id_epi = new ArrayList();
            this.id_nonepi =  new ArrayList();
        }
        public GraphCluster(String id, boolean epi){
            this.id = id;
            this.id_epi = new ArrayList();
            this.id_nonepi =  new ArrayList();
            if(epi){
                this.num_epi +=1;
            }else{
                this.num_ne +=1;
            }
        }

        public double[] getAntigen_features() {
            return antigen_features;
        }

        public void setAntigen_features(double[] antigen_features) {
            this.antigen_features = antigen_features;
        }
        
        public double[] getFeatures() {
            return features;
        }

        public void setFeatures(double[] features) {
            //System.out.println("set fitur:");
            this.features = features;
           // for(double d: features){
            //    System.out.print(d+ "\t");
            //}
        }
        
        public boolean isEpi_cluster() {
            return epi_cluster;
        }

        public void setEpi_cluster(boolean epi_cluster) {
            this.epi_cluster = epi_cluster;
        }
        
        
        
        public void addCE(String[]ce){
            if(ce[2].equals("e")){
                num_epi +=1;
                id_epi.add(ce[1]);
            }else{
                num_ne +=1;
                id_nonepi.add(ce[1]);
            }
        }
        public void calculate_pe(){
            this.persen_epi = (double)this.num_epi/(this.num_epi + this.num_ne);
        }
        public String getId() {
            return id;
        }

        public void setId(String id) {
            this.id = id;
        }
        public void addCount(boolean epi){
            if(epi){
                this.num_epi +=1;
            }else{
                this.num_ne +=1;
            }
        }
        public Integer getNum_epi() {
            return num_epi;
        }

        public void setNum_epi(Integer num_epi) {
            this.num_epi = num_epi;
        }

        public Integer getNum_ne() {
            return num_ne;
        }

        public void setNum_ne(Integer num_ne) {
            this.num_ne = num_ne;
        }
        public double getProcent(){
            return this.persen_epi;
        }
        
    }
      

