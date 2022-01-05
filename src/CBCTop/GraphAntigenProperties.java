/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CBCTop;

import graph.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;


/**
 *
 * @author dell
 */
public class GraphAntigenProperties {
    private ArrayList<String> vertexs;
    private ArrayList<Edge> edges;
    private double graph_density=0;
    private double degree_correlation=0;
    private ArrayList<Double> degree_centrality;
    private ArrayList<Double> clustering_coef;
    private ArrayList<Double> topologycal_coef;
    private ArrayList<Double> degree_dist;
    private ArrayList<Double> antigen_properties;
    
    public GraphAntigenProperties (ArrayList<Edge> edgelist ){
        this.vertexs = new ArrayList();
        this.edges = new ArrayList();
        HashSet<String> v_hash = new HashSet();
        for(Edge e: edgelist){
            edges.add(e);
            v_hash.add(e.v1);
            v_hash.add(e.v2);
                    
        }
        vertexs = new ArrayList<String>(v_hash);
        System.out.println("vertexs: "+ this.vertexs.size());
        this.graph_density = this.calc_Graph_Density();
        System.out.println("graph_density: "+ this.getGraph_density());
        this.degree_correlation = this.calc_degree_correlation();
        System.out.println("degree correlation: "+ this.degree_correlation);
        this.degree_centrality= calc_degree_centrality();
        /*
        System.out.print("degree centrality: ");
        for(Double d: this.degree_centrality){
            System.out.print(d+"\t");
        }
        System.out.println();
        */
        calc_clustering_coef();
        /*
        System.out.print("clustering coef: ");
        for(Double d: this.clustering_coef){
            System.out.print(d+"\t");
        }
        System.out.println();
        */
        calc_topologycal_coef();
        /*
        System.out.print("topologycal coef: ");
        for(Double d: this.topologycal_coef){
            System.out.print(d+"\t");
        }
        System.out.println();
        */
        this.calc_degree_distribution();
        /*
        System.out.println("degree dist: "+this.degree_dist.toString());
        for(Double d: this.degree_dist){
            System.out.print(d+"\t");
        }
        System.out.println();
        */
        
    }
    
    public double[] createFeatures(){
        double [] features = new double[12];
        for(double d: features){
            d =0;
        }
        
        features[0]= this.graph_density;
        if (this.degree_dist.size()>0){
            features[1]= this.mean(this.degree_dist);
        features[2]=this.variance(degree_dist);
        features[3]=this.median(degree_dist);
        features[4]=this.max(degree_dist);
        }
        if(this.clustering_coef.size()>0){
            features[5]=this.mean(this.clustering_coef);
        features[6]=this.variance(this.clustering_coef);
        features[7]=this.max(this.clustering_coef);
        
        }
        if(this.topologycal_coef.size()>0){
            features[8]=this.mean(this.topologycal_coef);
        features[9]=this.variance(this.topologycal_coef);
        features[10]=this.max(this.topologycal_coef);
        
        }
        
        features[11]=this.degree_correlation;
        return features;
         
    }
    private double variance(ArrayList<Double> feat){
        double ave = mean(feat);
        double sum=0;
        for(int i=0;i<feat.size();i++){
            sum += Math.pow(feat.get(i)-ave,2);
        }
        if(feat.size()==1){
            return 0;
        } else 
            return sum/(feat.size()-1);
        
    }
    private double median(ArrayList<Double> feat){
        int sz = feat.size();
        double med =0;
        if(sz == 1){
            med = feat.get(0);
        }
        else if(feat.size()%2 == 0){
            med = (feat.get((sz-2)/2)+feat.get(sz/2))/2;
        }else{
            med = feat.get((sz-1)/2);
        }
        return med;
    }
    private double mean(ArrayList<Double> feat){
        double sum=0;
        for(Double d: feat){
            sum += d;
        }
        if (feat.size()>0){
            return sum/feat.size();
        }else
            return 0;
        
    }
    private double max(ArrayList<Double> feat){
        double maxs = Double.MIN_VALUE;
        for(Double d: feat){
            if(d>maxs)
                maxs =d;
        }
        if (feat.size()>0){
            return maxs;
        }
        else
            return 0;
        
    }
    public ArrayList<String> getVertexs() {
        return vertexs;
    }

    public ArrayList<Edge> getEdges() {
        return edges;
    }

    public double getGraph_density() {
        return graph_density;
    }

    public double getDegree_correlation() {
        return degree_correlation;
    }

    public ArrayList<Double> getDegree_centrality() {
        return degree_centrality;
    }
    
    //graph density
    public double calc_Graph_Density(){
        //hitung jumlah edge
        int sz_edge = edges.size();
        int sz_v = vertexs.size();
        double pembagi = sz_v*(sz_v -1);
        double density = 0;
        if (pembagi != 0){
                density = 2* sz_edge/pembagi;
        }else
            density =0;
        return density;
    }
    //hasil degree correlation satu nilai
    public double calc_degree_correlation(){
    //seperti penjelasan newman 2000
        int M = this.edges.size();
        double pem_1 =0;
        double pem_2 =0; //pembilang 2 = penyebut 2
        double peny_1=0;
        for (Edge e: edges){
            String v1 = e.v1;
            String v2 = e.v2;
            int ji= this.degreeOfV(vertexs.indexOf(v1));
            int ki = this.degreeOfV(vertexs.indexOf(v2));
            double jiki = ji*ki/M;
            double ji_t_ki = (ji+ki)/2*M;
            double jip2_t_kip2 = (ji*ji + ki*ki)/2*M;
            pem_1+=jiki;
            pem_2+= ji_t_ki;
            peny_1 += jip2_t_kip2;
             
        }
        double peny = peny_1 - pem_2*pem_2;
        double r=0;
        if (peny !=0){
            r = (pem_1 - pem_2*pem_2)/(peny_1-pem_2*pem_2);
        }else
            r=0;
        return r;
    }
    public void calc_degree_distribution(){
    //jumlah node N(k) dengan k= 1,2,..link dibagi jumlah total dari nodes N
    //untuk setiap vertex hitung degreenya
    //dibagi jumlah node
        int v_size = this.vertexs.size();
        ArrayList<Integer> degree_v = new ArrayList();
        for (int i=0;i<v_size;i++){
            int d_v = degreeOfV(i);
            degree_v.add(d_v);            
        }
        ArrayList<Integer> degree_cum= new ArrayList();
        Collections.sort(degree_v, Collections.reverseOrder());
        //System.out.println("degree of v: "+degree_v.toString());
        for(int j=0;j< degree_v.get(0);j++){
            degree_cum.add(0);
        }
        //System.out.println("degree_cum tersedia: "+ degree_cum.size());
        for(int j=0;j<degree_v.size();j++){
            int idx = degree_v.get(j);
            //System.out.println("idx: "+ idx);
            int curr_val = degree_cum.get(idx-1);
            //System.out.println("curr val: "+ curr_val);
            degree_cum.set(idx-1, curr_val+1);
            //System.out.println("update degree cum: "+degree_cum.get(idx-1));
            
        }
        this.degree_dist= new ArrayList();
        for(int j=0;j<degree_cum.size();j++){
            double d = Double.valueOf(degree_cum.get(j))/v_size;
           // System.out.print(d+"\t");
            degree_dist.add(d);
        }
        //System.out.println("degree dist: "+ degree_dist.toString());
    
    }
    //
    public void calc_clustering_coef(){
        clustering_coef = new ArrayList();
        for(int i=0; i< this.vertexs.size();i++){
            double clus_coef = calc_local_clustering_coef(i);
            clustering_coef.add(clus_coef);
        }
    }
    public double calc_local_clustering_coef(Integer idx_v){
    /*
        Kv = degree (jumlah edge pada vertek v)
        Nv = jumlah edge antara neighboor vertek v
        Cc(v) = 2(Nv)/Kv*(Kv-1)
        */
        int K_v = degreeOfV(idx_v);
        String v = vertexs.get(idx_v);
        ArrayList<String> neigh_v = this.getVNeighbors(idx_v);
        //temukan edge antara neighbors
        int Nv=0;
        for(int i=0;i<neigh_v.size();i++){
            for(int j=0;j<neigh_v.size();j++){
                //buat pasangan vertek
                if (i!=j){
                    String v1 = neigh_v.get(i);
                    String v2 = neigh_v.get(j);
                    if(this.isEdgeExistBetween(v1, v2))
                        Nv +=1;
                }
            }                     
        }
        double cc_v=0;
        if(K_v >1)
            cc_v = 2*Nv/(K_v*(K_v -1));
        return cc_v;   
    }
    public boolean isEdgeExistBetween(String v1, String v2){
        for(Edge e: edges){
            if((v1.equals(e.v1)&& v2.equals(e.v2))|| (v1.equals(e.v2)&& v2.equals(e.v1))){
                return true;
            }
        }
        return false;
    }
    public ArrayList<String> getVNeighbors(int idx_v){
        //edge yang menghubungkan v dengan vertek lain
        
        HashSet<String> set_n = new HashSet();
        String v = vertexs.get(idx_v);
        for(Edge e:edges){
            String v1 = e.v1;
            String v2 = e.v2;
            if(v.equals(v1)){
                set_n.add(v2);
            }else if( v.equals(v2)){
                set_n.add(v1);
            }
        }
        ArrayList<String> neigh = new ArrayList(set_n);
        return neigh;
    }
    public int degreeOfV(int idx_v){
        String v = vertexs.get(idx_v);
        int count =0;
        for(Edge e:edges){
                String v1 = e.v1;
                String v2 = e.v2;
                if(v.equals(v1)|| v.equals(v2))    
                    count +=1;
            }
        return count;
    }
    public void calc_topologycal_coef(){
    //stelzl dkk, 2005
    //1. buat daftar node tetangga dan hitung jumlahnya
    ArrayList<ArrayList<String>> neighbors_of_vi = new ArrayList();
        for(int i=0;i< vertexs.size();i++){
            String vi = vertexs.get(i);
            ArrayList<String> vi_and_neigh_vi = this.getVNeighbors(i);
            //sisipkan vertek vi ke kedalam daftar pada posisi 0
            vi_and_neigh_vi.add(0, vi);
            neighbors_of_vi.add(vi_and_neigh_vi);            
        }
        
    //2. buat daftar node dan Arraylist untuk menyimpan vertex yanng dishare dengan tetangga
    //    untuk setiap node ambil tetangganya 
    //    cek node yang dishare dan simpan dalam daftar 
        this.topologycal_coef = new ArrayList();
        for(int j =0; j<neighbors_of_vi.size();j++){
            //untuk setiap tetangga cek node yang dishare dengan tetangga tersebut
            HashMap<String,Integer> shared_nei_list = new HashMap();
            String vi = neighbors_of_vi.get(j).get(0);
            ArrayList<String> tetangga_vi = neighbors_of_vi.get(j);
            tetangga_vi.remove(0);
            int jxy =0;
            for(String ni:tetangga_vi){
                //cek tetangganya
                int idx_ni = vertexs.indexOf(ni);
                ArrayList<String> tetangga_ni = this.getVNeighbors(idx_ni);
                //cek tetangga yang dishare
                ArrayList<String>shared_nei = new ArrayList();
                
                for(int ini=0; ini<tetangga_ni.size();ini++){
                    String nei_ni = tetangga_ni.get(ini);
                    if(tetangga_vi.contains(nei_ni))
                        shared_nei.add(nei_ni);
                }
                jxy += shared_nei.size()+1;
                shared_nei_list.put(ni, shared_nei.size()+1);
            }
            double ave_jxy = jxy/(tetangga_vi.size()*tetangga_vi.size());
            this.topologycal_coef.add(ave_jxy);
            
        }
    }
    public ArrayList<Double> calc_degree_centrality(){
    // hitung jumlah edge pada setiap vertek
    // dibagi (jumlah vertek -1)    
    //1. untuk setiap vertek hitung jumlah edge
            ArrayList<Double> degree = new ArrayList();
            for(int i=0;i<vertexs.size();i++){
                degree.add(0.0);
            }
            for(Edge e:edges){
                String v1 = e.v1;
                String v2 = e.v2;
                int idx_v1 = vertexs.indexOf(v1);
                int idx_v2 = vertexs.indexOf(v2);
                double d1 = degree.get(idx_v1)+1;
                degree.set(idx_v1, d1);
                double d2 = degree.get(idx_v2)+1;
                degree.set(idx_v2, d2);                 
            }
    //2. hitung jumlah vertek -1
            for(int ii=0;ii<degree.size();ii++){
                double n_d = degree.get(ii)/(vertexs.size()-1);
                degree.set(ii, n_d);
            }
        return degree;
    }
    
    
}
