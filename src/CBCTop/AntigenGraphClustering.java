/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CBCTop;




import static CBCTop.gpMCLCluSTree.classifyClusterWithAntigenicFeature;
import ComplexStructure.Complex;
import ExtractFromFile.IO;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import jsat.classifiers.CategoricalData;
import jsat.classifiers.ClassificationDataSet;
import jsat.classifiers.DataPoint;
import jsat.datatransform.Imputer;
import jsat.datatransform.PCA;
import jsat.io.ARFFLoader;
import jsat.linear.DenseVector;
import mcl.MCL;


/**
 *
 * @author dell
 */
public class AntigenGraphClustering {
    public static void main(String[] args) throws IOException{
      // cluster();
      //clusterbycn();
     // clusterbycnw();
      //analysisEpitopDistributionInClusterMCL();
      
      String path_epilist = "src/dataset/gao/epitopes/";
       String ending_epi = "";
        String path_gaocnw_cluster_bf075 = "src/dataset/gao/graph/cnw/mcl/bf0_75/";
        String path_gaocn_cluster_bf075 = "src/dataset/gao/graph/cn/mcl/bf0_75/";
       String path_gao_cluster_bf05 = "src/dataset/gao/graph/mcl/";
        String path_gao_cluster_bf025 = "src/dataset/gao/graph/mcl/bf0_25/";
        String path_gao_cluster_bf075 = "src/dataset/gao/graph/mcl/bf0_75/";
        String ending_gao_cluster_bf05 = ".assign";
     // ArrayList<String> epi_file_list = getEpiFilename(path_gao_cluster_bf05);
     // ArrayList<String> epi_file_list = getEpiFilename(path_gao_cluster_bf025);
      ArrayList<String> epi_file_list = getEpiFilename(path_gao_cluster_bf075);
     // ArrayList<String> epi_file_list = getEpiFilename(path_gaocn_cluster_bf075);
     //  ArrayList<String> epi_file_list = getEpiFilename(path_gaocnw_cluster_bf075);
      
     // ArrayList<String> epi_file_list = getEpiFilename(path_gao_cluster_bf025);
      //list_dist_cluster berisi residu_id + nomor cluster
      // ArrayList <ArrayList<String>> list_dist_cluster = cekDistribusiEpitopPdCluster(epi_file_list, path_gaocnw_cluster_bf075, path_epilist);
       //ArrayList <ArrayList<String>> list_dist_cluster = cekDistribusiEpitopPdCluster(epi_file_list, path_gaocn_cluster_bf075, path_epilist);
     // ArrayList <ArrayList<String>> list_dist_cluster = cekDistribusiEpitopPdCluster(epi_file_list, path_gao_cluster_bf05, path_epilist);
      //ArrayList <ArrayList<String>> list_dist_cluster = cekDistribusiEpitopPdCluster(epi_file_list, path_gao_cluster_bf025, path_epilist);
      ArrayList <ArrayList<String>> list_dist_cluster = cekDistribusiEpitopPdCluster(epi_file_list, path_gao_cluster_bf075, path_epilist);
      //cek nomor cluster dan jumlah anggota  epitop pada cluster
     //   resumeClusterEpi(list_dist_cluster);
        //cek distribusi epi non epi pada cluster
       // ArrayList<ArrayList<String>> epinonepidist= cekEpiNonEpiDistributionInCluster(epi_file_list, path_gao_cluster_bf05, path_epilist);
         //ArrayList<ArrayList<String>> epinonepidist= cekEpiNonEpiDistributionInCluster(epi_file_list, path_gao_cluster_bf025, path_epilist);
         ArrayList<ArrayList<String>> epinonepidist= cekEpiNonEpiDistributionInCluster(epi_file_list, path_gao_cluster_bf075, path_epilist);
        // ArrayList<ArrayList<String>> epinonepidist= cekEpiNonEpiDistributionInCluster(epi_file_list, path_gaocn_cluster_bf075, path_epilist);
       //  ArrayList<ArrayList<String>> epinonepidist= cekEpiNonEpiDistributionInCluster(epi_file_list, path_gaocnw_cluster_bf075, path_epilist);
        // printlist(epinonepidist);
         //System.out.println("jumlah file epilist:"+ epi_file_list.size()+"contoh: "+ epi_file_list.get(0));
         //System.out.println("jumlah epinonepidist:"+ epinonepidist.size());
        //cek prosentase epitop pada claster 
       // ArrayList<ArrayList<Cluster>> proc_epi_in_cluster = checkEpiPercentageInCluster(epinonepidist);
       //int kategori =0;
       ArrayList<Cluster> cluster_list= createListClusterElemen(epi_file_list, epinonepidist);
       //kategori =1;
       //ArrayList<Cluster> cluster_list_nonepi= createListClusterElemen(epi_file_list, epinonepidist,kategori);
       //System.out.println(cluster_list.get(0).pdbid);
       String out_filename = "";
       String graph_path = "src/dataset/gao/graph/";
       String path_ag = "src/dataset/gao/chain/";
       String ending = "_6.0.txt";
       ClassificationDataSet cds = GraphFromClusterList(out_filename, cluster_list, graph_path,ending, path_ag);
       System.out.println( "jumlah cluster"+ cluster_list.size());
       
       classifyClusterWithAntigenicFeature(cds);
       /*
       for(Cluster c: cluster_list){
           if (c.getFeatures()!=null){
                //String s= c.pdbid+ ","+c.getId()+",";
                String s="";
                //System.out.println("panjang fitur:"+c.getFeatures().length);
                for(int i=0;i<c.getFeatures().length;i++){
                    s += c.getFeatures()[i]+",";
                }
                s += c.isEpi_cluster();
                System.out.println( s);
               System.out.println( c.pdbid+ "\t"+c.getId()+"\t"+c.getFeatures().toString());
           }
           
       }
       
       //non_epi
       //GraphFromClusterList(out_filename, cluster_list_nonepi, graph_path,ending);
      */
    }
    
    public static ClassificationDataSet GraphFromClusterList(String out_filename, ArrayList<Cluster> cluster_list, String graph_path, String ending, String path_ag) throws IOException{
        //pca untuk reduksi dimensi antigenic properties
        String filetrain = "src/dataset/andersen_exc_3H42_A0.01.arff";
        File file = new File(filetrain);  
        ClassificationDataSet trainSet = ARFFLoader.loadArffFile(file).asClassificationDataSet(0);
        trainSet.applyTransform(new Imputer(trainSet));//impute missing values in the dataset
        PCA pca = new PCA(trainSet,25);
        
        HashMap<String, Integer> entry = new HashMap<String,Integer>();
            entry.put("true",0);
            entry.put("false",1);

            CategoricalData[] catData = new CategoricalData[1];
            catData[0] = new CategoricalData(2);
            catData[0].setCategoryName("class");
            for(Entry<String, Integer> map: entry.entrySet()){
                catData[0].setOptionName(map.getKey(), map.getValue());
            }
        //file graph
         System.out.println("mulai graph clustering");
         IO io = new IO();
         ArrayList<String> listData = new ArrayList();
         ArrayList<DataPoint> listDP = new ArrayList();
        for(Cluster c: cluster_list){
            //persiapan data categoric dari cluster            
            ArrayList<Double> numeric_feat = new ArrayList();
            ArrayList<Integer> cats = new ArrayList();
            
            String fitur ="";
            String graphfilename = graph_path + c.pdbid.substring(0,4).concat(c.pdbid.substring(c.pdbid.length()-2)).trim()+ending;
            System.out.println("graph file name: "+graphfilename);
            System.out.println("cpdbid chain "+ c.pdbid.substring(c.pdbid.length()-2));
            
            //extract edge from file
            ArrayList<Edge> edge_list = extractEdgeFromFile(graphfilename, c);
           // System.out.println("edge_list"+ edge_list.size());
            if(edge_list.size()>0){
               // System.out.println("ekstraksi fitur");
                GraphAntigenProperties gap = new GraphAntigenProperties (edge_list);
                double[]feature = gap.createFeatures();
                
                for(double d: feature){
                    numeric_feat.add(d);
                    fitur += d+",";
                    System.out.print(d+"\t");
                }
                
                c.setFeatures(feature);
               // System.out.println("fitur cluster: ");
              //  for(double d: c.getFeatures()){
              //      System.out.print(d+"\t");
              //  }
              //  System.out.println();
            }
            
            //double[]feature = gap.createFeatures();
            //c.setFeatures(feature);
            //for(Edge e: edge_list)
            //    System.out.println(e);
            ArrayList<String> id_residu = new ArrayList();
            id_residu.addAll(c.id_epi);
            id_residu.addAll(c.id_nonepi);
            String path_psai = "src/dataset/gao/psaia/";
            String pdbid = c.pdbid.substring(0, 4).trim();
            String chain = c.pdbid.substring(c.pdbid.length()-2).trim();
            String ending_psai = chain+"_202112081246_unbound";
            ComplexPreparation cp = new ComplexPreparation(path_ag,pdbid,chain);
            //data harus ditransformasi dengan fourier terlebih dahulu
           // ArrayList<Double> feat_cap = cp.calculateSelectedFeaturesInDouble(path_psai, ending_psai, id_residu, pca);
           
            ClusterAntigenicProperties cap = cp.calculateSelectedComplexFeatures(path_psai, ending_psai,id_residu, pca );
            
            if(!cap.getP().toString().contains("NaN")) {
                String prep = cap.getP().toString();
                int batas_awal = prep.indexOf("[");
                int batas_akhir = prep.indexOf("]");
                //System.out.println(cap.getVertexNum()+","+ prep);
                fitur += cap.getVertexNum()+","+prep.substring(batas_awal+1, batas_akhir);
                //System.out.println(cap.getVertexNum()+","+ fitur);
                //DataPoint p = cap.getP();
                
                double[] pca_feat = cap.getP().getNumericalValues().arrayCopy();
                for(double feat: pca_feat)
                    numeric_feat.add(feat);
                
                if (c.persen_epi>=0.7){
                    fitur +=","+true+"\n";
                    cats.add(entry.get("true"));
                }else{
                    fitur +=","+false+"\n";
                    cats.add(entry.get("false"));
                }
                System.out.println(fitur);
                listData.add(fitur);
                int[]a_cats = new int[cats.size()];
                for(int i=0;i<cats.size();i++){
                    a_cats[i]= cats.get(i);
                }
                DenseVector vec = new DenseVector(numeric_feat);
                DataPoint dp = new DataPoint(vec,a_cats,catData);
                System.out.println("data point: "+dp.toString());
                listDP.add(dp);
            }
                
            //System.out.println(cap.getP().toString());
            //c.setAntigen_features(cap.getfiturResume());
            //writeclusterfeaturetofile
            //String fitur = c.generateFiturForDataset(0.7);
            //if(fitur !=""){
            //    System.out.println(fitur);
            //    listData.add(fitur);
            //}
            
            
        }
        ClassificationDataSet cds = new ClassificationDataSet(listDP,0);
        String path_ca = "src/dataset/clustAntigenicProp.arff";
        io.writeArrayListStringToFile(path_ca, listData);
        return cds;
        //dari cluster baca propertiesnya
        
        //String graph_fn
        
    }
    public static void createDataSet(){
    
    }
    public static void calculateAntigenProperties(Cluster c, ComplexPreparation cp){
        //untuk setiap vertek ambil properties antigeniknya
        
    }
    public static ArrayList<Edge> extractEdgeFromFile(String graphfile, Cluster c){
        ArrayList <String> list_res_id = c.id_epi;
        list_res_id.addAll(c.id_nonepi);
        //for(String s: list_res_id)
        //    System.out.println(s);
        
        ArrayList<String> s_edge_list = getStringListFromFile(graphfile);
        ArrayList<Edge> edge_list = new ArrayList();
        for (String s: s_edge_list){
                //System.out.println("s s_edge_list: "+ s);
                //System.out.println("s s_edge_list: "+ s.split("\t")[0].trim());
                //System.out.println("s s_edge_list: "+ s.split("\t")[1].trim());
                //System.out.println("s s_edge_list: "+ s.split("\t")[2].trim());
                Edge e = new Edge(s);
                edge_list.add(e);
        }
        //ambil edge dengan salah satu vertex didalam cluster
         ArrayList<Edge> filtered_edge_list = new ArrayList();
        for(String s: list_res_id){
            for(Edge e : edge_list){
                if (s.equals(e.getV1())||s.equals(e.getV2())){
                    filtered_edge_list.add(e);                    
                }
            }
        }
        edge_list.clear();
        HashSet<Edge> edgelist = new HashSet();
        for( Edge e1: filtered_edge_list){
            //System.out.println(e1.toString());
            if(list_res_id.contains(e1.getV1())&& list_res_id.contains(e1.getV2())){
                //ambil 
                edgelist.add(e1);
                //System.out.println(e1.toString());
            }
        }
        ArrayList<Edge> edges = new ArrayList<Edge> (edgelist);
        return edges;
        
        
    
    }
    public static ArrayList<Cluster> createListClusterElemen(ArrayList<String> epi_file_list, ArrayList<ArrayList<String>> epinonepidist){
    //membuat object cluster
          ArrayList<Cluster> clus_list = new ArrayList();
          for(int i=0; i<epinonepidist.size();i++){
        //for(int i=0; i<epinonepidist.size();i++){
            System.out.println("epilist ke:"+i);
            String residu_name = epi_file_list.get(i).substring(0,10);
            System.out.println("epilist file: "+ epi_file_list.get(i));
            ArrayList<String> cluster = epinonepidist.get(i);
            System.out.println("jumlah baris: "+ cluster.size() );
            if(cluster.size()>0){
                //daftar cluster
                ArrayList<String[]> ec_list = new ArrayList();
                HashSet<String> c_no = new HashSet();                
                for(int j=0;j<cluster.size();j++){
                   String[] ec = extractCusterElemen(cluster.get(j));
                   ec_list.add(ec);             
                   c_no.add(ec[0]);
                }
               
                //create object Cluster                
                int count =0;
              //  ArrayList<Cluster> clus_list = new ArrayList();
              //daftar claster
                ArrayList<String> c_no_list = new ArrayList(c_no);
                // System.out.println("klaster: "+c_no_list);
                int count_cls =0;
                for(int k=0;k<c_no_list.size();k++){
                    String id_c = c_no_list.get(k);
                    Cluster cls = new Cluster(residu_name);
                    cls.id = id_c;
                    //menambah elemen cluster dari Arraylistlist
                    String[] s = ec_list.get(count_cls);                    
                    while (s[0].equals(id_c) && count_cls <cluster.size()){
                        cls.addCE(s);
                        count_cls +=1;
                        if(count_cls<cluster.size())
                            s = ec_list.get(count_cls);                        
                                                
                    }
                    cls.calculate_pe();
                    clus_list.add(cls);         
                    
                }  
                 
            }
        }
        //menampilkan jumlah cluster dengan persen epi
        ArrayList<Cluster> new_clus_list = new ArrayList();
        int num_c_epi_gt_70=0;
        int num_c_ne=0;
        //int num_c_epi_lt_70 =0;
        for (Cluster c: clus_list){
            if (c.getProcent()==0){
                num_c_ne +=1;
                c.setEpi_cluster(false);
                new_clus_list.add(c);               
                
            }
            if (c.getProcent()>0.7){
                num_c_epi_gt_70 +=1;
                c.setEpi_cluster(true);
                new_clus_list.add(c);
                }
                
            
            
        }
        System.out.println("num_c_epi_gt_70: "+ num_c_epi_gt_70);
        System.out.println("num_c_ne: "+num_c_ne);
        //System.out.println("num_c_epi_lt_70: "+ num_c_epi_lt_70);
        return new_clus_list;
    }
    private static String[] extractCusterElemen(String ce){
       // System.out.println(ce.split("\t")[0]);
       // System.out.println(ce.split("\t")[1]);
        String[] ec = new String[3];
        ec[0] = ce.split("\t")[0].substring(0,5).trim();
        ec[1] = ce.split("\t")[0].substring(5).trim();
        ec[2] = ce.split("\t")[1].trim();
        return ec;
    }
    public static void printlist(ArrayList<ArrayList<String>> epinonepidist){
        int count=0;
        for(ArrayList<String> dist: epinonepidist){
            System.out.println("id list: "+count);
            for(String s: dist){
                System.out.println(s);
            }
            count ++;
        }
    }
    public static ArrayList<ArrayList<Cluster>> checkEpiPercentageInCluster(ArrayList<ArrayList<String>> epinonepidist){
        ArrayList<ArrayList<Cluster>> cluster_list = new ArrayList();
        int count =0;
        int num_cluster=0;
        for(ArrayList<String> clus_dist :epinonepidist){
            System.out.println("clust_dist "+ count);
            ArrayList<Cluster> cluster_id = new ArrayList();            
            for(String s: clus_dist){
                // data : 0              80             	ne
                //jika displit menjadi:
                //0              80             	
                //ne
                String []splits = s.split("\t");
                String idc = splits[0].substring(0,12).trim();
                String idr = splits[0].substring(12).trim();
                Integer id_c = getIdClusterInArrayList(idc, cluster_id);
                //System.out.println("idc"+ idc+ " id_c: "+ id_c);
                if(id_c == -1){
                    //buat clusterbaru
                    Cluster c;
                    if(splits[1].equalsIgnoreCase("e")){
                        c = new Cluster(idc,true);
                    }else{
                        c = new Cluster(idc,false);
                    }
                    cluster_id.add(c);
                }else{
                    //cluster_id.get(id_c).
                    if(splits[1].equalsIgnoreCase("e")){
                        cluster_id.get(id_c).addCount(true);
                    }else{
                        cluster_id.get(id_c).addCount(false);
                    }
                        
                }
                //System.out.println(idc);
                //System.out.println(idr);
                //System.out.println(splits[1]);
            }
            
            cluster_list.add(cluster_id);
            
            for(Cluster c: cluster_id){
                c.calculate_pe();
                if(c.getProcent()>=0.5){
                    System.out.println(c.id+"\t"+c.getNum_epi()+"\t"+c.getNum_ne()+"\t"+c.getProcent());
                    num_cluster +=1;
                }                                    
            }
            count +=1;
            
        }
        System.out.println("jumlah cluster proc epi > 50%: "+num_cluster);
        return cluster_list;
    }
    private static Integer getIdClusterInArrayList(String idC, ArrayList<Cluster> clusters){
        Integer idcluster = -1;
        for(int i=0;i< clusters.size();i++){
            Cluster c = clusters.get(i);
            if (idC.equalsIgnoreCase(clusters.get(i).getId())){
                idcluster = i;
                break;
            }
        }
        return idcluster;
    }
    public static ArrayList<ArrayList<String>> cekEpiNonEpiDistributionInCluster(ArrayList<String> epi_file_list, String path_c, String path_epi){
        ArrayList<ArrayList<String>>list_epi_c = new ArrayList();
        System.out.println(path_c);
        for(String fn: epi_file_list){
            String[]splitfn = fn.split("\t");
            String cluster_file_name = path_c + splitfn[1];
           // System.out.println(splitfn[1]);
            String epi_file_name = path_epi + splitfn[0];
            //System.out.println(epi_file_name);
            //System.out.println(cluster_file_name);
             ArrayList<String> epitope = readEpiListFromFile(epi_file_name);
             ArrayList<String> cluster =  readClusterFromFile(cluster_file_name);
             
             ArrayList <String> epi_nonepi_dist_in_cluster = new ArrayList();
             if(epitope.size()>0 && cluster.size()>0){
                for(String c: cluster){
                    boolean found = false;
                 // System.out.println("c"+c.trim());
                    String c1 = c.trim().substring(0,10).trim();
                    String c2 = c.trim().substring(13).trim();
                 // System.out.println(c1 +","+ c2);
                    String concat_e_c="";
                    //cari apakah c2 ada didalam dafta epitop, kalau iya maka catat dan delete epi dari daftar
                    for(String s: epitope){
                        if (s.equals(c2)){
                            concat_e_c = c+"\t"+"e";
                            epitope.remove(s);
                            found = true;
                            break;                                
                        }
                    }
                    if (found == false){
                        concat_e_c = c+"\t"+"ne";
                    }
                    //System.out.println(concat_e_c);
                    epi_nonepi_dist_in_cluster.add(concat_e_c);
                    }
                }
                
            list_epi_c.add(epi_nonepi_dist_in_cluster);
            }
             
    
        return list_epi_c;
    }
    
    public static void resumeClusterEpi(ArrayList <ArrayList<String>> list_dist_cluster){
        
        for(int i= 0;i<list_dist_cluster.size();i++){
            HashSet noDupSet = new HashSet();
            ArrayList<String> list = list_dist_cluster.get(i);
            for(int j=0;j<list.size();j++){
                String s = list.get(j).split("\t")[1];
               // System.out.println(s.split("\t")[0]);
                noDupSet.add(s);
                
            }
            //System.out.println("jumlah cluster epi:"+ noDupSet.size());
            System.out.println( noDupSet.size());
        }
    }
    
    public static void cluster(){
      //  String input_path = "src/dataset/graph/exposed/weight/";
      
        String input_path = "src/dataset/gao/graph/";
	String outputPath= "src/dataset/gao/graph/mcl/bf0_75/";
	int coarseMode =1 ; //"-sc" atau "-hem"
	int coarseLevel = 2; //non negatif
	double bFactor = 0.75; //non negatif
	int mcl_mode=1; //"-reg" atau "-basic"
	int numThread= 4; //harus non negatif
	double epsilon=0.2; // lebih dari nol untuk menyetop MCL
	int rand_seed= 3; // 
	double skiprate =0.5; //bilangan antara 0 dan 1

        //baca file dari folder
        File folder = new File(input_path);
        FileReader fr = null;
        
        try{
            File[] listOfFiles = folder.listFiles();
            for (int i = 0; i < listOfFiles.length; i++) {
                if (listOfFiles[i].isFile()) {
                    String fn = listOfFiles[i].getName();
                    String fileName = input_path+fn;
                        System.out.println(fileName);
                        try{
                            
                            MCL.run(fileName, outputPath, coarseMode, coarseLevel, bFactor, mcl_mode, numThread, epsilon, rand_seed);

                        }catch(Exception e){
                            System.out.println("file can't read");
                        }
        }
        }
        }catch(Exception e){
        System.out.print("can't extract graph from file");
        }
    }
    public static void clusterbycn(){
      //  String input_path = "src/dataset/graph/exposed/weight/";
      
        String input_path = "src/dataset/gao/graph/cn/";
	String outputPath= "src/dataset/gao/graph/cn/mcl/bf0_75/";
	int coarseMode =1 ; //"-sc" atau "-hem"
	int coarseLevel = 2; //non negatif
	double bFactor = 0.75; //non negatif
	int mcl_mode=1; //"-reg" atau "-basic"
	int numThread= 4; //harus non negatif
	double epsilon=0.2; // lebih dari nol untuk menyetop MCL
	int rand_seed= 3; // 
	double skiprate =0.5; //bilangan antara 0 dan 1

        //baca file dari folder
        File folder = new File(input_path);
        FileReader fr = null;
        
        try{
            File[] listOfFiles = folder.listFiles();
            for (int i = 0; i < listOfFiles.length; i++) {
                if (listOfFiles[i].isFile()) {
                    String fn = listOfFiles[i].getName();
                    String fileName = input_path+fn;
                        System.out.println(fileName);
                        try{
                            
                            MCL.run(fileName, outputPath, coarseMode, coarseLevel, bFactor, mcl_mode, numThread, epsilon, rand_seed);

                        }catch(Exception e){
                            System.out.println("file can't read");
                        }
        }
        }
        }catch(Exception e){
        System.out.print("can't extract graph from file");
        }
    }
    public static void clusterbycnw(){
      //  String input_path = "src/dataset/graph/exposed/weight/";
      
        String input_path = "src/dataset/gao/graph/cnw/";
	String outputPath= "src/dataset/gao/graph/cnw/mcl/bf0_75/";
	int coarseMode =1 ; //"-sc" atau "-hem"
	int coarseLevel = 2; //non negatif
	double bFactor = 0.75; //non negatif
	int mcl_mode=1; //"-reg" atau "-basic"
	int numThread= 4; //harus non negatif
	double epsilon=0.2; // lebih dari nol untuk menyetop MCL
	int rand_seed= 3; // 
	double skiprate =0.5; //bilangan antara 0 dan 1

        //baca file dari folder
        File folder = new File(input_path);
        FileReader fr = null;
        
        try{
            File[] listOfFiles = folder.listFiles();
            for (int i = 0; i < listOfFiles.length; i++) {
                if (listOfFiles[i].isFile()) {
                    String fn = listOfFiles[i].getName();
                    String fileName = input_path+fn;
                        System.out.println(fileName);
                        try{
                            
                            MCL.run(fileName, outputPath, coarseMode, coarseLevel, bFactor, mcl_mode, numThread, epsilon, rand_seed);

                        }catch(Exception e){
                            System.out.println("file can't read");
                        }
        }
        }
        }catch(Exception e){
        System.out.print("can't extract graph from file");
        }
    }
    public static void analysisEpitopDistributionInClusterMCL(){
    //baca data epitop dari file
    String epitop_path = "src/dataset/gao/epitopes";
    ArrayList<ArrayList<String>> epi_list = new ArrayList();
     File folder = new File(epitop_path);
        FileReader fr = null;
        try{
            File[] listOfFiles = folder.listFiles();
             
            for (int i = 0; i < listOfFiles.length; i++) {
                if (listOfFiles[i].isFile()) {
                    String fileName = listOfFiles[i].getName();
                    //read data from file
                    ArrayList<String> epi = readStringListFromFile(epitop_path + fileName);
                    epi_list.add(epi);
                }
            }
        }catch(Exception e){
        
        }
    System.out.println("daftar epilist"+ epi_list.size());
    //baca data cluster setiap antigen
    String path_gao_cluster_bf05 = "src/dataset/gao/graph/mcl/";
    String path_gao_cluster_bf025 = "src/dataset/gao/graph/mcl/bf0_25/";
    String path_gao_cluster_bf075 = "src/dataset/gao/graph/mcl/bf0_75/";
    System.out.println("====");
    ArrayList<ArrayList<String>> gao_cluster_bf05= getClusterListFromFolder(path_gao_cluster_bf05);
    System.out.println("gao_cluster_bf05: "+gao_cluster_bf05.size());
    System.out.println("====");
    ArrayList<ArrayList<String>> gao_cluster_bf025= getClusterListFromFolder(path_gao_cluster_bf025);
    System.out.println("gao_cluster_bf025: "+gao_cluster_bf025.size());
   // for(gao_cluster_bf025.get(0).)
    System.out.println("====");
    ArrayList<ArrayList<String>> gao_cluster_bf075= getClusterListFromFolder(path_gao_cluster_bf075);
    System.out.println("gao_cluster_bf075: "+gao_cluster_bf075.size());    
    //rangkum sebaran cluster  #id_residu #residu #no_cluster
    
    }
    public static ArrayList<ArrayList<String>> cekDistribusiEpitopPdCluster(ArrayList<String> epi_file_list, String path_c, String path_epi){
        ArrayList<ArrayList<String>>list_epi_c = new ArrayList();
        for(String fn: epi_file_list){
            String[]splitfn = fn.split("\t");
            String cluster_file_name = path_c + splitfn[1];
            String epi_file_name = path_epi + splitfn[0];
            //System.out.println(epi_file_name);
            //System.out.println(cluster_file_name);
             ArrayList<String> epitope = readEpiListFromFile(epi_file_name);
             ArrayList<String> cluster =  readClusterFromFile(cluster_file_name);
             //cek epitop pada cluster
             ArrayList<String> epi_c_list = new ArrayList();
             if(epitope.size()>0 && cluster.size()>0){
                for(String s: epitope){
                    //System.out.println("s"+ s);
                    for(String c: cluster){
                       // System.out.println("c"+c.trim());
                        String c1 = c.trim().substring(0,10).trim();
                        String c2 = c.trim().substring(13).trim();
                       // System.out.println(c1 +","+ c2);
                        if (s.equals(c2)){
                                String concat_e_c = s +"\t"+ c1;
                                //System.out.println(concat_e_c);
                                epi_c_list.add(concat_e_c);
                        }
                        
                    }
                }
                list_epi_c.add(epi_c_list);
            }
        }
      
        return list_epi_c;
        
    }
    public static ArrayList<String> getEpiFilename(String path){
        //baca daftar file cluster
         String path_gao_cluster_bf05 = "src/dataset/gao/graph/mcl/";
        String path_gao_cluster_bf025 = "src/dataset/gao/graph/mcl/bf0_25/";
        String path_gao_cluster_bf075 = "src/dataset/gao/graph/mcl/bf0_75/";
        String ending_gao_cluster_bf05 = ".assign";
        //ambil daftar file pada path
        ArrayList<String> gao_cluster_bf05_file_list = getFileListFromFolder(path, ending_gao_cluster_bf05);
        Collections.sort(gao_cluster_bf05_file_list);   
       
        //baca daftar file epitop
        String path_epilist = "src/dataset/gao/epitopes/";
        String ending_epi = "";
        //ambil daftar file pada path
        ArrayList<String> epi_file_list = getFileListFromFolder(path_epilist, ending_epi);
        Collections.sort(epi_file_list);
        
        //jika sesuai maka buat string concat keduanya dan simpan pada daftar
        ArrayList<String> list_file = new ArrayList();
        for(String e: epi_file_list){
            String start = e.substring(0, 4)+ e.substring(8);
            //System.out.println(start);
            for(String f: gao_cluster_bf05_file_list){
                if(f.startsWith(start)){
                    String fn_concat = e+"\t"+f;
                    list_file.add(fn_concat);
                    break;
                }
            
            }
        }
        /*
        for(String s: list_file){
            System.out.println(s);
        }
        */
        return list_file;
    }
     public static ArrayList<String> getStringListFromFile(String filename){
         ArrayList<String> string_list = new ArrayList();
         try{
            if (new File(filename).exists()) {
                String dataline = null;
                FileReader fr = null;
                fr = new FileReader(new File(filename));
                BufferedReader br = new BufferedReader(fr);
                dataline = br.readLine(); 
                
                while (dataline != null)
                {
                    string_list.add(dataline);                        
                    dataline = br.readLine();
                }
                 br.close(); 
            }
            
        }catch(IOException e){
                e.printStackTrace();
            }
         return string_list;
     }
     public static ArrayList<String> getFileListFromFolder(String path, String ending){
        ArrayList<String> file_list = new ArrayList();
        File folder = new File(path);
        try{
            File[] listOfFiles = folder.listFiles();
            for (int i = 0; i < listOfFiles.length; i++) {
                if (listOfFiles[i].isFile()) {
                    String fileName = listOfFiles[i].getName();
                    //read data from file
                    if (fileName.endsWith(ending)){
                        file_list.add(fileName);
                    }
                    }
                }
        }catch(Exception e){

        }
        return file_list;
    }
    public static ArrayList<ArrayList<String>> getClusterListFromFolder(String path){
        ArrayList<ArrayList<String>> cluster_list = new ArrayList();
        File folder = new File(path);
        try{
            File[] listOfFiles = folder.listFiles();
            for (int i = 0; i < listOfFiles.length; i++) {
                if (listOfFiles[i].isFile()) {
                    String fileName = listOfFiles[i].getName();
                    //read data from file
                    if (fileName.endsWith("assign")){
                        ArrayList<String> clusters = readClusterFromFile(path + fileName);
                        cluster_list.add(clusters);
                    }


                    }
                }
        }catch(Exception e){

        }
        return cluster_list;
    }
    
    public static ArrayList<String> readClusterFromFile(String filename){
        ArrayList<String> cluster = new ArrayList();
        try{
            if (new File(filename).exists()) {
                String dataline = null;
                FileReader fr = null;
                fr = new FileReader(new File(filename));
                BufferedReader br = new BufferedReader(fr);
                dataline = br.readLine(); 
                dataline = br.readLine(); 
                if(dataline == null){
                        System.out.println("cluster null: "+filename);
                    }
                while (dataline != null)
                {
                    //System.out.println(dataline);
                    cluster.add(dataline);                        
                    dataline = br.readLine();
                }
                 br.close(); 
            }
            
        }catch(IOException e){
                e.printStackTrace();
            }
        //System.out.println("cluster size: "+ cluster.size());
        return cluster;
    }
    public static ArrayList<String> readEpiListFromFile(String filename){
        ArrayList<String> epilist = new ArrayList();
        try{
            if (new File(filename).exists()) {
                String dataline = null;
                FileReader fr = null;
                fr = new FileReader(new File(filename));
                    BufferedReader br = new BufferedReader(fr);
                    dataline = br.readLine();                
                    while (dataline != null)
                    {
                        //System.out.println(dataline);
                        String[] splitdl = dataline.split(" "); 
                        if(splitdl.length==2){                        
                            epilist.add(splitdl[1].trim());                        
                        }
                        dataline = br.readLine();
                    }
                    br.close();                
            }
        }catch(IOException e){
                e.printStackTrace();
            }
        //System.out.println("epilist size: "+ epilist.size());
        return epilist;
    }
    public static ArrayList<String> readStringListFromFile(String filename){
        ArrayList<String> epilist = new ArrayList();
        try{
            if (new File(filename).exists()) {
                String dataline = null;
                FileReader fr = null;
                fr = new FileReader(new File(filename));
                    BufferedReader br = new BufferedReader(fr);
                    dataline = br.readLine();                
                    while (dataline != null)
                    {
                        //System.out.println(dataline);
                        String[] splitdl = dataline.split(" "); 
                        if(splitdl.length==2){                        
                            epilist.add(splitdl[1].trim());                        
                        }
                        dataline = br.readLine();
                    }
                    br.close();                
            }
        }catch(IOException e){
                e.printStackTrace();
            }
        return epilist;
    }
    private static class Cluster{
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
        
        
        
        public String generateFiturForDataset(double batas_persen_epi){
            String fitur ="";
            if(features!=null){
                for(int i=0;i<this.features.length;i++){
                fitur = fitur + String.valueOf(this.features[i]).trim()+",";
                }}
            if(antigen_features!=null){
                for(int i=0;i<this.antigen_features.length;i++){
                fitur = fitur + String.valueOf(this.antigen_features[i]).trim()+",";
                }
                if(this.persen_epi >= batas_persen_epi){
                    fitur = fitur + "true"+"\n";
                }else{
                        fitur = fitur + "false"+"\n";
                }                   
            }
            System.out.println(fitur);
            return fitur;
            
        }
        public Cluster (String pdbid){
            this.pdbid=pdbid;
            this.id_epi = new ArrayList();
            this.id_nonepi =  new ArrayList();
        }
        public Cluster(String id, boolean epi){
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
    
}
