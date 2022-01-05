/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ExtractDataTraining;

import ComplexFeature.AAIndex;
import ComplexFeature.ExposedOnComplex;
import ComplexStructure.Aminoacid;
import ComplexStructure.Atom;
import ComplexStructure.Complex;
import ComplexStructure.ExposedResidu;
import ExtractData.CreateComplexStructure;
import ExtractFromFile.IO;
import SpatialAnalysis.PairIdResidu;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import jsat.classifiers.DataPoint;
import jsat.linear.DenseVector;

/**
 *
 * @author dell
 */
public class GraphAntigen {
    private String exceptelement = "O";
    private ArrayList<PairIdResidu> pirlist = new ArrayList();
    
    public ArrayList<PairIdResidu> loadGraphfromFile (String path, String filename){
        ArrayList<PairIdResidu> graph = new ArrayList();
        FileReader fr;
        BufferedReader br;
        String dataline;
        try{
            if(new File(path+ filename).exists()){
            fr = new FileReader(new File(path+ filename));
            br = new BufferedReader(fr);            
            dataline=br.readLine();
            while(dataline != null){}
                dataline=dataline.trim();
                String[] spl = dataline.split("\t");
                
            }
        }catch(Exception e){}
        return graph;
    }
   
    
    public ArrayList<Double>getAAIList(Aminoacid aa){
        IO io = new IO();
        AAIndex aai;
        aai= io.extractAAIndexFromFile();
        aai.composeAAI();
        double[] withoutNan =aai.getAaiWithoutNanData(aa.letterToNum());
        ArrayList<Double> aaid = new ArrayList();
        for(int i=0;i<withoutNan.length;i++){
            aaid.add(withoutNan[i]);
        }
        return aaid;
    }
    //proses develop graph
    // baca data complex didalam list
    //buat struktur asam aminonya
    //buat daftar berisi pasangan residu kalau jarak antara keduanya pada suatu nilai treshold
    //tambahkan bobot edgenya
    
    //daftar antigen dengan residuId string yang bisa dikonversi ke integer
    public ArrayList<String> cekIdResiduContainString(String path){
        File folder = new File(path);
        File[] listOfFiles = folder.listFiles();
        System.out.println("jumlah file: "+listOfFiles.length);
        ArrayList<String> list_file_with_id_string = new ArrayList();
        for (int i = 0; i < listOfFiles.length; i++) {
          if (listOfFiles[i].isFile()) {
            //System.out.println("File " + listOfFiles[i].getName());
            //lakukan proses baca dan jika mengandung string maka masukkan dalam daftar mengandung string tidk bisa dikonversi
            String filename = path + listOfFiles[i].getName();
            FileReader fr;
            BufferedReader br;
            String dataline;
           
            try{
                if(new File(filename).exists()){
                    fr = new FileReader(new File(filename));
                    br = new BufferedReader(fr);            
                    dataline=br.readLine();
                    while(dataline != null){
                        String[] spl = dataline.split("\t");
                        int id_1, id_2;
                        try {
                           id_1 = Integer.parseInt(spl[0]);
                           id_2 = Integer.parseInt(spl[1]);
                           
                        }
                        catch (NumberFormatException e)
                        {
                           // System.out.println(spl[0]+spl[1]);
                          //lakukan pemetaan ke id baru dengan hashmap
                            if (!list_file_with_id_string.contains(filename)){
                               list_file_with_id_string.add(listOfFiles[i].getName());
                               //System.out.println(listOfFiles[i].getName());
                            }   


                         }
                       dataline=br.readLine(); 
                    }}
            }catch(Exception e){
               e.printStackTrace();
            }
          } 
        }
        ArrayList<String> newlist = new ArrayList();
        for(int i= 0; i<list_file_with_id_string.size();i++){
                String s = list_file_with_id_string.get(i);
                if(!newlist.contains(s)){
                    newlist.add(s);
                }
        }
        
           
        return newlist;
    }
    
    public ArrayList<PairIdResidu> generateGraphFromAAList(ArrayList<Aminoacid> listAA, double treshold, double[][]weight){
        ArrayList<PairIdResidu> pirlist = new ArrayList(); 
        
         for (int i=0;i<listAA.size();i++){
             for(int j=1;j<listAA.size();j++){
                 if(i!= j){
                    Aminoacid aa1 = listAA.get(i);
                    Aminoacid aa2 = listAA.get(j);
                    double d = this.distance_amino_amino(aa1,aa2);
                    double w = weight[aa1.letterToNum()][aa2.letterToNum()];
                    if (d<= treshold){
                       PairIdResidu pir = new PairIdResidu(aa1.getresidueID(),aa2.getresidueID(),w);
                       pirlist.add(pir);
                 }
                 }
                 
             }
         }
         return pirlist;
    }
    public ArrayList<String> generateStringFromAAListLogWeight(ArrayList<Aminoacid> listAA, double treshold, double[][]weight){
         Set<String> pirlist = new HashSet<String>();
        
         for (int i=0;i<listAA.size();i++){
             for(int j=1;j<listAA.size();j++){
                 if(i!= j){
                    Aminoacid aa1 = listAA.get(i);
                    Aminoacid aa2 = listAA.get(j);
                    double d = this.distance_amino_amino(aa1,aa2);
                    double w = Math.log(weight[aa1.letterToNum()][aa2.letterToNum()]);
                    if (d<= treshold){     
                        String edge = null;
                        try{
                        int id1 = Integer.parseInt(aa1.getresidueID());
                        int id2 = Integer.parseInt(aa2.getresidueID());
                        
                        if(id1<=id2){
                            edge = aa1.getresidueID()+"\t"+aa2.getresidueID()+"\t"+ String.valueOf(w);
                        }else{
                                edge = aa2.getresidueID()+"\t"+ aa1.getresidueID()+"\t"+String.valueOf(w);
                        }
                        }catch(Exception e){
                            if(aa1.getresidueID().length()<aa2.getresidueID().length()){
                                 edge = aa1.getresidueID()+"\t"+aa2.getresidueID()+"\t"+ String.valueOf(w);
                            }else{
                                 edge = aa2.getresidueID()+"\t"+aa1.getresidueID()+"\t"+ String.valueOf(w);
                            }
                        }
                        pirlist.add(edge);
                    }
                 }
                 
             }
         }
         ArrayList<String> list = new ArrayList<>(pirlist);
         return list;
         /*
         ArrayList<PairIdResidu> newpirlist = new ArrayList(); 
         for(int i=0;i<pirlist.size();i++){
             PairIdResidu pir = pirlist.get(i);
             PairIdResidu pir2 = new PairIdResidu(pir.getIdResiduS(),pir.getIdResiduF(),pir.getWeight());
             if(!newpirlist.contains(pir)){
                 if( !newpirlist.contains(pir2))
                    newpirlist.add(pirlist.get(i));
             }
         }
         */
         
        
    } 
    public ArrayList<String> generateStringFromAAListWeight(ArrayList<Aminoacid> listAA, double treshold, double[][]weight){
         Set<String> pirlist = new HashSet<String>();
        
         for (int i=0;i<listAA.size();i++){
             for(int j=1;j<listAA.size();j++){
                 if(i!= j){
                    Aminoacid aa1 = listAA.get(i);
                    Aminoacid aa2 = listAA.get(j);
                    double d = this.distance_amino_amino(aa1,aa2);
                    //double w = Math.log(1000*weight[aa1.letterToNum()][aa2.letterToNum()]);
                    double w = 1000*weight[aa1.letterToNum()][aa2.letterToNum()];
                    if (d<= treshold){     
                        String edge = null;
                        try{
                        int id1 = Integer.parseInt(aa1.getresidueID());
                        int id2 = Integer.parseInt(aa2.getresidueID());
                        
                        if(id1<=id2){
                            edge = aa1.getresidueID()+"\t"+aa2.getresidueID()+"\t"+ String.valueOf(w);
                        }else{
                                edge = aa2.getresidueID()+"\t"+ aa1.getresidueID()+"\t"+String.valueOf(w);
                        }
                        }catch(Exception e){
                            if(aa1.getresidueID().length()<aa2.getresidueID().length()){
                                 edge = aa1.getresidueID()+"\t"+aa2.getresidueID()+"\t"+ String.valueOf(w);
                            }else{
                                 edge = aa2.getresidueID()+"\t"+aa1.getresidueID()+"\t"+ String.valueOf(w);
                            }
                        }
                        pirlist.add(edge);
                    }
                 }
                 
             }
         }
         ArrayList<String> list = new ArrayList<>(pirlist);
         return list;
         /*
         ArrayList<PairIdResidu> newpirlist = new ArrayList(); 
         for(int i=0;i<pirlist.size();i++){
             PairIdResidu pir = pirlist.get(i);
             PairIdResidu pir2 = new PairIdResidu(pir.getIdResiduS(),pir.getIdResiduF(),pir.getWeight());
             if(!newpirlist.contains(pir)){
                 if( !newpirlist.contains(pir2))
                    newpirlist.add(pirlist.get(i));
             }
         }
         */
         
        
    } 
    public ArrayList<PairIdResidu> generateGraphFromAAListWeight(ArrayList<Aminoacid> listAA, double treshold, double[][]weight){
         ArrayList<PairIdResidu> pirlist = new ArrayList();
        
         for (int i=0;i<listAA.size();i++){
             for(int j=1;j<listAA.size();j++){
                 if(i!= j){
                    Aminoacid aa1 = listAA.get(i);
                    Aminoacid aa2 = listAA.get(j);
                    double d = this.distance_amino_amino(aa1,aa2);
                    //double w = Math.log(weight[aa1.letterToNum()][aa2.letterToNum()]);
                    double w = weight[aa1.letterToNum()][aa2.letterToNum()]*1000;
                    if (d<= treshold){         
                        try{
                        int id1 = Integer.parseInt(aa1.getresidueID());
                        int id2 = Integer.parseInt(aa2.getresidueID());
                        if(id1<=id2){
                            PairIdResidu pir = new PairIdResidu(aa1.getresidueID(),aa2.getresidueID(),w);
                            if(!pirlist.contains(pir)){
                                pirlist.add(pir);
                            }
                        }else{
                            PairIdResidu pir = new PairIdResidu(aa2.getresidueID(),aa1.getresidueID(),w);
                            if(!pirlist.contains(pir)){
                                pirlist.add(pir);
                            }
                        }
                        }catch(Exception e){
                            if(aa1.getresidueID().length()<aa2.getresidueID().length()){
                                PairIdResidu pir = new PairIdResidu(aa1.getresidueID(),aa2.getresidueID(),w);
                                if(!pirlist.contains(pir)){
                                    pirlist.add(pir);
                                }
                            }else{
                                PairIdResidu pir = new PairIdResidu(aa2.getresidueID(),aa1.getresidueID(),w);
                                if(!pirlist.contains(pir)){
                                    pirlist.add(pir);
                                }
                            }
                        }
                       
                       
                 }
                 }
                 
             }
         }
         System.out.println();
         Set<PairIdResidu> pirs = new HashSet<PairIdResidu>(pirlist);
         pirlist.clear();
         pirlist.addAll(pirs);
         return pirlist;
         /*
         ArrayList<PairIdResidu> newpirlist = new ArrayList(); 
         for(int i=0;i<pirlist.size();i++){
             PairIdResidu pir = pirlist.get(i);
             PairIdResidu pir2 = new PairIdResidu(pir.getIdResiduS(),pir.getIdResiduF(),pir.getWeight());
             if(!newpirlist.contains(pir)){
                 if( !newpirlist.contains(pir2))
                    newpirlist.add(pirlist.get(i));
             }
         }
         */
         
        
    }
    private boolean addPIR(PairIdResidu pir, ArrayList<PairIdResidu> pirlist ){
        
        if (pirlist.contains(pir)){
            return false;
        }else{
            PairIdResidu newpir = new PairIdResidu(pir.getIdResiduS(),pir.getIdResiduF(),pir.getWeight());
            if (pirlist.contains(newpir)){
                return false;
            }else{
                pirlist.add(newpir);
                return true;
            }
        }
        
    }
    public ArrayList<PairIdResidu> generateGraph(Complex c, double treshold){
        ArrayList<PairIdResidu> pirlist = new ArrayList(); 
        ArrayList<Aminoacid> listAA= c.getAminoacidVector();
         for (int i=0;i<listAA.size();i++){
             for(int j=1;j<listAA.size();j++){
                 if(i!= j){
                    Aminoacid aa1 = listAA.get(i);
                    Aminoacid aa2 = listAA.get(j);
                    double d = this.distance_amino_amino(aa1,aa2);
                    if (d<= treshold){
                       PairIdResidu pir = new PairIdResidu(aa1.getresidueID(),aa2.getresidueID(),0);
                       pirlist.add(pir);
                 }
                 }
                 
             }
         }
         return pirlist;
    }
    
   
    
    private double distance_amino_amino(Aminoacid aa1, Aminoacid aa2)
  {
    int n = aa1.getatomNum();
    int m = aa2.getatomNum();
    double distance = 10000.0D;
    if ((aa1.getatomNum() == aa2.getatomNum()) && (aa1.getresidueID().equals(aa2.getresidueID())) && (aa1.getchainID() == aa2.getchainID())) {
      return distance;
    }
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        if ((!aa1.getAtombyOrder(i).equals(this.exceptelement)) && (!aa2.getAtombyOrder(j).equals(this.exceptelement)))
        {
          double new_distance = distanceBetweenPoint3D(aa1.getAtombyOrder(i), aa2.getAtombyOrder(j));
          if (new_distance < distance) {
            distance = new_distance;
          }
        }
      }
    }
    return distance;
  }
  
  private double distanceBetweenPoint3D(Atom pa, Atom pb)
  {
    double distance = 0.0D;
    double x = pa.getX() - pb.getX();
    double y = pa.getY() - pb.getY();
    double z = pa.getZ() - pb.getZ();
    distance = Math.sqrt(x * x + y * y + z * z);
    return distance;
  }
}
