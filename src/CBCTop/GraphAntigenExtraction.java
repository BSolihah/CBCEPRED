/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CBCTop;

import ComplexFeature.AminoAcidPairByFuncSubgrup;
import ComplexFeature.AminoacidPairByCategory;
import ComplexStructure.Complex;
import ExtractData.CreateComplexStructure;
import ExtractDataTraining.GraphAntigenExposedAtom;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author dell
 */
public class GraphAntigenExtraction {
    public static void main(String[] args) throws IOException {
    
    //String list_file ="list_antigen_test_seppa3.txt";
        String list_file = "list_antigen_train_seppa3.txt";
       //String list_part = "list_antigen_train_seppa3_2.txt";
       String path_aap = "src/dataset/graph/frekuensi/";
       String path_aap_fs = "src/dataset/graph/frekuensi/";
       String filename = "seppa3.txt";
       String filename_fs = "gao_fs.txt";
       String path_pdb ="src/dataset/gao/chain/";
       String path_exp_atom = "src/";
       //String path_pdb, String path_exp_atom, double t_asa, double t_dist, String path, String filename
       //AminoacidPairByCategory apc = new AminoacidPairByCategory(path_pdb,path_exp_atom, 0.01, 4,path_aap,filename);
       AminoAcidPairByFuncSubgrup apfs = new AminoAcidPairByFuncSubgrup(path_pdb,path_exp_atom, 0.01, 4,path_aap_fs,filename_fs);
       double[][] norm_weight = apfs.calculateWight();
       //double[][] weight = apc.calculateWight();
       //generate from exposed atom
       GraphAntigenExposedAtom ga = new GraphAntigenExposedAtom();
       //lakukan pembentukan graph dari file didalam folder
       String path_in_complex = "src/dataset/gao/chain/";
       String path_in_atom = "src/dataset/gao/antigens/"; //berisi pdbid+chain dan pdbid+chain
       String path_in_res = "src/dataset/gao/antigens/res/"; //berisi pdbid+chain dan pdbid+chain+res
       // graph dari bobot normal chisquare test
       //String path_out = "src/dataset/gao/graph/dist_8/";
       //  graph dari bobot cn
      // String path_out = "src/dataset/gao/graph/cn/";
       //String path_out = "src/dataset/gao/graph/cnw/";
       String path_out = "src/dataset/gao/graph/wfuncsubgrup/";
       //double t_distance = 6;
       double t_distance = 6.5;
       
        try{
            FileReader fr = null;
            File folder = new File(path_in_res);
            File[] listOfFiles = folder.listFiles();
            System.out.println("jumlah listoffile: "+ listOfFiles.length);
            for (int i = 0; i < listOfFiles.length; i++) {
                if (listOfFiles[i].isFile()) {   
                    String fileNameRes = listOfFiles[i].getName();
                    String pdbid = fileNameRes.substring(0, 4);
                    String chain = fileNameRes.substring(5,5);
                    String fileNameAtom = fileNameRes.substring(0, 5);
                    System.out.println("pdb: "+pdbid + ":"+chain );
                    try{
                        Complex c = new Complex(path_in_complex,pdbid,chain);
                    ga.generateGraph(c,path_in_res,path_in_atom, fileNameRes,fileNameAtom, t_distance, norm_weight, path_out);
                    }catch(Exception e){
                        System.out.println(e.getMessage());
                    }
                    
                    
                }
            }
        }catch(Exception e){
            System.out.println("exception in generate ga:" +e.getMessage());
        }
        // apc.concatAA();
     //  double [][] w = apc.calculateWight();
      /*
       String path_w = "src/dataset/graph/exposed/weight/";
        
        GraphAntigen ga = new GraphAntigen();
       
        //ga.createGraphAntigen(list_file, 4);
     //   ga.createGraphExposedAntigen(list_file, 4, 0.01,w,path_w);
     // apc.writeWeightAAPairToFile("seppa3_train",w);
      //String pathlogofweight = "src/dataset/graph/exposed/logofweight/";
       ArrayList<String> list_file_with_id_string= ga.cekIdResiduContainString(path_w);
     //for(int i=0;i< list_file_with_id_string.size();i++)
       //  System.out.println(list_file_with_id_string.get(i));
     // cekIdResiduContainString(path);
        removeFileFromDir(path_w,list_file_with_id_string);
        */
      //  createEpilist();
       }
    public static void extractASA(){
        CreateComplexStructure ccs = new CreateComplexStructure();
        
    }
    /*
    public static void createEpilist(){
        CreateComplexStructure ccs = new CreateComplexStructure();
        String path_epi ="src/dataset/";
        String list_file = "list_antigen_train_seppa3.txt";
        String fn_epilist = "epilist_seppa3.txt";              
        //ccs.loadStructureFromListFile(list_file,path_epi,fn_epilist);
        ccs.loadStructureFromContact(list_file, path_epi, fn_epilist);
        ArrayList<Complex> listc = ccs.getComplexlist();
        for(Complex c: listc)
            System.out.println(c.getComplexName()+ "\t"+ c.getMinAtomNo() + "\t"+ c.getMaxAtomNo());
    }
    */
    public static void removeFileFromDir( String path, ArrayList<String> list_file_with_id_string){
        for(String s: list_file_with_id_string){
            File currentFile = new File(path+s);
            currentFile.delete();
        }
    }
        
    }
    
   
    

