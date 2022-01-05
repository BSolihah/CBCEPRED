/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CBCTop;

import ComplexFeature.PSAI;
import ComplexStructure.*;
import ExtractDataTraining.ExtractProtusionIndex;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import jsat.classifiers.CategoricalResults;
import jsat.classifiers.ClassificationDataSet;
import jsat.classifiers.ClassificationModelEvaluation;
import jsat.classifiers.Classifier;
import jsat.classifiers.imbalance.SMOTE;
import jsat.classifiers.trees.DecisionTree;
import jsat.datatransform.PCA;
import jsat.datatransform.featureselection.ReliefF;
import jsat.io.ARFFLoader;
//import jsat.ARFFLoader;

import static jsattest.JsatTest.prepareDataSet;
//import jsat.preprocessing.HDBScanBasedUndersampling;
//import modified.NonSeedRandomBagging;

/**
 *
 * @author dell
 */
public class ConfBCTop {
    
    
    public static void main(String[] args) throws IOException {
    
//panggil method setpsaia
    //panggil method getdatasetfromexposed residue
   /*
    ComplexPreparation cp = new ComplexPreparation("1a2y","C",0.01);
    ArrayList<PairResiduIdFeature> pidf= cp.getPairIdNFeature();
    doClassificationOnResiduLevel(pidf);
    */
   String path ="F:\\S3\\disertasi\\program\\datasetuji\\";
   String pdbid="5GGR";
   String chain ="Z";
   // String path= args[0];
   // String pdbid = args[1];
   // String chain = args[2];
    ComplexPreparation cp = new ComplexPreparation(path,pdbid,chain,0.01);
    
    System.out.println(cp.getComplex().getComplexName());
    ArrayList<PairResiduIdFeature> pidf= cp.getPairIdNFeature(path,chain);
    System.out.println(pidf.size());
    
    //ArrayList<String> output= doClassificationOnResiduLevel(pidf);
    //System.out.println(output);
   // String pathout="src\\dataset\\";
    //writeOutputToFile(path,pdbid,chain,output);
    
    
    }
    
    public static void showExposedResidueAllVersion(){
        //read the pdbid chain antigen, file psaia,list pdbid chain
    ArrayList<String> pdbid_list= new ArrayList();
    ArrayList<String> chain_list=new ArrayList();
    loadListFromFile("list_antigen_andersen.txt",pdbid_list,chain_list);
    
    //untuk setiap pdbdi, chain create obyek kelas ComplexPreparation
    for(int i=0;i<pdbid_list.size();i++){
        System.out.println("pdbid:"+ pdbid_list.get(i));
        ComplexPreparation cp = new ComplexPreparation(pdbid_list.get(i), chain_list.get(i),0.01);
        System.out.println("residu:"+ cp.getComplex().getAminoacidVector().size());
        cp.setExposedResidueFromPSAIA();
        System.out.println(cp.getComplex().getChainByChianID(chain_list.get(i).charAt(0)).getDatasetFromExposedResidu());
       // cp.getComplex().getChainByChianID(chain_list.get(i).charAt(0)).getDatasetFromExposedResidu();
    }
    }
    public static void loadListFromFile(String list_file, ArrayList<String> pdbid_list,ArrayList<String> chain_list ){
    
    //baca file list
    //masukkan dalam list pdbid dan chain
    try{
        if (new File("src/dataset/"+list_file).exists()) {
            String dataline = null;
            FileReader fr = null;
            fr = new FileReader(new File("src/dataset/"+list_file));
                BufferedReader br = new BufferedReader(fr);
                dataline = br.readLine();                
                while (dataline != null)
                {
                   // System.out.println(dataline);
                    String[] splitdl = dataline.split(":"); 
                    if(splitdl.length==2){
                        pdbid_list.add(splitdl[0].trim());
                        chain_list.add(splitdl[1].trim());      
                    }
                    dataline = br.readLine();
                }
                br.close();
                
        }
    }catch(IOException e){
            e.printStackTrace();
        }
    
    }
       
    public static void getFeatureOfResidue(ArrayList<PairResiduIdFeature> featureVector, String pdbid,String chain)throws IOException{
        ComplexPreparation cp = new ComplexPreparation(pdbid, chain,0.01);
    //System.out.println("vektor aa:"+ cp.getComplex().getAminoacidVector().size());
    cp.setExposedResidueFromPSAIA();
    cp.setLogOdd();
    cp.setSphereExposure();
    //System.out.println(cp.getComplex().getChainByChianID(chain.charAt(0)).getDatasetFromExposedResidu());
    
    ArrayList<Aminoacid> listAA = cp.getComplex().getAminoacidVector();
    for(Aminoacid aa:listAA){
        if (aa.getRsaTien2013()>=0.01){
                     
            String aaidx =  cp.calculateAAI(aa);
            //String s = aa.getRsaTien2013()+","+aa.getpASA().printParamASA()+aa.getCaBfactor()+","+aa.getBfactor()+","+aa.getLo()+","+aa.getPsai().printPSAIParam()+aaidx+aa.isAsEpitope();
            System.out.println(aa.getresidueID()+";"+aa.getresidueName());
            System.out.println(aa.getRsaTien2013());
            System.out.println(aa.getpASA().printParamASA());
            System.out.println(aa.getCaBfactor());            
            System.out.println(aaidx);
        }
    }
    }
   
    //create complex and set the feature vector
   
    public static void writeOutputToFile(String path,String pdbid,String chain,ArrayList<String> output) throws IOException{
    BufferedWriter bw = null;
     FileWriter fw = null;
     String dataline;
     String output_path=path+pdbid+chain+"_output.txt";
     System.out.println(output_path);
     try{
         
         fw = new FileWriter(new File(output_path),true);
         bw = new BufferedWriter(fw);
         for(String s:output){
             System.out.println(s);
             bw.write(s+"\n");
         }
         
        // System.out.println("jumlah chain: "+ c.size());
         
     }catch(Exception e){
                    System.out.println("can't write");
     }finally{
         bw.close();
     }
     
    }
    
    public static SMOTE buildModel(ClassificationDataSet trainSet){
        //String filetrain = "src/dataset/andersen_exc_3H42_A0.01.arff";
        //System.out.println(filetrain +"\n");
        //File fileTrain = new File(filetrain);
        //proses ekstraksi data
        //ClassificationDataSet dataSet = ARFFLoader.loadArffFile(fileTrain).asClassificationDataSet(0);
        //dataSet.applyTransform(new Imputer(dataSet));
        Classifier classifier = new DecisionTree();
         //ClassificationDataSet trainSet = null;
       // ReliefF relieff = null;
       //prepareDataSet(String filepath, int mode, ClassificationDataSet new_cds, ReliefF relieff)
        //prepareDataSet(filetrain, 1, trainSet, relieff);
        //prepareDataSet(filetrain, 1, dataSet, relieff);
        System.out.println("trainset: "+ trainSet.size());
        //ClassificationDataSet trainSet = prepareDataSet(filetrain, 1);
        SMOTE smote = new SMOTE(classifier);
        //classifier.train(trainSet);
        smote.train(trainSet);
        return smote;
    }
    public static ArrayList<String> classifybySmote(ArrayList<PairResiduIdFeature> pridf, SMOTE smote, ReliefF relieff){
        ArrayList<String> output=new ArrayList();
        for(PairResiduIdFeature f:pridf){
            
            //relieff.fit(f.getFeature());
            CategoricalResults cr= smote.classify(relieff.transform(f.getFeature()));
            
         //   CategoricalResults cr= bag.classify(f.getFeature());
           //set hasil prediksi pada PairResiduIdFeature
           f.setPredictedCategory(cr);
            output.add((f.getResiduId()+":"+ cr.mostLikely()+":"+ cr.getProb(cr.mostLikely())));
            
        }
        return output;
    }
    public static ArrayList<String> classifybySmotePCA(ArrayList<PairResiduIdFeature> pridf, SMOTE smote, PCA relieff){
        ArrayList<String> output=new ArrayList();
        for(PairResiduIdFeature f:pridf){
            
            //relieff.fit(f.getFeature());
            CategoricalResults cr= smote.classify(relieff.transform(f.getFeature()));
            
         //   CategoricalResults cr= bag.classify(f.getFeature());
           //set hasil prediksi pada PairResiduIdFeature
           f.setPredictedCategory(cr);
            output.add((f.getResiduId()+":"+ cr.mostLikely()+":"+ cr.getProb(cr.mostLikely())));
            
        }
        return output;
    }
    public static ArrayList<String> classify(ArrayList<PairResiduIdFeature> pridf, Classifier classifier, PCA relieff){
        ArrayList<String> output=new ArrayList();
        for(PairResiduIdFeature f:pridf){
           
            //relieff.fit(f.getFeature());
            CategoricalResults cr= classifier.classify(relieff.transform(f.getFeature()));
            
         //   CategoricalResults cr= bag.classify(f.getFeature());
           //set hasil prediksi pada PairResiduIdFeature
           f.setPredictedCategory(cr);
            output.add((f.getResiduId()+":"+ cr.mostLikely()+":"+ cr.getProb(cr.mostLikely())));
            
        }
        return output;
    }
    public static ArrayList<String> doClassificationOnResiduLevel(ArrayList<PairResiduIdFeature> pridf, ReliefF relieff ){
        
        
                    
//String filetrain = "src/rescbc/andersen_all_0.01.arff";
        String filetrain = "src/dataset/andersen_exc_3H42_A0.01.arff";
        System.out.println(filetrain +"\n");
        File fileTrain = new File(filetrain);
        //proses ekstraksi data
        //ClassificationDataSet dataSet = ARFFLoader.loadArffFile(fileTrain).asClassificationDataSet(0);
        //dataSet.applyTransform(new Imputer(dataSet));
        Classifier classifier = new DecisionTree();
        ClassificationDataSet trainSet = null;
        //prepareDataSet(String filepath, int mode, ClassificationDataSet new_cds, ReliefF relieff)
        // prepareDataSet(filetrain, 1, trainSet, relieff);
        // ClassificationModelEvaluation modelEvaluation = new ClassificationModelEvaluation(classifier, trainSet);
        // modelEvaluation.evaluateTestSet(testSet);
        //modelEvaluation.prettyPrintConfusionMatrix(); 
        //System.out.println("training time: "+modelEvaluation.getTotalTrainingTime());
       // System.out.println("classification time: "+modelEvaluation.getTotalClassificationTime());
        //HDBScanBasedUndersampling hdbU = new HDBScanBasedUndersampling();
        // ClassificationDataSet trainSet= hdbU.clusterDBScanSMOTE(dataSet, 2);
        //CusBagging bag = new CusBagging(classifier);
       //bag.train(dataSet);
        classifier.train(trainSet);
       ArrayList<String> output=new ArrayList();
        for(PairResiduIdFeature f:pridf){
            CategoricalResults cr= classifier.classify(f.getFeature());
            
         //   CategoricalResults cr= bag.classify(f.getFeature());
           //set hasil prediksi pada PairResiduIdFeature
           f.setPredictedCategory(cr);
            output.add((f.getResiduId()+":"+ cr.mostLikely()+":"+ cr.getProb(cr.mostLikely())));
            
        }
        return output;
               
    }
    
    public static void extractFeatureOfComplex(Complex c, String chain){
        //siapkan file yang dibutuhkan: file pdb, file tbl, file psai
        
    }
}
