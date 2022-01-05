/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CBCTop;

import static CBCTop.AntigenGraphClustering.getFileListFromFolder;
import static CBCTop.AntigenGraphClustering.readClusterFromFile;
import static CBCTop.AntigenGraphClustering.readEpiListFromFile;
import static CBCTop.ConfBCTop.doClassificationOnResiduLevel;
import static CBCTop.TestGetDataSet.loadListContact;
import ComplexFeature.AAIndex;
import static ComplexFeature.AminoAcidPairByFuncSubgrup.findEpiFilenameinlist;
import ComplexFeature.ExposedOnComplex;
import ComplexStructure.Aminoacid;
import ComplexStructure.Atom;
import ComplexStructure.Complex;
import ComplexStructure.ExposedResidu;
import ExtractData.ExtractFolder;
import ExtractDataTraining.ExtractSphereExposure;
import ExtractDataTraining.GraphAntigenExposedAtom;
import static ExtractDataTraining.GraphAntigenExposedAtom.loadExposedAtom;
import static ExtractDataTraining.GraphAntigenExposedAtom.loadExposedAtomResiduLevel;
import ExtractFromFile.IO;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import jsat.classifiers.CategoricalData;
import jsat.classifiers.ClassificationDataSet;
import jsat.classifiers.ClassificationModelEvaluation;
import jsat.classifiers.Classifier;
import jsat.classifiers.DataPoint;
import jsat.classifiers.DataPointPair;
import jsat.classifiers.bayesian.NaiveBayes;
import jsat.classifiers.imbalance.SMOTE;
import jsat.classifiers.svm.DCSVM;
import jsat.classifiers.svm.LSSVM;
import jsat.classifiers.trees.DecisionTree;
import jsat.classifiers.trees.TreePruner.PruningMethod;
import static jsat.classifiers.trees.TreePruner.PruningMethod.REDUCED_ERROR;
import jsat.datatransform.Imputer;
import jsat.datatransform.PCA;
import jsat.datatransform.featureselection.ReliefF;
import jsat.distributions.kernels.PolynomialKernel;
import jsat.distributions.kernels.RBFKernel;
import jsat.io.ARFFLoader;
import jsat.linear.DenseVector;
import jsattest.ClassificationDataSetUpdate;
import jsattest.CluSMOTE;
import jsattest.CluSMOTEBagDTModel;
import jsattest.ClusterbasedSampling;
import static jsattest.JsatTest.prepareDataSet;
import mcl.MCL;




/**
 *
 * @author dell
 */
public class gpMCLCluSTree {
    public static void main(String[]args) throws IOException{
        //generateClusterByMCL();
        //CreateDataset();        
        //classifyClusterMCL();
        List<Integer> indices = new ArrayList();
         List<Integer> sumber = new ArrayList();
         for(int i=0;i<100;i++){
             sumber.add(i);
         }
         
         int persen = (int) Math.floor(indices.size()*0.3);
         Random r = new Random();
         for(int i=0;i< 30;i++){
             int idx = r.nextInt(100);
             System.out.print(idx+"\t");
         }
         for(Integer i: sumber){
            // System.out.print(i+", ");
         }
         
        
    }
    public static void classifyClusterWithAntigenicFeature(ClassificationDataSet dataSet){
        dataSet.applyTransform(new Imputer(dataSet));
         //System.out.println(new_cds.getDataPoint(0).numNumericalValues());
         List<Integer> indices = new ArrayList();
         List<String> sumber = new ArrayList();
         for(int i=0;i<dataSet.size();i++){
             sumber.add(String.valueOf(i));
         }
         
         int persen = (int) Math.floor(indices.size()*0.3);
         Random r = new Random();
         for(int i=0;i< persen;i++){
             int idx = r.nextInt(dataSet.size());
             indices.add(idx);
             if(sumber.contains(String.valueOf(idx))){
                 sumber.remove(String.valueOf(idx));
             }
         }
         List<Integer> sumberInt = new ArrayList();
         for(String i: sumber){
             sumberInt.add(Integer.valueOf(i));
         }
         
        
        CategoricalData categori = dataSet.getPredicting();
        List<DataPointPair<Integer>> dppl = dataSet.getAsDPPList();
        ClassificationDataSetUpdate cds_update = new ClassificationDataSetUpdate(dppl,categori);
        ClassificationDataSet test = cds_update.getSamples(indices);
        ClassificationDataSet train = cds_update.getSamples(sumberInt);
        
        
        DecisionTree classifier = new DecisionTree();
        classifier.setMaxDepth(10);
        //classifier.train(train);
        
        
        ClassificationModelEvaluation cme = new ClassificationModelEvaluation(classifier, train);
        //Random r = new Random();
        //cme.evaluateCrossValidation(5, r);
        //cme.prettyPrintConfusionMatrix();
        cme.evaluateTestSet(test);
        cme.prettyPrintConfusionMatrix();
        
    }
    public static void classifyClusterMCL() throws IOException{
    String filetrain = "src/dataset/andersen_exc_3H42_A0.01.arff";
   // String filetrain = "src/dataset/datasetgao_v2.arff";
    int mode = 1;
    //System.out.println(filetrain +"\n");
    //File fileTrain = new File(filetrain);
    //melakukan hdbscan setelah ditransformasi
    File file = new File(filetrain);  
    ClassificationDataSet trainSet = ARFFLoader.loadArffFile(file).asClassificationDataSet(0);
    trainSet.applyTransform(new Imputer(trainSet));//impute missing values in the dataset
    //ClassificationDataSet trainSet = prepareDataSet(filetrain,1);
    //ReliefF relieff = new ReliefF(25);
    PCA relieff = new PCA(trainSet,25);
    //relieff.fit(trainSet);
    for(int i=0; i<trainSet.size();i++){
         DataPoint p = relieff.transform(trainSet.getDataPoint(i));
         trainSet.setDataPoint(i, p);
    }
    
    //lakukan clustring disini
    int num_fitur = trainSet.getNumFeatures();
        for (int i =0;i<1;i++)
            System.out.println(trainSet.getNumericName(i));
        int num_class = trainSet.getClassSize();
        System.out.println("num_fitur:"+ num_fitur + " num_class: "+ num_class);
          
        List<DataPoint> listdata0 = trainSet.getDataPoints();
          /*
          for(int i=0; i<100;i++)
            System.out.println(listdata0.get(i).toString());
          */
        CategoricalData categori = trainSet.getPredicting();
        System.out.println("categori: "+categori.getCategoryName());
        List<DataPointPair<Integer>> dppl = trainSet.getAsDPPList();
          /*
          // pasangan data point dan datapair
          for (int i=0;i<dppl.size();i++){
              if (dppl.get(i).getPair()==0){
                  System.out.println("idx: "+ i);
                  System.out.println("dp" + dppl.get(i).getDataPoint());
                  System.out.println("pair" + dppl.get(i).getPair());
              }
          }
          */
        ClassificationDataSetUpdate cds_update = new ClassificationDataSetUpdate(dppl,categori);
        ClassificationDataSet dataset_minority = cds_update.getMinorityClassSample();
        ClassificationDataSet dataset_majority = cds_update.getMajorityClassSample();
          
          //System.out.println("sampel minority: "+cds_update.getIdxFromMinorityClass().size());
          //System.out.println("sampel majority: "+cds_update.getIdxFromMajorityClass().size());
          
          System.out.println("dataset minority: "+dataset_minority.size());
          System.out.println("dataset majority: "+dataset_majority.size());
          ClassificationDataSet new_cds= null;
          ClusterbasedSampling cbs = new ClusterbasedSampling();
          new_cds = cbs.pamSampling(cds_update, 2);
              System.out.println("datasetawal:"+ cds_update.size());
              System.out.println("datasetbaru:"+ new_cds.size());
          
          //System.out.println(new_cds.getDataPoint(0).numNumericalValues());
    DecisionTree classifier = new DecisionTree();
    classifier.setMaxDepth(5);
    
    //classifier.setPruningMethod(PruningMethod.ERROR_BASED);
    //Classifier classifier = new NaiveBayes();
    //PolynomialKernel poly = new PolynomialKernel(3);
    //RBFKernel rbf = new RBFKernel(0.01);
    //LSSVM lssvm = new LSSVM(rbf);
    //DCSVM dcsvm = new DCSVM(rbf);
    
   SMOTE smote = new SMOTE(classifier);
    smote.train(new_cds);
    //CluSMOTEBagDTModel cusbag = new CluSMOTEBagDTModel(classifier);
    //ClassificationModelEvaluation cme = new ClassificationModelEvaluation(classifier, new_cds);
   // ClassificationModelEvaluation cme = new ClassificationModelEvaluation(smote, new_cds);
   // ClassificationModelEvaluation cme = new ClassificationModelEvaluation(dcsvm, new_cds);
   //classifier.train(new_cds);
    Random r = new Random();
   // cme.evaluateCrossValidation(5, r);
    //cme.prettyPrintConfusionMatrix();
    
    //buat structure complex
    String path = "src/dataset/gao/chain/";
    String fileseppa = "list_antigen_train_seppa3.txt";
    String path_mcl = "src/dataset/gao/graph/mcl/bf0_75/";
    ArrayList<String> listusedfile = getFileToCreateDataSet(fileseppa, path);
    System.out.println("jumlah file digunakan: "+listusedfile.size());
    for (int i=0;i<listusedfile.size();i++){
        String fn = listusedfile.get(i);
        String pdb = fn.substring(0, 4);
        String chain = fn.substring(4);
        System.out.println(i+": "+pdb+chain);
        fn = fn +"_6.0.txt_B-MCL-0.75_SC_Level-2_numT-4_conv-norm.assign";
        ArrayList<String> clusters = readClusterFromFile(path_mcl + fn);
        if(clusters.size()>0){
            //set residu terekspose pada complex
            ComplexPreparation cp = new ComplexPreparation(path,pdb,chain);
            for(String s: clusters){
                String id = s.substring(13).trim();
                cp.getComplex().getAminoById(id).setAsExposed(true);
                //System.out.println(id+ ": "+cp.getComplex().getAminoById(id).isAsExposed());
            }
            String path_epi = "src/dataset/gao/epitopes/";
            String fn_epi = getFilename(path_epi,pdb,chain);
            
            if(fn_epi.length()>0){
                fn_epi = path_epi.concat(fn_epi);
                System.out.println("fn_epi: "+fn_epi);
                ArrayList<String> epitope = readEpiListFromFile(fn_epi);
                System.out.println("epi: "+epitope.size());
                if(epitope !=null){
                    for(String epi:epitope){
                    System.out.println(epi);
                    cp.getComplex().getAminoById(epi).setAsEpitope(true);
                    System.out.println(epi+ ": "+cp.getComplex().getAminoById(epi).isAsEpitope());
                    }
                }
                //klasifikasi menggunakan level residu pada dataset didalam cluster
                //ekstraksi fitur dari complex berdasarkan identifikasi sebagai residu terekspose
                String path_psaia = "src/dataset/gao/psaia/";
                String akhiran = chain+"_202112081246_unbound";
                ArrayList<PairResiduIdFeature> pidf= cp.setComplexFeatures(path_psaia,akhiran);
                System.out.println("pidf: "+pidf.size());
                //ArrayList<String> output= ConfBCTop.classify(pidf,classifier);
                //ArrayList<String> output = ConfBCTop.classify(pidf, classifier, relieff);
               // ArrayList<String> output= ConfBCTop.classifybySmote(pidf,smote, relieff);classifybySmotePCA
               ArrayList<String> output= ConfBCTop.classifybySmotePCA(pidf,smote, relieff);
                //System.out.println("output: "+ output.size());
                //for(String s: output)
                //    System.out.println(s);
                calculateConfMatrix(output, cp.getComplex());
            }          
        
        }
        
    }
    
    }
    public static void generateClusterByMCL() throws IOException{
    String filetrain = "src/dataset/andersen_exc_3H42_A0.01.arff";
    //System.out.println(filetrain +"\n");
    //File fileTrain = new File(filetrain);
    
    ClassificationDataSet trainSet = prepareDataSet(filetrain,1);
    ReliefF relieff = new ReliefF(30);
    relieff.fit(trainSet);
    for(int i=0; i<trainSet.size();i++){
         DataPoint p = relieff.transform(trainSet.getDataPoint(i));
         trainSet.setDataPoint(i, p);
    }
          //System.out.println(new_cds.getDataPoint(0).numNumericalValues());
    Classifier classifier = new DecisionTree();
    SMOTE smote = new SMOTE(classifier);
    smote.train(trainSet);
    
    //buat structure complex
    String path = "src/dataset/gao/chain/";
    String fileseppa = "list_antigen_train_seppa3.txt";
    String path_mcl = "src/dataset/gao/graph/mcl/bf0_75/";
    ArrayList<String> listusedfile = getFileToCreateDataSet(fileseppa, path);
    System.out.println("jumlah file digunakan: "+listusedfile.size());
    for (int i=0;i<listusedfile.size();i++){
        String fn = listusedfile.get(i);
        String pdb = fn.substring(0, 4);
        String chain = fn.substring(4);
        System.out.println(i+": "+pdb+chain);
        fn = fn +"_6.0.txt_B-MCL-0.75_SC_Level-2_numT-4_conv-norm.assign";
        ArrayList<String> clusters = readClusterFromFile(path_mcl + fn);
        if(clusters.size()>0){
            //set residu terekspose pada complex
            ComplexPreparation cp = new ComplexPreparation(path,pdb,chain);
            for(String s: clusters){
                String id = s.substring(13).trim();
                cp.getComplex().getAminoById(id).setAsExposed(true);
                //System.out.println(id+ ": "+cp.getComplex().getAminoById(id).isAsExposed());
            }
            String path_epi = "src/dataset/gao/epitopes/";
            String fn_epi = getFilename(path_epi,pdb,chain);
            
            if(fn_epi.length()>0){
                fn_epi = path_epi.concat(fn_epi);
                System.out.println("fn_epi: "+fn_epi);
                ArrayList<String> epitope = readEpiListFromFile(fn_epi);
                System.out.println("epi: "+epitope.size());
                if(epitope !=null){
                    for(String epi:epitope){
                    System.out.println(epi);
                    cp.getComplex().getAminoById(epi).setAsEpitope(true);
                    System.out.println(epi+ ": "+cp.getComplex().getAminoById(epi).isAsEpitope());
                    }
                }
                //klasifikasi menggunakan level residu pada dataset didalam cluster
                //ekstraksi fitur dari complex berdasarkan identifikasi sebagai residu terekspose
                String path_psaia = "src/dataset/gao/psaia/";
                String akhiran = chain+"_202112081246_unbound";
                ArrayList<PairResiduIdFeature> pidf= cp.setComplexFeatures(path_psaia,akhiran);
                System.out.println("pidf: "+pidf.size());
                //ArrayList<String> output= ConfBCTop.classify(pidf,classifier);
                //ArrayList<String> classifybySmote(ArrayList<PairResiduIdFeature> pridf, SMOTE smote, ReliefF relieff)
                ArrayList<String> output= ConfBCTop.classifybySmote(pidf,smote, relieff);
                //System.out.println("output: "+ output.size());
                //for(String s: output)
                //    System.out.println(s);
                calculateConfMatrix(output, cp.getComplex());
            }          
        
        }
        
    }
    
    }
    public static ArrayList<ExposedOnComplex> getECOListfromDataset(){
        GraphAntigenExposedAtom ga = new GraphAntigenExposedAtom();
        String path_in_complex = "src/dataset/gao/chain/";
        String path_in_atom = "src/dataset/gao/antigens/"; //berisi pdbid+chain dan pdbid+chain
        String path_in_res = "src/dataset/gao/antigens/res/"; //berisi pdbid+chain dan pdbid+chain+res
        String path_epi = "src/dataset/gao/epitopes/";
        ArrayList<String> epi_file_list = getFileListFromFolder(path_epi, "");
        System.out.println("contoh epi file: "+epi_file_list.get(0));
        ArrayList<ExposedOnComplex> eoclist = new ArrayList();
        try{
            FileReader fr = null;
            File folder = new File(path_in_res);
            File[] listOfFiles = folder.listFiles();
            System.out.println("jumlah listoffile: "+ listOfFiles.length);
            for (int i = 0; i < listOfFiles.length; i++) {
                if (listOfFiles[i].isFile()) {   
                    String fileNameRes = listOfFiles[i].getName();
                    
                    String pdbid = fileNameRes.substring(0,4);
                    String chain = fileNameRes.substring(4,5);
                    System.out.println("chain: "+chain);
                    String fileNameAtom = fileNameRes.substring(0, 5);
                    IO io = new IO();
                    AAIndex aai;
                    aai= io.extractAAIndexFromFile();
                    aai.composeAAI();
                    System.out.println("pdb: "+pdbid + ":"+chain );
                    try{
                        Complex c = new Complex(path_in_complex,pdbid,chain);
                        ExposedOnComplex eoc = new ExposedOnComplex();
                        eoc.setPdbid(pdbid);
                        eoc.setChain(chain);
                        
                        System.out.println("complex: "+c.getComplexName()+ c.getAminoacidVector().size());
                        String epi_file_name = path_epi + findEpiFilenameinlist(epi_file_list,pdbid,chain);
                        System.out.println("epi_file_name"+epi_file_name);
                        ArrayList<String> epitope = readEpiListFromFile(epi_file_name);
                        //baca dari file exposed atom dan exposed residu
                        ArrayList<Atom> list_atom_res = loadExposedAtomResiduLevel(path_in_res, fileNameRes);
                        System.out.println("list_atom_res: "+ list_atom_res.get(0).getX()+" "+list_atom_res.get(0).getY()+" "+list_atom_res.get(0).getZ());
                        ArrayList<Atom> list_atom_atom = loadExposedAtom(path_in_atom, fileNameAtom);
                        System.out.println("list_atom_atom: "+ list_atom_atom.get(0).getX()+" "+list_atom_atom.get(0).getY()+" "+list_atom_atom.get(0).getZ()+" "+list_atom_atom.get(0).getAtomName() );
                        ga.concatResiduInfo(list_atom_res,list_atom_atom);
                        //mulai melakukan extract exposedoncomplex
                        ArrayList<ExposedResidu> listExposedR= new ArrayList();
                        ArrayList<Aminoacid> listAminoacid= new ArrayList();
                        boolean aaResIsSet = false;
                        for(Atom atom : list_atom_res){
                            String res_id = atom.getResidueID();
                            if(aaResIsSet == false){
                                Aminoacid aa = c.getAminoById(res_id);
                                c.getAminoById(res_id).setAsExposed(true);
                                if(epitope.contains(res_id))
                                    c.getAminoById(res_id).setAsEpitope(true);
                                aaResIsSet = true;
                                ExtractSphereExposure ese = new ExtractSphereExposure();
                                ese.calculateSE(c.getAminoById(res_id), c);
                                List<Double> atrib = new ArrayList();
                                atrib.add(aa.getRsaTien2013());            
                                atrib.addAll(aa.getpASA().getParamASAD());
                                atrib.add(aa.getCaBfactor());
                                atrib.add(aa.getBfactor());
                                atrib.add(aa.getLo());
                                atrib.addAll(aa.getPsai().getListPSAI());
                                double[] withoutNan =aai.getAaiWithoutNanData(aa.letterToNum());
                                ArrayList<Double> aaid = new ArrayList();
                                for(int ii=0;ii<withoutNan.length;ii++){
                                    aaid.add(withoutNan[i]);
                                }
                                atrib.addAll(aaid);
                                DenseVector vec = new DenseVector(atrib);
                                DataPoint p= new DataPoint(vec);
                                ExposedResidu expR = new ExposedResidu(aa,p);
                                listExposedR.add(expR);
                                listAminoacid.add(aa);
                                
                            }                           
                          System.out.println("eoc"+ eoc.getPdbid()+eoc.getChain()+eoc.getListExposedR().size());
                        }
                        eoc.setListExposedR(listExposedR);
                        eoc.setListAA(listAminoacid);
                        eoclist.add(eoc);
        
                    }catch(Exception e){
                        System.out.println(e.getMessage());
                    }
                    
                    
                }
            }
        }catch(Exception e){
            System.out.println("exception in generate ga:" +e.getMessage());
        }

               
        return eoclist;
        
    }
    public static boolean lookUpOnArrayList(ArrayList<String> str, String s){
        boolean found = false;
        for(String ss: str){
            if (ss.equals(s)){
                found = true;
                return found;
            }
        }
        return found;
    }
    public static void CreateDataset() throws IOException{
       IO io = new IO();
        GraphAntigenExposedAtom ga = new GraphAntigenExposedAtom();
        String path_in_complex = "src/dataset/gao/chain/";
        String path_in_atom = "src/dataset/gao/antigens/"; //berisi pdbid+chain dan pdbid+chain
        String path_in_res = "src/dataset/gao/antigens/res/"; //berisi pdbid+chain dan pdbid+chain+res
        String path_epi = "src/dataset/gao/epitopes/";
        String path_fitur = "src/dataset/datasetgao_v2.arff";
        ArrayList<String> epi_file_list = getFileListFromFolder(path_epi, "");
        //baca file pdb
        String path = "src/dataset/gao/chain/";
        String path_psai = "src/dataset/gao/psaia/";
        String end_psai = "_202112081246_unbound";
        ExtractFolder ef = new ExtractFolder();
        File f = new File(path_epi);
        ArrayList<String> fns = ef.listAllFiles(f);
        ArrayList<String>fns_part = new ArrayList();
        fns_part.addAll(fns.subList(0, 3));
        ArrayList<String> fiturs = new ArrayList();
        for(String fn: fns_part){
            String fileNameRes = fn;
            
            String pdb = fn.substring(0,4);
            String chain = fn.substring(fn.length()-1);
            String fileNameAtom = pdb+chain;
            ComplexPreparation cp = new ComplexPreparation(path, pdb,chain);
           // String epi_file_name = path_epi + findEpiFilenameinlist(epi_file_list,pdb,chain);
            System.out.println(pdb +" "+chain);
             ArrayList<String> epitope = readEpiListFromFile(path_epi+fn);
                        //baca dari file exposed atom dan exposed residu
            ArrayList<Atom> list_atom_res = loadExposedAtomResiduLevel(path_in_res, pdb+chain+"res");
            System.out.println("list_atom_res: "+ list_atom_res.get(0).getX()+" "+list_atom_res.get(0).getY()+" "+list_atom_res.get(0).getZ());
            System.out.println("fileNameAtom: "+fileNameAtom);
            ArrayList<Atom> list_atom_atom = loadExposedAtom(path_in_atom, fileNameAtom);
            System.out.println("list_atom_atom: "+ list_atom_atom.get(0).getX()+" "+list_atom_atom.get(0).getY()+" "+list_atom_atom.get(0).getZ()+" "+list_atom_atom.get(0).getAtomName() );
            ga.concatResiduInfo(list_atom_res,list_atom_atom);
            
            for(Atom atom : list_atom_res){
                String res_id = atom.getResidueID();
                //System.out.println("res_id: "+res_id+ "\n"+"epi:");
                Aminoacid aa = cp.getComplex().getAminoById(res_id);
                if(aa!=null){
                    cp.getComplex().getAminoById(res_id).setAsExposed(true);  
                    cp.getComplex().getAminoById(res_id).getAtomById(atom.getAtomNo()).setIsexposed(true);
                    if(lookUpOnArrayList(epitope,res_id)){
                           cp.getComplex().getAminoById(res_id).setAsEpitope(true);
                           //System.out.println("res_id: "+res_id);
                           
                    }                                                 
                }
                
           
            //System.out.println(fitur);
            }
             ArrayList<String> fitur = cp.calculateComplexFeatures(path_psai, chain+end_psai);
             
             io.writeArrayListStringToFile(path_fitur, fitur);
             fiturs.addAll(fitur);
        }
        System.out.println("jumlah baris: "+ fiturs.size());
   
    
    }
    public static void calculateConfMatrix(ArrayList<String> output,Complex c){
        HashMap<String,Aminoacid> hashset = new HashMap();
        ArrayList<Aminoacid> aa = c.getAminoacidVector();
        int countTP=0;
        int countFP=0;
        int countTN=0;
        int countFN=0;
        // System.out.println("ukuran s: "+ output.size());
        for(String s: output){
            
            String[]splits = s.split(":");
            String id = splits[0];
            String label = splits[1];
            String predictedState =null;
            Aminoacid a = getAminoacidFromList(id,aa);
           // System.out.println("epi= "+a.isAsEpitope());
            if(a.isAsEpitope()==true && label.equals("0")){
                predictedState= "TP";
                countTP+=1;
            }else if(a.isAsEpitope()==true && label.equals("1") ){
                predictedState= "FP";
                countFP+=1;
            }else if(a.isAsEpitope()==false && label.equals("1")){
                predictedState= "TN";
                countTN+=1;
            }else if(a.isAsEpitope()==false && label.equals("0")){
                predictedState= "FN";
                countFN+=1;
            }
            if (a.isAsEpitope()){
                System.out.println("s: "+ s);
                System.out.println("epi= "+a.isAsEpitope()+ ", label: "+ label+ "predictedState "+ predictedState);
            }
            
        }
        System.out.println("TP FP "+countTP + "\t"+ countFP);
        System.out.println("FN TN: "+countFN + "\t"+ countTN);
        
    }
    public static Aminoacid getAminoacidFromList(String id, ArrayList<Aminoacid> aa){
        Aminoacid faa= null;
        for(Aminoacid a: aa){
            if(a.getresidueID().equals(id)){
                faa = a;
                return faa;
            }
        }
        return faa;
    }
    public static ArrayList<String> getFileToCreateDataSet(String list_file, String path_gao){
        //String list_file = "list_antigen_train_seppa3.txt";
         ArrayList<Contact> contacts = loadListContact(list_file);
         //String path_gao = "src/dataset/gao/chain/";
         File folder = new File(path_gao);
         ExtractFolder ef = new ExtractFolder();
         ArrayList<String> filenames = ef.listAllFiles(folder);
         ArrayList<String>filetouse = new ArrayList();
         for(Contact c: contacts){
             String pdb= c.pdbid;
             String ag = c.Ag;
            // System.out.println(pdb+"_"+ag.toUpperCase());
             if(filenames.contains(pdb.toLowerCase()+"_"+ag.toUpperCase()))
                 filetouse.add(pdb.toLowerCase()+ag.toUpperCase());
             
         }
         return filetouse;
    }
    public static String getFilename(String path, String pdb, String c){
        IO io = new IO();
        File folder = new File(path);
        String fn = "";
        ArrayList<String> files = io.listFilesForFolder(folder);
        for(String f: files){
            String start = f.substring(0,4);
            if(start.equals(pdb)&&f.endsWith(c)){
                fn=f;
               // System.out.println("fn:"+ fn);
                break;
            }
        }
        return fn;
    }

    
}
