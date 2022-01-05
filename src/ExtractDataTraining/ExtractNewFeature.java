/*
 * 1. create complex
 * 2. 
 */
package ExtractDataTraining;

import ComplexFeature.AAIndex;
import ComplexFeature.ExposedOnComplex;
import ComplexStructure.Aminoacid;
import ComplexStructure.Complex;
import ComplexStructure.ExposedResidu;
import ExtractData.CreateComplexStructure;
import ExtractFromFile.IO;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import jsat.classifiers.DataPoint;
import jsat.linear.DenseVector;

/**
 *
 * @author dell
 */
public class ExtractNewFeature {
    public void extractFeatureOfComplex(Complex c) throws IOException{
        ExtractLogOdds lo = new ExtractLogOdds();
        lo.calculateLogOdds();
        IO io = new IO();
        AAIndex aai;
        aai= io.extractAAIndexFromFile();
        aai.composeAAI();
        lo.calculateLOTunggal(c, c.getAntigenChainIDs().charAt(0));
        ExtractSphereExposure ese = new ExtractSphereExposure();
        for(Aminoacid aa: c.getAminoacidVector()){
            
            
        }
    }
    public void extractFeatureOfDataTest(String filename, double threshold){
        CreateComplexStructure ccs = new CreateComplexStructure(filename);
        ExtractSphereExposure ese;
        ArrayList<Complex> complexList= ccs.getComplexlist();
     ExtractLogOdds lo = new ExtractLogOdds();
     lo.calculateLogOdds();
     IO io = new IO();
     AAIndex aai;
     aai= io.extractAAIndexFromFile();
     aai.composeAAI();
     for(int i=0;i< complexList.size();i++){
         lo.calculateLOTunggal(complexList.get(i), ccs.getChainlist().get(i).charAt(0));
         ArrayList<Aminoacid> listAA= complexList.get(i).getAminoacidVector();
         BufferedWriter bw;
            FileWriter fw;
         try{
             fw = new FileWriter(new File("src/datatest/NewFeatures_Epilist"+complexList.get(i).getComplexName()+".arff"),true);
             bw = new BufferedWriter(fw);
             ese = new ExtractSphereExposure();//pastikan sphere exposure dari file kumpulan data train atau data set
             for(Aminoacid aa: listAA){
                double  rsa = aa.getRsaTien2013();
              //  System.out.println(aa.getresidueID()+": "+rsa);
                
                if (rsa > threshold){
                    ese.calculateSE(aa, complexList.get(i));
                  //  System.out.println(aa.getpASA().printParamASA());
                    String aaidx = extractAAIndexWithoutNaN(aai,aa,threshold);
                    String s = aa.getRsaTien2013()+","+aa.getCaBfactor()+","+aa.getBfactor()+","+aa.getLo()+","+aa.getPsai().printPSAIParam()+aaidx+aa.isAsEpitope();
                  //  String s = aa.getRsaTien2013()+","+aa.getpASA().printParamASA()+aa.getCaBfactor()+","+aa.getBfactor()+","+aa.getLo()+","+aa.getPsai().printPSAIParam()+aaidx+aa.isAsEpitope();
                 //   String s = aa.getRsaTien2013()+","+aa.getpASA().printCNHSE()+aa.getCaBfactor()+","+aa.getBfactor()+","+aa.getLo()+","+aa.getPsai().printPSAIParam()+aaidx+aa.isAsEpitope();
                   if(!s.equals("")){bw.write(s+"\n");}
                }
             
         
            }
             bw.close();
         }     
         catch(Exception except){
                System.out.println(except.getMessage());
            }
     }
        
    }
    public ArrayList<ExposedOnComplex> extractFeatureExposedResidues(double threshold,String list_file){
    
    //create structure
    // String list_file ="list_antigen_andersen_1.txt";
    // CreateComplexStructure ccs = new CreateComplexStructure(list_file);
     CreateComplexStructure ccs = new CreateComplexStructure();
     ccs.loadStructureFromContact(list_file);
     //System.out.println("sampai sini");
     ExtractSphereExposure ese;
     ArrayList<Complex> complexList= ccs.getComplexlist();
     //System.out.println("struktur terbentuk: "+complexList.size());
     ExtractLogOdds lo = new ExtractLogOdds();
     lo.calculateLogOdds();
     IO io = new IO();
     AAIndex aai;
     aai= io.extractAAIndexFromFile();
     aai.composeAAI();
     ArrayList<ExposedOnComplex> eocList = new ArrayList();
     for(int i=0;i< complexList.size();i++){
         lo.calculateLOTunggal(complexList.get(i), complexList.get(i).getAntigenChainIDs().charAt(0));
         ArrayList<Aminoacid> listAA= complexList.get(i).getAminoacidVector();
         ExposedOnComplex eoc = new ExposedOnComplex();
         eoc.setPdbid(complexList.get(i).getComplexName());
         eoc.setChain(complexList.get(i).getChainIDs());
         System.out.println(eoc.getPdbid()+"_"+eoc.getChain());
         try{
             ese = new ExtractSphereExposure();
             ArrayList<ExposedResidu> listExposedR= new ArrayList();
             ArrayList<Aminoacid> listAminoacid= new ArrayList();
             for(Aminoacid aa: listAA){
                double  rsa = aa.getRsaTien2013();
               // System.out.println("rsa: "+aa.getresidueID()+": "+rsa);
                if (rsa > threshold){
                    ese.calculateSE(aa, complexList.get(i));
                    List<Double> atrib = new ArrayList();
                    atrib.add(aa.getRsaTien2013());            
                    atrib.addAll(aa.getpASA().getParamASAD());
                    atrib.add(aa.getCaBfactor());
                    atrib.add(aa.getBfactor());
                    atrib.add(aa.getLo());
                    atrib.addAll(aa.getPsai().getListPSAI());
                    atrib.addAll(this.getAAIList(aa));
                    DenseVector vec = new DenseVector(atrib);
                    DataPoint p= new DataPoint(vec);
                    ExposedResidu expR = new ExposedResidu(aa,p);
                    listExposedR.add(expR);
                    listAminoacid.add(aa);
                }
            }
             eoc.setListExposedR(listExposedR);
            // eoc.setListAA(listAminoacid);
         }     
         catch(Exception except){
                System.out.println(except.getMessage());
            }
       eocList.add(eoc);
     }
    return eocList;
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
    public void extractFeatureFromSelectedID(ArrayList<Complex> complexList,ArrayList<Integer> listIdx, String filename, double threshold){
        ExtractLogOdds lo = new ExtractLogOdds();
        ExtractSphereExposure ese;
     lo.calculateLogOdds();
     IO io = new IO();
     AAIndex aai;
     aai= io.extractAAIndexFromFile();
     aai.composeAAI();
        for(int i=0;i< listIdx.size();i++){
         lo.calculateLOTunggal(complexList.get(listIdx.get(i)), complexList.get(listIdx.get(i)).getAntigenChainIDs().charAt(0));
         ArrayList<Aminoacid> listAA= complexList.get(listIdx.get(i)).getAminoacidVector();
         BufferedWriter bw;
            FileWriter fw;
         try{
            //   fw = new FileWriter(new File("src/datatraining/seppa3/tes/"+filename+"_"+threshold+".arff"),true);
            //  fw = new FileWriter(new File("src/datatraining/ren2015/"+filename+"_"+threshold+".arff"),true);
              fw = new FileWriter(new File("src/dataset/andersen/train/"+filename+"_"+threshold+".arff"),true);
            // fw = new FileWriter(new File("src/datatraining/testset/"+complexList.get(i).getComplexName()+complexList.get(i).getAntigenChainIDs().charAt(0)+threshold+".arff"),true);
             //fw = new FileWriter(new File("src/datatraining/experiment_loocv/"+"andersen_exc_3H42_BO"+threshold+".arff"),true);            
             bw = new BufferedWriter(fw);
             ese = new ExtractSphereExposure();
             int exposedCount=0;
             for(Aminoacid aa: listAA){
                double  rsa = aa.getRsaTien2013();
               // System.out.println("rsa: "+aa.getresidueID()+": "+rsa);
                
                if (rsa > threshold){
                    exposedCount+=1;
                    ese.calculateSE(aa, complexList.get(listIdx.get(i)));
                  //  System.out.println(aa.getpASA().printParamASA());
                    String aaidx = extractAAIndexWithoutNaN(aai,aa,threshold);
                  //  String s = aa.getRsaTien2013()+","+aa.getCaBfactor()+","+aa.getBfactor()+","+aa.getLo()+","+aa.getPsai().printPSAIParam()+aaidx+aa.isAsEpitope();
                    String s = aa.getRsaTien2013()+","+aa.getpASA().printParamASA()+aa.getCaBfactor()+","+aa.getBfactor()+","+aa.getLo()+","+aa.getPsai().printPSAIParam()+aaidx+aa.isAsEpitope();
                 //   String s = aa.getRsaTien2013()+","+aa.getpASA().printCNHSE()+aa.getCaBfactor()+","+aa.getBfactor()+","+aa.getLo()+","+aa.getPsai().printPSAIParam()+aaidx+aa.isAsEpitope();
                   if(!s.equals("")){bw.write(s+"\n");}
                 // System.out.println()
                }
             
         
            }
             bw.close();
         }     
         catch(Exception except){
                System.out.println(except.getMessage());
            }

         
     }
     
    
    }
    public void extractFeatureOfDS10CV(double threshold, String list_file,int part){
        //
    CreateComplexStructure ccs = new CreateComplexStructure();
    ccs.loadStructureFromListFile(list_file);
    
     ArrayList<Complex> complexList= ccs.getComplexlist();
     //System.out.println("struktur terbentuk: "+complexList.size());
     
     //buat untuk membentuk dataset training dan testing dengan pembagian 9:1
     //bagi dataset menjadi datatrai dan datatest
     ArrayList<Integer> trainIdx= new ArrayList();
     ArrayList<Integer> testIdx= new ArrayList();
     
     this.createListDataSetCV(trainIdx, testIdx, complexList.size(), part);
     System.out.println("complex: "+ complexList.size()+ "train: "+ trainIdx.size()+"test: "+ testIdx.size());
     for(Integer i:testIdx){
         System.out.print(i+"\t");
     }
     String filenameTrain = "train_ren2015"+"exc"+part;
     String filenameTest ="test_ren2015"+part;
     this.extractFeatureFromSelectedID(complexList,trainIdx, filenameTrain, threshold);
    this.extractFeatureFromSelectedID(complexList,testIdx, filenameTest, threshold);
     
    }
    //metode untuk membentuk 10-cv
    public void createListDataSetCV(ArrayList<Integer> trainIdx,ArrayList<Integer>testIdx, int dtSize, int part){
        
        
        int set = dtSize/10;//90/10=9 0;9;18;27
        for(int i=0;i<=set;i++){
            int start =i*set;
            int stop = (i+1)*set -1;
            System.out.println(start+"-->"+stop);
            if(i==part){
            //masukkan ke testset
            for(int j=start;j<=stop;j++){
                testIdx.add(j);
                
            }
            }else{
            //masukkan ke trainset
                for(int j=start;j<=stop;j++){
                trainIdx.add(j);
                }
            }
        }
        
                
        
    }
    public void extractFeatureOfWholeDS(double threshold,String list_file){
    //create structure
    // String list_file ="list_antigen_andersen_1.txt";
    // CreateComplexStructure ccs = new CreateComplexStructure(list_file);
     CreateComplexStructure ccs = new CreateComplexStructure();
     //ccs.loadStructureFromFileTest(list_file);
    // ccs.loadStructureFromFile(list_file);
    ccs.loadStructureFromContact(list_file);
    //System.out.println("sampai sini");
     ExtractSphereExposure ese;
     ArrayList<Complex> complexList= ccs.getComplexlist();
     System.out.println("struktur terbentuk: "+complexList.size());
     ExtractLogOdds lo = new ExtractLogOdds();
     lo.calculateLogOdds();
     IO io = new IO();
     AAIndex aai;
     aai= io.extractAAIndexFromFile();
     aai.composeAAI();
     for(int i=0;i< complexList.size();i++){
         lo.calculateLOTunggal(complexList.get(i), complexList.get(i).getAntigenChainIDs().charAt(0));
         ArrayList<Aminoacid> listAA= complexList.get(i).getAminoacidVector();
         BufferedWriter bw;
            FileWriter fw;
         try{
        fw = new FileWriter(new File("src/dataset/andersen/train/andersen_all_exc_LO"+"_"+threshold+".arff"),true);             
    //  fw = new FileWriter(new File("src/dataset/andersen/train/"+complexList.get(i).getComplexName()+complexList.get(i).getAntigenChainIDs().charAt(0)+"_"+threshold+".arff"),true);
//  fw = new FileWriter(new File("src/dataset/seppa3/test/"+complexList.get(i).getComplexName()+complexList.get(i).getAntigenChainIDs().charAt(0)+"_"+threshold+".arff"),true);
            //  fw = new FileWriter(new File("src/dataset/seppa3/"+"alltrain"+threshold+".arff"),true);
            //  fw = new FileWriter(new File("src/datatraining/testset_zhang2011/"+complexList.get(i).getComplexName()+complexList.get(i).getAntigenChainIDs().charAt(0)+threshold+".arff"),true);
            // fw = new FileWriter(new File("src/datatraining/testset/"+complexList.get(i).getComplexName()+complexList.get(i).getAntigenChainIDs().charAt(0)+threshold+".arff"),true);
             //fw = new FileWriter(new File("src/datatraining/experiment_loocv/"+"andersen_exc_3H42_BO"+threshold+".arff"),true);            
             bw = new BufferedWriter(fw);
             ese = new ExtractSphereExposure();
             for(Aminoacid aa: listAA){
                double  rsa = aa.getRsaTien2013();
               // System.out.println("rsa: "+aa.getresidueID()+": "+rsa);
                
                if (rsa > threshold){
                    ese.calculateSE(aa, complexList.get(i));
                  //  System.out.println(aa.getpASA().printParamASA());
                    String aaidx = extractAAIndexWithoutNaN(aai,aa,threshold);
                  //  String s = aa.getRsaTien2013()+","+aa.getCaBfactor()+","+aa.getBfactor()+","+aa.getLo()+","+aa.getPsai().printPSAIParam()+aaidx+aa.isAsEpitope();
                 //   String s = aa.getRsaTien2013()+","+aa.getpASA().printParamASA()+aa.getCaBfactor()+","+aa.getBfactor()+","+aa.getLo()+","+aa.getPsai().printPSAIParam()+aaidx+aa.isAsEpitope();
                    String s = aa.getRsaTien2013()+","+aa.getpASA().printParamASA()+aa.getCaBfactor()+","+aa.getBfactor()+","+aa.getPsai().printPSAIParam()+aaidx+aa.isAsEpitope();
                 //   String s = aa.getRsaTien2013()+","+aa.getpASA().printCNHSE()+aa.getCaBfactor()+","+aa.getBfactor()+","+aa.getLo()+","+aa.getPsai().printPSAIParam()+aaidx+aa.isAsEpitope();
                   if(!s.equals("")){bw.write(s+"\n");}
                }
             
         
            }
             bw.close();
         }     
         catch(Exception except){
                System.out.println(except.getMessage());
            }

         
     }
     
    //
    }
     public String extractAAIndexWithoutNaN(AAIndex aai,Aminoacid aa, double threshold){
        String aaiwn="";
        boolean exposedStatus = aa.getRsaTien2013()>=threshold ;
        boolean exposedStatusOfEpi = aa.isAsEpitope()&& exposedStatus;
        if(exposedStatus||exposedStatusOfEpi){
            String aaId = aa.getPdbid()+"_"+aa.getchainID()+"_"+aa.getresidueID()+"_"+aa.letterToNum();
            double[] withoutNan =aai.getAaiWithoutNanData(aa.letterToNum());
            aaiwn = aai.printArray1D(withoutNan);
           // aaisun = aaId+","+aaiwn+aa.isAsEpitope();
            
        }
        return aaiwn;
    }
}
