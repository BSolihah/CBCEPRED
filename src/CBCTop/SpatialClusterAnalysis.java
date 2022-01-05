/*
 * extract exposedResidu pada complex
 * klasifikasi dengan model dan assign pada exposedResidu
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CBCTop;

import ComplexFeature.ExposedOnComplex;
import ComplexStructure.Aminoacid;
import ComplexStructure.Atom;
import SpatialAnalysis.GraphPredictedEpitope;
import SpatialAnalysis.PairOfPredictedEpitope;
import SpatialAnalysis.PairIdResidu;
import SpatialAnalysis.Cluster;
import ComplexStructure.ExposedResidu;
import ComplexStructure.PredictedEpitope;
import ExtractDataTraining.ExtractNewFeature;
import java.io.File;
import jsat.utils.IntList;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import jsat.io.ARFFLoader;
import jsat.DataSet;
import jsat.SimpleDataSet;
import jsat.classifiers.CategoricalResults;
import jsat.classifiers.ClassificationDataSet;
import jsat.classifiers.Classifier;
import jsat.classifiers.DataPoint;

import jsat.classifiers.trees.DecisionTree;
import jsat.clustering.HDBSCAN;
import jsat.datatransform.Imputer;
import jsat.linear.DenseVector;

import jsat.utils.random.RandomUtil;
import static jsattest.JsatTest.prepareDataSet;


/**
 *
 * @author dell
 */
public class SpatialClusterAnalysis {
    private String exceptelement = "O";
    private String tp="TP";
    private String fp="FP";
    private String tn="TN";
    private String fn="FN";
    private ArrayList<ExposedOnComplex> eocList;
    private ArrayList<PredictedEpitope> predEpiList;

    
   
    public SpatialClusterAnalysis(String fileList){
        ExtractNewFeature enf = new ExtractNewFeature();
        
        eocList = enf.extractFeatureExposedResidues(0.01, fileList);
    }
    public ArrayList<ExposedOnComplex> getEocList() {
        return eocList;
    }
    //analysis asam amino aromatik pada tp dan tn pada cluster
    
    //setelah melakukan clustering spatial maka langkah selanjutnya yang harus dilakukan adalah
    //memilih sub cluster mana yang akan dipilih. 
    //beberapa cara memilih sub cluster:
    //1. dengan menambahkan parameter Log odd ratio pasangan residu 
    //2. dengan melihat karakteristik asam amino 
    
    //melakukan scoring terhadap subcluster
    //parameter apa yang dapat digunakan?
    //pertama cek banyaknya anggota subcluster dan jumlah tp dan fp pada masing-masing subcluster
    public ArrayList<GraphPredictedEpitope> removeSubClusterWithRareEdge(ArrayList<GraphPredictedEpitope> gpeList){
        ArrayList<GraphPredictedEpitope> newGPEList = new ArrayList();
       // boolean[] idx= new boolean[gpeList.size()];
        
        for(int i=0;i<gpeList.size();i++){
            GraphPredictedEpitope gpe= gpeList.get(i);
            if(gpeList.get(i).getEdgeList().size()<gpeList.get(i).getPeList().size()){
              //  gpeListTemp.remove(i);
            //  idx[i]=true;
            }else{
             //   idx[i]=false;
                newGPEList.add(gpe);
            }
           
        }
        
        return newGPEList;
    }
    public String cekSubClusterCondition(ArrayList<GraphPredictedEpitope> gpeList){
        String condition="";
        for(GraphPredictedEpitope gpe:gpeList){
            condition+="jumlah node: "+gpe.getPeList().size()+ " jumlah edge: "+gpe.getEdgeList().size()+"\n";
            //countTNandFPonGpe(gpe);
            gpe.printGraphPE();
            
            
        }
        return condition;
    }
    public void printSubCluster(ArrayList<GraphPredictedEpitope> gpeList){
        for(GraphPredictedEpitope gpe:gpeList){
            
            
        }
    }
    private void countTNandFPonGpe(GraphPredictedEpitope gpe){
        int counttp=0;
        int countfp = 0;
        ArrayList<PredictedEpitope> pe = gpe.getPeList();
        for(PredictedEpitope p:pe){
            if(p.getPredictedState().equalsIgnoreCase(this.tp)){
                counttp+=1;
            }else if(p.getPredictedState().equalsIgnoreCase(this.fp)){
                countfp+=1;
            }
        }
        System.out.println("TP:"+counttp+" FP:"+countfp);
    }
    private void countConfussionMetricOnGpe(GraphPredictedEpitope gpe){
        int counttp=0;
        int countfp = 0;
        int counttn=0;
        int countfn=0;
        
        ArrayList<PredictedEpitope> pe = gpe.getPeList();
        for(PredictedEpitope p:pe){
            if(p.getPredictedState().equalsIgnoreCase(this.tp)){
                counttp+=1;
            }else if(p.getPredictedState().equalsIgnoreCase(this.fp)){
                countfp+=1;
            }else if(p.getPredictedState().equalsIgnoreCase(this.tn)){
                counttn+=1;
            }else if(p.getPredictedState().equalsIgnoreCase(this.fn)){
                countfn+=1;
            }
        }
        System.out.println("TP:"+counttp+" FP:"+countfp);
    }
    public ppeOnGraphList checkPpeOnGpEList(PairOfPredictedEpitope ppe, ArrayList<GraphPredictedEpitope> graphList){
    //kondisi: 
    //1. pe1 dan pe2 ditemukan disalah satu graph
    //2. pe1 dan pe2 ditemukan di graph berbeda
    //3. salah satu dari pe1 atau pe2 ditemukan pada salah satu graph
    //4. tidak ditemukan dimanapun
    //output:
    //1. ditemukan/tidak ditemukan/disalah
    PredictedEpitope pe1= ppe.getPe1();
    PredictedEpitope pe2= ppe.getPe2();
    PairIdResidu pir = new PairIdResidu(pe1.getIdResidue(),pe2.getIdResidue(),ppe.getWeight());
    int idxpe1=-1;
    for (int i=0;i<graphList.size();i++){
        if(graphList.get(i).getPeList().contains(pe1)){
            idxpe1=i;
         //   System.out.println("idxpe1"+idxpe1);
            break;
        }
        
    }
    int idxpe2=-1;
    for (int i=0;i<graphList.size();i++){
        if( graphList.get(i).getPeList().contains(pe2)){
            idxpe2=i;
         //   System.out.println("idxpe2"+idxpe2);
            break;
        }
        
    }
    int idxpir = -1;
    for (int i=0;i<graphList.size();i++){
        if( graphList.get(i).getEdgeList().contains(pir)){
        idxpir=i;
       // System.out.println("idxpir"+idxpir);
        break;
        }
        
    }
    ppeOnGraphList pog = new ppeOnGraphList(idxpe1,idxpe2,idxpir);
    return pog;
    
    }
    public ArrayList<GraphPredictedEpitope> doSpatialClusterOnPairOfPE(ArrayList<PairOfPredictedEpitope> pelist){
        // kelompokkan pasangan pe dalam sub sub cluster dan temukan subcluster yang mengandung tp
        //buat klaster spasial untuk pasangan pe dengan jarak <=4
        ArrayList<GraphPredictedEpitope> gpeList = new ArrayList();
        
        for(int j=0;j<pelist.size();j++){
            System.out.println("j:"+j);
            if(gpeList.isEmpty()){
                //buat sub klaster baru
                GraphPredictedEpitope gpe = new GraphPredictedEpitope(pelist.get(j));
                gpeList.add(gpe);
              //  System.out.println("gpe baru ditambahkan"+ gpeList.size());
            }else{
                ppeOnGraphList pog= checkPpeOnGpEList(pelist.get(j),gpeList);
               // System.out.println(pog.idxpe1+","+pog.idxpe2+","+pog.idxpir);
                if(pog.idxpe1 == -1 &&pog.idxpe2==-1){
                    GraphPredictedEpitope gpe = new GraphPredictedEpitope(pelist.get(j));
                    gpeList.add(gpe);
                  //  System.out.println("gpe baru ditambahkan"+ gpeList.size());
                }else if(pog.idxpe1==pog.idxpe2 && pog.idxpir ==-1){
                    PairIdResidu pir = new PairIdResidu(pelist.get(j).getPe1().getIdResidue(),pelist.get(j).getPe2().getIdResidue(),pelist.get(j).getWeight());
                    gpeList.get(pog.idxpe1).addEdge(pir);
                  //  System.out.println("node baru ditambahkan ke gpeList"+pog.idxpe1);
                }else if((pog.idxpe1==-1||pog.idxpe2==-1)&& (pog.idxpe1!=pog.idxpe2)){
                    if(pog.idxpe1!=-1){
                        //tambahkan pe2 ke graph pe1
                        gpeList.get(pog.idxpe1).addPe(pelist.get(j).getPe2());
                        PairIdResidu pir = new PairIdResidu(pelist.get(j).getPe1().getIdResidue(),pelist.get(j).getPe2().getIdResidue(),pelist.get(j).getWeight());
                        gpeList.get(pog.idxpe1).addEdge(pir);
                      //  System.out.println("node dan edge baru ditambahkan");
                    }else if(pog.idxpe2!=-1){
                        //tambahkan pe1 ke graph pe2
                        gpeList.get(pog.idxpe2).addPe(pelist.get(j).getPe1());
                        PairIdResidu pir = new PairIdResidu(pelist.get(j).getPe2().getIdResidue(),pelist.get(j).getPe1().getIdResidue(),pelist.get(j).getWeight());
                        gpeList.get(pog.idxpe2).addEdge(pir);
                       // System.out.println("node dan edge baru ditambahkan");
                        
                    }
                }else if(pog.idxpe1!=pog.idxpe2 && pog.idxpe1!= -1 &&pog.idxpe2!=-1){
                        //
                       // System.out.println("jumlah subgraph"+gpeList.size());
                        gpeList.get(pog.idxpe1).addSubGraph(gpeList.get(pog.idxpe2));
                        PairIdResidu pir = new PairIdResidu(pelist.get(j).getPe1().getIdResidue(),pelist.get(j).getPe2().getIdResidue(),pelist.get(j).getWeight());
                        gpeList.get(pog.idxpe1).addEdge(pir);
                        gpeList.remove(pog.idxpe2);
                       // System.out.println("gabung graph" +", jumlah subgraph"+gpeList.size());
                }
            }
        }
        for(GraphPredictedEpitope g:gpeList)
            this.countTNandFPonGpe(g);
        return gpeList;
        
    }
    
    public ArrayList<PairOfPredictedEpitope> distanceAnalysis(ArrayList<PredictedEpitope> pe){
        //hitung jarak setiap asam amino
        //buat pasangan 
        ArrayList<Integer> p1List = new ArrayList();
        ArrayList<Integer> p2List = new ArrayList();
        for(int i=0;i<pe.size();i++){
            for(int j=i+1;j<pe.size();j++){
                if(i!=j){
                    p1List.add(i);
                    p2List.add(j);
                   // System.out.println(i+ "vs" + j);
                }
            }
        }
       
        //hitung jarak antar pasangan dan simpan dlm arrayList
        ArrayList<PairOfPredictedEpitope> pairList= new ArrayList();
        if(p1List.size()==p2List.size()){
            for(int i= 0;i<p1List.size();i++){                
                double d= distance_amino_amino(pe.get(p1List.get(i)).getAmino(),pe.get(p2List.get(i)).getAmino());
                if(d<=4){
                  //  System.out.println(pe.get(p1List.get(i)).getIdResidue()+","+pe.get(p2List.get(i)).getIdResidue()+","+pe.get(p1List.get(i)).getPredictedState()+","+pe.get(p2List.get(i)).getPredictedState()+","+ d);
                   PairOfPredictedEpitope p = new PairOfPredictedEpitope(pe.get(p1List.get(i)),pe.get(p2List.get(i)),d);
                   pairList.add(p);
                   
                }                
            }            
        }
        
        return pairList;
        
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
    public int numberCluster(int[] cluster){
        int num = -1;
        for(int i=0;i<cluster.length;i++){
            if(cluster[i]>num)
                num = cluster[i];
        }
        return num;
            
        
    }
     public DataSet createDataSetFromPredictedEpitopeList(ArrayList<PredictedEpitope> pelist){
         ArrayList<DataPoint> dpList = new ArrayList();
         for(PredictedEpitope pe: pelist){
             dpList.add(pe.getSpatialDP());
         }
         SimpleDataSet dataset = new SimpleDataSet(dpList);
         return dataset;
    }
     
    public void clusterOnPredictedEpitope(ArrayList<PredictedEpitope> pe){
        
        DataSet dataset = createDataSetFromPredictedEpitopeList(pe);
        HDBSCAN dbscan = new HDBSCAN(3);     
      // DBSCAN dbscan = new DBSCAN();
        int[] clusterResult= new int[pe.size()];
      //  clusterResult = dbscan.cluster(majClass, r, 0, clusterResult)
       clusterResult = dbscan.cluster(dataset, clusterResult);
       int num =numberCluster(clusterResult);
       System.out.println("jumlah cluster"+num);
       List<Cluster> listC = new ArrayList();
       for (int i=0;i<clusterResult.length;i++){
            Cluster c = new Cluster(i,clusterResult[i]);
            listC.add(c);
        }
        Collections.sort(listC, new Comparator<Cluster>() {
			@Override
			public int compare(Cluster c1, Cluster c2) {
				return c1.getNumber()-c2.getNumber();
			}
		});
        List<List<Integer>> clusterElemen = new ArrayList();
        for(int i=0;i<num+1;i++){
            List<Integer> listkei = new ArrayList();
            clusterElemen.add(listkei);
        }
        for(int j=0;j<listC.size();j++){
            if(listC.get(j).getNumber()>=0)
            clusterElemen.get(listC.get(j).getNumber()).add(listC.get(j).getIdx());
        }
        /*
        List<Integer> newCdsIdx= new ArrayList();
        for(int k=0;k<clusterElemen.size();k++){
            List<Integer> listkek = clusterElemen.get(k);
            int sci = listkek.size();
            int sizeMACi = Math.round((float)r*(float)minClass.getSampleSize()*((float)sci/(float)majClass.getSampleSize()));
            for(int i=0;i<sizeMACi;i++){
                int nextRand = RandomUtil.getRandom().nextInt(listkek.size()-1);
                newCdsIdx.add(listkek.get(nextRand));
            }
            
        }*/
    }
            
    public  ArrayList<PredictedEpitope> extractSpatialFeatureOnResiduTpOrFp(ArrayList<ExposedResidu> erList){
        ArrayList<PredictedEpitope> predEpiList= new ArrayList();
        System.out.println("Residu terekspose: "+erList.size());
        int countTP=0;
            int countFP=0;
            int countTN=0;
            int countFN=0;
        for(ExposedResidu e:erList){
            
            String predictedState =null;
            if(e.isIsEpitope()==true && e.getPredictedAs().mostLikely()==0){
                predictedState= "TP";
                countTP+=1;
            }else if(e.isIsEpitope()==false && e.getPredictedAs().mostLikely()==0){
                predictedState= "FP";
                countFP+=1;
            }else if(e.isIsEpitope()==false && e.getPredictedAs().mostLikely()==1){
                //predictedState= "TN";
                countTN+=1;
            }else if(e.isIsEpitope()==true && e.getPredictedAs().mostLikely()==1){
                //predictedState= "FN";
                countFN+=1;
            }
            if(predictedState== "TP" ||predictedState== "FP" ){
                PredictedEpitope predEpi = new PredictedEpitope();
                predEpi.setIdResidue(e.getIdResidue());
                predEpi.setAmino(e.getAmino());
                predEpi.setPredictedState(predictedState);
                List<Double> atrib = new ArrayList();
                atrib.add(e.getcAlpha().getX());
                atrib.add(e.getcAlpha().getY());
                atrib.add(e.getcAlpha().getZ());
                DenseVector vec = new DenseVector(atrib);
                DataPoint p= new DataPoint(vec);
                predEpi.setSpatialDP(p);
                predEpiList.add(predEpi);
            }
            
             
            
        }
        System.out.println("TP:"+ countTP+","+"TN:"+ countTN+","+"FP: "+countFP+"FN:"+ countFN);
        return predEpiList;
    }
    public void classifyAllEOC(){
        String filetrain = "src/dataset/andersen_all_0.01.arff";
       // File fileTrain = new File(filetrain);
        //proses ekstraksi data
        //ClassificationDataSet dataSet = ARFFLoader.loadArffFile(fileTrain).asClassificationDataSet(0);
        ClassificationDataSet trainSet = prepareDataSet(filetrain, 1);
        Classifier classifier = new DecisionTree();
        classifier.train(trainSet);
        for(ExposedOnComplex eoc:eocList){
            ArrayList<ExposedResidu> erList =eoc.getListExposedR();
            for(ExposedResidu e:erList){
            CategoricalResults cr= classifier.classify(e.getFeature());
           // CategoricalResults cr= bag.classify(e.getFeature());
           //set hasil prediksi pada PairResiduIdFeature
           e.setPredictedAs(cr);     
           
            }
        }
    }
    public void clasifyOnResiduLevel(ExposedOnComplex eoc){
        String filetrain = "src/rescbc/andersen_all_0.01.arff";
        System.out.println(filetrain +"\n");
        File fileTrain = new File(filetrain);
        //proses ekstraksi data
        ClassificationDataSet dataSet = ARFFLoader.loadArffFile(fileTrain).asClassificationDataSet(0);
        dataSet.applyTransform(new Imputer(dataSet));
        Classifier classifier = new DecisionTree();
        //HDBScanBasedUndersampling hdbU = new HDBScanBasedUndersampling();
        //ClassificationDataSet trainSet= hdbU.clusterDBScanSMOTE(dataSet, 2);
       // CusBagging bag = new CusBagging(classifier);
        //bag.train(dataSet);
        classifier.train(dataSet);
        ArrayList<ExposedResidu> erList =eoc.getListExposedR();
        for(ExposedResidu e:erList){
            CategoricalResults cr= classifier.classify(e.getFeature());
           // CategoricalResults cr= bag.classify(f.getFeature());
           //set hasil prediksi pada PairResiduIdFeature
           e.setPredictedAs(cr);
            
            
        }
    }
   

    
private class ppeOnGraphList{
    int idxpe1;
    int idxpe2;
    int idxpir;

        public ppeOnGraphList(int idxpe1, int idxpe2, int idxpir) {
            this.idxpe1 = idxpe1;
            this.idxpe2 = idxpe2;
            this.idxpir = idxpir;
        }
    
}
    
    
}

