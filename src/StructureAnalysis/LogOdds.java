/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package StructureAnalysis;

import ComplexCalculate.UnboundDSSP;
import ComplexStructure.Aminoacid;
import ComplexStructure.Complex;
import ComplexStructure.Epitope;
import ExtractDataTraining.Benchmark;
import ExtractDataTraining.ExtractBenchmark;
import ExtractFromFile.IO;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;

/**
 *
 * @author OpenGress
 */
public final class LogOdds {
    
    private Epitope epi;
    private HashMap<String,String> sequence;
    private HashMap<Integer,String> idxSequence;
    private ArrayList<String> idxExposedResidue;
    private ArrayList<String> overlapmer;
    private ArrayList<String> overlapmerAsEpi;
    private ArrayList<String> overlapmerAsNonEpi;
    private double[][] residueFrequency ;
    private int[]backgroundFrequency;
    private int[]BF;
    double [] lo;

    public double[] getLo() {
        return lo;
    }

    public void setLo(double[] lo) {
        this.lo = lo;
    }

    
    /*
    public LogOdds(){
        this.sequence= new HashMap<String,String>();
        this.idxSequence = new HashMap<Integer,String>();
        this.doLogOddsAnalysis();
    }
    */
     //membuat daftar epitope dan daftar non epitope dari overlap mer 
    
    public int[] getBF() {
        return BF;
    }

    public void setBF(int[] BF) {
        this.BF = BF;
    }

    public int[] getBackgroundFrequency() {
        return backgroundFrequency;
    }

    public void setBackgroundFrequency(int[] backgroundFrequency) {
        this.backgroundFrequency = backgroundFrequency;
    }
    public int totalCount(int[]f){
        int sum=0;
        for(int i=0;i<f.length;i++){
               sum = sum+f[i];
        }
         return sum;
        
    }
    public void updateFrequency(int[] f, int count, double []uf){
        
        for(int i=0;i<f.length;i++){
            uf[i]= (double)f[i]/(double)count;
        }
        
    }
    public void printArray2D(int[]ar){
        for(int i=0;i<ar.length;i++){
            System.out.print(ar[i]+"\t");
        }
        System.out.println();
    }
    public void printArray2D(double[]ar){
        for(int i=0;i<ar.length;i++){
            System.out.print(ar[i]+"\t");
        }
        System.out.println();
    }
    public void printString(ArrayList<String> s){
        for(String m:s){
            System.out.println(m);
        }
    }
    public void LogOddsAnalysisWithDiffBg(double[] lo, double []loce,double []lononce){
        ExtractBenchmark b = new ExtractBenchmark();
        ArrayList<Benchmark> listb= b.extractBenchmarkFromFile();
        
        ArrayList <String> merEpi = new ArrayList();
        ArrayList<String> merexposed = new ArrayList();
        
        b.createMerEpiNonEpi(listb, merEpi, merexposed);
        
        this.countBFFromBenchmark(listb);
        
        int[]crfepi= new int[20];
        int[]bfepi = new int[20];
        this.countFrequency(merEpi, bfepi, crfepi);
        
       // this.countFrequencyCR(merEpi, crfepi);
      //  System.out.println("Frequensi epi");
       // this.printArray2D(crfepi);
       // System.out.println("Frequensi background epi");
       // this.printArray2D(bfepi);
        int[]crfnonepi= new int[20];
        int[]bfnonepi = new int[20];
        this.countFrequency(merexposed, bfnonepi, crfnonepi);
       // this.countFrequencyCR(merexposed, crfnonepi);
      //  System.out.println("Frequensi nonepi");
      //  this.printArray2D(crfnonepi);
      //  System.out.println("Frequensi background nonepi");
      ///  this.printArray2D(bfnonepi);
      //  System.out.println("frekuensi background:");
        
        this.calculateLogOdd(crfepi, bfepi, crfnonepi, bfnonepi,lo,loce,lononce);
       
        
    }
    public double[]LogOddsAnalysisWithSeqAsBg(){
        ExtractBenchmark b = new ExtractBenchmark();
        ArrayList<Benchmark> listb= b.extractBenchmarkFromFile();
        ArrayList <String> merEpi = new ArrayList();
        ArrayList<String> merexposed = new ArrayList();
        b.createMerEpiNonEpi(listb, merEpi, merexposed);
        this.countBFFromBenchmark(listb);
        
        int[]crfepi= new int[20];
        
        this.countFrequencyCR(merEpi, crfepi);
       // System.out.println("Frequensi epi");
      //  this.printArray2D(crfepi);
        int[]crfnonepi= new int[20];
        this.countFrequencyCR(merexposed, crfnonepi);
       // System.out.println("Frequensi nonepi");
       // this.printArray2D(crfnonepi);
      //  System.out.println("frekuensi background:");
      //  this.printArray2D(this.backgroundFrequency);
        double[] lo= this.calculateLogOdd(crfepi, crfnonepi, this.backgroundFrequency);
        return lo;
        
    }
    public void countBFFromBenchmark(ArrayList<Benchmark> blist){
        this.backgroundFrequency = new int[20];
        for(int i=0;i<20;i++){
            this.backgroundFrequency[i]=0;
        }
        for(Benchmark b:blist){
            this.addBackgroundFreq(b.getSequence());
        }
    }
    /*
    public double[] LogOddsAnalysis(){
    //bentuk mer epi dan non epi
        ExtractBenchmark b = new ExtractBenchmark();
        ArrayList<Benchmark> listb= b.extractBenchmarkFromFile();
        ArrayList <String> merEpi = new ArrayList();
        ArrayList<String> merexposed = new ArrayList();
        b.createMerEpiNonEpi(listb, merEpi, merexposed);
      //  System.out.println("merepi:");
      //  printString(merEpi);
      //  System.out.println("mernonepi:");
       // printString(merexposed);
    //hitung frekuensi
        int[] bfepi = new int[20];
        int []rfepi = new int[20];
        countFrequency(merEpi,bfepi,rfepi);
        System.out.println("frekuensi background pada list epi:");
        this.printArray2D(bfepi);
        System.out.println("frekuensi aa pada list epi:");
        this.printArray2D(rfepi);
        //hitung total aa pada background dan pada posisi central
        int aainbfepi = totalCount(bfepi);
        int aainrfepi= totalCount(rfepi);
     //   System.out.println("total bf:"+ aainbfepi);
     //   System.out.println("total rf:"+ aainrfepi);
        // update frekuensi pada bf dan pada rf
    // fungsi update urf = rf/aainrf; ubf = bf/aainbf
        double[] ubfepi= new double[20];
        updateFrequency(bfepi, aainbfepi,ubfepi);
    //    System.out.println("update background frequency:");
      //  this.printArray2D(ubfepi);
        double[] urfepi= new double[20];
        updateFrequency(rfepi, aainrfepi,urfepi);
     //   System.out.println("update frequency aa:");
    //    this.printArray2D(urfepi);
        //hitung log odd
       
        
        int[] bfnonepi = new int[20];
        int []rfnonepi = new int[20];
        countFrequency(merexposed,bfnonepi,rfnonepi);
        int aainbfnonepi = totalCount(bfnonepi);
        int aainrfnonepi= totalCount(rfnonepi);
        double[] ubfnonepi= new double[20];
        updateFrequency(bfnonepi, aainbfnonepi,ubfnonepi);
        double[] urfnonepi= new double[20];
        updateFrequency(rfnonepi, aainrfnonepi,urfnonepi);
        double[] lo= this.calculateLogOdd(urfepi, ubfepi,urfnonepi, ubfnonepi);
        return lo;
        //hitung log odd epi aa sebagai log odd epi - log odd non epi
        
        
    }
    */
   public void countSequenceBasedBackgroundFrequency(int []bf){
        IO io = new IO();
        ArrayList<Epitope> allepilist = io.getEpitopeListFromFile();
        this.backgroundFrequency= new int[20];
        
        for(int i=0;i<20;i++){
            this.backgroundFrequency[i]=0;
            
        }
        
        for (Epitope e: allepilist){
           String sequence = e.getSequence();
              //System.out.println(sequence.values().toString());
            this.addBackgroundFreq(sequence);
           
            
        }
   }
   public void countFrequencyCR(ArrayList<String> mer, int[]crf){
       
        for (int j=0;j<20;j++){
            crf[j]=1;
        }
       for(String m: mer){
           int idxcenter = letterToNum(String.valueOf(m.charAt(4)));
           crf[idxcenter]+=1;
       }
       
   }
    public void countFrequency(ArrayList<String> mer, int[]bf, int[] rf){
        
        
        //hitung background frekuensi setiap asam amino pada mer epi 
        //hitung frekuensi aa pd residu pusat
      //  int[] bf = new int[20];
      //  int[] rf = new int[20];
        //pseudocunt setiap aa ditambah 1 dengan beta =20;
        for (int j=0;j<20;j++){
            bf[j]=1;
            rf[j]=1;
        }
        
        for(String m: mer){
            for(int i=0;i<m.length();i++){
                int idx = letterToNum(String.valueOf(m.charAt(i)));
                bf[idx]= bf[idx]+1;
            }
            int idxcenter = letterToNum(String.valueOf(m.charAt(4)));
            bf[idxcenter]= bf[idxcenter]-1;
            rf[idxcenter]= rf[idxcenter]+1;
        }
          
    }
     public int letterToNum(String aa){
        
        int pos;
        switch (aa){
            case "A":
                    pos=0;
                    break;
                case "R":
                    pos=1;
                    break;
                case "N":
                    pos=2;
                    break;
                case "D":
                   pos=3;
                    break;
                case "C":
                    pos=4;
                    break;
                case "Q":
                    pos=5;
                    break;
                case "E":
                    pos=6;
                    break;
                case "G":
                    pos=7;
                    break;
                case "H":
                    pos=8;
                    break;
                case "I":
                    pos=9;
                    break;
                case "L":
                    pos=10;
                    break;
                case "K":
                    pos=11;
                    break;
                case "M":
                    pos=12;
                    break;
                case "F":
                    pos=13;
                    break;
                case "P":
                    pos=14;
                    break;
                case "S":
                    pos=15;
                    break;
                case "T":
                    pos=16;
                    break;
                case "W":
                    pos=17;
                    break;
                case "Y":
                    pos=18;
                    break;
                case "V":
                    pos=19;
                    break;
                default:
                    pos=-1;
        }
        return pos;
    
    }
     public void printLogOddsOfAA(double[]lo){
         for(int i=0;i<20;i++){
             System.out.println(numToletter(i)+": "+ lo[i]);
         }
         
     }
     //bobot diambil dari nilai log odds yang dinormalkan menjadi bernilai antara 0 sd 1
     //rumus new value = (val - min)/(max - min)
     public double[] LogOddsAsWeight(double[]lo){
        double[] newlo= new double[20];
         double min = Double.MAX_VALUE;
         double max = Double.MIN_VALUE;
         for(int i=0;i<lo.length;i++){
             if(min > lo[i]){
                 min = lo[i];
             }
             if(max < lo[i]){
                 max = lo[i];
             }
         }
         for(int j=0;j<lo.length;j++){
             newlo[j]= (lo[j]-min)/(max - min);
         }
         return newlo;
         

     }
     public String numToletter(int aa){
        
        String letter="";
        switch (aa){
            case 0:
                    letter="A";
                    break;
                case 1:
                    letter="R";
                    break;
                case 2:
                    letter="N";
                    break;
                case 3:
                   letter="D";
                    break;
                case 4:
                    letter="C";
                    break;
                case 5:
                    letter="Q";
                    break;
                case 6:
                    letter="E";
                    break;
                case 7:
                    letter="G";
                    break;
                case 8:
                    letter="H";
                    break;
                case 9:
                    letter="I";
                    break;
                case 10:
                    letter="L";
                    break;
                case 11:
                    letter="K";
                    break;
                case 12:
                    letter="M";
                    break;
                case 13:
                    letter="F";
                    break;
                case 14:
                    letter="P";
                    break;
                case 15:
                    letter="S";
                    break;
                case 16:
                    letter="T";
                    break;
                case 17:
                    letter="W";
                    break;
                case 18:
                    letter="Y";
                    break;
                case 19:
                    letter="V";
                    break;
                
        }
        return letter;
    
    }
    public void doLogOdds(){
        ExtractBenchmark b = new ExtractBenchmark();
        ArrayList<Benchmark> listb= b.extractBenchmarkFromFile();
        ArrayList <String> merEpi = new ArrayList();
        ArrayList<String> merexposed = new ArrayList();
        int[]bf = new int[20];
        //membuat overlap mer epitope dan exposed residu non epitope
        //menghitung distribusi asam amino atau background frequency
        b.createOverlap(listb, merEpi, merexposed, bf);
        //sequence weighting dengan menghitung frekuensi residu setiap posisi
        double[][] frepi=this.calculateResidueFrequency(merEpi);
        
        int[]drepi = this.getDistinctInPos(frepi);
        IO ioobjek = new IO();
        double [][] blosum= ioobjek.getBlosumMatrixFromFile();
        int drepiInCenter = drepi[4];
        double []psa=new double[20];
        if(drepiInCenter == 20){
            //hitung pseudocount untuk semua aa
            for(int i=0;i<20;i++){
                psa[i]=this.calculatePseudoCount(i, 5, frepi, blosum);
            }
        }
        double[] rf = this.getRFAtPosisionX(4, frepi);
        double alpha = this.calculateAlpha(rf, drepiInCenter);
        double beta = 200;
        double []PA= this.calculatePA(alpha, beta,rf , psa);
         System.out.println("PA:");
        this.printArraySatuD(PA);
        //hitung logodd (PA dan backgroun)
        double[] lo= this.calculateLogOdd(PA, bf);
        
        System.out.println("Log Odds");
        this.printArraySatuD(lo);
        
        
        
        
        
    }
    public void doLogOddsAnalysisInBlosum(){
        IO io = new IO();
        ArrayList<Epitope> allepilist = io.getEpitopeListFromFile();
        ArrayList<String> merEpi= new ArrayList();
        ArrayList<String> merNonEpi= new ArrayList();
        this.backgroundFrequency= new int[20];
        
        for(int i=0;i<20;i++){
            this.backgroundFrequency[i]=0;
            
        }
        
        for (int i=0; i<allepilist.size();i++){
           
            //set sequence, idxexposedresidue and epilist of chain
            System.out.println(allepilist.get(i).getPdbid() +"_" + allepilist.get(i).getChain());
            this.setSequenceAndEpilistOfPdbidChain(allepilist.get(i).getPdbid() , allepilist.get(i).getChain(),1);
            //hitung frekuensi background
           
            //this.addBackgroundFrequency(this.sequence,idxepi);
            System.out.println(sequence.values().toString());
            this.addBackgroundFreq(sequence.values().toString());
           this.setOverlapMer();
           this.setOverlapmerEpiandNonEpi(this.epi);
                 
            merEpi.addAll(this.getOverlapmerAsEpi());
            merNonEpi.addAll(this.getOverlapmerAsNonEpi());
            
        }
      //  this.printArraySatuD(backgroundFrequency);
        //this.writeLogOddAnalysisResult(merEpi, "daftarOverlapMerEpi.txt");
        //this.writeLogOddAnalysisResult(merNonEpi, "daftarOverlapMerNonEpi.txt");
        double[][] rfepi=this.calculateResidueFrequency(merEpi);
        System.out.println("rf epitope:");
        this.printResidueFrequency(rfepi);
        
        int[]drepi = this.getDistinctInPos(rfepi);
        this.printDistinctPos(drepi);
        //position weight tanpa pseudo-count frequency
        //this.getPositionWeigth(rfepi, drepi);
        //perbaikan 
        
        IO ioobjek = new IO();
         double [][] blosum= ioobjek.getBlosumMatrixFromFile();
         // this.printBLOSUM(blosum);
         //hitung pseudocount
        int drepiInCenter = drepi[4];
        double []psa=new double[20];
        if(drepiInCenter == 20){
            //hitung pseudocount untuk semua aa
            for(int i=0;i<20;i++){
                psa[i]=this.calculatePseudoCount(i, 5, rfepi, blosum);
            }
        }
        double[] rf = this.getRFAtPosisionX(4, rfepi);
        double alpha = this.calculateAlpha(rf, drepiInCenter);
        /*
        //definisi alpha 1
        //hitung alpha = rerata jumlah asam amino berbeda pada posisi alignment
        double alpha=0;
        for(int i=0;i<drepi.length;i++)
            alpha = alpha+drepi[i];
        alpha = alpha/20;
        */
        double beta = 200;
        //hitung frekuensi efektif aa
        
         double []PA= this.calculatePA(alpha, beta,rf , psa);
         System.out.println("PA:");
        this.printArraySatuD(PA);
        //hitung logodd (PA dan backgroun)
        double[] lo= this.calculateLogOdd(PA, this.backgroundFrequency);
        
        System.out.println("Log Odds");
        this.printArraySatuD(lo);
        double [][] rfnonepi= this.calculateResidueFrequency(merNonEpi);
         //System.out.println("rf non epitope:");
        //this.printResidueFrequency(rfnonepi);
        int[]drnonepi = this.getDistinctInPos(rfnonepi);
        //this.printDistinctPos(drnonepi);
        //this.getPositionWeigth(rfnonepi, drnonepi);
        this.printResidueFrequency(rfnonepi);
        //hitung sequence weight
        /*
        for(int i=0;i<this.overlapmer.size();i++){
            this.calculateSequenceWeight(rfnonepi, sequence);
        }*/
        
    }
    public void acumulateBF(int[] distAA){
        for(int i=0;i<20;i++){
            this.backgroundFrequency[i]+=distAA[i];
        }
    }
    public double calculateAlpha(double[]rfpi,int dr){
        double alpha=0;
        for(int i=0;i<rfpi.length;i++){
            alpha += rfpi[i];
        }
        return alpha/dr;
    }
    public void calculateLogOdd(int[] fepi, int[] bfepi, int[]rfnonepi, int[]bfnonepi, double[]loce,double[]lononce,double[]lo){
        double sumbf = 0;
        double sumepif=0;
        double sumbfnonepi=0;
        double sumnonepif=0;
        for(int i=0;i<bfepi.length;i++)
            sumbf += (double)bfepi[i];
     //   System.out.println("sumbf: "+sumbf);
        for(int k=0;k<bfnonepi.length;k++)
            sumbfnonepi += (double)bfnonepi[k];
      //  System.out.println("sumbfnonepi: "+sumbfnonepi);
        for(int j=0;j<fepi.length;j++)
            sumepif += (double)fepi[j];
     //   System.out.println("sumfepi: "+sumepif);
         for(int j=0;j<rfnonepi.length;j++)
            sumnonepif += (double)rfnonepi[j];
      //   System.out.println("sumfnonepi: "+sumnonepif);
        double[] cfepi = new double[20];
        for(int m=0;m<fepi.length;m++)
            cfepi[m]=(double)fepi[m]/sumepif;
      //  System.out.println("probabilitas epi: ");
      //   this.printArray2D(cfepi);
        double[] cbfepi = new double[20];
        for(int n=0;n<bfepi.length;n++)
            cbfepi[n]=(double)bfepi[n]/sumbf;
     //   System.out.println("probabilitas bfepi: ");
     //    this.printArray2D(cbfepi);
        double[] crfnonepi = new double[20];
        for(int p=0;p<rfnonepi.length;p++)
            crfnonepi[p]=(double)rfnonepi[p]/sumnonepif;
    //    System.out.println("probabilitas rfnonepi: ");
    //     this.printArray2D(crfnonepi);
       double[] cbfnonepi = new double[20];
        for(int q=0;q<cbfnonepi.length;q++)
            cbfnonepi[q]=(double)bfnonepi[q]/sumbfnonepi;
      //   System.out.println("probabilitas bfnonepi: ");
       //  this.printArray2D(cbfnonepi);
        for(int r=0;r<cfepi.length;r++)
            cfepi[r]=cfepi[r]/cbfepi[r];
    //    System.out.println("odd epi: ");
      //   this.printArray2D(cfepi);
        for(int s=0;s<crfnonepi.length;s++)
            crfnonepi[s]=crfnonepi[s]/cbfnonepi[s];
     //   System.out.println("odd nonepi: ");
      //   this.printArray2D(crfnonepi);
        
       // double log2 = Math.log(2);
        for(int i=0;i<lo.length;i++){
            loce[i]= Math.log(cfepi[i])/Math.log(2);
            lononce[i]=Math.log(crfnonepi[i])/Math.log(2);
            lo[i]= loce[i]-lononce[i];
           
        }

    }
    public double[]calculateLogOdd(int[] crfepi, int[] crfnonepi, int[]bf){
        double[] lo = new double[20];
        double log2 = Math.log(2);
        for(int i=0;i<lo.length;i++){
            lo[i]= (Math.log((double)crfepi[i]/(double)bf[i])/log2)- (Math.log((double)crfnonepi[i]/(double)bf[i])/log2);
        }
        return lo;
    }
    public double[]calculateLogOdd(double[] pa, int[] qa){
        double[] lo = new double[20];
        double log2 = Math.log(2);
        for(int i=0;i<pa.length;i++){
            lo[i]= Math.log(pa[i]/qa[i])/log2;
        }
        return lo;
    }
   public void addBackgroundFreq(String newsequence){
       for(int i=0;i<newsequence.length();i++){
           char c= newsequence.charAt(i);
           int pos = this.getResiduePosition(String.valueOf(c));
           this.backgroundFrequency[pos]+=1;
           
       }
   }
    public void addBackgroundFrequency(String newsequence, ArrayList<Integer> idxepi){
        for(int i=0;i<newsequence.length();i++){
            if(!idxepi.contains(i)){
                char c= newsequence.charAt(i);
                int pos = this.getResiduePosition(String.valueOf(c));
                this.backgroundFrequency[pos]+=1;
            }
            
        }
    }
    public void printArraySatuD(int[] arr){
        for(int i=0;i<arr.length;i++){
            System.out.print(arr[i]+"\t");
        }
        System.out.print("\n");
    }
    public void printArraySatuD(double[]arr){
        for(int i=0;i<arr.length;i++){
            System.out.print(arr[i]+"\t");
        }
        System.out.print("\n");
    }
    public double[]getRFAtPosisionX(int x, double[][]rfmer){
    double []rf = new double[20];
    for(int i=0;i<20;i++){
        rf[i]=rfmer[i][5];
    }
    return rf;
    }
    //frekuensi efektif aa
    public double[] calculatePA(double alpha, double beta,double []fa, double[]ga){
        double []PA= new double[20];
        for(int i=0;i<fa.length;i++){
            PA[i]=(fa[i]*alpha+ga[i]*beta)/(alpha+beta);
        }
        return PA;
    }
    //print blosum
    public void printBLOSUM(double[][]blos){
     System.out.print("dimensi: "+ blos.length +"x"+blos[0].length+"\n");
    for (int i=0;i<blos.length;i++){
        for(int j=0;j<blos[0].length;j++){
            System.out.print(blos[i][j]+"\t");}
        System.out.print("\n");
    }
    }
    //Pca=nca + bca/Nc+Bc
    //Pca=probabilitas residu pada kolom c
    //nca= jumlah a pada kolom c
    //bca=pseudocount a pada kolom c
    //Nc=total count pada kolom c
    //total pseudocount pada kolom c
    
    //pseudofrekuensi aminoacid b
    //sum of fa.qb|a
    //fa: kemunculan asam amino a
    //qb|a: kemungkinan b dalam ikatan pada posisi t jika a dalam ikatan pada posisi t.
    //dari matrik blosum qb|a: nilai pada kolom b dari baris a
    public double calculatePseudoCount(String aa, int posInAlignment,double[][]rf,double[][]blosum){
        //posisi pada blosum matrik
        int pos= this.getResiduePosition(aa);
        double pseudocount=0;
        System.out.println("rf length: "+ rf.length);
        
        for(int i=0;i<rf.length;i++){
            if(i!=pos){
                pseudocount = pseudocount + rf[i][posInAlignment-1]*blosum[i][pos];
                System.out.print("rf*b: "+rf[i][posInAlignment-1]+"\t"+blosum[i][pos]+"\n");
            }
            
        }
        return pseudocount;
    }
    public double calculatePseudoCount(String aa, int posInAlignment,double[]rf,double[][]blosum){
        //posisi pada blosum matrik
        int pos= this.getResiduePosition(aa);
        double pseudocount=0;
        System.out.println("rf length: "+ rf.length);
        
         for(int i=0;i<rf.length;i++){
            if(i!=pos){
                pseudocount = pseudocount + rf[i]*blosum[i][pos];
                System.out.print("rf*b: "+rf[i]+"\t"+blosum[i][pos]+"\n");
            }
         }
        return pseudocount;
    }
    public double calculatePseudoCount(int pos, int posInAlignment,double[][]rf,double[][]blosum){
        //posisi pada blosum matrik
        
        double pseudocount=0;
        System.out.println("rf length: "+ rf.length);
        
        for(int i=0;i<rf.length;i++){
            if(i!=pos){
                pseudocount = pseudocount + rf[i][posInAlignment-1]*blosum[i][pos];
                System.out.print("rf*b: "+rf[i][posInAlignment-1]+"\t"+blosum[i][pos]+"\n");
            }
            
        }
        return pseudocount;
    }
    public double calculateSequenceWeight(double [][]positionweight,String mer){
        double weight=0;
        for(int i=0;i<mer.length();i++){
            int posAA = this.getResiduePosition(mer.substring(i, i+1));
            weight=weight + positionweight[posAA][i];
        }
        return weight;
    }
    public void getPositionWeigth(double[][] rf, int[]dr){
        for (int i=0;i<rf.length;i++){
            for(int j=0;j<dr.length;j++){
                rf[i][j]= 1/(rf[i][j]*dr[j]);
            }
        }
    }
    public void printDistinctPos(int[] dp){
        for(int i=0;i<dp.length;i++)
            System.out.print(dp[i]+"\t");
        System.out.println();
    }
    public int[]getDistinctInPos(double[][] pbrw){
        int [] distinctResidu= new int[9];
        for(int k=0;k<distinctResidu.length;k++)
            distinctResidu[k]=0;
        if(pbrw.length==20 && pbrw[0].length==9){
            
             for(int j=0;j<pbrw[0].length;j++)
                for(int i=0;i<pbrw.length;i++){
                    if(pbrw[i][j]!=0)
                        distinctResidu[j]=distinctResidu[j]+1;
                }
        }
        return distinctResidu;
            
    }
    public void printResidueFrequency(double[][] rf){
        for (int i=0;i<rf.length;i++){
                for(int j=0;j<rf[0].length;j++){
                System.out.print(rf[i][j]+"\t");
            }
                System.out.println();
        }
        
            
    }
    public int getResiduePosition(String s){
        int pos;
        switch (s){
            case "A":
                    pos=0;
                    break;
                case "R":
                    pos=1;
                    break;
                case "N":
                    pos=2;
                    break;
                case "D":
                   pos=3;
                    break;
                case "C":
                    pos=4;
                    break;
                case "Q":
                    pos=5;
                    break;
                case "E":
                    pos=6;
                    break;
                case "G":
                    pos=7;
                    break;
                case "H":
                    pos=8;
                    break;
                case "I":
                    pos=9;
                    break;
                case "L":
                    pos=10;
                    break;
                case "K":
                    pos=11;
                    break;
                case "M":
                    pos=12;
                    break;
                case "F":
                    pos=13;
                    break;
                case "P":
                    pos=14;
                    break;
                case "S":
                    pos=15;
                    break;
                case "T":
                    pos=16;
                    break;
                case "W":
                    pos=17;
                    break;
                case "Y":
                    pos=18;
                    break;
                case "V":
                    pos=19;
                    break;
                default:
                    pos=-1;
        }
        return pos;
    }
    public double[][] calculateResidueFrequency(ArrayList<String> mer){
       
        double [][] pbrw = new double[20][9];
        for(int i=0;i<20;i++)
            for(int j=0;j<9;j++)
                pbrw[i][j]=0;
        for(int k=0;k<9;k++)
            for(int el=0;el<mer.size();el++){
            switch(mer.get(el).substring(k,k+1)){
                case "A":
                    pbrw[0][k]=pbrw[0][k]+1;
                    break;
                case "R":
                    pbrw[1][k]=pbrw[1][k]+1;
                    break;
                case "N":
                    pbrw[2][k]=pbrw[2][k]+1;
                    break;
                case "D":
                    pbrw[3][k]=pbrw[3][k]+1;
                    break;
                case "C":
                    pbrw[4][k]=pbrw[4][k]+1;
                    break;
                case "Q":
                    pbrw[5][k]=pbrw[5][k]+1;
                    break;
                case "E":
                    pbrw[6][k]=pbrw[6][k]+1;
                    break;
                case "G":
                    pbrw[7][k]=pbrw[7][k]+1;
                    break;
                case "H":
                    pbrw[8][k]=pbrw[8][k]+1;
                    break;
                case "I":
                    pbrw[9][k]=pbrw[9][k]+1;
                    break;
                case "L":
                    pbrw[10][k]=pbrw[10][k]+1;
                    break;
                case "K":
                    pbrw[11][k]=pbrw[11][k]+1;
                    break;
                case "M":
                    pbrw[12][k]=pbrw[12][k]+1;
                    break;
                case "F":
                    pbrw[13][k]=pbrw[13][k]+1;
                    break;
                case "P":
                    pbrw[14][k]=pbrw[14][k]+1;
                    break;
                case "S":
                    pbrw[15][k]=pbrw[15][k]+1;
                    break;
                case "T":
                    pbrw[16][k]=pbrw[16][k]+1;
                    break;
                case "W":
                    pbrw[17][k]=pbrw[17][k]+1;
                    break;
                case "Y":
                    pbrw[18][k]=pbrw[18][k]+1;
                    break;
                case "V":
                    pbrw[19][k]=pbrw[19][k]+1;
                    break;
                }
            }
        return pbrw;
        
    }
    public String sequenceHashValueToString(String s){
        int end = s.length();
       String[] as= s.substring(1, end-1).split(", ");
       String finals="";
       for(int i=0;i<as.length;i++){
           //if(as[i].)    
           finals += as[i];
               
       } 
       return finals;
    }
    public LogOdds(){
      // this.lo = this.LogOddsAnalysis();
    }
    public LogOdds(int mode){
       // this.lo = this.LogOddsAnalysisWithSeqAsBg();
    }
   public LogOdds(String pdbid, String antigenchain){
       this.sequence= new HashMap<String,String>();
        this.idxSequence = new HashMap<Integer,String>();
       this.setSequenceAndEpilistOfPdbidChain(pdbid, antigenchain,1);
       //String s= sequenceHashValueToString(this.sequence.values().toString());
       this.setOverlapMer();
       this.setOverlapmerEpiandNonEpi(this.epi);
   }
    //menandai epitope atau bukan dari setiap pdbid
    //dari epitope list diperoleh id dari asam amino sebagai epitope
    //jika pusat overlapmer adalah epitope maka masukkan kedalam daftar
    
    public void setSequenceAndEpilistOfPdbidChain(String pdbid, String antigenChain, int mode){
        //dari pdbid dan pdbchain
        //bentuk kompleksnya
        Complex complex = new Complex(pdbid, antigenChain);
        //ambil daftar epitope dari sequence
//        acumulateBF(complex.getChainByChianID(antigenChain.charAt(0)).getDistributionOfAA());
        complex.getChainByChianID(antigenChain.charAt(0)).getAASequenceOfChain(this.sequence, this.idxSequence);
        this.epi = this.getEpilistOfPdbid(pdbid, antigenChain);
        
       //set index of exposed residu
        this.idxExposedResidue= this.getIdxExposedAA(complex, antigenChain.charAt(0), mode);
            
        
    }
    public ArrayList<String> getIdxExposedAA(Complex c, char antigenChain, int mode){
         UnboundDSSP udssp = new UnboundDSSP(c,String.valueOf(antigenChain));
        ArrayList<Aminoacid> exposedaa= udssp.getExposedAA(c,antigenChain,mode );
        ArrayList<String> idxEAA = new ArrayList();
        if(exposedaa.size()>0){
            for (int i=0;i<exposedaa.size();i++){
                //System.out.println(exposedaa.get(i).getresidueID());
                try{
                    idxEAA.add(exposedaa.get(i).getresidueID());
                }catch(Exception e){
                    System.out.println(e.getMessage());
                }
                
            }
        }
        return idxEAA;
        
        
    }
    
    //buat hashmap sequence dan index
    //ambil mer yang pusatnya dengan RSA lebih dari 0.05
    public void setOverlapMer(){
        //kalau idx sama dengan 
        this.overlapmer = new ArrayList();
        int seqLength = this.sequence.size();
        int sumofexposedRes = this.idxExposedResidue.size();
        //ambil idxexposed residu
        int count=0;
        String seqInString = this.sequenceHashValueToString(this.sequence.values().toString());
        System.out.println(seqInString);
        for(int i=0;i<seqInString.length()-9;i++){
           String key = this.idxSequence.get(i);
           if(count<sumofexposedRes){
               if(key.equalsIgnoreCase(this.idxExposedResidue.get(count))){
                   System.out.println(seqInString.substring(i, i+9));
                this.overlapmer.add(seqInString.substring(i, i+9));
                count = count+1;
            }
           }
            
            
        }
    }
    public void setOverlapMer(String primarysequence){
        this.overlapmer = new ArrayList();
        int seqLength = primarysequence.length();
        for(int i=0;i<seqLength-9;i++){
            this.overlapmer.add(primarysequence.substring(i, i+9));
        }
    }
    public void setMerEpiandNonEpi(Epitope e){
        this.overlapmerAsEpi= new ArrayList();
        this.overlapmerAsNonEpi=new ArrayList();
        //ambil index sequence sebagai epitope
        ArrayList<Integer> idxEpi =e.getIdxsequenceAsEpitope();
        if(this.overlapmer.size()>0){
            int count =0;
            int omcount=0;
            while(count<idxEpi.size()&&omcount < this.overlapmer.size()){
                
                if (omcount+4 ==idxEpi.get(count).intValue()){
                this.overlapmerAsEpi.add(this.overlapmer.get(omcount));
              // this.overlapmerAsNonEpi.remove(omcount);
                int idx = this.overlapmerAsNonEpi.indexOf(this.overlapmer.get(omcount));
                this.overlapmerAsNonEpi.remove(idx);
                count = count+1;
            }
            omcount= omcount+1;
        }
        }
    }
    public void setOverlapmerEpiandNonEpi(Epitope e){
        this.overlapmerAsEpi= new ArrayList();
        this.overlapmerAsNonEpi=new ArrayList();
        //ambil index sequence sebagai epitope
        ArrayList<Integer> idxEpi =e.getIdxsequenceAsEpitope();
        if(this.overlapmer.size()>0){
            this.overlapmerAsNonEpi = new ArrayList<> (this.overlapmer);
        }else{
            this.setOverlapMer();
            this.overlapmerAsNonEpi = new ArrayList<> (this.overlapmer);
        }
        int count =0;
        int omcount=0;
        while(count<idxEpi.size()&&omcount < this.overlapmer.size()){
            if (omcount+4 ==idxEpi.get(count).intValue()){
                this.overlapmerAsEpi.add(this.overlapmer.get(omcount));
              // this.overlapmerAsNonEpi.remove(omcount);
                int idx = this.overlapmerAsNonEpi.indexOf(this.overlapmer.get(omcount));
                this.overlapmerAsNonEpi.remove(idx);
                count = count+1;
            }
            omcount= omcount+1;
        }
        
    }

    public Epitope getEpi() {
        return epi;
    }

    public ArrayList<String> getOverlapmerAsEpi() {
        return overlapmerAsEpi;
    }

    public ArrayList<String> getOverlapmerAsNonEpi() {
        return overlapmerAsNonEpi;
    }
    
    public Epitope getEpilistOfPdbid(String pdbid, String pdbchain){
        IO io = new IO();
        Epitope epi = new Epitope();
        ArrayList<Epitope> epilist = io.getEpitopeListFromFile();
        
        for(int i=0;i<epilist.size();i++){
            if(epilist.get(i).getPdbid().equalsIgnoreCase(pdbid)){
                epi =epilist.get(i);
                break;
                
            }
            
        }
        return epi;
    }

    public ArrayList<String> getOverlapmer() {
        return overlapmer;
    }
    //menuliskan daftar epi dan non epi
    public void writeLogOddAnalysisResult(ArrayList<String> mer, String filename){
       
        BufferedWriter bw = null;
	FileWriter fw = null;
	try {
            fw = new FileWriter(new File("src/datatraining/" + filename));
            bw = new BufferedWriter(fw);
            //tulis disini 
            for(int i=0;i<mer.size();i++){
                String content = mer.get(i)+ "\n";
                bw.write(content);
            }
  //         
		
	} catch (IOException e) {
		e.printStackTrace();
	} finally {
		try {
                    if (bw != null)
			bw.close();

                    if (fw != null)
			fw.close();
		} catch (IOException ex) {
			ex.printStackTrace();

		}

	}
    }
   
    
}
