/*
 *extract benchmark from file
 * buat fungsi pembentukan overlapmer epi dan non epi
 */
package ExtractDataTraining;

import ExtractData.CreateComplexStructure;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author OpenGress
 */
public class ExtractBenchmark {
    /*
     baca string sequence
     */
    //menghitung bf pada mer epi dan non epi
    //parameter mer
    public int[]accumulateBF(ArrayList<String> s){
        int []idx = new int[20];
        for(int i=0;i<20;i++){
            idx[i]=0;
        }
        for(String sub : s){
            for(int j=0;j<sub.length();j++){
                int num = this.letterToNum(String.valueOf(sub.charAt(j)));
                idx[num]= idx[num]+1;
            }
        }
        return idx;
    }
    public void createMerExposedResidu(Benchmark b, ArrayList<String>merExpAAP,ArrayList<Integer> idxfourthpos){
        
            ArrayList<Integer> idxexposedres = b.getIdxexp();
            
            for(int i=0;i<idxexposedres.size()-1;i++){
                int idx1= idxexposedres.get(i);
                if(idx1<4){
                    i++;
                }else if(idx1<b.getSequence().length()-6){
                    int idx2= idxexposedres.get(i+1);
                    String mer;
                    if(idx2==idx1+1){
                    //ambil mer sebagai exposed residu
                        mer = b.getSequence().substring(idx1-4, idx1+6);
                        merExpAAP.add(mer);
                        idxfourthpos.add(idx1);
                        
                    }
                }
           
            }
        
    }
    public void getMerEpiNonEpiFromExposedMer(ArrayList<String> expmerAAP, ArrayList<Integer> pos,ArrayList<Integer> idxepi, ArrayList<String> epimerAAP,ArrayList<String> nonepimerAAP ){
        
        for (String s:expmerAAP){
            nonepimerAAP.add(s);
        }
       // System.out.println("pos epi: ");
        for(int i=0;i<idxepi.size()-1;i++){
            Integer posEpi1 = idxepi.get(i);
            Integer nextEpi = posEpi1+1;
            Integer posEpi2 = idxepi.get(i+1);
            if(posEpi2.equals(nextEpi)){
               // posEpi.add(idxepi.indexOf(posEpi1));
                epimerAAP.add(expmerAAP.get(idxepi.indexOf(posEpi1)));
                nonepimerAAP.remove(idxepi.indexOf(posEpi1));
            }
        }
        
    }
    
    public void createMerEpiNonEpiAAP(ArrayList<Benchmark> blist, ArrayList<String>merEpiAAP, ArrayList<String> mernonEpiAAP){
        for(Benchmark b : blist){
            try{
                ArrayList<Integer> idxfourthpos= new ArrayList();
                ArrayList<String>merExpAAP = new ArrayList();
                createMerExposedResidu( b, merExpAAP,idxfourthpos);
                ArrayList<Integer> idxepi = b.getIdxEpi();
                ArrayList<String> epimerAAP= new ArrayList();
                ArrayList<String> nonepimerAAP = new ArrayList();
                getMerEpiNonEpiFromExposedMer(merExpAAP, idxfourthpos,idxepi, epimerAAP,nonepimerAAP );
                merEpiAAP.addAll(epimerAAP);
                mernonEpiAAP.addAll(nonepimerAAP);
            }catch(Exception e){
                System.out.println(e.getMessage());
            }
            
            
        }
    }
    public void calculateFrequencyAAP(ArrayList <String>mer, int [][]pos4pos5, int [][]background){
        
        for(int i=0;i<20;i++)
            for(int j=0;j<20;j++){
                pos4pos5[i][j]=0;
                background[i][j]=0;
            }
        //harus dihitung sebagai kemunculan kejadian pasangan aa secara bersama-sama
        int lengthMer = mer.get(0).length();
        for(String s:mer){
            for(int k=0;k<lengthMer -1;k++){
            int pos1 = letterToNum(String.valueOf(s.charAt(k)));
            int pos2 = letterToNum(String.valueOf(s.charAt(k+1)));
            if(k==4){
                pos4pos5[pos1][pos2]+=1;
            }else{
                background[pos1][pos2]+=1;
            }
            
        }
        }
        
    }
    public void calculateBFAAP(ArrayList <String>mer, int[]bg, int[]pos4, int[]pos5){
        for(String s: mer){
            for(int i=0;i<s.length();i++){
                int idx = letterToNum(String.valueOf(s.charAt(i)));
            if(i==4){
                pos4[idx]= pos4[idx]+1;
            }else if(i==5){
                pos5[idx]= pos5[idx]+1;
            }else{
                bg[idx]=bg[idx]+1;
            }
            }
        }
    }
    public void updateWithLaplaceCount(int[][]pos4pos5){
        for(int i=0;i<20;i++)
        for (int j=0;j<20;j++){
            pos4pos5[i][j] += 1;
            
        }
    }
   public void calculateLogOddRatio(double[][]epi, double[][]nonepi, double[][]ratio){
       for(int i=0;i<20;i++){
       for (int j=0;j<20;j++){
           ratio[i][j]=epi[i][j]-nonepi[i][j];
       }
       }
    }
    public void calculateLogOdd(int[][]bg, int[][]pos4pos5, double[][]logodd){
     //membagi frekuensi motif dengan frekuensi background
       
            for(int i=0;i<20;i++)
            for(int j=0;j<20;j++){
                logodd[i][j]=Math.log((double) pos4pos5[i][j]/(double)bg[i][j]);
               
            }
        
    }
    public void createMerEpiNonEpi(ArrayList<Benchmark> blist, ArrayList<String>merEpi, ArrayList<String> expmer){
        for(Benchmark b : blist){
            try{
                String s= b.getSequence();
                ArrayList<Integer> idxepitope = b.getIdxEpi();
                ArrayList<Integer> idxexposedres = b.getIdxexp();
                //System.out.println("epi: "+ idxepitope.size()+"; exp: "+ idxexposedres.size());
                int countepi =0;
                int countexp =4;
                for(int i=0;i<s.length()-9;i++){
                    String mer = s.substring(i,i+9);
                    String mer2 = s.substring(i,i+9);
                    int indexepi = idxepitope.get(countepi);
                    int indexexp = idxexposedres.get(countexp);
                    //System.out.println("idx seq:"+ (i+4) + "idx:"+idx.get(count));
                    if(indexepi ==i+4){
                        merEpi.add(mer);
                        countepi= countepi+1;
                    }
                    if(indexexp ==i+4){
                        expmer.add(mer2);
                        countexp= countexp+1;
                    }
                
                }
            }catch(Exception e){
               // System.out.println(e.getMessage());
            }
            
            
        }
    
    }
   
    public void createOverlap(ArrayList<Benchmark> blist, ArrayList<String>merEpi, ArrayList<String> expmer, int[]bf){
       
        for(int j=0;j<20;j++){
            bf[j]=0;
        }
        
        for(Benchmark b : blist){
            try{
                String s= b.getSequence();
                //hitung backgroundfrequency
                this.calculateBF(s, bf);                
            ArrayList<Integer> idxepitope = b.getIdxEpi();
            ArrayList<Integer> idxexposedres = b.getIdxexp();
          //  System.out.println("epi: "+ idxepitope.size()+"; exp: "+ idxexposedres.size());
            int countepi =0;
            int countexp =4;
            for(int i=0;i<s.length()-9;i++){
                String mer = s.substring(i,i+9);
                String mer2 = s.substring(i,i+9);
                int indexepi = idxepitope.get(countepi);
                int indexexp = idxexposedres.get(countexp);
                //System.out.println("idx seq:"+ (i+4) + "idx:"+idx.get(count));
                if(indexepi ==i+4){
                    merEpi.add(mer);
                    countepi= countepi+1;
                }
                if(indexexp ==i+4){
                    expmer.add(mer2);
                    countexp= countexp+1;
                }
                
            }
            }catch(Exception e){
                System.out.println(e.getMessage());
            }
            
            
        }
    }
    public void calculateBF(String sequence, int[]distAA){
        
        for (int i=0;i<sequence.length();i++){
            int idx = letterToNum(String.valueOf(sequence.charAt(i)));
            distAA[idx]= distAA[idx]+1;
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
    
    /*
        * baca file banchmark
        * ambil pdbid dan chain
        * ambil data sequence sebagai sequence
        * ambil data status epitope
     */
    public ArrayList<Benchmark> extractBenchmarkFromFile(){
        ArrayList<Benchmark> listB= new ArrayList();
        String filename = "seqepiexposed.txt";
      //  String filename = "seqepiexposedRen2015.txt";
        //file seqepiexposed berisi sequence, posisi epi dan posisi exposed
        String dataline = null;
        FileReader fr = null;
    try
    {
      if (new File("src/dataset/" + filename).exists()) {
        fr = new FileReader(new File("src/dataset/" + filename));
        BufferedReader br = new BufferedReader(fr);
        dataline = br.readLine();
       // System.out.println(dataline);
        while (dataline != null)
        {
           // System.out.println(dataline);
             String[] linestring= dataline.split("\t");
            // String[] linestring= dataline.split(" ");
             Benchmark b = new Benchmark();
             //ambil 3 bagian data
            String seq = linestring[0];
            String epi = linestring[1];
            String exposed = linestring[2];
            
            ArrayList<Integer> idxepi = new ArrayList();
           
          
            for(int i=0;i<epi.length();i++){
                String statusepi = String.valueOf(epi.charAt(i));
                if(statusepi.equalsIgnoreCase("1")){
                   idxepi.add(Integer.valueOf(i));
                }                
            }
             ArrayList<Integer>idxepiaap = new ArrayList();
             ArrayList<Integer>idxoneepiaap = new ArrayList();
            for(int i=0;i<epi.length()-1;i++){
                String statusepi1 = String.valueOf(epi.charAt(i));
                String statusepi2 = String.valueOf(epi.charAt(i+1));
                if(statusepi1.equalsIgnoreCase("1")&& statusepi2.equalsIgnoreCase("1")){
                   idxepiaap.add(Integer.valueOf(i));
                } else{
                    idxoneepiaap.add(Integer.valueOf(i));
                }
                    
            }
             ArrayList<Integer> idxexp = new ArrayList();
           for(int i=0;i<exposed.length();i++){
                String statusexp = String.valueOf(exposed.charAt(i));
                if(statusexp.equalsIgnoreCase("1")){
                   idxexp.add(Integer.valueOf(i));
                }                
            }
            ArrayList<Integer> idxexpaap = new ArrayList();
            for(int i=0;i<exposed.length()-1;i++){
                String statusexp1 = String.valueOf(exposed.charAt(i));
                String statusexp2 = String.valueOf(exposed.charAt(i+1));
                if(statusexp1.equalsIgnoreCase("1")&& statusexp2.equalsIgnoreCase("1")){
                   idxexpaap.add(Integer.valueOf(i));
                }                
            }
           b.setSequence(seq);           
           b.setEpipos(epi);
           b.setExposedpos(exposed);
           b.setIdxEpi(idxepi);//index epitope
           b.setIdxexp(idxexp);//index exposed residu
           b.setIdxEpiAAP(idxepiaap);//awal index pasangan epitope
           b.setIdxOneEpiAAP(idxoneepiaap);//awal index pasangan berisi satu epitope
            listB.add(b);
             dataline = br.readLine();
        }
        br.close();
      }
      
    }
    catch (IOException e)
    {
      e.printStackTrace();
     
    }finally{
        return listB;
    }
   
    }
}
