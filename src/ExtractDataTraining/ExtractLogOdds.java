/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ExtractDataTraining;

import ComplexCalculate.UnboundDSSP;
import ComplexStructure.*;
import StructureAnalysis.LogOdds;
import StructureAnalysis.StructureBasedAnalysis;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author OpenGress
 */
public class ExtractLogOdds {
    private double []lo = new double[20];
    private double []loce =new double[20];;
    private double []lononce=new double[20];
    
    public ExtractLogOdds(){
        this.calculateLogOdds();
    }

    public double[] getLo() {
        return lo;
    }

    public void setLo(double[] lo) {
        this.lo = lo;
    }

    public double[] getLoce() {
        return loce;
    }

    public void setLoce(double[] loce) {
        this.loce = loce;
    }

    public double[] getLononce() {
        return lononce;
    }

    public void setLononce(double[] lononce) {
        this.lononce = lononce;
    }
    
    //hitung log odd ratio masing-masing asam amino
    //plot nilainya ke masing masing residu
    //tuliskan kedalam file LogOdd.arff
    public void calculateLogOdds(){
        LogOdds lod = new LogOdds();
        lod.LogOddsAnalysisWithDiffBg(this.lo,this.loce,this.lononce);
    }
    
    public void normalizeValue(double[]val,double[]newVal){
        double minValue =Double.MAX_VALUE;
        double maxValue = Double.MIN_VALUE;
        for(int i=0;i<20;i++){
            double temp = val[i];
            if(minValue >val[i]){
                minValue = val[i];
            }
            if(maxValue < val[i]){
                maxValue = val[i];
            }
        }
        for(int j=0;j<20;j++){
            newVal[j]= (val[j]-minValue)/(maxValue-minValue);
        }
        
        
    }
    
    public void ExtractLogOddsFromBenchmark(double threshold){
        ExtractDataBenchmark edb = new ExtractDataBenchmark();
        ArrayList<Epitope> epilist = edb.loadBenchmarkData("epitopelist.txt");
    //buat obyek UnboundDSSP untuk setiap pdbid
       // System.out.println("amino acid vector length");
        for(Epitope e:epilist){
            String pdbid = e.getPdbid();
            String chain = e.getChain();
            String chain_num = e.getNo_chain();
            Complex complex = new Complex(pdbid, chain);
            //complex.getChainByChianID(chain.charAt(0)).setEpitopeStatus(e, chain);
            complex.getChainByChianID(chain.charAt(0)).setAAAsEpitope(e);
            //System.out.println (pdbid+"\t"+complex.getChainByChianID(chain.charAt(0)).getAminoacidVector().size());
           // complex.getChainByChianID(chainid.charAt(0)).getAminoacidVector().size();
          //  System.out.println("panjang sequence"+ e.getSequence().length());
           // complex.setEpitopeStatus(e, chainid.charAt(0));
            UnboundDSSP unbound_dssp = new UnboundDSSP(complex, chain);
            //
            this.writeWindowedLO(complex, chain.charAt(0), 9, threshold);
            
            
        }
    }
    public void normalizeLO(double[]newlo,double[]newloce, double[]newlononce){
         this.normalizeValue(loce, newloce);
         this.normalizeValue(lo, newlononce);
         for(int k=0;k<20;k++){
                newlo[k]= newloce[k] - newlononce[k];
          }
    }
    public void calculateStructWindSmoothLO(Complex c, char chainId,double distc){
        StructureBasedAnalysis sba = new StructureBasedAnalysis();
        ArrayList<Aminoacid> listaa = c.getAminoacidVector();
         ArrayList<NeighborResidu> neighbor = new ArrayList();
         int j=0;
         for(Aminoacid aa: listaa){
             neighbor = sba.NeighboorAminoAcid(aa, c, distc);
             double cumaa =0;
             for(int i=0;i<neighbor.size();i++)
                  cumaa += this.lo[neighbor.get(i).getAa().letterToNum()];
             cumaa = cumaa/neighbor.size();
             c.getAminoacidVector().get(j).setLostructwinsmooth(cumaa);
             j+=1;
         }
    }
    public void calculateSeqWindSmoothLO(Complex c, char chainId,int winsz){
         ArrayList<Aminoacid> listaa = c.getAminoacidVector();
         ArrayList<Integer> lolist = new ArrayList();
         int j=0;
         for(Aminoacid aa: listaa){
             lolist = c.getNeighborhoodAA(aa, winsz);
             double cumaa =0;
             for(int i=0;i<lolist.size();i++)
                  cumaa += this.lo[lolist.get(i)];
             cumaa = cumaa/lolist.size();
             c.getAminoacidVector().get(j).setLoseqwinsmooth(cumaa);
             j+=1;
         }
    }
    public void writeWindowedLO(Complex c, char chainId,int winsz, double threshold){
        ArrayList<Aminoacid> listaa = c.getAminoacidVector();
        //double []windowedLo = new double[listaa.size()];
        int idx=0;
        BufferedWriter bw = null;
	FileWriter fw = null;
	try {
            fw = new FileWriter(new File("src/datatraining/Features/LO"+winsz+".arff"),true);
            bw = new BufferedWriter(fw);
            ArrayList<Integer> lolist = new ArrayList();
            for(Aminoacid aa: listaa){
                boolean exposedStatus = aa.getRsaTien2013()>threshold ;
                boolean exposedStatusOfEpi = aa.isAsEpitope()&& exposedStatus;
        
                if(exposedStatus||exposedStatusOfEpi){
                    lolist = c.getNeighborhoodAA(aa, winsz);
                    double cumaa =0;
                    for(int i=0;i<lolist.size();i++)
                        cumaa += this.lo[lolist.get(i)];
                    cumaa = cumaa/lolist.size();
                    String aaId = aa.getPdbid()+"_"+aa.getchainID()+"_"+aa.getresidueID()+"_"+aa.letterToNum();
                    String extendedSeqLO = aaId + "," +cumaa+","+aa.isAsEpitope();
                     bw.write(extendedSeqLO);
                }  
                     
            }
             bw.close();
        }catch(Exception e){
            e.printStackTrace();
        }
    }
    public void calculateLOTunggal(Complex c, char chainId){
       Chain chain = c.getChainByChianID(chainId);
       for (int i=0; i<chain.getaminoacidNum();i++){
          Aminoacid aa= chain.getAminoacidbyorder(i);
          int idxaa = aa.letterToNum();
          c.getChainByChianID(chainId).getAminoacidbyorder(i).setLo(lo[idxaa]);
          
      }
        
    }
    public void writeLOTunggalToFile(Complex c, char chainId,int winsz){
       Chain chain = c.getChainByChianID(chainId);
      BufferedWriter bw = null;
	FileWriter fw = null;
	try {
            String createdfile;
            fw = new FileWriter(new File("src/datatraining/LO/normalizelogoddstunggal.arff"),true);
            bw = new BufferedWriter(fw);
            /*//tulis disini 
            String content="%1. Title: Amino acid Log Odds Database"+"\n";
            content += "% 2. Sources:"+"\n";
            content += "%" + "\t"+"(a) Creator: B. Solihah"+"\n";
            content += "%" + "\t"+"(b) Date: March, 2018"+"\n";
            content += "%"+"\n";
            content += "@ Relation LogOddsTunggal"+"\n";
            
            
            content +="@ATTRIBUTE logodds 	REAL"+"\n";
            content +="@ATTRIBUTE logodds_ce	REAL"+"\n";
            content +="@ATTRIBUTE logodds_cnone	REAL"+"\n";
            content += "@ATTRIBUTE class 	{true, false}"+"\n";
            content +="@DATA"+"\n";
            bw.write(content);
            */
      for (int i=0; i<chain.getaminoacidNum();i++){
          Aminoacid aa= chain.getAminoacidbyorder(i);
          if(aa.getRsaRose1985()>0.05 || aa.getRsaMiller1987()>0.05 ||aa.getRsaRostSander2005()>0.05 ||aa.getRsaTien2013()>0.05||aa.isAsEpitope()){
            int idxaa = aa.letterToNum();
            
              
             
              String dataline= lo[idxaa] +","+aa.isAsEpitope()+"\n";
             // String dataline= c.getComplexName()+"_"+chainId+"_"+aa.getresidueID()+"_"+aa.getresidueName1()+"\t"+aa.isAsEpitope()+"\t"+aa.getpASA().getCn()+"\n";
          bw.write(dataline);
          }
          
          
      }
     
      bw.close();
		
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
