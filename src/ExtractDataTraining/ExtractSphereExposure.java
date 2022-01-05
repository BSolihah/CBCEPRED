/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ExtractDataTraining;

import ComplexCalculate.UnboundDSSP;
import ComplexStructure.*;
import StructureAnalysis.ParameterASA;
import StructureAnalysis.StructureBasedAnalysis;
import java.io.*;
import java.util.ArrayList;

/**
 *
 * @author OpenGress
 */

public class ExtractSphereExposure {
    public ExtractSphereExposure() throws IOException{
        /*
     ExtractDataBenchmark edb = new ExtractDataBenchmark();
//epitope list data training        
    ArrayList<Epitope> epilist = edb.loadBenchmarkData("epitopelist.txt");
//ArrayList<Epitope> epilist = edb.loadBenchmarkData("epitopelisttestset.txt");
    //buat obyek UnboundDSSP untuk setiap pdbid
       // System.out.println("amino acid vector length");
        
        for(Epitope e:epilist){
            String pdbid = e.getPdbid();
            String chain = e.getChain();
            String chain_num = e.getNo_chain();
            Complex complex = new Complex(pdbid, chain);
            //complex.getChainByChianID(chain.charAt(0)).setEpitopeStatus(e, chain);
            complex.getChainByChianID(chain.charAt(0)).setAAAsEpitope(e);
            int numaa = complex.getChainByChianID(chain.charAt(0)).getaminoacidNum();
            //System.out.println (pdbid+"\t"+complex.getChainByChianID(chain.charAt(0)).getAminoacidVector().size());
           // complex.getChainByChianID(chainid.charAt(0)).getAminoacidVector().size();
          //  System.out.println("panjang sequence"+ e.getSequence().length());
           // complex.setEpitopeStatus(e, chainid.charAt(0));
            UnboundDSSP unbound_dssp = new UnboundDSSP(complex, chain);
            //lakukan perhitungan sphere exposure pada residu terekspose dari complex
            //this.writeSphereExposureToFile(complex, chain.charAt(0), pdbid);
           // this.writeHSEToFile(complex, chain.charAt(0));
           // this.writeQSEToFile(complex, chain.charAt(0));
           // this.writeEFSEToFile(complex, chain.charAt(0));
           // this.writeSFSEToFile(complex, chain.charAt(0));
           //  this.writeFSEToFile(complex, chain.charAt(0));
            //     this.writeCNToFile(complex, chain.charAt(0));
            
            
        }*/
    }
    
    
    public ExtractSphereExposure(int noepitope) throws IOException{
    
         ExtractDataBenchmark edb = new ExtractDataBenchmark();
        ArrayList<Epitope> epilist = edb.loadBenchmarkData("epitopelist.txt");
        Epitope e = epilist.get(noepitope);
         String pdbid = e.getPdbid();
            String chain = e.getChain();
            String chain_num = e.getNo_chain();
        Complex complex = new Complex(pdbid, chain);
            //complex.getChainByChianID(chain.charAt(0)).setEpitopeStatus(e, chain);
            complex.getChainByChianID(chain.charAt(0)).setAAAsEpitope(e);
            int numaa = complex.getChainByChianID(chain.charAt(0)).getaminoacidNum();
            //System.out.println (pdbid+"\t"+complex.getChainByChianID(chain.charAt(0)).getAminoacidVector().size());
           // complex.getChainByChianID(chainid.charAt(0)).getAminoacidVector().size();
          //  System.out.println("panjang sequence"+ e.getSequence().length());
           // complex.setEpitopeStatus(e, chainid.charAt(0));
            UnboundDSSP unbound_dssp = new UnboundDSSP(complex, chain);
            //lakukan perhitungan sphere exposure pada residu terekspose dari complex
           // this.writeSphereExposureToFile(complex, chain.charAt(0), pdbid);
             
    
    }
    
   public String cekNeighboorAtom(Aminoacid aa, Complex c){
       StructureBasedAnalysis sba = new StructureBasedAnalysis();
       int cn= sba.ContactNumber(aa, c);
       int na= sba.listNeighboorAtom(aa, c).size();
       ArrayList<AtomCoordinate> frameref = sba.calculateFrameRef(aa);
        ArrayList<AtomCoordinate> newneighbooratomcoord= sba.defineCoordOnFrameRef(frameref, sba.listNeighboorAtom(aa, c));
       String cek ="\t "+ cn+"\t"+ na;
       return cek;
   }
   
    public void calculateSE(Aminoacid aa, Complex complex){
         StructureBasedAnalysis sba = new StructureBasedAnalysis();
                //cek disini untuk print koordinat
                int[]efse = sba.calculateEightSubAreaFSE(aa, complex);
                
                int[]sfse = sba.calculateSixSubAreaFSE(aa, complex);
                int[] fse= sba.calculateFSE(aa, complex);
                int[] qse = sba.calculateQSE(aa, complex);
                int[] hse = sba.calculateHSE(qse);
                int cnfromhse = sba.calculateCNFromHSE(hse);
                int cn = sba.ContactNumber(aa, complex);
                ParameterASA pAsa = new ParameterASA();
                pAsa.setCn(cn);
                pAsa.setCnfromHSE(cnfromhse);
                pAsa.setFse(fse);
                pAsa.setHse(hse);
                pAsa.setQse(qse);
                pAsa.setEfse(efse);
                pAsa.setSfse(sfse);
                aa.setpASA(pAsa);
    }
     public void  writeSphereExposureToFile(Complex c, char chainId,String filename) throws IOException {
      
      Chain chain = c.getChainByChianID(chainId);
      
     // File fout = new File("src/datatraining/SE/"+"SE_"+filename+".txt");
      File fout = new File("src/datatraining/SE/"+"SE_"+"all"+".txt");
      
  //    FileOutputStream fos = new FileOutputStream(fout);
      FileWriter fw = new FileWriter(fout,true);
      //BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fos));
      BufferedWriter bw = new BufferedWriter(fw);
      for (int i=0; i<chain.getaminoacidNum();i++){
          Aminoacid aa= chain.getAminoacidbyorder(i);
          if(aa.getRsaRose1985()>0.05 || aa.getRsaMiller1987()>0.05 ||aa.getRsaRostSander2005()>0.05 ||aa.getRsaTien2013()>0.05||aa.isAsEpitope()){
              this.calculateSE(aa, c);
              if(aa.getpASA().getCnfromHSE()==0){
                  System.out.println(c.getComplexName()+"_"+chainId+"_"+aa.getresidueID()+"_"+aa.getresidueName1()+ this.cekNeighboorAtom(aa, c));
              }
              String dataline= c.getComplexName()+"_"+chainId+"_"+aa.getresidueID()+"_"+aa.getresidueName1()+"\t"+aa.isAsEpitope()+"\t"+aa.getpASA().printParamASA()+"\n";
             // String dataline= c.getComplexName()+"_"+chainId+"_"+aa.getresidueID()+"_"+aa.getresidueName1()+"\t"+aa.isAsEpitope()+"\t"+aa.getpASA().getCn()+"\n";
          bw.write(dataline);
          }
          
          
      }
     
      bw.close();
  }  
  public String printHSE(Aminoacid aa, Complex c){
       String dataline=null;
      if(aa.getRsaRose1985()>0.05 || aa.getRsaMiller1987()>0.05 ||aa.getRsaRostSander2005()>0.05 ||aa.getRsaTien2013()>0.05||aa.isAsEpitope()){
              this.calculateSE(aa, c);
             dataline= aa.getpASA().printHSE()+aa.isAsEpitope()+"\n";
             // String dataline= c.getComplexName()+"_"+chainId+"_"+aa.getresidueID()+"_"+aa.getresidueName1()+"\t"+aa.isAsEpitope()+"\t"+aa.getpASA().getCn()+"\n";
          
          }
      return dataline;
      
  }
     
  public void writeHSEToFile(Complex c, char chainId){
       Chain chain = c.getChainByChianID(chainId);
      BufferedWriter bw = null;
	FileWriter fw = null;
	try {
            String createdfile;
            fw = new FileWriter(new File("src/datatraining/SE/HSE.arff"),true);
            bw = new BufferedWriter(fw);
            /*//tulis disini 
            String content="%1. Title: Amino acid HSE Database"+"\n";
            content += "% 2. Sources:"+"\n";
            content += "%" + "\t"+"(a) Creator: B. Solihah"+"\n";
            content += "%" + "\t"+"(b) Date: March, 2018"+"\n";
            content += "%"+"\n";
            content += "@ Relation HSE"+"\n";
            
            content +="@ATTRIBUTE residuName  string"+"\n";
            content +="@ATTRIBUTE HSEup	REAL"+"\n";
            content +="@ATTRIBUTE HSEdown	REAL"+"\n";
            content += "@ATTRIBUTE class 	{true, false}"+"\n";
            content +="@DATA"+"\n";
            bw.write(content);
            */
      for (int i=0; i<chain.getaminoacidNum();i++){
          Aminoacid aa= chain.getAminoacidbyorder(i);
          if(aa.getRsaRose1985()>0.05 || aa.getRsaMiller1987()>0.05 ||aa.getRsaRostSander2005()>0.05 ||aa.getRsaTien2013()>0.05||aa.isAsEpitope()){
              this.calculateSE(aa, c);
              
             // String dataline= aa.getresidueName1()+","+aa.getpASA().printHSE()+aa.isAsEpitope()+"\n";
               String dataline= aa.letterToNum()+","+aa.getpASA().printHSE()+aa.isAsEpitope()+"\n";
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
  public void writeQSEToFile(Complex c, char chainId){
       Chain chain = c.getChainByChianID(chainId);
      BufferedWriter bw = null;
	FileWriter fw = null;
	try {
            String createdfile;
            fw = new FileWriter(new File("src/datatraining/SE/QSE.arff"),true);
            bw = new BufferedWriter(fw);
           /* //tulis disini 
            String content="%1. Title: Amino acid QSE Database"+"\n";
            content += "% 2. Sources:"+"\n";
            content += "%" + "\t"+"(a) Creator: B. Solihah"+"\n";
            content += "%" + "\t"+"(b) Date: March, 2018"+"\n";
            content += "%"+"\n";
            content += "@ Relation QSE"+"\n";
            
            content +="@ATTRIBUTE residuName  String"+"\n";
            content +="@ATTRIBUTE QSE1        REAL"+"\n";
            content +="@ATTRIBUTE QSE2	      REAL"+"\n";
            content +="@ATTRIBUTE QSE3        REAL"+"\n";
            content +="@ATTRIBUTE QSE4	      REAL"+"\n";
            content +="@ATTRIBUTE QSE5        REAL"+"\n";
            content +="@ATTRIBUTE QSE6	      REAL"+"\n";
            content +="@ATTRIBUTE QSE7        REAL"+"\n";
            content +="@ATTRIBUTE QSE8	      REAL"+"\n";
            
            content += "@ATTRIBUTE class 	{true, false}"+"\n";
            content +="@DATA"+"\n";
            bw.write(content);
            */
      for (int i=0; i<chain.getaminoacidNum();i++){
          Aminoacid aa= chain.getAminoacidbyorder(i);
          if(aa.getRsaRose1985()>0.05 || aa.getRsaMiller1987()>0.05 ||aa.getRsaRostSander2005()>0.05 ||aa.getRsaTien2013()>0.05||aa.isAsEpitope()){
              this.calculateSE(aa, c);
              
             // String dataline= aa.getresidueName1()+","+aa.getpASA().printQSE()+aa.isAsEpitope()+"\n";
               String dataline= aa.getpASA().printQSE()+aa.isAsEpitope()+"\n";
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
  public void writeFSEToFile(Complex c, char chainId){
       Chain chain = c.getChainByChianID(chainId);
      BufferedWriter bw = null;
	FileWriter fw = null;
	try {
            String createdfile;
            fw = new FileWriter(new File("src/datatraining/SE/FSE.arff"),true);
            bw = new BufferedWriter(fw);
            /*//tulis disini 
            String content="%1. Title: Amino acid FSE Database"+"\n";
            content += "% 2. Sources:"+"\n";
            content += "%" + "\t"+"(a) Creator: B. Solihah"+"\n";
            content += "%" + "\t"+"(b) Date: March, 2018"+"\n";
            content += "%"+"\n";
            content += "@ Relation FSE"+"\n";
            
            
            content +="@ATTRIBUTE FSE1        REAL"+"\n";
            content +="@ATTRIBUTE FSE2	      REAL"+"\n";
            content +="@ATTRIBUTE FSE3        REAL"+"\n";
            content +="@ATTRIBUTE FSE4	      REAL"+"\n";
            content +="@ATTRIBUTE FSE5        REAL"+"\n";
            content +="@ATTRIBUTE FSE6	      REAL"+"\n";
            content +="@ATTRIBUTE FSE7        REAL"+"\n";
            content +="@ATTRIBUTE FSE8	      REAL"+"\n";
            
            content += "@ATTRIBUTE class 	{true, false}"+"\n";
            content +="@DATA"+"\n";
            bw.write(content);
            */
      for (int i=0; i<chain.getaminoacidNum();i++){
          Aminoacid aa= chain.getAminoacidbyorder(i);
          if(aa.getRsaRose1985()>0.05 || aa.getRsaMiller1987()>0.05 ||aa.getRsaRostSander2005()>0.05 ||aa.getRsaTien2013()>0.05||aa.isAsEpitope()){
              this.calculateSE(aa, c);
              
              String dataline= aa.getpASA().printFSE()+aa.isAsEpitope()+"\n";
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
  public void writeCNToFile(Complex c, char chainId){
       Chain chain = c.getChainByChianID(chainId);
      BufferedWriter bw = null;
	FileWriter fw = null;
	try {
            String createdfile;
            fw = new FileWriter(new File("src/datatraining/CN/CN.arff"),true);
            bw = new BufferedWriter(fw);
            
      for (int i=0; i<chain.getaminoacidNum();i++){
          Aminoacid aa= chain.getAminoacidbyorder(i);
          if(aa.getRsaRose1985()>0.05 || aa.getRsaMiller1987()>0.05 ||aa.getRsaRostSander2005()>0.05 ||aa.getRsaTien2013()>0.05||aa.isAsEpitope()){
              this.calculateSE(aa, c);
              
              String dataline= aa.getpASA().getCn()+","+aa.isAsEpitope()+"\n";
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
  public void writeEFSEToFile(Complex c, char chainId){
       Chain chain = c.getChainByChianID(chainId);
      BufferedWriter bw = null;
	FileWriter fw = null;
	try {
            String createdfile;
            fw = new FileWriter(new File("src/datatraining/SE/EFSE.arff"),true);
            bw = new BufferedWriter(fw);
            /*//tulis disini 
            String content="%1. Title: Amino acid EFSE Database"+"\n";
            content += "% 2. Sources:"+"\n";
            content += "%" + "\t"+"(a) Creator: B. Solihah"+"\n";
            content += "%" + "\t"+"(b) Date: March, 2018"+"\n";
            content += "%"+"\n";
            content += "@ Relation EFSE"+"\n";
            
            content +="@ATTRIBUTE residuName  String"+"\n";
            content +="@ATTRIBUTE EFSE1        REAL"+"\n";
            content +="@ATTRIBUTE EFSE2	      REAL"+"\n";
            content +="@ATTRIBUTE EFSE3        REAL"+"\n";
            content +="@ATTRIBUTE EFSE4	      REAL"+"\n";
            content +="@ATTRIBUTE EFSE5        REAL"+"\n";
            content +="@ATTRIBUTE EFSE6	      REAL"+"\n";
            content +="@ATTRIBUTE EFSE7        REAL"+"\n";
            content +="@ATTRIBUTE EFSE8	      REAL"+"\n";
            
            content += "@ATTRIBUTE class 	{true, false}"+"\n";
            content +="@DATA"+"\n";
            bw.write(content);
            */
      for (int i=0; i<chain.getaminoacidNum();i++){
          Aminoacid aa= chain.getAminoacidbyorder(i);
          if(aa.getRsaRose1985()>0.05 || aa.getRsaMiller1987()>0.05 ||aa.getRsaRostSander2005()>0.05 ||aa.getRsaTien2013()>0.05||aa.isAsEpitope()){
              this.calculateSE(aa, c);
              
              String dataline= aa.getpASA().printEFSE()+aa.isAsEpitope()+"\n";
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
  public void writeSFSEToFile(Complex c, char chainId){
       Chain chain = c.getChainByChianID(chainId);
      BufferedWriter bw = null;
	FileWriter fw = null;
	try {
            String createdfile;
            fw = new FileWriter(new File("src/datatraining/SE/SFSE.arff"),true);
            bw = new BufferedWriter(fw);
            /*//tulis disini 
            String content="%1. Title: Amino acid SFSE Database"+"\n";
            content += "% 2. Sources:"+"\n";
            content += "%" + "\t"+"(a) Creator: B. Solihah"+"\n";
            content += "%" + "\t"+"(b) Date: March, 2018"+"\n";
            content += "%"+"\n";
            content += "@ Relation SFSE"+"\n";
            
            content +="@ATTRIBUTE residuName  string"+"\n";
            content +="@ATTRIBUTE SFSE1        REAL"+"\n";
            content +="@ATTRIBUTE SFSE2	      REAL"+"\n";
            content +="@ATTRIBUTE SFSE3        REAL"+"\n";
            content +="@ATTRIBUTE SFSE4	      REAL"+"\n";
            content +="@ATTRIBUTE SFSE5        REAL"+"\n";
            content +="@ATTRIBUTE SFSE6	      REAL"+"\n";
            content += "@ATTRIBUTE class 	{true, false}"+"\n";
           
            content +="@DATA"+"\n";
            bw.write(content);
            */
      for (int i=0; i<chain.getaminoacidNum();i++){
          Aminoacid aa= chain.getAminoacidbyorder(i);
          if(aa.getRsaRose1985()>0.05 || aa.getRsaMiller1987()>0.05 ||aa.getRsaRostSander2005()>0.05 ||aa.getRsaTien2013()>0.05||aa.isAsEpitope()){
              this.calculateSE(aa, c);
              
              String dataline= aa.getpASA().printSFSE()+aa.isAsEpitope()+"\n";
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
