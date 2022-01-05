package ComplexCalculate;

import ComplexStructure.Aminoacid;
import ComplexStructure.Atom;
import ComplexStructure.Chain;
import ComplexStructure.Complex;
import java.io.*;
import java.util.ArrayList;
import java.util.Vector;

public class UnboundDSSP
{
  public UnboundDSSP(Vector<Complex> unboundcomplexVector)
  {
    for (int i = 0; i < unboundcomplexVector.size(); i++) {
      DSSPfile((Complex)unboundcomplexVector.get(i));
    }
  }
  
  public UnboundDSSP(Complex complex, String chain)
  {
    String dsspfile = complex.getcomplexName()+"_"+chain;
    
    String dataline = null;
    FileReader fr = null;
    
    try
    {
      //if (new File("src/dssp/" + complexName + "_" + complex.getAntigenChainIDs().charAt(0) + ".dssp").exists()) {
        if (new File("src/data/dssp/" + dsspfile + ".dssp").exists()) {
        //fr = new FileReader(new File("src/dssp/" + complexName + "_" + complex.getAntigenChainIDs().charAt(0) + ".dssp"));
           fr = new FileReader(new File("src/data/dssp/" + dsspfile + ".dssp"));
      } 
        /*
        else {
        fr = new FileReader(new File("src/data/dssp/" + complexName + "_" + "save.txt"));
      }*/
        
      BufferedReader br = new BufferedReader(fr);
      dataline = br.readLine();
      while (!dataline.trim().startsWith("#")) {
        dataline = br.readLine();
      }
      dataline = br.readLine();
      //System.out.println(dataline.substring(5, 10).trim());
      
      while (dataline != null)
      {
        if ((!dataline.substring(5, 10).trim().equals("")) || (dataline.substring(5, 10).trim() == null))
        {
          String residueID = dataline.substring(5, 11).trim();
          char chainID;
        
          if (dataline.substring(11, 12).trim().equals("")) {
            chainID = 'X';
          } else {
            chainID = dataline.substring(11, 12).trim().charAt(0);
          }
          String residueName = dataline.substring(12, 14).trim();
          
          char temp = chainID;
          char secondaryStructure;

          if (dataline.substring(16, 17).trim().equals("")) {
            secondaryStructure = '\000';
          } else {
            secondaryStructure = dataline.substring(16, 17).trim().charAt(0);
          }
          try{
            //  String sarea = dataline.substring(35, 38).trim();
          //System.out.print(sarea+"\t");
          double area = Double.valueOf(dataline.substring(35, 38).trim()).doubleValue();
         
          if (complex.getChainByChianID(chainID).getAminoacidbyIDs(residueID) != null)
          {
            complex.getChainByChianID(chainID).getAminoacidbyIDs(residueID).setsecondaryStructure(calculateSecondaryStructure(secondaryStructure));
            complex.getChainByChianID(chainID).getAminoacidbyIDs(residueID).setsolventAccessibility(area);
            //rsa diset berdasarkan sejumlah metode yang dijadikan rujukan
            //complex.getChainByChianID(chainID).getAminoacidbyIDs(residueID).setRsa(this.calculateRSAByTienMSA(residueName, area));
            complex.getChainByChianID(chainID).getAminoacidbyIDs(residueID).setRsaMiller1987(this.calculateRSAMiller(residueName, area));
            complex.getChainByChianID(chainID).getAminoacidbyIDs(residueID).setRsaRose1985(this.calculateRSARose(residueName, area));
            complex.getChainByChianID(chainID).getAminoacidbyIDs(residueID).setRsaRostSander2005(this.calculateRSARostAndSander(residueName, area));
            complex.getChainByChianID(chainID).getAminoacidbyIDs(residueID).setRsaTien2013(this.calculateRSAByTienMSA(residueName, area));
          }
          }catch(Exception e){
              System.out.println( chainID +";"+ residueID);
              //e.printStackTrace();
          }
          
          
          
        }
        dataline = br.readLine();
      }
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
  }
  
  
  public UnboundDSSP(Complex complex)
  {
    String complexName = complex.getcomplexName();
    
    String dataline = null;
    FileReader fr = null;
    try
    {
      //if (new File("src/dssp/" + complexName + "_" + complex.getAntigenChainIDs().charAt(0) + ".dssp").exists()) {
        if (new File("src/data/dssp/" + complexName + ".dssp").exists()) {
        //fr = new FileReader(new File("src/dssp/" + complexName + "_" + complex.getAntigenChainIDs().charAt(0) + ".dssp"));
           fr = new FileReader(new File("src/data/dssp/" + complexName + ".dssp"));
      } 
        /*
        else {
        fr = new FileReader(new File("src/data/dssp/" + complexName + "_" + "save.txt"));
      }*/
      BufferedReader br = new BufferedReader(fr);
      dataline = br.readLine();
      while (!dataline.trim().startsWith("#")) {
        dataline = br.readLine();
      }
      dataline = br.readLine();
      
      while (dataline != null)
      {
        if ((!dataline.substring(5, 10).trim().equals("")) || (dataline.substring(5, 10).trim() == null))
        {
          String residueID = dataline.substring(5, 11).trim();
          char chainID;

          if (dataline.substring(11, 12).trim().equals("")) {
            chainID = 'X';
          } else {
            chainID = dataline.substring(11, 12).trim().charAt(0);
          }
          String residueName = dataline.substring(12, 14).trim();
          
          char temp = chainID;
          char secondaryStructure;

          if (dataline.substring(16, 17).trim().equals("")) {
            secondaryStructure = '\000';
          } else {
            secondaryStructure = dataline.substring(16, 17).trim().charAt(0);
          }
          try{
            //  String sarea = dataline.substring(35, 38).trim();
          //System.out.print(sarea+"\t");
          double area = Double.valueOf(dataline.substring(35, 38).trim()).doubleValue();
          if (complex.getChainByChianID(chainID).getAminoacidbyIDs(residueID) != null)
          {
            complex.getChainByChianID(chainID).getAminoacidbyIDs(residueID).setsecondaryStructure(calculateSecondaryStructure(secondaryStructure));
            complex.getChainByChianID(chainID).getAminoacidbyIDs(residueID).setsolventAccessibility(area);
            complex.getChainByChianID(chainID).getAminoacidbyIDs(residueID).setRsa(this.calculateRSAByTienMSA(residueName, area));
          }
          }catch(Exception e){
          }
          
          
          
        }
        dataline = br.readLine();
      }
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
  }
  
  private void DSSPfile(Complex unboundcomplex)
  {
    String complexName = unboundcomplex.getcomplexName();
    String dataline = null;
    FileReader fr = null;
    try
    {
      fr = new FileReader(new File("train/??dssp??/" + complexName + "_save.txt"));
      BufferedReader br = new BufferedReader(fr);
      dataline = br.readLine();
      while (!dataline.trim().startsWith("#")) {
        dataline = br.readLine();
      }
      dataline = br.readLine();
      while (dataline != null)
      {
        if ((!dataline.substring(5, 10).trim().equals("")) || (dataline.substring(5, 10).trim() == null))
        {
          String residueID = dataline.substring(5, 11).trim();
          char chainID;

          if (dataline.substring(11, 12).trim().equals("")) {
            chainID = 'X';
          } else {
            chainID = dataline.substring(11, 12).trim().charAt(0);
          }
          String residueName = dataline.substring(12, 14).trim();
          
          char temp = chainID;
          char secondaryStructure;

          if (dataline.substring(16, 17).trim().equals("")) {
            secondaryStructure = '\000';
          } else {
            secondaryStructure = dataline.substring(16, 17).trim().charAt(0);
          }
          double area = Double.valueOf(dataline.substring(35, 38).trim()).doubleValue();
          if (unboundcomplex.getChainByChianID(chainID).getAminoacidbyIDs(residueID) != null)
          {
            unboundcomplex.getChainByChianID(chainID).getAminoacidbyIDs(residueID).setsecondaryStructure(calculateSecondaryStructure(secondaryStructure));
            //asa
            unboundcomplex.getChainByChianID(chainID).getAminoacidbyIDs(residueID).setsolventAccessibility(area);
            //rsa
            
          }
        }
        dataline = br.readLine();
      }
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
  }
  
  private double calculateRSARostAndSander(String residueName, double area)
  {
    if ((residueName.equals("ALA")) || (residueName.equals("A"))) {
      return area / 106.0D;
    }
    if ((residueName.equals("ARG")) || (residueName.equals("R"))) {
      return area / 248.0D;
    }
    if ((residueName.equals("ASN")) || (residueName.equals("N"))) {
      return area / 157.0D;
    }
    if ((residueName.equals("ASP")) || (residueName.equals("D"))) {
      return area / 163.0D;
    }
    if ((residueName.equals("CYS")) || (residueName.equals("C"))) {
      return area / 135.0D;
    }
    if ((residueName.equals("GLN")) || (residueName.equals("Q"))) {
      return area / 198.0D;
    }
    if ((residueName.equals("GLU")) || (residueName.equals("E"))) {
      return area / 194.0D;
    }
    if ((residueName.equals("GLY")) || (residueName.equals("G"))) {
      return area / 84.0D;
    }
    if ((residueName.equals("HIS")) || (residueName.equals("H"))) {
      return area / 184.0D;
    }
    if ((residueName.equals("ILE")) || (residueName.equals("I"))) {
      return area / 169.0D;
    }
    if ((residueName.equals("LEU")) || (residueName.equals("L"))) {
      return area / 164.0D;
    }
    if ((residueName.equals("LYS")) || (residueName.equals("K"))) {
      return area / 205.0D;
    }
    if ((residueName.equals("MET")) || (residueName.equals("M"))) {
      return area / 188.0D;
    }
    if ((residueName.equals("PHE")) || (residueName.equals("F"))) {
      return area / 197.0D;
    }
    if ((residueName.equals("PRO")) || (residueName.equals("P"))) {
      return area / 136.0D;
    }
    if ((residueName.equals("SER")) || (residueName.equals("S"))) {
      return area / 130.0D;
    }
    if ((residueName.equals("THR")) || (residueName.equals("T"))) {
      return area / 142.0D;
    }
    if ((residueName.equals("TRP")) || (residueName.equals("W"))) {
      return area / 227.0D;
    }
    if ((residueName.equals("TYR")) || (residueName.equals("Y"))) {
      return area / 222.0D;
    }
    if ((residueName.equals("VAL")) || (residueName.equals("V"))) {
      return area / 142.0D;
    }
    return -1.0D;
  }
  private double calculateRSAByTienMSA(String residueName, double area)
  {
    if ((residueName.equals("ALA")) || (residueName.equals("A"))) {
      return area / 129.0D;
    }
    if ((residueName.equals("ARG")) || (residueName.equals("R"))) {
      return area / 274.0D;
    }
    if ((residueName.equals("ASN")) || (residueName.equals("N"))) {
      return area / 195.0D;
    }
    if ((residueName.equals("ASP")) || (residueName.equals("D"))) {
      return area / 193.0D;
    }
    if ((residueName.equals("CYS")) || (residueName.equals("C"))) {
      return area / 167.0D;
    }
    if ((residueName.equals("GLN")) || (residueName.equals("Q"))) {
      return area / 223.0D;
    }
    if ((residueName.equals("GLU")) || (residueName.equals("E"))) {
      return area / 225.0D;
    }
    if ((residueName.equals("GLY")) || (residueName.equals("G"))) {
      return area / 104.0D;
    }
    if ((residueName.equals("HIS")) || (residueName.equals("H"))) {
      return area / 224.0D;
    }
    if ((residueName.equals("ILE")) || (residueName.equals("I"))) {
      return area / 197.0D;
    }
    if ((residueName.equals("LEU")) || (residueName.equals("L"))) {
      return area / 201.0D;
    }
    if ((residueName.equals("LYS")) || (residueName.equals("K"))) {
      return area / 236.0D;
    }
    if ((residueName.equals("MET")) || (residueName.equals("M"))) {
      return area / 224.0D;
    }
    if ((residueName.equals("PHE")) || (residueName.equals("F"))) {
      return area / 240.0D;
    }
    if ((residueName.equals("PRO")) || (residueName.equals("P"))) {
      return area / 159.0D;
    }
    if ((residueName.equals("SER")) || (residueName.equals("S"))) {
      return area / 155.0D;
    }
    if ((residueName.equals("THR")) || (residueName.equals("T"))) {
      return area / 172.0D;
    }
    if ((residueName.equals("TRP")) || (residueName.equals("W"))) {
      return area / 285.0D;
    }
    if ((residueName.equals("TYR")) || (residueName.equals("Y"))) {
      return area / 263.0D;
    }
    if ((residueName.equals("VAL")) || (residueName.equals("V"))) {
      return area / 174.0D;
    }
    return -1.0D;
  }
  private double calculateRSAMiller(String residueName, double area)
  {
    if ((residueName.equals("ALA")) || (residueName.equals("A"))) {
      return area / 113.0D;
    }
    if ((residueName.equals("ARG")) || (residueName.equals("R"))) {
      return area / 241.0D;
    }
    if ((residueName.equals("ASN")) || (residueName.equals("N"))) {
      return area / 151.0D;
    }
    if ((residueName.equals("ASP")) || (residueName.equals("D"))) {
      return area / 151.0D;
    }
    if ((residueName.equals("CYS")) || (residueName.equals("C"))) {
      return area / 140.0D;
    }
    if ((residueName.equals("GLN")) || (residueName.equals("Q"))) {
      return area / 183.0D;
    }
    if ((residueName.equals("GLU")) || (residueName.equals("E"))) {
      return area / 189.0D;
    }
    if ((residueName.equals("GLY")) || (residueName.equals("G"))) {
      return area / 85.0D;
    }
    if ((residueName.equals("HIS")) || (residueName.equals("H"))) {
      return area / 194.0D;
    }
    if ((residueName.equals("ILE")) || (residueName.equals("I"))) {
      return area / 182.0D;
    }
    if ((residueName.equals("LEU")) || (residueName.equals("L"))) {
      return area / 180.0D;
    }
    if ((residueName.equals("LYS")) || (residueName.equals("K"))) {
      return area / 211.0D;
    }
    if ((residueName.equals("MET")) || (residueName.equals("M"))) {
      return area / 204.0D;
    }
    if ((residueName.equals("PHE")) || (residueName.equals("F"))) {
      return area / 218.0D;
    }
    if ((residueName.equals("PRO")) || (residueName.equals("P"))) {
      return area / 143.0D;
    }
    if ((residueName.equals("SER")) || (residueName.equals("S"))) {
      return area / 122.0D;
    }
    if ((residueName.equals("THR")) || (residueName.equals("T"))) {
      return area / 146.0D;
    }
    if ((residueName.equals("TRP")) || (residueName.equals("W"))) {
      return area / 259.0D;
    }
    if ((residueName.equals("TYR")) || (residueName.equals("Y"))) {
      return area / 229.0D;
    }
    if ((residueName.equals("VAL")) || (residueName.equals("V"))) {
      return area / 160.0D;
    }
    return -1.0D;
  }
  private double calculateRSARose(String residueName, double area)
  {
    if ((residueName.equals("ALA")) || (residueName.equals("A"))) {
      return area / 118.1D;
    }
    if ((residueName.equals("ARG")) || (residueName.equals("R"))) {
      return area / 256.0D;
    }
    if ((residueName.equals("ASN")) || (residueName.equals("N"))) {
      return area / 158.7D;
    }
    if ((residueName.equals("ASP")) || (residueName.equals("D"))) {
      return area / 158.0D;
    }
    if ((residueName.equals("CYS")) || (residueName.equals("C"))) {
      return area / 146.1D;
    }
    if ((residueName.equals("GLN")) || (residueName.equals("Q"))) {
      return area / 186.2D;
    }
    if ((residueName.equals("GLU")) || (residueName.equals("E"))) {
      return area / 193.0D;
    }
    if ((residueName.equals("GLY")) || (residueName.equals("G"))) {
      return area / 88.1D;
    }
    if ((residueName.equals("HIS")) || (residueName.equals("H"))) {
      return area / 202.5D;
    }
    if ((residueName.equals("ILE")) || (residueName.equals("I"))) {
      return area / 181.0D;
    }
    if ((residueName.equals("LEU")) || (residueName.equals("L"))) {
      return area / 193.0D;
    }
    if ((residueName.equals("LYS")) || (residueName.equals("K"))) {
      return area / 225.8D;
    }
    if ((residueName.equals("MET")) || (residueName.equals("M"))) {
      return area / 203.4D;
    }
    if ((residueName.equals("PHE")) || (residueName.equals("F"))) {
      return area / 222.8D;
    }
    if ((residueName.equals("PRO")) || (residueName.equals("P"))) {
      return area / 146.8D;
    }
    if ((residueName.equals("SER")) || (residueName.equals("S"))) {
      return area / 129.8D;
    }
    if ((residueName.equals("THR")) || (residueName.equals("T"))) {
      return area / 152.5D;
    }
    if ((residueName.equals("TRP")) || (residueName.equals("W"))) {
      return area / 266.0D;
    }
    if ((residueName.equals("TYR")) || (residueName.equals("Y"))) {
      return area / 236.0D;
    }
    if ((residueName.equals("VAL")) || (residueName.equals("V"))) {
      return area / 164.5D;
    }
    return -1.0D;
  }
  private char calculateSecondaryStructure(char ss)
  {
    if ((ss == 'H')|| (ss == 'G')) {
      return 'h';
    }
    if ((ss == 'E') || (ss == 'B')) {
      return 's';
    }
    if ((ss == 'I')  || (ss == 'T') || (ss == 'S') || (ss == 0)) {
      return 'c';
    }
    return 'u';
  }
  
  public ArrayList<Aminoacid> getExposedAA(Complex complex, char chainID, int mode){
      // dari komplek yang sudah dipanggil kontruksotnya
      ArrayList<Aminoacid> exposedAA = new ArrayList();
      Chain c = complex.getChainByChianID(chainID);
      for (int i=0; i<c.getaminoacidNum();i++){
           switch(mode){
              case 1:
                  if(c.getAminoacidbyorder(i).getRsaRose1985()<0.05){
                      c.getAminoacidbyorder(i).setPdbid(c.getPdbid());
                      exposedAA.add(c.getAminoacidbyorder(i));
                    
                    }
                  break;
              case 2:
                  if(c.getAminoacidbyorder(i).getRsaMiller1987()<0.05){
                      c.getAminoacidbyorder(i).setPdbid(c.getPdbid());
                      exposedAA.add(c.getAminoacidbyorder(i));
                    } 
                  break;
              case 3:
                  if(c.getAminoacidbyorder(i).getRsaRostSander2005()<0.05){
                      c.getAminoacidbyorder(i).setPdbid(c.getPdbid());
                      exposedAA.add(c.getAminoacidbyorder(i));
                    } 
                  break;
              case 4:
                  if(c.getAminoacidbyorder(i).getRsaTien2013()<0.05){
                      c.getAminoacidbyorder(i).setPdbid(c.getPdbid());
                      exposedAA.add(c.getAminoacidbyorder(i));
                    } 
                  break;
          }
         
          //Ambil residu aa (syarat RSA > 0.05
          
      }
     // System.out.print("jumlah exposed aa"+ exposedAA.size()+"\n");
      return exposedAA;
  }
  public ArrayList<Aminoacid> getEpitopeFromNonExposedAA(Complex complex, String chain, int mode){
      ArrayList<Aminoacid> nonexposedEpi = new ArrayList();
      Chain c = complex.getChainByChianID(chain.charAt(0));
      for (int i=0; i<c.getaminoacidNum();i++){
          //Ambil residu aa (syarat RSA > 0.05
          switch(mode){
              case 1:
                  if(c.getAminoacidbyorder(i).getRsaRose1985()<0.05&& c.getAminoacidbyorder(i).isAsEpitope()){
                      c.getAminoacidbyorder(i).setPdbid(c.getPdbid());
                      nonexposedEpi.add(c.getAminoacidbyorder(i));
                    
                    }
                  break;
              case 2:
                  if(c.getAminoacidbyorder(i).getRsaMiller1987()<0.05&& c.getAminoacidbyorder(i).isAsEpitope()){
                      c.getAminoacidbyorder(i).setPdbid(c.getPdbid());
                      nonexposedEpi.add(c.getAminoacidbyorder(i));
                    } 
                  break;
              case 3:
                  if(c.getAminoacidbyorder(i).getRsaRostSander2005()<0.05&& c.getAminoacidbyorder(i).isAsEpitope()){
                      c.getAminoacidbyorder(i).setPdbid(c.getPdbid());
                      nonexposedEpi.add(c.getAminoacidbyorder(i));
                    } 
                  break;
              case 4:
                  if(c.getAminoacidbyorder(i).getRsaTien2013()<0.05&& c.getAminoacidbyorder(i).isAsEpitope()){
                      c.getAminoacidbyorder(i).setPdbid(c.getPdbid());
                      nonexposedEpi.add(c.getAminoacidbyorder(i));
                    } 
                  break;
          }
         
          
      }
       return nonexposedEpi;
     // System.out.print("jumlah exposed aa"+ exposedAA.size()+"\n");
      
  }
  public void printExposedAA(Complex complex, char chainID, int mode){
      ArrayList<Aminoacid> aa = getExposedAA(complex,chainID, mode);
      System.out.print("jumlah  AA pada chain: "+ complex.getChainByChianID(chainID).getaminoacidNum()+"\n");
      //System.out.print("jumlah exposes AA: "+ aa.size());
      for (int i=0;i<aa.size();i++){
          System.out.print(aa.get(i).getresidueID()+"\t "+ aa.get(i).getresidueName()+"\t"+aa.get(i).getRsa()+"\n");
      }
  }
  public double[] writeRSAToFile(Complex c, char chainId, String filename)throws IOException{
      double[] minRSAOfEpitope= new double [4];
      
      double minRSARose= Double.MAX_VALUE;
      double minRSARostSander = Double.MAX_VALUE;
      double minRSAMiller = Double.MAX_VALUE;
      double minRSATien = Double.MAX_VALUE;
      Chain chain = c.getChainByChianID(chainId);
      
      File fout = new File("src/datatraining/RSA/"+"RSA_"+filename+".txt");
      FileOutputStream fos = new FileOutputStream(fout);
      BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fos));
      for (int i=0; i<chain.getaminoacidNum();i++){
          Aminoacid aa= chain.getAminoacidbyorder(i);
          if(aa.getRsaRose1985()>0.05 || aa.getRsaMiller1987()>0.05 ||aa.getRsaRostSander2005()>0.05 ||aa.getRsaTien2013()>0.05){
          String dataline= c.getComplexName()+"_"+chainId+"_"+aa.getresidueID()+"_"+aa.getresidueName1()+"\t"+aa.isAsEpitope()+"\t"+aa.getRsaRose1985()+"\t"+aa.getRsaMiller1987()+"\t"+aa.getRsaRostSander2005()+"\t"+aa.getRsaTien2013()+"\n";
          bw.write(dataline);
          }
          
          
           if(aa.isAsEpitope()){
              if(minRSARose > aa.getRsaRose1985()){
                  minRSARose = aa.getRsaRose1985();
              }
              
              if(minRSAMiller > aa.getRsaMiller1987()){
                  minRSAMiller = aa.getRsaMiller1987();
              }
              if(minRSARostSander > aa.getRsaRostSander2005()){
                  minRSARostSander = aa.getRsaRostSander2005();
              }
              if(minRSATien > aa.getRsaTien2013()){
                  minRSATien = aa.getRsaTien2013();
              }
          }
      }
      bw.close();
       minRSAOfEpitope[0]= minRSARose;
      minRSAOfEpitope[1]= minRSAMiller;
      minRSAOfEpitope[2]= minRSARostSander;
      minRSAOfEpitope[3]= minRSATien;
      return minRSAOfEpitope;
  
  }
  public void  writeASAToFile(Complex c, char chainId,String filename) throws IOException {
      
      Chain chain = c.getChainByChianID(chainId);
      
      File fout = new File("src/datatraining/ASA/"+"ASA_"+filename+".txt");
      FileOutputStream fos = new FileOutputStream(fout);
      BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fos));
      for (int i=0; i<chain.getaminoacidNum();i++){
          Aminoacid aa= chain.getAminoacidbyorder(i);
          String dataline= c.getComplexName()+"_"+chainId+"_"+aa.getresidueID()+"_"+aa.getresidueName1()+"\t"+aa.isAsEpitope()+"\t"+aa.getsolventAccessibility()+"\n";
          bw.write(dataline);
          
          
      }
     
      bw.close();
  }
  public String simpleSecStruct(char c){
      String secStruc="";
      switch(c){
          case 'h':
              secStruc ="100";
              break;
          case 's':
              secStruc ="010";
              break;
          case 'c':
              secStruc ="001";
              break;
      }
      return secStruc;
  }
  public void  writeSecondaryStructureToFile(Complex c, char chainId,double threshold) throws IOException {
      
      Chain chain = c.getChainByChianID(chainId);
      
      File fout = new File("src/datatraining/SS/"+"SecondarySructure_"+threshold+".txt");
      FileOutputStream fos = new FileOutputStream(fout);
      BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fos));
      for (int i=0; i<chain.getaminoacidNum();i++){
          Aminoacid aa= chain.getAminoacidbyorder(i);
          boolean exposedStatus = aa.getRsaTien2013()>threshold ;
          boolean exposedStatusOfEpi = aa.isAsEpitope()&& exposedStatus;
          String secondaryStructure="";
          if(exposedStatus||exposedStatusOfEpi){
             
            String aaId = aa.getPdbid()+"_"+aa.getchainID()+"_"+aa.getresidueID()+"_"+aa.letterToNum();
            secondaryStructure = aaId + "," +this.simpleSecStruct(aa.getsecondaryStructure())+","+aa.isAsEpitope();
            bw.write(secondaryStructure +"\n");
        }  
          
          
          
          
          
      }
     
      bw.close();
  }
  public void writeExposedAAAllModeToFile(Complex c, char chainid, String filename)throws IOException {
      File fout = new File(filename+".txt");
      FileOutputStream fos = new FileOutputStream(fout);
      BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fos));
      ArrayList<Aminoacid> aa1 = this.getExposedAA(c, chainid, 1);
      ArrayList<Aminoacid> aa2 = this.getExposedAA(c, chainid, 2);
      ArrayList<Aminoacid> aa3 = this.getExposedAA(c, chainid, 3);
      for (int i=1 ; i<5;i++){
          
      }
  }
  public void writeExposedAAToFile(Complex c, char chainId,String filename, int mode) throws IOException {
    ArrayList<Aminoacid> aa = getExposedAA(c,chainId, mode);
      File fout = new File(filename+".txt");
	FileOutputStream fos = new FileOutputStream(fout);
 
	BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fos));
        for (int i=0;i<aa.size();i++){
            //cetak perdasarkan pilihan mode
          String exposed = aa.get(i).getresidueID()+"\t "+ aa.get(i).getresidueName()+"\t"+aa.get(i).getRsa()+"\n";
          bw.write(exposed);
          //bw.newLine();
      }
	 
	bw.close();
}
  
  public void output(Complex complex)
  {
    for (int i = 0; i < complex.getAntigenChainIDs().length(); i++) {
      for (int j = 0; j < complex.getChainByChianID(complex.getAntigenChainIDs().charAt(i)).getaminoacidNum(); j++)
      {
        Aminoacid aminoacid = complex.getChainByChianID(complex.getAntigenChainIDs().charAt(i)).getAminoacidbyorder(j);
        System.out.println(complex.getcomplexName() + " " + aminoacid.getchainID() + " " + aminoacid.getresidueID() + " " + aminoacid.getresidueName() + " " + aminoacid.getBinding() + " " + aminoacid.getsecondaryStructure() + " " + aminoacid.getsolventAccessibility() +" "+aminoacid.getRsa()+ " " + aminoacid.getAtombyOrder(0).getbFactor());
      }
    }
  }
}
