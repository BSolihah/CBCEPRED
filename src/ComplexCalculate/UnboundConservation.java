package ComplexCalculate;

import ComplexStructure.Aminoacid;
import ComplexStructure.Chain;
import ComplexStructure.Complex;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;


public class UnboundConservation
{
  private double[] blosum62;
  public UnboundConservation(Complex unboundComplex){
      this.blosum62 = new double[] { 4.0D, 5.0D, 6.0D, 6.0D, 9.0D, 5.0D, 5.0D, 6.0D, 8.0D, 4.0D, 4.0D, 5.0D, 5.0D, 6.0D, 7.0D, 4.0D, 5.0D, 11.0D, 7.0D, 4.0D };
      //analysisFile(unboundComplex);
      calculateConservation(unboundComplex);
  }
  
  public void calculateConservation(Complex unboundComplex){
      //definisikan blosum62 : komponen diagonal dari blosum62
      this.blosum62 = new double[] { 4.0D, 5.0D, 6.0D, 6.0D, 9.0D, 5.0D, 5.0D, 6.0D, 8.0D, 4.0D, 4.0D, 5.0D, 5.0D, 6.0D, 7.0D, 4.0D, 5.0D, 11.0D, 7.0D, 4.0D };
      
      String  antigenchain = unboundComplex.getAntigenChainIDs();
      Chain chain =unboundComplex.getChainByChianID(antigenchain.charAt(0));
      
      double[][] PSSM;
      PSSM = getPSSM(unboundComplex, chain.getChainID());
      //untuk setiap chain, hitung skor conservation
          for (int j=0;j<chain.getaminoacidNum();j++){
              //ambil pssm
              //ambil brr
              String rn = chain.getAminoacidbyorder(j).getresidueName();
              int r = letterTonum(chain.getAminoacidbyorder(j).getresidueName());
              double br = this.blosum62[r] ;
              double pssm = PSSM[j][r];
              double conserv;
              if(pssm<br){
                  conserv = br-pssm;
              }else
                  conserv =0;
              //set conserv
              System.out.print(rn+" "+r +":"+br + ","+ pssm+":" +conserv +"\t");
              unboundComplex.getChainByChianID(antigenchain.charAt(0)).getAminoacidbyorder(j).setConservation(conserv);
          }
             
      }
               
  

  public UnboundConservation(ArrayList<Complex> unboundcomplexVector)
  {
    System.out.println("***********beign to get complex constructures' conservation**********");
    this.blosum62 = new double[] { 4.0D, 5.0D, 6.0D, 6.0D, 9.0D, 5.0D, 5.0D, 6.0D, 8.0D, 4.0D, 4.0D, 5.0D, 5.0D, 6.0D, 7.0D, 4.0D, 5.0D, 11.0D, 7.0D, 4.0D };
    for (int i = 0; i < unboundcomplexVector.size(); i++) {
      analysisFile((Complex)unboundcomplexVector.get(i));
    }
    System.out.println("***********finish to get complex constructures' conservation**********");
  }
  
  private void analysisFile(Complex unboundcomplex)
  {
    for (int i = 0; i < unboundcomplex.getAntigenChainIDs().length(); i++)
    {
      //Chain chain = unboundcomplex.getChainByChianID(unboundcomplex.getAntibodyChainIDs().charAt(i));
        Chain chain = unboundcomplex.getChainByChianID(unboundcomplex.getAntigenChainIDs().charAt(i));
      double[][] PSSM = new double[unboundcomplex.getChainByChianID(unboundcomplex.getAntigenChainIDs().charAt(i)).getaminoacidNum()][20];
      PSSM = getPSSM(unboundcomplex, chain.getChainID());
      this.printPSSM(PSSM);
      
      for (int j = 0; j < chain.getaminoacidNum(); j++)
      {
          for (int k=0;k< PSSM.length;k++){
              double conservation;
              if (PSSM[j][letterTonum(chain.getAminoacidbyorder(k).getresidueName())] < this.blosum62[letterTonum(chain.getAminoacidbyorder(k).getresidueName())]) {
              String rn = chain.getAminoacidbyorder(k).getresidueName();
              int r = letterTonum(chain.getAminoacidbyorder(k).getresidueName());
              double blosumidx = this.blosum62[letterTonum(chain.getAminoacidbyorder(k).getresidueName())] ;
              double pssm = PSSM[j][r];
              conservation = this.blosum62[letterTonum(chain.getAminoacidbyorder(k).getresidueName())] - PSSM[j][letterTonum(chain.getAminoacidbyorder(k).getresidueName())];
          System.out.print ("j: "+ j+" r: "+ r+"\t");
          System.out.print(rn+" "+r +":"+blosumidx + ","+ pssm+":" +conservation +"\t");
        } else {
          conservation = 0.0D;
        }
        chain.getAminoacidbyorder(j).setConservation(conservation);
        
      }
          }
        
    }
  }
  
  private double[][] getPSSM(Complex unboundcomplex, char chainID)
  {
    Chain chain = unboundcomplex.getChainByChianID(chainID);
    double[][] PSSM = new double[chain.getaminoacidNum()][20];
    try
    {
      FileReader fr = new FileReader("src/pssm/" + unboundcomplex.getcomplexName() + "_" + chainID + ".PSSM");
      BufferedReader br = new BufferedReader(fr);
      String dataline = br.readLine();
      dataline = br.readLine();
      dataline = br.readLine();
      dataline = br.readLine();
      for (int i = 0; i < chain.getaminoacidNum(); i++)
      {
        String[] temp = dataline.trim().split("\\s+");
        for (int j = 0; j < 20; j++) {
          PSSM[i][j] = Double.valueOf(temp[(j + 2)]).doubleValue();
        }
        dataline = br.readLine();
      }
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
    return PSSM;
  }
  
  private int letterTonum(String letter)
  {
    if ((letter.equals("ALA")) || (letter.equals("A"))) {
      return 0;
    }
    if ((letter.equals("ARG")) || (letter.equals("R"))) {
      return 1;
    }
    if ((letter.equals("ASN")) || (letter.equals("N"))) {
      return 2;
    }
    if ((letter.equals("ASP")) || (letter.equals("D"))) {
      return 3;
    }
    if ((letter.equals("CYS")) || (letter.equals("C"))) {
      return 4;
    }
    if ((letter.equals("GLN")) || (letter.equals("Q"))) {
      return 5;
    }
    if ((letter.equals("GLU")) || (letter.equals("E"))) {
      return 6;
    }
    if ((letter.equals("GLY")) || (letter.equals("G"))) {
      return 7;
    }
    if ((letter.equals("HIS")) || (letter.equals("H"))) {
      return 8;
    }
    if ((letter.equals("ILE")) || (letter.equals("I"))) {
      return 9;
    }
    if ((letter.equals("LEU")) || (letter.equals("L"))) {
      return 10;
    }
    if ((letter.equals("LYS")) || (letter.equals("K"))) {
      return 11;
    }
    if ((letter.equals("MET")) || (letter.equals("M"))) {
      return 12;
    }
    if ((letter.equals("PHE")) || (letter.equals("F"))) {
      return 13;
    }
    if ((letter.equals("PRO")) || (letter.equals("P"))) {
      return 14;
    }
    if ((letter.equals("SER")) || (letter.equals("S"))) {
      return 15;
    }
    if ((letter.equals("THR")) || (letter.equals("T"))) {
      return 16;
    }
    if ((letter.equals("TRP")) || (letter.equals("W"))) {
      return 17;
    }
    if ((letter.equals("TYR")) || (letter.equals("Y"))) {
      return 18;
    }
    if ((letter.equals("VAL")) || (letter.equals("V"))) {
      return 19;
    }
    return -1;
  }
  public void printPSSM(double[][] pssm){
      for (int i=0; i<pssm.length;i++){
          for (int j=0;j<pssm[0].length;j++){
                  System.out.print(pssm[i][j]+"\t");
          }
          System.out.print("\n");
      }
          
  
  }
  
  
}
