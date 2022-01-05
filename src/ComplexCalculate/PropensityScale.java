/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ComplexCalculate;

import ComplexStructure.Aminoacid;
import java.util.ArrayList;

/**
 *
 * @author OpenGress
 */
public class PropensityScale {
    private double [] AndersonLogOddArr;
    private double []ParkerHArr;
    
    public PropensityScale(){
    setAndersonLogOddArr();
    setParkerHArr();
    }
    //untuk menghemat memory sebaiknya hanya diset log ogg dari asam amino terekspose
    public void setAndersonLogOddOfAA(ArrayList<Aminoacid> AA){
         for(int i=0;i<AA.size();i++){
               
               AA.get(i).setAndersonlogodd(this.getAndersonLogOdd(AA.get(i) .getresidueName()));
            }
     
    }
    public void setParkerHofAA(ArrayList<Aminoacid> AA){
        for (int i=0;i<AA.size();i++){
             AA.get(i).setParkerh(this.getParkerH(AA.get(i) .getresidueName()));
        }
    }
        
        
    
    private void setParkerHArr(){
    ParkerHArr = new double[20];
    ParkerHArr[1]= 0.870; 
    ParkerHArr[2]= 1.640 ;
    ParkerHArr[3]=  2.460; 
    ParkerHArr[4]= 0.110; 
    ParkerHArr[5]= 1.370 ;
    ParkerHArr[6]=  1.860 ;
    ParkerHArr[7]= 1.280;
    ParkerHArr[8]= 0.300;
    ParkerHArr[9]= -2.450 ;
    ParkerHArr[10]= -2.870 ; 
    ParkerHArr[11]= 1.260; 
    ParkerHArr[12]= -1.410; 
    ParkerHArr[13]= -2.780 ;
    ParkerHArr[14]= 0.300; 
    ParkerHArr[15]= 1.500 ;
    ParkerHArr[16]= 1.150 ;
    ParkerHArr[17]= -3.000;
    ParkerHArr[18]= -0.780 ;
    ParkerHArr[19]= -1.270 ;
    }
    
    private void setAndersonLogOddArr(){
        AndersonLogOddArr = new double[20];
        AndersonLogOddArr[0]= -1.522;

        AndersonLogOddArr[1]= 1.180;
        AndersonLogOddArr[2]= 1.242;
        AndersonLogOddArr[3]= 0.691;
        AndersonLogOddArr[4]= -3.519;
        AndersonLogOddArr[5]= 1.082;
        AndersonLogOddArr[6]= 0.346;
        AndersonLogOddArr[7]= 0.189;
        AndersonLogOddArr[8]= 1.098;
        AndersonLogOddArr[9]= -0.713;
        AndersonLogOddArr[10]= -1.836;
        AndersonLogOddArr[11]= 1.136;
        AndersonLogOddArr[12]= 0.273;
        AndersonLogOddArr[13]= -1.147;
        AndersonLogOddArr[14]= 1.164 ;
        AndersonLogOddArr[15]= -0.145;
        AndersonLogOddArr[16]= -0.233 ;
        AndersonLogOddArr[17]= -0.064;
        AndersonLogOddArr[18]= 0.030;
        AndersonLogOddArr[19]= -1.474;
        
    }
    public double getAndersonLogOdd(String letter)
  { 
      if ((letter.equals("ALA")) || (letter.equals("A"))) {
      return AndersonLogOddArr[0];
    }
    if ((letter.equals("ARG")) || (letter.equals("R"))) {
      return AndersonLogOddArr[1];
    }
    if ((letter.equals("ASN")) || (letter.equals("N"))) {
      return AndersonLogOddArr[2];
    }
    if ((letter.equals("ASP")) || (letter.equals("D"))) {
      return AndersonLogOddArr[3];
    }
    if ((letter.equals("CYS")) || (letter.equals("C"))) {
      return AndersonLogOddArr[4];
    }
    if ((letter.equals("GLN")) || (letter.equals("Q"))) {
      return AndersonLogOddArr[5];
    }
    if ((letter.equals("GLU")) || (letter.equals("E"))) {
      return AndersonLogOddArr[6];
    }
    if ((letter.equals("GLY")) || (letter.equals("G"))) {
      return AndersonLogOddArr[7];
    }
    if ((letter.equals("HIS")) || (letter.equals("H"))) {
      return AndersonLogOddArr[8];
    }
    if ((letter.equals("ILE")) || (letter.equals("I"))) {
      return AndersonLogOddArr[9];
    }
    if ((letter.equals("LEU")) || (letter.equals("L"))) {
      return AndersonLogOddArr[10];
    }
    if ((letter.equals("LYS")) || (letter.equals("K"))) {
      return AndersonLogOddArr[11];
    }
    if ((letter.equals("MET")) || (letter.equals("M"))) {
      return AndersonLogOddArr[12];
    }
    if ((letter.equals("PHE")) || (letter.equals("F"))) {
      return AndersonLogOddArr[13];
    }
    if ((letter.equals("PRO")) || (letter.equals("P"))) {
      return AndersonLogOddArr[14];
    }
    if ((letter.equals("SER")) || (letter.equals("S"))) {
      return AndersonLogOddArr[15];
    }
    if ((letter.equals("THR")) || (letter.equals("T"))) {
      return AndersonLogOddArr[16];
    }
    if ((letter.equals("TRP")) || (letter.equals("W"))) {
      return AndersonLogOddArr[17];
    }
    if ((letter.equals("TYR")) || (letter.equals("Y"))) {
      return AndersonLogOddArr[18];
    }
    if ((letter.equals("VAL")) || (letter.equals("V"))) {
      return AndersonLogOddArr[19];
    }
    
    
    return 0;
  }



public double getParkerH(String letter){
    
    if ((letter.equals("ALA")) || (letter.equals("A"))) {
      return ParkerHArr[0];
    }
    if ((letter.equals("ARG")) || (letter.equals("R"))) {
      return ParkerHArr[1];
    }
    if ((letter.equals("ASN")) || (letter.equals("N"))) {
      return ParkerHArr[2];
    }
    if ((letter.equals("ASP")) || (letter.equals("D"))) {
      return ParkerHArr[3];
    }
    if ((letter.equals("CYS")) || (letter.equals("C"))) {
      return ParkerHArr[4];
    }
    if ((letter.equals("GLN")) || (letter.equals("Q"))) {
      return ParkerHArr[5];
    }
    if ((letter.equals("GLU")) || (letter.equals("E"))) {
      return ParkerHArr[6];
    }
    if ((letter.equals("GLY")) || (letter.equals("G"))) {
      return ParkerHArr[7];
    }
    if ((letter.equals("HIS")) || (letter.equals("H"))) {
      return ParkerHArr[8];
    }
    if ((letter.equals("ILE")) || (letter.equals("I"))) {
      return ParkerHArr[9];
    }
    if ((letter.equals("LEU")) || (letter.equals("L"))) {
      return ParkerHArr[10];
    }
    if ((letter.equals("LYS")) || (letter.equals("K"))) {
      return ParkerHArr[11];
    }
    if ((letter.equals("MET")) || (letter.equals("M"))) {
      return ParkerHArr[12];
    }
    if ((letter.equals("PHE")) || (letter.equals("F"))) {
      return ParkerHArr[13];
    }
    if ((letter.equals("PRO")) || (letter.equals("P"))) {
      return ParkerHArr[14];
    }
    if ((letter.equals("SER")) || (letter.equals("S"))) {
      return ParkerHArr[15];
    }
    if ((letter.equals("THR")) || (letter.equals("T"))) {
      return ParkerHArr[16];
    }
    if ((letter.equals("TRP")) || (letter.equals("W"))) {
      return ParkerHArr[17];
    }
    if ((letter.equals("TYR")) || (letter.equals("Y"))) {
      return ParkerHArr[18];
    }
    if ((letter.equals("VAL")) || (letter.equals("V"))) {
      return ParkerHArr[19];
    }
    return 0.0;
}
}