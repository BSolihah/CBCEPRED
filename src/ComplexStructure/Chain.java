package ComplexStructure;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;


public class Chain
{
  private char chainID;
  private ArrayList<Aminoacid> aminoacidVector;
  private int aminoacidNum = 0;
  private String type = null;
  private String pdbid;
  private int numExposed;
  private int numEpi;
  
  public Chain()
  {
    this.aminoacidVector = null;
  }

    public int getNumExposed() {
        return numExposed;
    }

    public int getNumEpi() {
        return numEpi;
    }
  
  public String getSecondaryStructure(){
    String sequence="";
    for(int i=0;i<aminoacidVector.size();i++){
          sequence = sequence + aminoacidVector.get(i).getsecondaryStructure();
      }
      return sequence;
  }
  public String getSequenceOfChain(){
      String sequence="";
      for(int i=0;i<aminoacidVector.size();i++){
          sequence = sequence + aminoacidVector.get(i).getresidueName1();
      }
      return sequence;
  }
  public String getExposesResidue(){
      String sequence="";
      
      for(int i=0;i<aminoacidVector.size();i++){
          if (aminoacidVector.get(i).getRsa()>0.05){
              sequence = sequence + "1";
          }else{
              sequence = sequence + "0";
          }
      }
      return sequence;
  }
  public String getAsabasedExposesResidue(double th){
      this.numExposed =0;
      String sequence="";
      for(int i=0;i<aminoacidVector.size();i++){
          if (aminoacidVector.get(i).getsolventAccessibility()>th){
              sequence = sequence + "1";
              this.numExposed +=1;
          }else{
              sequence = sequence + "0";
          }
      }
      return sequence;
  }
  public String getExposedResidue(double th){
      String sequence="";
      this.numExposed =0;
      for(int i=0;i<aminoacidVector.size();i++){
          if (aminoacidVector.get(i).getRsaTien2013()>th){
              sequence = sequence + "1";
              this.numExposed +=1;
          }else{
              sequence = sequence + "0";
          }
      }
      return sequence;
  }
  public String getSequenceOfEpitopeStatus(){
      String sequence="";
      this.numEpi=0;
      for(int i=0;i<aminoacidVector.size();i++){
          if (aminoacidVector.get(i).isAsEpitope()==true){
              sequence = sequence + "1";
              this.numEpi+=1;
          }else{
              sequence = sequence + "0";
          }
      }
      return sequence;
  }
  public String getDatasetFromExposedResidu(){
  String smiller=    "miller     : ";
  String srose=      "rose       : ";
  String srostsender="rose sender: ";
  String stien=      "tien       : ";
  
  for(int i=0;i<aminoacidVector.size();i++){
      if(aminoacidVector.get(i).getRsaMiller1987()>0.01){
          smiller = smiller + aminoacidVector.get(i).getresidueName1();
      }else{
              smiller = smiller + "-";
          }
      if(aminoacidVector.get(i).getRsaRose1985()>0.01){
          srose = srose + aminoacidVector.get(i).getresidueName1();
      }else{
          srose = srose + "-";
      }
      if(aminoacidVector.get(i).getRsaRostSander2005()>0.01){
          srostsender = srostsender + aminoacidVector.get(i).getresidueName1();
      }else{
        srostsender = srostsender +"-";
      }
      if(aminoacidVector.get(i).getRsaTien2013()>0.01){
          stien = stien + aminoacidVector.get(i).getresidueName1();
      }else{
          stien = stien +"-";
      }
      
  }
  String sAll = this.getPdbid()+"_"+this.getChainID()+"\n"+smiller+"\n"+srose+"\n"+srostsender+"\n"+stien+"\n";
  return sAll;
  
  }
  public String getSequenceOfChainAsEpitopeOrExposedResidue(){
      String sequence="";
      for(int i=0;i<aminoacidVector.size();i++){
          if(aminoacidVector.get(i).isAsEpitope()){
              sequence = sequence + aminoacidVector.get(i).getresidueName1();
          }else if (aminoacidVector.get(i).getRsa()>0.05){
              sequence = sequence + aminoacidVector.get(i).getresidueName1();
          }else{
              sequence = sequence + "-";
          }
          
      }
      return sequence;
  }
  public ArrayList<Aminoacid> getAAAsEpitope(){
      ArrayList<Aminoacid> epilist = new ArrayList();
      for(int i=0;i<aminoacidVector.size();i++){
          if (aminoacidVector.get(i).isAsEpitope()){
              epilist.add(aminoacidVector.get(i));
          }
      }
      return epilist;
  }
  public void setAAAsEpitope(Epitope e){
      if(aminoacidVector.size()>0){
          for(int i=0;i<e.getEpipos().length();i++ ){
          if(aminoacidVector.get(i).getresidueName1().equalsIgnoreCase(e.getSequence().substring(i, i+1))){
              if(e.getEpipos().substring(i, i+1).equalsIgnoreCase("1")){
                  aminoacidVector.get(i).setAsEpitope(true);
              }
          }
          }
          
      }
      
  }
  public void setAAAsEpitope(ArrayList<Integer> idx){
      for(int i=0;i<idx.size();i++){
          int index = idx.get(i).intValue();
         // System.out.println(index);
          aminoacidVector.get(index).setAsEpitope(true);
      }
  }
  public Chain(ArrayList<Aminoacid> aminoacidVector)
  {
    this.aminoacidVector = aminoacidVector;
    this.aminoacidNum = aminoacidVector.size();
    this.chainID = ((Aminoacid)aminoacidVector.get(0)).getchainID();
    this.pdbid = ((Aminoacid)aminoacidVector.get(0)).getPdbid();
  }
 
  public void setEpitopeStatus(Epitope e, String chainID){
      
      int idxlength = e.getEpipos().length();
      int numAA = this.getaminoacidNum();
      //System.out.println("set epi"+this.pdbid+": "+idxlength +"; "+numAA);
      for(int i=0;i<idxlength;i++){
          boolean epi = Boolean.valueOf(e.getEpipos().substring(i,i+1 )).booleanValue();
          System.out.println("aa: "+this.aminoacidVector.get(i).toString()+", status epi: "+ epi);
          this.aminoacidVector.get(i).setAsEpitope(epi);
      }
     
      
  }
  public static String aainOneLetter(String letter)
  {
    if (letter.equals("ALA")) {
      return "A";
    }
    if (letter.equals("ARG")) {
      return "R";
    }
    if (letter.equals("ASN")) {
      return "N";
    }
    if (letter.equals("ASP")) {
      return "D";
    }
    if (letter.equals("CYS")) {
      return "C";
    }
    if (letter.equals("GLN")) {
      return "Q";
    }
    if (letter.equals("GLU")) {
      return "E";
    }
    if (letter.equals("GLY")) {
      return "G";
    }
    if (letter.equals("HIS")) {
      return "H";
    }
    if (letter.equals("ILE")) {
      return "I";
    }
    if (letter.equals("LEU")) {
      return "L";
    }
    if (letter.equals("LYS")) {
      return "K";
    }
    if (letter.equals("MET")) {
      return "M";
    }
    if (letter.equals("PHE")) {
      return "F";
    }
    if (letter.equals("PRO")) {
      return "P";
    }
    if (letter.equals("SER")) {
      return "S";
    }
    if (letter.equals("THR")) {
      return "T";
    }
    if (letter.equals("TRP")) {
      return "W";
    }
    if (letter.equals("TYR")) {
      return "Y";
    }
    if (letter.equals("VAL")) {
      return "V";
    }
    System.out.println("abnormal value");
    return "X";
  }
    
  
  public void getAASequenceOfChain(HashMap<String,String> hashseq, HashMap<Integer,String> hashIdx){
      int count=0;
      for(Aminoacid aa: aminoacidVector){
          hashseq.put(aa.getresidueID(), aa.getresidueName1());
          hashIdx.put(count,aa.getresidueID());
          count++;
      }
      
  }
  public String getAminoacidSequenceofChain(){
      String s = "";
      for (int i=0;i<aminoacidVector.size();i++)
          s = s+ this.aainOneLetter(this.aminoacidVector.get(i).getresidueName());
      return s;
  }
  public int getaminoacidNum()
  {
    return this.aminoacidVector.size();
  }
  
  public ArrayList<Aminoacid> getAminoacidVector()
  {
    return this.aminoacidVector;
  }
  
  public Aminoacid getAminoacidbyorder(int n)
  {
   // return (Aminoacid)this.aminoacidVector.elementAt(n);
    return (Aminoacid)this.aminoacidVector.get(n);
  }
  
  public Aminoacid getAminoacidbyIDs(String residueID)
  {
    int n = this.aminoacidVector.size();
    for (int i = 0; i < n; i++) {
      if (((Aminoacid)this.aminoacidVector.get(i)).getresidueID().equals(residueID)) {
        return (Aminoacid)this.aminoacidVector.get(i);
      }
    }
    System.out.println("no such residue"+ residueID);
    return null;
  }
  public ArrayList<Integer> getIdxAminoacidbyIDs(ArrayList<String> residueID)
  {
     ArrayList<Integer> idx=new ArrayList();
    int n = this.aminoacidVector.size();
    for(int j=0;j<residueID.size();j++){
        for (int i = 0; i < n; i++) {
            if (((Aminoacid)this.aminoacidVector.get(i)).getresidueID().equals(residueID)) {
                idx.add(i);          
            }
        }
    }
    
    
    
    return idx;
  }
  
  public void setChain(ArrayList<Aminoacid> aminoacidVector)
  {
    this.aminoacidVector = ((ArrayList)aminoacidVector.clone());
    this.aminoacidNum = aminoacidVector.size();
  }
  
  public char getChainID()
  {
    return this.chainID;
  }

    public String getPdbid() {
        return pdbid;
    }

    public void setPdbid(String pdbid) {
        this.pdbid = pdbid;
    }
  
  public void setchainID(char ChainID)
  {
    this.chainID = ChainID;
    for (int i = 0; i < this.aminoacidVector.size(); i++) {
      ((Aminoacid)this.aminoacidVector.get(i)).setChainID(ChainID);
    }
  }
  public int[] getDistributionOfAA(){
      int[] dist = new int [20];
      for(int i=0;i<20;i++)
          dist[i]=0;
      int num= this.getAminoacidVector().size();
      for(int i=0;i<num;i++){
          int idx = this.getAminoacidVector().get(i).letterToNum();
          dist[idx]= dist[idx]+1;
      }
      return dist;
  }
  public void deleteAminoacidByOrder(int n)
  {
    this.aminoacidVector.remove(n);
  }
}
