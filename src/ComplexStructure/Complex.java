package ComplexStructure;

import ExtractDataTraining.ExtractDataBenchmark;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.Vector;

public class Complex
{
  private String complexName;
  public static int neigborNum = 20;
  private ArrayList<Chain> chainVector;
  private int chainNum;
  private String chainIDs = "";
  private String antigenChainIDs;
  private String antibodyChainIDs;
  private ArrayList<Atom> atomVector;
  private Boolean rsaisset=false;
  private Boolean asepitopeisset = false;
  private ArrayList<Aminoacid> aminoacidVector;
  
  
  public Complex(String complexName, ArrayList<Chain> ChainVector)
  {
      
    this.chainVector = new ArrayList();
    this.complexName = complexName;
    this.chainNum = ChainVector.size();
    char[] chainID = new char[this.chainNum];
    for (int i = 0; i < this.chainNum; i++) {
      chainID[i] = ((Chain)ChainVector.get(i)).getChainID();
    }
    this.chainIDs = new String(chainID);
    for (int i = 0; i < this.chainNum; i++) {
      this.chainVector.add((Chain)ChainVector.get(i));
    }
    this.antigenChainIDs = "";
    this.antibodyChainIDs = "";
  }

    public String getComplexName() {
        return complexName;
    }
    
    public void setComplexName(String complexName) {
        this.complexName = complexName;
    }
   public Aminoacid getAminoById(String id){
        Aminoacid aa=null;
        for(Aminoacid a: aminoacidVector){
            if (a.getresidueID().equals(id)){
               // System.out.println("resid: "+a.getresidueID());
                return a;
            }           
            
        }
        return aa;
    }
  public void setComplex(ArrayList<Chain> ChainVector)
  {
    this.chainVector = ((ArrayList)ChainVector.clone());
    this.chainNum = this.chainVector.size();
    char[] chainID = new char[this.chainNum];
    for (int i = 0; i < this.chainNum; i++) {
      chainID[i] = ((Chain)this.chainVector.get(i)).getChainID();
    }
    this.chainIDs = new String(chainID);
  }
 
  
  public Complex(String path, String pdbid, String chain){
      this.complexName = pdbid;
    this.chainIDs = chain;
    this.antigenChainIDs = chain;
    this.antibodyChainIDs = "not set";
    
    this.chainVector = new ArrayList();
    this.atomVector = new ArrayList();
    this.aminoacidVector= new ArrayList();
    
    String dataline = null;
    FileReader fr = null;
    //String filename = "src/rescbc/"+pdb_file+"_"+selectedChain+".pdb";
     String filename = path +pdbid+"_"+chain;
    System.out.println("path pdbid:"+ filename);
    File file = new File(filename);
    /*
    if (!file.exists()) {
      file = new File("/src/data/" + pdb_file + ".pdb");
    }*/
    try
    {
      fr = new FileReader(file);
      BufferedReader br = new BufferedReader(fr);
      dataline = br.readLine();
      while (dataline != null)
      {
        if (dataline.startsWith("ATOM"))
        {
          Atom atom = new Atom();
          String atomName = new String(dataline.substring(12, 16).trim());
          String residueName = new String(dataline.substring(17, 20).trim());
          char chainName;
          
          if (dataline.substring(21, 22).trim().equals("")) {
            chainName = 'X';
          } else {
            chainName = dataline.substring(21, 22).trim().charAt(0);
          }
          String residueID = new String(dataline.substring(22, 27).trim());
          //System.out.println(residueID);
          double x = Double.parseDouble(dataline.substring(30, 38).trim());
          double y = Double.parseDouble(dataline.substring(38, 46).trim());
          double z = Double.parseDouble(dataline.substring(46, 54).trim());
          atom.setAtomName(atomName);
          atom.setChainName(chainName);
          if (dataline.length() > 60)
          {
            double bFactor = Double.valueOf(dataline.substring(60, 66).trim()).doubleValue();
           // System.out.println(bFactor);
            atom.setbFactor(bFactor);
          }
          else
          {
            int state = (int)Double.parseDouble(dataline.substring(55, dataline.length()));
            atom.setType(state);
          }
          atom.setResidueName(residueName);
          atom.setResidueID(residueID);
          
          atom.setX(x);
          atom.setY(y);
          atom.setZ(z);
          
          this.atomVector.add(atom);
        }
        dataline = br.readLine();
      }
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
    anlysisPDB();
  }
  public Complex(String pdb_file, String selectedChain)
  {
    this.complexName = pdb_file;
    this.chainIDs = selectedChain;
    this.antigenChainIDs = selectedChain;
    this.antibodyChainIDs = "not set";
    
    this.chainVector = new ArrayList();
    this.atomVector = new ArrayList();
    this.aminoacidVector= new ArrayList();
    
    String dataline = null;
    FileReader fr = null;
    //String filename = "src/rescbc/"+pdb_file+"_"+selectedChain+".pdb";
    // String filename = "src/dataset/andersen/chain/"+pdb_file+"_"+selectedChain+".pdb";
     String filename = "src/dataset/seppa3/test/chain/"+pdb_file+"_"+selectedChain+".pdb";
    //System.out.println(filename);
    File file = new File(filename);
    /*
    if (!file.exists()) {
      file = new File("/src/data/" + pdb_file + ".pdb");
    }*/
    try
    {
      fr = new FileReader(file);
      BufferedReader br = new BufferedReader(fr);
      dataline = br.readLine();
      while (dataline != null)
      {
        if (dataline.startsWith("ATOM"))
        {
          Atom atom = new Atom();
          String atomName = new String(dataline.substring(12, 16).trim());
          String residueName = new String(dataline.substring(17, 20).trim());
          char chainName;
          
          if (dataline.substring(21, 22).trim().equals("")) {
            chainName = 'X';
          } else {
            chainName = dataline.substring(21, 22).trim().charAt(0);
          }
          String residueID = new String(dataline.substring(22, 27).trim());
          //System.out.println(residueID);
          double x = Double.parseDouble(dataline.substring(30, 38).trim());
          double y = Double.parseDouble(dataline.substring(38, 46).trim());
          double z = Double.parseDouble(dataline.substring(46, 54).trim());
          atom.setAtomName(atomName);
          atom.setChainName(chainName);
          if (dataline.length() > 60)
          {
            double bFactor = Double.valueOf(dataline.substring(60, 66).trim()).doubleValue();
           // System.out.println(bFactor);
            atom.setbFactor(bFactor);
          }
          else
          {
            int state = (int)Double.parseDouble(dataline.substring(55, dataline.length()));
            atom.setType(state);
          }
          atom.setResidueName(residueName);
          atom.setResidueID(residueID);
          
          atom.setX(x);
          atom.setY(y);
          atom.setZ(z);
          
          this.atomVector.add(atom);
        }
        dataline = br.readLine();
      }
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
    anlysisPDB();
  }

  
    public Boolean getAsepitopeisset() {
        return asepitopeisset;
    }

    public void setAsepitopeisset(Boolean asepitopeisset) {
        this.asepitopeisset = asepitopeisset;
    }

    public Boolean getRsaisset() {
        return rsaisset;
    }

    public void setRsaisset(Boolean rsaisset) {
        this.rsaisset = rsaisset;
    }
  
  public String printEpiStatusOfExposedResidue(){
      
      String chainid = this.getAntigenChainIDs();
      Chain chain = this.getChainByChianID(chainid.charAt(0));
      String s= chain.getSequenceOfChain();
      System.out.println(s+";");
      s.concat("\n"+chain.getExposesResidue());
      System.out.println(chain.getExposesResidue()+";");
      s.concat("\n"+ chain.getSequenceOfEpitopeStatus()+";");
      System.out.println("\n"+chain.getSequenceOfEpitopeStatus()+";");
      return s;
      
  }

    public ArrayList<Aminoacid> getAminoacidVector() {
        return aminoacidVector;
    }
     //return value daftar residu dalam format num
  
  
   //return value daftar residu dalam format num
  public ArrayList<Integer> getNeighborhoodAA(Aminoacid a,int oddwinsize){
//      int residuId = Integer.valueOf(a.getresidueID());
      int c = (oddwinsize-1)/2;
      int idxaa = aminoacidVector.indexOf(a);
     
      //System.out.println("residuId: "+ residuId + "idxaa: "+ idxaa);
      ArrayList<Integer> listnaa = new ArrayList();
      int i=-c;
      while(i<=c){
          if(i+idxaa >=0&& (i+idxaa)<aminoacidVector.size()){
          listnaa.add(aminoacidVector.get(i+idxaa).letterToNum());
          }
          i++;
      }
     return listnaa;
      
  }
  
  public ArrayList<Aminoacid> getListNeighborhoodAA(Aminoacid a,int oddwinsize){
//      int residuId = Integer.valueOf(a.getresidueID());
      int c = (oddwinsize-1)/2;
      int idxaa = aminoacidVector.indexOf(a);
     
      //System.out.println("residuId: "+ residuId + "idxaa: "+ idxaa);
      ArrayList<Aminoacid> listnaa = new ArrayList();
      int i=-c;
      while(i<=c){
          if(i+idxaa >=0&& (i+idxaa)<aminoacidVector.size()){
          listnaa.add(aminoacidVector.get(i+idxaa));
          }
          i++;
      }
     return listnaa;
      
  }
  public ArrayList<Integer> getListIdxNeighborhoodAA(Aminoacid a,int oddwinsize){
//      int residuId = Integer.valueOf(a.getresidueID());
      int c = (oddwinsize-1)/2;
      int idxaa = aminoacidVector.indexOf(a);
     
      //System.out.println("residuId: "+ residuId + "idxaa: "+ idxaa);
      ArrayList<Integer> listnaa = new ArrayList();
      int i=-c;
      while(i<=c){
          if(i+idxaa >=0&& (i+idxaa)<aminoacidVector.size()){
          listnaa.add(i+idxaa);
          }
          i++;
      }
     return listnaa;
      
  }
  private void anlysisPDB()
  {
    
    int len = this.atomVector.size();
    Atom firstAtom = new Atom();
    Atom nextAtom = new Atom();
    firstAtom = (Atom)this.atomVector.get(0);
    
//    Vector<Aminoacid> aminoacidVector = new Vector(100, 100);
   //aminoacidVector = new ArrayList();
    
   // Vector<Atom> atomVector_temp = new Vector(100, 100);
    
    ArrayList <Atom> atomVector_temp = new ArrayList();
    atomVector_temp.add(firstAtom);
    for (int i = 1; i < len; i++)
    {
      nextAtom = (Atom)this.atomVector.get(i);
      if ((firstAtom.getChainName() == nextAtom.getChainName()) && (firstAtom.getResidueID().equals(nextAtom.getResidueID())))
      {
        atomVector_temp.add(nextAtom);
      }
      else if ((firstAtom.getChainName()== nextAtom.getChainName()) && (!firstAtom.getResidueID().equals(nextAtom.getResidueID())))
      {
        Aminoacid aminoacid = new Aminoacid(atomVector_temp);
        aminoacid.setPdbid(this.complexName);
        
        if (aminoacid.getatomNum() > 1) {
          aminoacidVector.add(aminoacid);
        }
        atomVector_temp = new ArrayList();
        firstAtom = nextAtom;
        atomVector_temp.add(firstAtom);
      }
      else if (firstAtom.getChainName() != nextAtom.getChainName())
      {
        Aminoacid aminoacid = new Aminoacid(atomVector_temp);
        if (aminoacid.getatomNum() > 1) {
          aminoacidVector.add(aminoacid);
        }
        Chain chain = new Chain(aminoacidVector);
        this.chainVector.add(chain);
        aminoacidVector = new ArrayList();
        atomVector_temp = new ArrayList();
        firstAtom = nextAtom;
        atomVector_temp.add(firstAtom);
      }
    }
    Aminoacid aminoacid = new Aminoacid(atomVector_temp);
    if (aminoacid.getatomNum() > 1) {
        aminoacid.setPdbid(this.complexName);
        aminoacidVector.add(aminoacid);
    }
    Chain chain = new Chain(aminoacidVector);
    chain.setPdbid(this.complexName);
    this.chainVector.add(chain);
    
    this.chainNum = this.chainVector.size();
    
    char[] chainID = new char[this.chainNum];
    for (int i = 0; i < this.chainNum; i++) {
      //chainID[i] = ((Chain)this.chainVector.elementAt(i)).getChainID();
      chainID[i] = ((Chain) this.chainVector.get(i)).getChainID();
    }
    this.chainIDs = new String(chainID);
  }
  
    public int getChainNum() {
        return chainNum;
    }

    public void setChainNum(int chainNum) {
        this.chainNum = chainNum;
    }

    public ArrayList<Chain> getChainVector() {
        return chainVector;
    }

    public void setChainVector(ArrayList<Chain> chainVector) {
        this.chainVector = chainVector;
    }
  
  public Chain getChainByChianID(char chainID)
  {
    for (int i = 0; i < this.chainVector.size(); i++) {
      if (chainID == ((Chain)this.chainVector.get(i)).getChainID()) {
        return (Chain)this.chainVector.get(i);
      }
    }
    System.out.println("no scuh chain");
    return null;
  }
  public ArrayList<Aminoacid> getAminoAcidList(ArrayList<Integer> idxN, char chainID){
      ArrayList<Aminoacid> listAa = null;
      for (int i = 0; i < this.chainVector.size(); i++) {
      if (chainID == ((Chain)this.chainVector.get(i)).getChainID()) {
          ArrayList<Aminoacid> vectorAA= this.chainVector.get(i).getAminoacidVector();
          listAa = new ArrayList();
          for(int j=0;j<idxN.size();j++){
              listAa.add(vectorAA.get(idxN.get(j)));
          }
          
      }
    }
    return listAa;
  }
  
  public int getAminoacidNumByChainID(char chainID)
  {
    for (int i = 0; i < this.chainVector.size(); i++) {
      if (chainID == ((Chain)this.chainVector.get(i)).getChainID()) {
        return ((Chain)this.chainVector.get(i)).getaminoacidNum();
      }
    }
    System.out.println("no scuh chain");
    return -1;
  }
  
  public String getChainIDs()
  {
    return this.chainIDs;
  }
  
  public String getAntigenChainIDs()
  {
    return this.antigenChainIDs;
  }
  
  public String getAntibodyChainIDs()
  {
    return this.antibodyChainIDs;
  }
  
  public String getcomplexName()
  {
    return this.complexName;
  }
  
  public void setchainID(char[] ChainID)
  {
    this.chainIDs = new String(ChainID);
  }
  public void printConserv(){
      System.out.print("\n");
      for (int i=0; i< this.getChainByChianID(("C").charAt(0)).getaminoacidNum();i++){
          System.out.print(this.getChainByChianID(("C").charAt(0)).getAminoacidbyorder(i).getresidueName()+": "+getChainByChianID(("C").charAt(0)).getAminoacidbyorder(i).getConservation()+"\n");
      }
  }
  public ArrayList<Atom> getCAlphaOfComplex(){
      ArrayList listca= new ArrayList();
          ArrayList <Chain> cv = this.getChainVector();
          for (int i= 0; i<cv.size();i++){
              for(int j=0;j<cv.get(i).getaminoacidNum();j++){
                  
                 listca.add((Atom)cv.get(i).getAminoacidbyorder(j).getCAlpha());
              }
              
          }
          return listca;
  }
  //Calpha, Cbeta, N, of residu
  public ArrayList<AtomCoordinate> getFrameReference(Aminoacid aa){
      ArrayList <AtomCoordinate> frameRef=new ArrayList();
      Atom cbeta = aa.getAtomofAA("CB");
      if(cbeta==null){
          //buat cbetapseudo
          // rotasikan N sebesar -120 terhadap sumbu CAC
          
      }
      frameRef.add(new AtomCoordinate(aa.getAtomCoordinateofAA("CA").getX(), aa.getAtomCoordinateofAA("CA").getY(),aa.getAtomCoordinateofAA("CA").getZ()));
      frameRef.add(new AtomCoordinate(aa.getAtomCoordinateofAA("N").getX(), aa.getAtomCoordinateofAA("N").getY(),aa.getAtomCoordinateofAA("N").getZ()));
      //frameRef.add(cbeta);
      return frameRef;
  }
  
  public void printKoordCA(){
      ArrayList<Atom> lca = this.getCAlphaOfComplex();
      for(int i=0;i<lca.size();i++){
          System.out.print(lca.get(i).getX()+","+ lca.get(i).getY()+", "+lca.get(i).getZ()+"\n");
      }
  }

    public void setEpitopeState(ArrayList<Double> deltaASA) {
        for(int i=0;i<this.aminoacidVector.size();i++){
            Boolean state = deltaASA.get(i)>5;
            this.getAminoacidVector().get(i).setAsEpitope(state);
        }
    }
    public String printAA(){
        String s="";
        for(Aminoacid aa: this.aminoacidVector){
            s+=aa.getresidueName1();
        }
        return s;
        
    }
}
