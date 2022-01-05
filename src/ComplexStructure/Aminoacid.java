package ComplexStructure;
import ComplexFeature.PSAI;
import StructureAnalysis.ParameterASA;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Vector;

public class Aminoacid
{
  //private Vector<Atom> atomvector;
    ArrayList<Atom> atomvector;
  private String residueName;
  private String pdbid;
  private char chainID;
  private String residueID;
  private boolean asEpitope;
  private boolean asExposed;
  public int binding = 0;
  private int alphaNum=0;
  //public Vector<Aminoacid> neighborAminoVector;
   
  public ArrayList<Aminoacid> neighborAminoVector;
  public ArrayList neighborDistance;
  public Aminoacid AbAmino;
  private ParameterASA pASA;
  private PSAI psai;
  public double antibodyneighborDistance;
  private char secondaryStructure;
  private double solventAccessibility;
  private double rsa;
  private double rsaTien2013;
  private double rsaRostSander2005;
  private double rsaMiller1987;
  private double rsaRose1985;
  private double conservation;
  private int exposed;
  private double andersonlogodd;
  private double parkerh;
  private double bfactor;
  private double lo;
  private double loseqwinsmooth;
  private double lostructwinsmooth;

    public Atom getAtomById (int no_atom){
        for( Atom a: this.atomvector){
            if(no_atom == a.getAtomNo()){
                return a;
            }
        }
        return null;
    }
  
  public PSAI getPsai() {
        return psai;
    }

    public void setPsai(PSAI psai) {
        this.psai = psai;
    }

    public boolean isAsExposed() {
        return asExposed;
    }

    public void setAsExposed(boolean asExposed) {
        this.asExposed = asExposed;
    }

    
  public double getLo() {
        return lo;
    }

    public void setLo(double lo) {
        this.lo = lo;
    }

    public double getLoseqwinsmooth() {
        return loseqwinsmooth;
    }

    public void setLoseqwinsmooth(double loseqwinsmooth) {
        this.loseqwinsmooth = loseqwinsmooth;
    }

    public double getLostructwinsmooth() {
        return lostructwinsmooth;
    }

    public void setLostructwinsmooth(double lostructwinsmooth) {
        this.lostructwinsmooth = lostructwinsmooth;
    }
  

    
  public ArrayList<Atom> getAtomvector() {
        return atomvector;
    }
  

    public double getCaBfactor() {
        return this.getCAlpha().getbFactor();
    }
    public double getBfactor() {
        double sumofBF=0;
        for(Atom a: this.atomvector){
            sumofBF += a.getbFactor();
        }
        return sumofBF/this.atomvector.size();
    }
    public double getRsaMiller1987() {
        return rsaMiller1987;
    }

    public void setRsaMiller1987(double rsaMiller1987) {
        this.rsaMiller1987 = rsaMiller1987;
    }

    public double getRsaRose1985() {
        return rsaRose1985;
    }

    public void setRsaRose1985(double rsaRose1985) {
        this.rsaRose1985 = rsaRose1985;
    }

    public double getRsaRostSander2005() {
        return rsaRostSander2005;
    }

    public void setRsaRostSander2005(double rsaRostSander2005) {
        this.rsaRostSander2005 = rsaRostSander2005;
    }

    public double getRsaTien2013() {
        return rsaTien2013;
    }

    public void setRsaTien2013(double rsaTien2013) {
        this.rsaTien2013 = rsaTien2013;
    }

   
  public String getPdbid() {
        return pdbid;
    }

    public void setPdbid(String pdbid) {
        this.pdbid = pdbid;
    }
  
  
  public Aminoacid()
  {
    this.atomvector = null;
  }

    public String getresidueName1() {
        return this.LetterToletter(this.getresidueName());
        
    }
  
  
  public Aminoacid(ArrayList<Atom> atomVector)
  {
    this.neighborAminoVector = new ArrayList();
    this.neighborDistance = new ArrayList();
    this.antibodyneighborDistance = -1.0D;
    this.AbAmino = new Aminoacid();
    
    this.atomvector = atomVector;
    this.residueName = ((Atom)atomVector.get(0)).getResidueName();
    this.chainID = ((Atom)atomVector.get(0)).getChainName();
    this.residueID = ((Atom)atomVector.get(0)).getResidueID();
  }

    public boolean isAsEpitope() {
        return asEpitope;
    }

    public ParameterASA getpASA() {
        return pASA;
    }
    public String intToString(int[] data){
        String s="";
        for(int i=0;i<data.length;i++){
            s= s+ String.valueOf(data[i])+"\t";
        }
        return s;
        
    }
    public String printParamASA(){
        int cn =pASA.getCn();
        String s = String.valueOf(cn);
        int[]hse =pASA.getHse();
        s= s+ intToString(hse);
        int[]qse = pASA.getQse();
        s= s+ intToString(qse);
        int[] fse = pASA.getFse();
        s= s+ intToString(fse);
        int[] efse = pASA.getEfse();
        s= s+ intToString(efse);
        int[] sfse= pASA.getSfse();
        s=s+intToString(sfse);
        return s;
    }

    public void setpASA(ParameterASA pASA) {
        this.pASA = pASA;
    }
    
    public void setAsEpitope(boolean asEpitope) {
        this.asEpitope = asEpitope;
    }
  
    public double getAndersonlogodd() {
        return andersonlogodd;
    }

    public void setAndersonlogodd(double andersonlogodd) {
        this.andersonlogodd = andersonlogodd;
    }

    public double getParkerh() {
        return parkerh;
    }

    public void setParkerh(double parkerh) {
        this.parkerh = parkerh;
    }

    public double getRsa() {
        return rsa;
    }

    public void setRsa(double rsa) {
        
        this.rsa = rsa;
    }
  
  public void setsecondaryStructure(char secondaryStructure)
  {
    this.secondaryStructure = secondaryStructure;
  }
  
  public void setsolventAccessibility(double solventAccessibility)
  {
      
    this.solventAccessibility = solventAccessibility;
  }
  
  public void setConservation(double conservation)
  {
    this.conservation = conservation;
  }
  
  public char getsecondaryStructure()
  {
    return this.secondaryStructure;
  }
  
  public double getsolventAccessibility()
  {
    return this.solventAccessibility;
  }
  
  public int getatomNum()
  {
    return this.atomvector.size();
  }
  
  public String getresidueName()
  {
    return this.residueName;
  }
  
  public char getchainID()
  {
    return this.chainID;
  }
  
  public String getresidueID()
  {
    return this.residueID;
  }
  
  public Atom getAtombyOrder(int i)
  {
    return (Atom)this.atomvector.get(i);
  }
  
  public double getConservation()
  {
    return this.conservation;
  }
  
  public int getalphaNum()
  {
    return this.alphaNum;
  }
  
  public void setalphaNum(int alphaNum)
  {
    this.alphaNum = alphaNum;
  }
  public Atom getCAlpha(){
      Atom a =null;
      for (int i=0;i<this.getatomNum();i++){          
          if (this.getAtombyOrder(i).getAtomName().equals("CA")){
            a = this.getAtombyOrder(i);
            break;
          }
          
      }
      return a;
  }
  public Atom getAtomofAA(String atomname){
      Atom a = new Atom();
      for (int i=0;i<this.getatomNum();i++){          
          if (this.getAtombyOrder(i).getAtomName().equals(atomname)){
            a= this.getAtombyOrder(i);
            break;
          }
          
      }
      return a;
  }
   public AtomCoordinate getAtomCoordinateofAA(String atomname){
      AtomCoordinate a = new AtomCoordinate();
      for (int i=0;i<this.getatomNum();i++){          
          if (this.getAtombyOrder(i).getAtomName().equals(atomname)){
            a.setX(this.getAtombyOrder(i).getX());
            a.setY(this.getAtombyOrder(i).getY());
            a.setZ(this.getAtombyOrder(i).getZ());
            break;
          }
          
      }
      return a;
  }
  //jarak dihitung menggunakan euclidean distance dan dipilih jarak terpendek antar atomVector
  public static double distanceCalculate(Aminoacid a1, Aminoacid a2)
  {
    String exceptelement = "O";
    int n = a1.getatomNum();
    int m = a2.getatomNum();
    double distance = 10000.0D;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        if ((!a1.getAtombyOrder(i).equals(exceptelement)) && (!a2.getAtombyOrder(j).equals(exceptelement)))
        {
          double new_distance = distanceBetweenPoint3D(a1.getAtombyOrder(i), a2.getAtombyOrder(j));
          if (new_distance < distance) {
            distance = new_distance;
          }
        }
      }
    }
    return distance;
  }
  
  public static double distanceBetweenPoint3D(Atom pa, Atom pb)
  {
    double distance = 0.0D;
    double x = pa.getX() - pb.getX();
    double y = pa.getY() - pb.getY();
    double z = pa.getZ() - pb.getZ();
    distance = Math.sqrt(x * x + y * y + z * z);
    return distance;
  }
  
  public int getBinding()
  {
    return this.binding;
  }
  
  public void setBinding(int binding)
  {
    if ((binding != 0) && (binding != 1) && (binding != -1)) {
      System.out.println("no proper value for binding");
    }
    this.binding = binding;
  }
  
  public void setChainID(char ChainID)
  {
    this.chainID = ChainID;
  }
  public  int letterToNum(){
  String letter = this.residueName;
      int num=-1;
  if (letter.equals("A")||letter.equals("ALA")) {
      num= 0;
    }
  if(letter.equals("R")||letter.equals("ARG")){
      num= 1;
  } 
  if(letter.equals("N")||letter.equals("ASN")){
      num= 2;
  }
  if(letter.equals("D")||letter.equals("ASP")){
     num= 3;
  }
  if(letter.equals("C")||letter.equals("CYS")){
      num= 4;
  }
  if(letter.equals("Q")||letter.equals("GLN")){
      num= 5;
  }
  if(letter.equals("E")||letter.equals("GLU")){
      num= 6;
  }
  if(letter.equals("G")||letter.equals("GLY")){
      num= 7;
  }
  if(letter.equals("H")||letter.equals("HIS")){
      num= 8;
  }
  if(letter.equals("I")||letter.equals("ILE")){
      num= 9;
  }
  if(letter.equals("L")||letter.equals("LEU")){
      num= 10;
  }
  if(letter.equals("K")||letter.equals("LYS")){
      num= 11;
  }
  if(letter.equals("M")||letter.equals("MET")){
      num= 12;
  }
  if(letter.equals("F")||letter.equals("PHE")){
      num= 13 ;
  }
  if(letter.equals("P")||letter.equals("PRO")){
      num= 14;
  }
  if(letter.equals("S")||letter.equals("SER")){
      num= 15;
  }
  if(letter.equals("T")||letter.equals("THR")){
      num= 16 ;
  }
  if(letter.equals("W")||letter.equals("TRP")){
      num= 17;
  }
  if(letter.equals("Y")||letter.equals("TYR")){
      num= 18;
  }
  if(letter.equals("V")||letter.equals("VAL")){
      num= 19;
  }
  return num;
     
  
  }
  public String LetterToletter(String letter)
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
}