package ComplexStructure;
public class Atom
{
  private String atomName;
  private String residueName;
  private char chainName;
  private String residueID;
  private double bFactor;
  private double x;
  private double y;
  private double z;
  private int type;
  private int atomNo;
  private boolean isexposed = false;

    public String getAtomName() {
        return atomName;
    }

    public void setAtomName(String atomName) {
        this.atomName = atomName;
    }

    public String getResidueName() {
        return residueName;
    }

    public void setResidueName(String residueName) {
        this.residueName = residueName;
    }

    public char getChainName() {
        return chainName;
    }

    public void setChainName(char chainName) {
        this.chainName = chainName;
    }

    public String getResidueID() {
        return residueID;
    }

    public void setResidueID(String residueID) {
        this.residueID = residueID;
    }

    public double getbFactor() {
        return bFactor;
    }

    public void setbFactor(double bFactor) {
        this.bFactor = bFactor;
    }

    public double getX() {
        return x;
    }

    public void setX(double x) {
        this.x = x;
    }

    public double getY() {
        return y;
    }

    public void setY(double y) {
        this.y = y;
    }

    public double getZ() {
        return z;
    }

    public void setZ(double z) {
        this.z = z;
    }

    public int getType() {
        return type;
    }

    public void setType(int type) {
        this.type = type;
    }

    public int getAtomNo() {
        return atomNo;
    }

    public void setAtomNo(int atomNo) {
        this.atomNo = atomNo;
    }

    public boolean isIsexposed() {
        return isexposed;
    }

    public void setIsexposed(boolean isexposed) {
        this.isexposed = isexposed;
    }
  
  
}