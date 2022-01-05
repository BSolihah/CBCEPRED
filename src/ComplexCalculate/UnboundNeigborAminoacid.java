package ComplexCalculate;

import ComplexStructure.Aminoacid;
import ComplexStructure.Atom;
import ComplexStructure.Chain;
import ComplexStructure.Complex;
import java.util.ArrayList;


public class UnboundNeigborAminoacid
{
  private String exceptelement = "O";
  double threshold;
  
  public UnboundNeigborAminoacid(ArrayList<Complex> complexVector, double threshold)
  {
    this.threshold = threshold;
    for (int i = 0; i < complexVector.size(); i++) {
      distance_antigen_antigen((Complex)complexVector.get(i));
    }
    setAlphaTom(complexVector);
    System.out.print("Alpha tom sudah diset");
  }
  
  private void distance_antigen_antigen(Complex complex)
  {
    for (int i = 0; i < complex.getAntigenChainIDs().length(); i++) {
      for (int j = 0; j < complex.getChainByChianID(complex.getAntigenChainIDs().charAt(i)).getaminoacidNum(); j++) {
        distance_antigen_antigen(complex, complex.getChainByChianID(complex.getAntigenChainIDs().charAt(i)).getAminoacidbyorder(j));
      }
    }
  }
  
  private void distance_antigen_antigen(Complex complex, Aminoacid antigen_amino)
  {
    ArrayList<Aminoacid> antigen_amino_all = new ArrayList();
    for (int i = 0; i < complex.getAntigenChainIDs().length(); i++) {
      for (int j = 0; j < complex.getChainByChianID(complex.getAntigenChainIDs().charAt(i)).getaminoacidNum(); j++) {
        antigen_amino_all.add(complex.getChainByChianID(complex.getAntigenChainIDs().charAt(i)).getAminoacidbyorder(j));
      }
    }
    double[] distance = new double[antigen_amino_all.size()];
    for (int i = 0; i < antigen_amino_all.size(); i++)
    {
      distance[i] = distance_amino_amino(antigen_amino, (Aminoacid)antigen_amino_all.get(i));
      if (distance[i] < 1.0E-13D) {
        distance[i] = 10000.0D;
      }
    }
    int[] index = new int[antigen_amino_all.size()];
    for (int i = 0; i < antigen_amino_all.size(); i++) {
      index[i] = i;
    }
    for (int i = 0; i < antigen_amino_all.size(); i++) {
      for (int j = i + 1; j < antigen_amino_all.size(); j++) {
        if (distance[i] > distance[j])
        {
          int temp1 = index[i];
          index[i] = index[j];
          index[j] = temp1;
          double temp2 = distance[i];
          distance[i] = distance[j];
          distance[j] = temp2;
        }
      }
    }
    int k = 0;
    for (int i = 0;; i++)
    {
      if (k == Complex.neigborNum) {
        return;
      }
      if (((Aminoacid)antigen_amino_all.get(index[i])).getsolventAccessibility() >= this.threshold)
      {
        antigen_amino.neighborAminoVector.add((Aminoacid)antigen_amino_all.get(index[i]));
        antigen_amino.neighborDistance.add(Double.valueOf(distance[i]));
        k++;
      }
    }
  }
  
  private double distance_amino_amino(Aminoacid aa1, Aminoacid aa2)
  {
    int n = aa1.getatomNum();
    int m = aa2.getatomNum();
    double distance = 10000.0D;
    if ((aa1.getatomNum() == aa2.getatomNum()) && (aa1.getresidueID().equals(aa2.getresidueID())) && (aa1.getchainID() == aa2.getchainID())) {
      return distance;
    }
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < m; j++) {
        if ((!aa1.getAtombyOrder(i).equals(this.exceptelement)) && (!aa2.getAtombyOrder(j).equals(this.exceptelement)))
        {
          double new_distance = distanceBetweenPoint3D(aa1.getAtombyOrder(i), aa2.getAtombyOrder(j));
          if (new_distance < distance) {
            distance = new_distance;
          }
        }
      }
    }
    return distance;
  }
  
  private double distanceBetweenPoint3D(Atom pa, Atom pb)
  {
    double distance = 0.0D;
    double x = pa.getX() - pb.getX();
    double y = pa.getY() - pb.getY();
    double z = pa.getZ() - pb.getZ();
    distance = Math.sqrt(x * x + y * y + z * z);
    return distance;
  }
  
  private void setAlphaTom(ArrayList<Complex> complexVector)
  {
    for (int i = 0; i < complexVector.size(); i++) {
      setAlphaTom((Complex)complexVector.get(i));
    }
  }
  
  private void setAlphaTom(Complex complex)
  {
    ArrayList<Atom> allatomVector = new ArrayList();
    for (int i = 0; i < complex.getAntigenChainIDs().length(); i++) {
      for (int j = 0; j < complex.getChainByChianID(complex.getAntigenChainIDs().charAt(i)).getaminoacidNum(); j++) {
        for (int k = 0; k < complex.getChainByChianID(complex.getAntigenChainIDs().charAt(i)).getAminoacidbyorder(j).getatomNum(); k++) {
          allatomVector.add(complex.getChainByChianID(complex.getAntigenChainIDs().charAt(i)).getAminoacidbyorder(j).getAtombyOrder(k));
        }
      }
    }
    for (int i = 0; i < complex.getAntigenChainIDs().length(); i++) {
      for (int j = 0; j < complex.getChainByChianID(complex.getAntigenChainIDs().charAt(i)).getaminoacidNum(); j++)
      {
        ArrayList<Atom> targetatomVector = new ArrayList();
        for (int k = 0; k < complex.getChainByChianID(complex.getAntigenChainIDs().charAt(i)).getAminoacidbyorder(j).getatomNum(); k++) {
          targetatomVector.add(complex.getChainByChianID(complex.getAntigenChainIDs().charAt(i)).getAminoacidbyorder(j).getAtombyOrder(k));
        }
        int num = countAlphaNum(targetatomVector, allatomVector);
        complex.getChainByChianID(complex.getAntigenChainIDs().charAt(i)).getAminoacidbyorder(j).setalphaNum(num);
      }
    }
  }
  
  
  private int countAlphaNum(ArrayList<Atom> targetatomVector, ArrayList<Atom> allatomVector)
  {
    int num = 0;
    for (int i = 0; i < allatomVector.size(); i++)
    {
      double distance = minDistanceBetweenAtomandAtomset((Atom)allatomVector.get(i), targetatomVector);
      if ((distance > 1.0E-9D) && (distance < 10.0D) && (((Atom)allatomVector.get(i)).getAtomName().equals("CA"))) {
        num++;
      }
    }
    return num;
  }
  
  private double minDistanceBetweenAtomandAtomset(Atom atom, ArrayList<Atom> atomVector)
  {
    double minDistance = 10000.0D;
    for (int i = 0; i < atomVector.size(); i++)
    {
      double temp = distanceBetweenPoint3D(atom, (Atom)atomVector.get(i));
      if (temp < minDistance) {
        minDistance = temp;
      }
    }
    return minDistance;
  }
 
 
  public void printCN (Complex c, char chain){
    
      for (int i=0;i< c.getAminoacidNumByChainID(chain); i++){
          System.out.print(c.getChainByChianID(chain).getAminoacidbyorder(i).getresidueID()+"\t");
          System.out.print(c.getChainByChianID(chain).getAminoacidbyorder(i).getresidueName()+"\t");
          System.out.print(c.getChainByChianID(chain).getAminoacidbyorder(i).getalphaNum()+"\n");
      }
}
}

