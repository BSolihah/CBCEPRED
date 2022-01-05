package ComplexCalculate;

import ComplexStructure.Aminoacid;
import ComplexStructure.Atom;
import ComplexStructure.Chain;
import ComplexStructure.Complex;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Vector;

public class BindingSite
{
  private String exceptelement = "O";
  private ArrayList<Complex> ComplexVector;
  
  public BindingSite(ArrayList<Complex> complexVector)
  {
    this.ComplexVector = ((ArrayList)complexVector.clone());
    for (int i = 0; i < complexVector.size(); i++) {
      distance_antigen_antibody((Complex)complexVector.get(i));
    }
  }
  
  public BindingSite(Complex complex)
  {
    this.ComplexVector = new ArrayList();
    distance_antigen_antibody(complex);
  }
  
  private void distance_antigen_antibody(Complex complex)
  {
    for (int i = 0; i < complex.getAntigenChainIDs().length(); i++) {
      for (int j = 0; j < complex.getChainByChianID(complex.getAntigenChainIDs().charAt(i)).getaminoacidNum(); j++) {
        distance_amino_antibody(complex, complex.getChainByChianID(complex.getAntigenChainIDs().charAt(i)).getAminoacidbyorder(j));
      }
    }
  }
  
  private void distance_amino_antibody(Complex complex, Aminoacid antigen_amino)
  {
    double distance = 10000.0D;
    
    Aminoacid nearestAmino = null;
    for (int i = 0; i < complex.getAntibodyChainIDs().length(); i++) {
      for (int j = 0; j < complex.getChainByChianID(complex.getAntibodyChainIDs().charAt(i)).getaminoacidNum(); j++)
      {
        double new_distance = distance_amino_amino(antigen_amino, complex.getChainByChianID(complex.getAntibodyChainIDs().charAt(i)).getAminoacidbyorder(j));
        if (new_distance < distance)
        {
          distance = new_distance;
          nearestAmino = complex.getChainByChianID(complex.getAntibodyChainIDs().charAt(i)).getAminoacidbyorder(j);
        }
      }
    }
    antigen_amino.AbAmino = nearestAmino;
    antigen_amino.antibodyneighborDistance = distance;
    if (antigen_amino.antibodyneighborDistance <= 4.0D) {
      antigen_amino.setBinding(1);
    } else {
      antigen_amino.setBinding(-1);
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
  
  public void output(Complex complex)
  {
    int bindingNum = 0;
    int nonbindingNum = 0;
    for (int i = 0; i < complex.getAntigenChainIDs().length(); i++) {
      for (int j = 0; j < complex.getChainByChianID(complex.getAntigenChainIDs().charAt(i)).getaminoacidNum(); j++)
      {
        Aminoacid aminoacid = complex.getChainByChianID(complex.getAntigenChainIDs().charAt(i)).getAminoacidbyorder(j);
        System.out.println(complex.getcomplexName() + " " + aminoacid.getchainID() + " " + aminoacid.getresidueID() + " " + aminoacid.getresidueName() + " " + aminoacid.getBinding());
        if (aminoacid.getBinding() == 1) {
          bindingNum++;
        }
        if (aminoacid.getBinding() == -1) {
          nonbindingNum++;
        }
      }
    }
    System.out.println("binding site: " + bindingNum + " " + "nonbinding site: " + nonbindingNum);
  }
  
  public void output()
  {
    int bindingNum = 0;
    int nonbindingNum = 0;
    for (int k = 0; k < this.ComplexVector.size(); k++) {
      for (int i = 0; i < ((Complex)this.ComplexVector.get(k)).getAntigenChainIDs().length(); i++) {
        for (int j = 0; j < ((Complex)this.ComplexVector.get(k)).getChainByChianID(((Complex)this.ComplexVector.get(k)).getAntigenChainIDs().charAt(i)).getaminoacidNum(); j++)
        {
          Aminoacid aminoacid = ((Complex)this.ComplexVector.get(k)).getChainByChianID(((Complex)this.ComplexVector.get(k)).getAntigenChainIDs().charAt(i)).getAminoacidbyorder(j);
          System.out.println(((Complex)this.ComplexVector.get(k)).getcomplexName() + " " + aminoacid.getchainID() + " " + aminoacid.getresidueID() + " " + aminoacid.getresidueName() + " " + aminoacid.getBinding());
          if (aminoacid.getBinding() == 1) {
            bindingNum++;
          }
          if (aminoacid.getBinding() == -1) {
            nonbindingNum++;
          }
        }
      }
    }
    System.out.println("binding site: " + bindingNum + " " + "nonbinding site: " + nonbindingNum);
  }
}
