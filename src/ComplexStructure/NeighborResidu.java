/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ComplexStructure;

/**
 *
 * @author OpenGress
 */
public class NeighborResidu {
    private Aminoacid aa;
    private int idxOnComplex;
    private double distance;
    private double beta;
    private double lo;

    public Aminoacid getAa() {
        return aa;
    }

    public void setAa(Aminoacid aa) {
        this.aa = aa;
    }

    
    public double getDistance() {
        return distance;
    }

    public void setDistance(double distance) {
        this.distance = distance;
    }

    public int getIdxOnComplex() {
        return idxOnComplex;
    }

    public void setIdxOnComplex(int idxOnComplex) {
        this.idxOnComplex = idxOnComplex;
    }

    public double getBeta() {
        return beta;
    }

    public void setBeta(double beta) {
        this.beta = beta;
    }

    public double getLo() {
        return lo;
    }

    public void setLo(double lo) {
        this.lo = lo;
    }
    
}
