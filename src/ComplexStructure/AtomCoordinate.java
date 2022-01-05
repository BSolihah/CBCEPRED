/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ComplexStructure;

/**
 *
 * @author OpenGress
 */
public class AtomCoordinate {
    private double x;
    private double y;
    private double z;

    public AtomCoordinate(){
        this.x=0;
        this.y=0;
        this.z=0;
    }
    public AtomCoordinate(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }
    public AtomCoordinate(AtomCoordinate ac){
         this.x = ac.getX();
        this.y = ac.getY();
        this.z = ac.getZ();
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
    
}
