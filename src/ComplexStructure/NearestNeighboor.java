/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ComplexStructure;

/**
 *
 * @author OpenGress
 */
public class NearestNeighboor {
    private boolean asepitope;
    private int neighborSum;
    private String residuId;
    private String neighboor;

    public NearestNeighboor(boolean epi, int sum, String id,String s){
        this.asepitope=epi;
        this.neighborSum=sum;
        this.residuId =id;
        this.neighboor=s;
    }
    
    public String printNN(){
    String dt= this.residuId+ ";"+ String.valueOf(this.neighborSum)+";"+this.neighboor+ String.valueOf( this.asepitope)+"\n";
    return dt;
    }
}
