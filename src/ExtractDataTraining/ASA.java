/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ExtractDataTraining;

/**
 *
 * @author OpenGress
 */
public class ASA {
    private String id;
    private boolean epitopeStatus;
    private double asaValue;

    public ASA(String id, boolean status, double value){
        this.asaValue=value;
        this.epitopeStatus = status;
        this.id=id;
    }
    public double getAsaValue() {
        return asaValue;
    }

    public void setAsaValue(double asaValue) {
        this.asaValue = asaValue;
    }

    public boolean isEpitopeStatus() {
        return epitopeStatus;
    }

    public void setEpitopeStatus(boolean epitopeStatus) {
        this.epitopeStatus = epitopeStatus;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }
    
    
}
