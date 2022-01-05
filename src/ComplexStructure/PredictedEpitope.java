/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ComplexStructure;

import jsat.classifiers.DataPoint;

/**
 *
 * @author dell
 */
public class PredictedEpitope {
    private String idResidue;
    private String predictedState;//TF or FP
    private DataPoint spatialDP;
    private Aminoacid amino;

    public Aminoacid getAmino() {
        return amino;
    }

    public void setAmino(Aminoacid amino) {
        this.amino = amino;
    }
    
    public String getIdResidue() {
        return idResidue;
    }

    public void setIdResidue(String idResidue) {
        this.idResidue = idResidue;
    }

    public String getPredictedState() {
        return predictedState;
    }

    public void setPredictedState(String predictedState) {
        this.predictedState = predictedState;
    }

    public DataPoint getSpatialDP() {
        return spatialDP;
    }

    public void setSpatialDP(DataPoint spatialDP) {
        this.spatialDP = spatialDP;
    }
    
    
    
}
