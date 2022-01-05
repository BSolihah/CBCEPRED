/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ComplexStructure;

import jsat.classifiers.CategoricalResults;
import jsat.classifiers.DataPoint;

/**
 *
 * @author dell
 */
public class ExposedResidu {
    private String idResidue;
    private Atom cAlpha;
    private boolean isEpitope;
    private CategoricalResults predictedAs;
    private DataPoint feature;
    private Aminoacid amino;
    
    public ExposedResidu(Aminoacid aa, DataPoint f){
     this.idResidue = aa.getresidueID();
     this.cAlpha = aa.getCAlpha();
     this.isEpitope = aa.isAsEpitope();
     this.feature = f;
     this.amino = aa;
    }

    public Aminoacid getAmino() {
        return amino;
    }

    public void setAmino(Aminoacid amino) {
        this.amino = amino;
    }

    public CategoricalResults getPredictedAs() {
        return predictedAs;
    }

    public void setPredictedAs(CategoricalResults predictedAs) {
        this.predictedAs = predictedAs;
    }

    public String getIdResidue() {
        return idResidue;
    }

    public Atom getcAlpha() {
        return cAlpha;
    }

    public boolean isIsEpitope() {
        return isEpitope;
    }

    public DataPoint getFeature() {
        return feature;
    }
    
}
