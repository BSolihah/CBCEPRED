
package ExtractDataTraining;

import java.util.ArrayList;

/**
 *
 * @author OpenGress
 */
public class Benchmark {
    private String sequence;
    private String epipos;
    private String exposedpos;
    private ArrayList<Integer> idxEpi;
    private ArrayList<Integer> idxexp;
    private ArrayList<Integer>idxEpiAAP;
    private ArrayList<Integer>idxexpAAP;
    private ArrayList<Integer>idxOneEpiAAP;
    

    public String getEpipos() {
        return epipos;
    }

    public void setEpipos(String epipos) {
        this.epipos = epipos;
    }

    public String getExposedpos() {
        return exposedpos;
    }

    public void setExposedpos(String exposedpos) {
        this.exposedpos = exposedpos;
    }

    public ArrayList<Integer> getIdxEpi() {
        return idxEpi;
    }

    public ArrayList<Integer> getIdxEpiAAP() {
        return idxEpiAAP;
    }

    public void setIdxEpiAAP(ArrayList<Integer> idxEpiAAP) {
        this.idxEpiAAP = idxEpiAAP;
    }

    public ArrayList<Integer> getIdxOneEpiAAP() {
        return idxOneEpiAAP;
    }

    public void setIdxOneEpiAAP(ArrayList<Integer> idxOneEpiAAP) {
        this.idxOneEpiAAP = idxOneEpiAAP;
    }

    public ArrayList<Integer> getIdxexpAAP() {
        return idxexpAAP;
    }

    public void setIdxexpAAP(ArrayList<Integer> idxexpAAP) {
        this.idxexpAAP = idxexpAAP;
    }

    public void setIdxEpi(ArrayList<Integer> idxEpi) {
        this.idxEpi = idxEpi;
    }

    public ArrayList<Integer> getIdxexp() {
        return idxexp;
    }

    public void setIdxexp(ArrayList<Integer> idxexp) {
        this.idxexp = idxexp;
    }

    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }
    
}
