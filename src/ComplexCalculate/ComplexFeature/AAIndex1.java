/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ComplexFeature;

/**
 *
 * @author OpenGress
 */
public class AAIndex1 {
    private String H_AccessionNumber;                                                   
    private String D_DataDescription;                                                 
    private String R_LITDBEntryNumber;                                                 
    private String A_Authors;                                                          
    private String T_TitleOfTheArticle;                                               
    private String J_JournalReference;
    private double[] I_AminoAcidIndex;

    public double[] getI_AminoAcidIndex() {
        return I_AminoAcidIndex;
    }

    public void setI_AminoAcidIndex(double[] I_AminoAcidIndex) {
        this.I_AminoAcidIndex = I_AminoAcidIndex;
    }

    
    public String getA_Authors() {
        return A_Authors;
    }

    public void setA_Authors(String A_Authors) {
        this.A_Authors = A_Authors;
    }

    public String getD_DataDescription() {
        return D_DataDescription;
    }

    public void setD_DataDescription(String D_DataDescription) {
        this.D_DataDescription = D_DataDescription;
    }

    public String getH_AccessionNumber() {
        return H_AccessionNumber;
    }

    public void setH_AccessionNumber(String H_AccessionNumber) {
        this.H_AccessionNumber = H_AccessionNumber;
    }

    public String getJ_JournalReference() {
        return J_JournalReference;
    }

    public void setJ_JournalReference(String J_JournalReference) {
        this.J_JournalReference = J_JournalReference;
    }

    public String getR_LITDBEntryNumber() {
        return R_LITDBEntryNumber;
    }

    public void setR_LITDBEntryNumber(String R_LITDBEntryNumber) {
        this.R_LITDBEntryNumber = R_LITDBEntryNumber;
    }

    public String getT_TitleOfTheArticle() {
        return T_TitleOfTheArticle;
    }

    public void setT_TitleOfTheArticle(String T_TitleOfTheArticle) {
        this.T_TitleOfTheArticle = T_TitleOfTheArticle;
    }
    
    
}
