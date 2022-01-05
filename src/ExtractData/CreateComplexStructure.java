/*
 * membaca data struktur dari pdb
 * menandai setiap residu binding
 * menandai residu terekspose
 */
package ExtractData;

import ComplexFeature.PSAI;
import ComplexStructure.Aminoacid;
import ComplexStructure.Complex;
import ComplexStructure.Epitope;
import ExtractDataTraining.ExtractDataBenchmark;
import ExtractDataTraining.ExtractProtusionIndex;
import java.util.ArrayList;

/**
 *
 * @author dell
 */
public class CreateComplexStructure {
    ArrayList<Complex> complexlist;
    ExtractDataBinding edb;
    ArrayList<String> pdblist;
    ArrayList<String> chainlist;
    Boolean complexIsLoad = false;
    
    
    public CreateComplexStructure(String list_file){
        complexlist = new ArrayList();
        edb = new ExtractDataBinding();//id dan residuId dari binding site
        //loadStructureFromFile(list_file);//membentuk struktur dari file
       // loadStructureFromFileTest(list_file);
    }

    public CreateComplexStructure() {
        //throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
         complexlist = new ArrayList();
        edb = new ExtractDataBinding();
    }

    public void setPdblist(ArrayList<String> pdblist) {
        this.pdblist = pdblist;
    }

    public void setChainlist(ArrayList<String> chainlist) {
        this.chainlist = chainlist;
    }

    
    public ArrayList<Complex> getComplexlist() {
        return complexlist;
    }

    public ArrayList<String> getPdblist() {
        return pdblist;
    }

    public ArrayList<String> getChainlist() {
        return chainlist;
    }

    public Boolean getComplexIsLoad() {
        return complexIsLoad;
    }
    public void loadStructureFromContact(String listContact){
         ExtractDataBinding edb = new ExtractDataBinding();
         ArrayList<Epitope> epilist= edb.loadBindingFromContact(listContact);
         System.out.println("epilist"+ epilist.get(0).getPdbid() );
         
         for(int i=0;i<epilist.size();i++){
            String pdbid = epilist.get(i).getPdbid();
            String chain = epilist.get(i).getChain();
         //   System.out.println(pdbid + ": "+ chain);
            Complex complex = new Complex(pdbid, chain);   
           setExposedResidueFromPSAIA(complex,pdbid,chain);
           setBindingFromEpilist(complex,epilist.get(i));
            complexlist.add(complex);
         }
         
    }
    public void loadStructureFromListFile(String list_file){
        ExtractDataBinding edb = new ExtractDataBinding();
         ArrayList<Epitope> epilist= edb.loadEpiFromBinding(list_file);
         
         for(int i=0;i<epilist.size();i++){
            String pdbid = epilist.get(i).getPdbid();
            String chain = epilist.get(i).getChain();
         //   System.out.println(pdbid + ": "+ chain);
            Complex complex = new Complex(pdbid, chain);   
           setExposedResidueFromPSAIA(complex,pdbid,chain);
           setBindingFromEpilist(complex,epilist.get(i));
            complexlist.add(complex);
         }
    }
    public void loadStructureFromPdbChList(String filelist, int part){
    
    }
   
    public void loadStructureFromFile(String list_file){
        complexlist = new ArrayList();
        //load list pdb dari file
        //create complex chain
        System.out.println("list_file: "+ list_file);
        loadchainlist(list_file);
        System.out.println("jumlah pdbid: "+ pdblist.size());
        ExtractDataBenchmark extractDB = new ExtractDataBenchmark();
       // ArrayList<Epitope> epilist = extractDB.loadBenchmarkDataTest("epitopelisttestset.txt");
       // ArrayList<Epitope> epilist = extractDB.loadBenchmarkData("benchmarksequencedataedit1.txt");
        /*
        for(Epitope e: epilist){
        System.out.println(e.getPdbid()+":"+e.getChain());
        }*/
        for(int i=0;i<pdblist.size();i++){
            String pdbid = pdblist.get(i).trim();
            String chain = chainlist.get(i).trim();
            Complex complex = new Complex(pdbid, chain);   
           // System.out.println("complex: "+complex.getComplexName()+complex.getAntigenChainIDs());
            setExposedResidueFromPSAIA(complex,pdbid,chain);
            // System.out.println("pdbid: "+ pdbid + ", chain: "+ chain );
          //   System.out.println("from epilist:"+ epilist.get(i).getPdbid()+":"+epilist.get(i).getChain()+":"+epilist.get(i).getIdxsequenceAsEpitope().size());
         //   complex.getChainByChianID(chain.charAt(0)).setAAAsEpitope(epilist.get(i));
            ArrayList<Binding> bindlist = edb.loadDataBindingFromFile(pdbid,chain);
            ArrayList<String> bindlistS = edb.getResidueIDBinding(bindlist);
            setBindingState(complex,bindlistS);
            
           
           // System.out.println("bindlist: "+ bindlist.size());
            
//                System.out.println("id residu binding :"+ bindlist.get(bindlist.size()-1).residueId);
            
            // System.out.println("id residu terakhir complex: "+complex.getAminoacidVector().get(complex.getAminoacidVector().size()-1).getresidueID());
            
            complexlist.add(complex);
        }
        this.complexIsLoad=true;
        
    }
    public void loadStructureFromFileTest(String list_file){
        
        //load list pdb dari file
        //create complex chain
      //  System.out.println("list_file: "+ list_file);
       // loadchainlist(list_file);
//        System.out.println("jumlah pdbid: "+ pdblist.size());
        ExtractDataBenchmark extractDB = new ExtractDataBenchmark();
        ArrayList<Epitope> epilist = extractDB.loadBenchmarkDataTest("epitopelisttestset.txt");
        // ArrayList<Epitope> epilist = extractDB.loadBenchmarkDataTest("epitopelisttestset.txt");
       // ArrayList<Epitope> epilist = extractDB.loadBenchmarkData("benchmarksequencedataedit1.txt");
        /*
        for(Epitope e: epilist){
        System.out.println("epilist: "+e.getPdbid()+":"+e.getChain());
        }*/
        for(int i=0;i<epilist.size();i++){
            String pdbid = epilist.get(i).getPdbid();
            String chain = epilist.get(i).getChain();
            Complex complex = new Complex(pdbid, chain);   
           // System.out.println("complex: "+complex.getComplexName()+complex.getAntigenChainIDs());
            setExposedResidueFromPSAIA(complex,pdbid,chain);
            setBindingFromEpilist(complex,epilist.get(i));
            // System.out.println("pdbid: "+ pdbid + ", chain: "+ chain );
          //   System.out.println("from epilist:"+ epilist.get(i).getPdbid()+":"+epilist.get(i).getChain()+":"+epilist.get(i).getIdxsequenceAsEpitope().size());
         //   complex.getChainByChianID(chain.charAt(0)).setAAAsEpitope(epilist.get(i));
           // ArrayList<Binding> bindlist = edb.loadDataBindingFromFile(pdbid,chain);
           // ArrayList<String> bindlistS = edb.getResidueIDBinding(bindlist);
           // setBindingState(complex,bindlistS);
            
           
           // System.out.println("bindlist: "+ bindlist.size());
            
//                System.out.println("id residu binding :"+ bindlist.get(bindlist.size()-1).residueId);
            
            // System.out.println("id residu terakhir complex: "+complex.getAminoacidVector().get(complex.getAminoacidVector().size()-1).getresidueID());
            
            complexlist.add(complex);
        }
        this.complexIsLoad=true;
        
    }
     public void loadchainlist(String list){
        
        edb.loadListFile(list);
        this.setPdblist(edb.getPdbid_list());
        this.setChainlist( edb.getChain_list());
    }
     public void setBindingFromEpilist(Complex c, Epitope e){
         String pdbid = e.getPdbid();
//         String state = e.getEpipos();
 //        System.out.println("jumlah state: "+ state.length());
         ArrayList<String> idxEpi = e.getIdxEpi();
         
         //if(c.getAminoacidVector().size()==state.length()){
       //  System.out.println("jumlah vektor"+ c.getAminoacidVector().size());
                     
         for(int i=0;i<c.getAminoacidVector().size();i++){
             //System.out.println(c.getAminoacidVector().get(i).getresidueID());
             String idAA = c.getAminoacidVector().get(i).getresidueID();
             if(idxEpi.contains(idAA)){
                 c.getAminoacidVector().get(i).setAsEpitope(true);
             }
             //c.getAminoacidVector().get(idxepi.get(i)).setAsEpitope(true);
             }
         //}
     }
    public void setBindingResidu(String list){
    //untuk setiap complex dlm list
    //set binding residunya
        loadchainlist(list);
        if (this.complexIsLoad==false){
            this.loadStructureFromFile(list);}
        
        for(int i=0;i<pdblist.size();i++){
            String pdbid = pdblist.get(i);
            String chain = chainlist.get(i);
            ArrayList<Binding> bindlist = edb.loadDataBindingFromFile(pdbid,chain);
           ArrayList<String> bindlistS = edb.getResidueIDBinding(bindlist);
            Complex c= complexlist.get(i);
            setBindingState(c,bindlistS);
            
        }
    }
    private void setBindingStateFromFile(Complex c,String chain,Epitope e ){
        
        
    }
    private void setBindingState(Complex c, ArrayList<String> bindlist){
        int count=0;
        for(int i=0;i<c.getAminoacidVector().size();i++){
            String id=  c.getAminoacidVector().get(i).getresidueID();
            if(count<bindlist.size()){
                String bid = bindlist.get(count);
                if(bid.equalsIgnoreCase(id)&&count<bindlist.size()){
                count +=1;
               // System.out.print(id+";");
                c.getAminoacidVector().get(i).setAsEpitope(true);
            }
            }
            
           // System.out.println(id+";"+bid);
            
        }
        //System.out.println();
        
        
    }
    public void setExposedResidueFromPSAIA(Complex c, String pdbid,String chain){
        ExtractProtusionIndex ePI = new ExtractProtusionIndex();
        ArrayList<PSAI> psailist =ePI.readPIFromFile(pdbid);
        System.out.println("psailist: "+psailist.size());
        ArrayList<PSAI> list =ePI.copyPIOfChain(chain, psailist);
        //untuk setiap aa pada c set rsabytien
        for(int i=0;i<list.size();i++){
            double[]psaiparam =list.get(i).getPsaiParam();
               String residu = c.getAminoacidVector().get(i).getresidueName();
           //    System.out.println("psaiparam"+ psaiparam[0]);
            c.getAminoacidVector().get(i).setRsaTien2013(calculateRSAByTienMSA(residu,psaiparam[0]));
            c.getAminoacidVector().get(i).setsolventAccessibility(psaiparam[0]);
            //set PSAI
            c.getAminoacidVector().get(i).setPsai(list.get(i));
           //  System.out.println("rsa tien:"+ c.getAminoacidVector().get(i).getRsaTien2013());
        }
                
    }
     
    private double calculateRSAByTienMSA(String residueName, double area)
  {
    if ((residueName.equals("ALA")) || (residueName.equals("A"))) {
      return area / 129.0D;
    }
    if ((residueName.equals("ARG")) || (residueName.equals("R"))) {
      return area / 274.0D;
    }
    if ((residueName.equals("ASN")) || (residueName.equals("N"))) {
      return area / 195.0D;
    }
    if ((residueName.equals("ASP")) || (residueName.equals("D"))) {
      return area / 193.0D;
    }
    if ((residueName.equals("CYS")) || (residueName.equals("C"))) {
      return area / 167.0D;
    }
    if ((residueName.equals("GLN")) || (residueName.equals("Q"))) {
      return area / 223.0D;
    }
    if ((residueName.equals("GLU")) || (residueName.equals("E"))) {
      return area / 225.0D;
    }
    if ((residueName.equals("GLY")) || (residueName.equals("G"))) {
      return area / 104.0D;
    }
    if ((residueName.equals("HIS")) || (residueName.equals("H"))) {
      return area / 224.0D;
    }
    if ((residueName.equals("ILE")) || (residueName.equals("I"))) {
      return area / 197.0D;
    }
    if ((residueName.equals("LEU")) || (residueName.equals("L"))) {
      return area / 201.0D;
    }
    if ((residueName.equals("LYS")) || (residueName.equals("K"))) {
      return area / 236.0D;
    }
    if ((residueName.equals("MET")) || (residueName.equals("M"))) {
      return area / 224.0D;
    }
    if ((residueName.equals("PHE")) || (residueName.equals("F"))) {
      return area / 240.0D;
    }
    if ((residueName.equals("PRO")) || (residueName.equals("P"))) {
      return area / 159.0D;
    }
    if ((residueName.equals("SER")) || (residueName.equals("S"))) {
      return area / 155.0D;
    }
    if ((residueName.equals("THR")) || (residueName.equals("T"))) {
      return area / 172.0D;
    }
    if ((residueName.equals("TRP")) || (residueName.equals("W"))) {
      return area / 285.0D;
    }
    if ((residueName.equals("TYR")) || (residueName.equals("Y"))) {
      return area / 263.0D;
    }
    if ((residueName.equals("VAL")) || (residueName.equals("V"))) {
      return area / 174.0D;
    }
    return -1.0D;
  }
    
}
