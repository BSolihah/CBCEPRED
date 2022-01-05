/*
 * dihitung berdasarkan pasangan yang ada didalam dataset
 * 1. buat graph berdasarkan jarak dengan treshold tertentu
 * 2. identifikasi jumlah pasangan epi-epi, non-epi non epi, epi- non epi, dan non epi- epi 
* (cek untuk beberapa dataset, bgm perbandingannya dan pikirkan bagaimana adjustmennya sehingga bisa 
* digunakan pada sembarang data set
 *  
 */
package ComplexFeature;

import static CBCTop.AntigenGraphClustering.getEpiFilename;
import static CBCTop.AntigenGraphClustering.getFileListFromFolder;
import static CBCTop.AntigenGraphClustering.readEpiListFromFile;
import ComplexStructure.Aminoacid;
import ComplexStructure.Atom;
import ComplexStructure.Complex;
import ComplexStructure.ExposedResidu;
import ExtractDataTraining.ExtractSphereExposure;
import ExtractDataTraining.GraphAntigenExposedAtom;
import static ExtractDataTraining.GraphAntigenExposedAtom.loadExposedAtom;
import static ExtractDataTraining.GraphAntigenExposedAtom.loadExposedAtomResiduLevel;
import ExtractFromFile.IO;
import java.util.ArrayList;

import SpatialAnalysis.PairIdResidu;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import jsat.classifiers.DataPoint;
import jsat.linear.DenseVector;

/**
 *
 * @author dell
 */
public class AminoAcidPairByFuncSubgrup {
    private double[][]expected_value_epi_epi;
    private double [][] epi_epi;
    private double [][] epi_nonepi;
    private double [][] nonepi_nonepi;
    private String listfile;
    private GraphAntigenExposedAtom ga = null;
    
    
    
    
    public AminoAcidPairByFuncSubgrup(String path_pdb, String path_exp_atom, double t_asa, double t_dist, String path, String filename){
       
        //this.listfile = file_list;
        epi_epi = new double[13][13];
        epi_nonepi = new double[13][13];
        nonepi_nonepi = new double[13][13];
        expected_value_epi_epi = new double[13][13];
        for(int i=0;i<13;i++){
            for(int j=0;j<13;j++){
                epi_epi[i][j]=1;
                epi_nonepi[i][j]=1;
                nonepi_nonepi[i][j]=1;
                expected_value_epi_epi[j][i]=0;
                
            }
        }
        if(new File(path+ filename).exists()){
            readFromFileAAPByCategory(path,filename);
           // printAAPByCategory();
                      
        }else{
            ArrayList<ExposedOnComplex> eoc = getECOListfromDataset(path_pdb, path_exp_atom);
            System.out.println("eoc: "+ eoc.size());
            calculateFrequenciesAminoacidPair(eoc, t_dist);
            savetofileAAPByCategory(path, filename);
        }
        calculateExpectedValue();
        printExpectedValue();
    }
    public static String findEpiFilenameinlist(ArrayList<String> epifiles, String pdbid, String chain){
        String fn="";
        for(String s: epifiles){
            if(s.startsWith(pdbid.toLowerCase()) && s.endsWith(chain.toUpperCase())){
                fn = s;
                epifiles.remove(s);
                break;
            }
        }
        return fn;
    }
    public ArrayList<ExposedOnComplex> getECOListfromDataset(String path_pdb, String path_exp_atom){
        this.ga = new GraphAntigenExposedAtom();
        String path_in_complex = "src/dataset/gao/chain/";
        String path_in_atom = "src/dataset/gao/antigens/"; //berisi pdbid+chain dan pdbid+chain
        String path_in_res = "src/dataset/gao/antigens/res/"; //berisi pdbid+chain dan pdbid+chain+res
        String path_epi = "src/dataset/gao/epitopes/";
        ArrayList<String> epi_file_list = getFileListFromFolder(path_epi, "");
        System.out.println("contoh epi file: "+epi_file_list.get(0));
        ArrayList<ExposedOnComplex> eoclist = new ArrayList();
        try{
            FileReader fr = null;
            File folder = new File(path_in_res);
            File[] listOfFiles = folder.listFiles();
            System.out.println("jumlah listoffile: "+ listOfFiles.length);
            for (int i = 0; i < listOfFiles.length; i++) {
                if (listOfFiles[i].isFile()) {   
                    String fileNameRes = listOfFiles[i].getName();
                    
                    String pdbid = fileNameRes.substring(0,4);
                    String chain = fileNameRes.substring(4,5);
                    System.out.println("chain: "+chain);
                    String fileNameAtom = fileNameRes.substring(0, 5);
                    IO io = new IO();
                    AAIndex aai;
                    aai= io.extractAAIndexFromFile();
                    aai.composeAAI();
                    System.out.println("pdb: "+pdbid + ":"+chain );
                    try{
                        Complex c = new Complex(path_in_complex,pdbid,chain);
                        ExposedOnComplex eoc = new ExposedOnComplex();
                        eoc.setPdbid(pdbid);
                        eoc.setChain(chain);
                        
                        System.out.println("complex: "+c.getComplexName()+ c.getAminoacidVector().size());
                        String epi_file_name = path_epi + findEpiFilenameinlist(epi_file_list,pdbid,chain);
                        System.out.println("epi_file_name"+epi_file_name);
                        ArrayList<String> epitope = readEpiListFromFile(epi_file_name);
                        //baca dari file exposed atom dan exposed residu
                        ArrayList<Atom> list_atom_res = loadExposedAtomResiduLevel(path_in_res, fileNameRes);
                        System.out.println("list_atom_res: "+ list_atom_res.get(0).getX()+" "+list_atom_res.get(0).getY()+" "+list_atom_res.get(0).getZ());
                        ArrayList<Atom> list_atom_atom = loadExposedAtom(path_in_atom, fileNameAtom);
                        System.out.println("list_atom_atom: "+ list_atom_atom.get(0).getX()+" "+list_atom_atom.get(0).getY()+" "+list_atom_atom.get(0).getZ()+" "+list_atom_atom.get(0).getAtomName() );
                        ga.concatResiduInfo(list_atom_res,list_atom_atom);
                        //mulai melakukan extract exposedoncomplex
                        ArrayList<ExposedResidu> listExposedR= new ArrayList();
                        ArrayList<Aminoacid> listAminoacid= new ArrayList();
                        boolean aaResIsSet = false;
                        for(Atom atom : list_atom_res){
                            String res_id = atom.getResidueID();
                            if(aaResIsSet == false){
                                Aminoacid aa = c.getAminoById(res_id);
                                c.getAminoById(res_id).setAsExposed(true);
                                if(epitope.contains(res_id))
                                    c.getAminoById(res_id).setAsEpitope(true);
                                aaResIsSet = true;
                                ExtractSphereExposure ese = new ExtractSphereExposure();
                                ese.calculateSE(c.getAminoById(res_id), c);
                                List<Double> atrib = new ArrayList();
                                atrib.add(aa.getRsaTien2013());            
                                atrib.addAll(aa.getpASA().getParamASAD());
                                atrib.add(aa.getCaBfactor());
                                atrib.add(aa.getBfactor());
                                atrib.add(aa.getLo());
                                atrib.addAll(aa.getPsai().getListPSAI());
                                double[] withoutNan =aai.getAaiWithoutNanData(aa.letterToNum());
                                ArrayList<Double> aaid = new ArrayList();
                                for(int ii=0;ii<withoutNan.length;ii++){
                                    aaid.add(withoutNan[i]);
                                }
                                atrib.addAll(aaid);
                                DenseVector vec = new DenseVector(atrib);
                                DataPoint p= new DataPoint(vec);
                                ExposedResidu expR = new ExposedResidu(aa,p);
                                listExposedR.add(expR);
                                listAminoacid.add(aa);
                                
                            }                           
                          System.out.println("eoc"+ eoc.getPdbid()+eoc.getChain()+eoc.getListExposedR().size());
                        }
                        eoc.setListExposedR(listExposedR);
                        eoc.setListAA(listAminoacid);
                        eoclist.add(eoc);
        
                    }catch(Exception e){
                        System.out.println(e.getMessage());
                    }
                    
                    
                }
            }
        }catch(Exception e){
            System.out.println("exception in generate ga:" +e.getMessage());
        }

               
        return eoclist;
        
    }
    //fungsi menghitung bobot berdasarkan expected value
    //
    public double[][] calculateNormalizeWeight(){
        // nw = sum((nxy - expxy)^2)/expex
        double[][] norm_weight = new double[13][13];
        for(int i=0;i<13;i++){
            for(int j=0;j<13;j++){
                if(i!=j){
                    norm_weight[i][j]= Math.pow(this.epi_epi[i][j]+this.epi_epi[j][i] - this.expected_value_epi_epi[i][j],2)/this.expected_value_epi_epi[i][j];
                }else{
                    norm_weight[i][j]= Math.pow(this.epi_epi[i][j] - this.expected_value_epi_epi[i][j],2)/this.expected_value_epi_epi[i][j];
                }
                
            }
        }
        return norm_weight;
    }
    //fungsi menghitung expected value dari setiap pasang asam amino
    //untuk setiap pasang asam amino dihitung jumlah total epi-epi nonepi-epi dan nonepi-nonepi
    //hitung frekuensi seluruh epi epi, nonepi-epi, dan nonepi-nonepi
    
    // expected value pasangan [i][j] = marginal column frequency * marginal row frequency / total sampel
    //setiap pasang yang mungkin
    private void calculateExpectedValue(){
        //hitung frequensi epi-epi pada setiap pasang asam amino
        double[][] sumofe_ne_b = new double[13][13];
        double sumofepi = 0;
        double sumofepinonepi = 0;
        double sumofnonepinonepi = 0;
        for(int i=0;i<13;i++){
            for(int j=0; j<13;j++){
                // hitung per baris
                sumofe_ne_b[i][j] = epi_epi[i][j]+epi_nonepi[i][j] +nonepi_nonepi[i][j];
                //total kolom e-e
                sumofepi += epi_epi[i][j];
                //total kolom e-ne
                sumofepinonepi += epi_nonepi[i][j];
                //total kolom ne-ne
                sumofnonepinonepi += nonepi_nonepi[i][j];
            }
        }
        double total = sumofepi+  sumofepinonepi+ sumofnonepinonepi;
         
        expected_value_epi_epi = new double[13][13];
        double sum =0;
        for (int i=0;i<13;i++){
            for(int j=0;j<13;j++){
                expected_value_epi_epi[i][j]= (sumofe_ne_b[i][j]* sumofepi)/total;
                
                //System.out.print( expected_value_epi_epi[i][j]+"\t");
            }
            //System.out.print("\n");
        }
        //gabungkan pasangan AB dengan BA 
        double[][] n_e_v = new double[13][13];
        for (int i=0;i<13;i++){
            for(int j=0;j<13;j++){
                if(i!=j){
                    n_e_v[i][j]= expected_value_epi_epi[i][j]+ expected_value_epi_epi[j][i];
                }else{
                    n_e_v[i][j]= expected_value_epi_epi[i][j];
                }
                    
                
                
                //System.out.print( expected_value_epi_epi[i][j]+"\t");
            }
            //System.out.print("\n");
        }
        expected_value_epi_epi = n_e_v;
        
        
    }
    // berdasarkan paper 
    
    
    //bobot dihitung berdasarkan probabilitas epi/semua pada setiap pasangan asam amino
    public double[][] calculateWight(){
        System.out.println("calculate weight");
        if (this.getEpi_epi().length != 0){
        double[][]w_epi_epi = new double[13][13];
        for (int i=0;i<13;i++){
            for(int j=0;j<13;j++){
                w_epi_epi[i][j]= (epi_epi[i][j]+epi_epi[j][i])/(epi_epi[i][j]+nonepi_nonepi[i][j]+epi_epi[j][i]+nonepi_nonepi[j][i]);
                System.out.print( w_epi_epi[i][j]+"\t");
            }
            System.out.print("\n");
        }
        return w_epi_epi;
        }else
            return null;
        
    }
    public void concatAA(){
        String []aa_letter = {"A","B","X","J","C","G","H","U","M","Z","P","O","Y"};
        for(int i=0;i<aa_letter.length;i++)
            for(int j=0;j<aa_letter.length;j++)
                System.out.println(aa_letter[i]+aa_letter[j]);
    }
    public void writeWeightAAPairToFile(String file, double[][] weight){
     BufferedWriter bw = null;
	FileWriter fw = null;
	try {
            fw = new FileWriter(new File("src/dataset/graph/exposedaa/weight_"+file+".txt"));
            bw = new BufferedWriter(fw);
            //tulis disini 
          
            for(int i=0;i<weight.length;i++){
                for(int j=0;j<weight[0].length;j++){
                
                String content = weight[i][j]+"\t";
                bw.write(content);
                }
             String newline = "\n" ;
             bw.write(newline);
            }
		
	} catch (IOException e) {
		e.printStackTrace();
	}
 }
    public void readFromFileAAPByCategory(String path, String filename){
        FileReader fr;
        BufferedReader br;
        String dataline;
        try{
            if(new File(path+ filename).exists()){
                fr = new FileReader(new File(path+ filename));
                br = new BufferedReader(fr);            
                dataline=br.readLine();
                int count =0;
                while(dataline != null){
                    if (count<13){     
                        System.out.println("read epi epi");
                        for(int i=0;i<epi_epi.length;i++){
                            String[]splitstr = dataline.split("\t");
                            for (int j=0;j<epi_epi[0].length;j++)
                                epi_epi[i][j] = Double.parseDouble(splitstr[j]);
                            dataline = br.readLine();
                            count +=1;
                        }
                    } else if(count>=13 & count<26){
                        System.out.println("read epi nonepi");
                         for(int i=0;i<epi_nonepi.length;i++){
                            String[]splitstr = dataline.split("\t");
                            for (int j=0;j<epi_nonepi[0].length;j++)
                                epi_nonepi[i][j] = Double.parseDouble(splitstr[j]);
                            dataline = br.readLine();
                            count +=1;
                        }
                    }else{
                        System.out.println("read nonepi nonepi");
                         for(int i=0;i<nonepi_nonepi.length;i++){
                            String[]splitstr = dataline.split("\t");
                            for (int j=0;j<nonepi_nonepi[0].length;j++)
                                nonepi_nonepi[i][j] = Double.parseDouble(splitstr[j]);
                            dataline = br.readLine();
                            count +=1;
                        }
                    }
                        
                    
                    
                }
            }
        } catch(Exception e){
        
        }
    }
    public void savetofileAAPByCategory(String path, String filename){
        BufferedWriter bw;
        FileWriter fw;
            try{
              //  fw = new FileWriter(new File("src/dataset/graphexposedaa/"+eoc.getPdbid()+eoc.getChain()+"_"+t_distance+"_"+t_asa+".txt"),true);
              fw = new FileWriter(new File(path+filename+".txt"),true);
                bw = new BufferedWriter(fw);
                String s = null;
                for(int i=0;i<epi_epi.length;i++){
                    for (int j=0;j<epi_epi[0].length;j++)
                        s += epi_epi[i][j]+ "\t";
                    s+="\n";
                    }
                bw.write(s+"\n");
                s= null;
                for(int i=0;i<epi_nonepi.length;i++){
                    for (int j=0;j<epi_nonepi[0].length;j++)
                        s += epi_nonepi[i][j]+ "\t";
                    s+="\n";
                    }
                    bw.write(s+"\n");
                
                s= null;
                for(int i=0;i<nonepi_nonepi.length;i++){
                    for (int j=0;j<nonepi_nonepi[0].length;j++)
                        s += nonepi_nonepi[i][j]+ "\t";
                    s+="\n";
                    }
                    bw.write(s+"\n");
                
                bw.close();
                fw.close();    
            }
                
                
            catch(Exception except){
                   System.out.println(except.getMessage());
               }
            
        
        
        
    }
     public void savetofileExpectedValue(String path, String filename){
        BufferedWriter bw;
        FileWriter fw;
            try{
              //  fw = new FileWriter(new File("src/dataset/graphexposedaa/"+eoc.getPdbid()+eoc.getChain()+"_"+t_distance+"_"+t_asa+".txt"),true);
              fw = new FileWriter(new File(path+filename+".txt"),true);
                bw = new BufferedWriter(fw);
                String s = null;
                for(int i=0;i<this.expected_value_epi_epi.length;i++){
                    for (int j=0;j<expected_value_epi_epi[0].length;j++)
                        s += expected_value_epi_epi [i][j]+ "\t";
                    s+="\n";
                    }
                bw.write(s+"\n");
                
                
                bw.close();
                fw.close();    
            }
                
                
            catch(Exception except){
                   System.out.println(except.getMessage());
               }
            
        
        
        
    }
    public void printAAPByCategory(){
        System.out.println("epi-epi");
        for(int i=0;i<epi_epi.length;i++){
            for (int j=0;j<epi_epi[0].length;j++)
                System.out.print(epi_epi[i][j]+ "\t");
            System.out.println();
        }
        System.out.println("epi-nonepi");
        for(int i=0;i<epi_nonepi.length;i++){
            for (int j=0;j<epi_nonepi[0].length;j++)
                System.out.print(epi_nonepi[i][j]+ "\t");
            System.out.println();
        }
        System.out.println("nonepi-nonepi");
        for(int i=0;i<nonepi_nonepi.length;i++){
            for (int j=0;j<nonepi_nonepi[0].length;j++)
                System.out.print(nonepi_nonepi[i][j]+ "\t");
            System.out.println();
        }
    }
    public void printExpectedValue(){
        System.out.println("expected value");
        for(int i=0;i<expected_value_epi_epi.length;i++){
            for (int j=0;j<expected_value_epi_epi[0].length;j++)
                System.out.print(expected_value_epi_epi[i][j]+ "\t");
            System.out.println();
        }
        
    }

    public double[][] getEpi_epi() {
        return epi_epi;
    }

    public double[][] getEpi_nonepi() {
        return epi_nonepi;
    }

    public double[][] getNonepi_nonepi() {
        return nonepi_nonepi;
    }
    
    
    
   
    
    public void calculateFrequenciesAminoacidPair(ArrayList<ExposedOnComplex> eocList, double t_distance){
        //mapping aa to fs
        HashMap<String,Integer> map_aa = new HashMap();
        map_aa.put("A",0);
        map_aa.put("R",1);
        map_aa.put("K", 1);
        map_aa.put("N",2);
        map_aa.put("Q",2);
        map_aa.put("D",3);
        map_aa.put("E",3);
        map_aa.put("C",4);
        map_aa.put("G",5);
        map_aa.put("H",6);
        map_aa.put("I",7);
        map_aa.put("L",7);
        map_aa.put("V",7);
        map_aa.put("M",8);
        map_aa.put("F",9);
        map_aa.put("W",9);
        map_aa.put("P",10);
        map_aa.put("S",11);
        map_aa.put("T",11);
        map_aa.put("Y",12);
        
        
        //buat graph dari residu terekspose
        //dari setiap node yang terbentuk berapa jumlah pasangan e-e, e-ne, ne-ne
        for (int i =0;i<eocList.size();i++){
            ExposedOnComplex eoc = eocList.get(i);
            ArrayList<Aminoacid> listAA = eoc.getListAA();
            //System.out.println("jumlah residu terekspose: "+e_aa.size());
            for (int j=0;j<listAA.size();j++){
                for (int k =1;k<listAA.size()-1;k++ ){
                    Aminoacid aa_1 = listAA.get(j);
                    Aminoacid aa_2 = listAA.get(k);
                    String res_1 = aa_1.getresidueName1();
                    String res_2 = aa_2.getresidueName1();
                    if (aa_1.isAsEpitope()& aa_2.isAsEpitope())
                        epi_epi[map_aa.get(res_1)][map_aa.get(res_2)]+=1;
                    else if (aa_1.isAsEpitope()== false & aa_2.isAsEpitope()== false)
                        nonepi_nonepi[map_aa.get(res_1)][map_aa.get(res_2)] +=1;
                    else if ((aa_1.isAsEpitope()== false & aa_2.isAsEpitope()== true)|(aa_1.isAsEpitope()== true & aa_2.isAsEpitope()== false))
                        epi_nonepi[map_aa.get(res_1)][map_aa.get(res_2)] +=1;
                    
                }
            }
            
            
        }
    }
    /*
     1 2 3 4 ==> 1-2, 1-3,1-4, 2-3, 2-4, 3-4
    */
    
     private int letterTonum(String letter)
      {
    if ((letter.equals("ALA")) || (letter.equals("A"))) {
      return 0;
    }
    if ((letter.equals("ARG")) || (letter.equals("R"))) {
      return 1;
    }
    if ((letter.equals("ASN")) || (letter.equals("N"))) {
      return 2;
    }
    if ((letter.equals("ASP")) || (letter.equals("D"))) {
      return 3;
    }
    if ((letter.equals("CYS")) || (letter.equals("C"))) {
      return 4;
    }
    if ((letter.equals("GLN")) || (letter.equals("Q"))) {
      return 5;
    }
    if ((letter.equals("GLU")) || (letter.equals("E"))) {
      return 6;
    }
    if ((letter.equals("GLY")) || (letter.equals("G"))) {
      return 7;
    }
    if ((letter.equals("HIS")) || (letter.equals("H"))) {
      return 8;
    }
    if ((letter.equals("ILE")) || (letter.equals("I"))) {
      return 9;
    }
    if ((letter.equals("LEU")) || (letter.equals("L"))) {
      return 10;
    }
    if ((letter.equals("LYS")) || (letter.equals("K"))) {
      return 11;
    }
    if ((letter.equals("MET")) || (letter.equals("M"))) {
      return 12;
    }
    if ((letter.equals("PHE")) || (letter.equals("F"))) {
      return 13;
    }
    if ((letter.equals("PRO")) || (letter.equals("P"))) {
      return 14;
    }
    if ((letter.equals("SER")) || (letter.equals("S"))) {
      return 15;
    }
    if ((letter.equals("THR")) || (letter.equals("T"))) {
      return 16;
    }
    if ((letter.equals("TRP")) || (letter.equals("W"))) {
      return 17;
    }
    if ((letter.equals("TYR")) || (letter.equals("Y"))) {
      return 18;
    }
    if ((letter.equals("VAL")) || (letter.equals("V"))) {
      return 19;
    }
    return -1;
  }

    private ArrayList<ExposedOnComplex> extractFeatureOfExposedAntigen(String file_list, double t_asa) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
