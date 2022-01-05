/*Protusion Index dihitung menggunakan PSAIA
 * dibaca dari file
 * disimpan dalam parameter PI
 * 
 */
package ExtractDataTraining;

import ComplexCalculate.UnboundDSSP;
import ComplexFeature.PSAI;
import ComplexStructure.Aminoacid;
import ComplexStructure.Complex;
import ComplexStructure.Epitope;
import ExtractData.CreateComplexStructure;
import java.io.*;
import java.util.ArrayList;

/**
 *
 * @author OpenGress
 */
public class ExtractProtusionIndex {
    private ArrayList<PSAI> psailist;
    
    //read from file
    public ExtractProtusionIndex(){
        psailist= new ArrayList();
       
    }

   
   //gabungkan dengan fungsi pembentukan complex 
    
    public void doExtraction(double threshold){
        ExtractDataBenchmark edb = new ExtractDataBenchmark();
        ArrayList<Epitope> epilist = edb.loadBenchmarkData("epitopelist.txt");
      for(Epitope e:epilist){
            String pdbid = e.getPdbid();
            String chain = e.getChain();
            String chain_num = e.getNo_chain();
            Complex complex = new Complex(pdbid, chain);
            complex.getChainByChianID(chain.charAt(0)).setAAAsEpitope(e);
            UnboundDSSP unbound_dssp = new UnboundDSSP(complex, chain);
            ArrayList<Aminoacid> listAA=complex.getAminoacidVector();
            //read pi from file
            ArrayList<PSAI> psailist =this.readPIFromFile(pdbid);
            ArrayList<PSAI> chainpsailist =this.copyPIOfChain(chain, psailist);
         //   System.out.println("jumlah aa:"+listAA.size()+" jumlah psai"+ chainpsailist.size());
            BufferedWriter bw;
            FileWriter fw;
            try {
                String createdfile;
                fw = new FileWriter(new File("src/datatraining/FEATURES/PSAI_"+threshold+".arff"),true);
                bw = new BufferedWriter(fw);
               
                
                //untuk setiap asam amino yang memenuhi syarat maka temukan psainya
               // System.out.println("proses pencocokan");
                for(int count=0;count<listAA.size();count++){
                    Aminoacid a = listAA.get(count);
                    boolean exposedStatus = a.getRsaTien2013()>threshold ;
                    boolean exposedStatusOfEpi = a.isAsEpitope()&& exposedStatus;
                    String aaId = a.getPdbid()+"_"+a.getchainID()+"_"+a.getresidueID()+"_"+a.letterToNum();
                    String psaiid = chainpsailist.get(count).getPsaiId();
                   // System.out.println((exposedStatus && aaId.equals(psaiid))||(exposedStatusOfEpi&& aaId.equals(psaiid)));
                    String dataline="";
                    if((exposedStatus ||exposedStatusOfEpi)){
                        dataline = psaiid + "," + chainpsailist.get(count).printPSAIParam() +a.isAsEpitope();
                    }
                    bw.write(dataline+"\n");
                    
                }
                bw.close();
            }catch(Exception except){
               // System.out.println(except.getMessage());
                except.printStackTrace();
            }
            
      }
        
    }
    public void printArray2D(double[] array){
        for(int i=0;i<array.length;i++){
            System.out.print(array[i]+",");
        }
        System.out.println();
            
    }
    public ArrayList<PSAI> copyPIOfChain(String chain, ArrayList<PSAI>psailist){
        ArrayList<PSAI> chainPSAI = new ArrayList();
        for(PSAI psai :psailist){
            if(psai.getPsaiId().charAt(5)==chain.charAt(0)){
             chainPSAI.add(psai);
            }
        }
        return chainPSAI;
    }
    //
    public String printEpitopeState(Complex c){
        String sequence = c.printAA()+" ";
        for(int i=0;i<c.getAminoacidVector().size();i++){
            boolean epi = c.getAminoacidVector().get(i).isAsEpitope();
            if (epi==true){
                sequence += "1";
            }else{
                sequence += "0";
            }
            
        }
        return sequence;
    }
    public void setEpitopeFromPSAI(){
        //String list_file ="list_antigen_andersen.txt";
        String list_file ="list_test_zhang2011.txt";
        CreateComplexStructure cs = new CreateComplexStructure(list_file);
        cs.loadchainlist(list_file);
        ArrayList<String> pdblist= cs.getPdblist();
        for(int i=0;i<pdblist.size();i++){
            //baca psai bound
            String pathBound="src/data/struktur_psai/bound/";
            String akhiranBound ="_201804140420_bound.tbl";
            ArrayList<PSAI> listBound = readPIFromFile(pdblist.get(i),pathBound,akhiranBound);
            //baca psai unbound
            String pathUnbound="src/data/struktur_psai/unbound/";
            String akhiranUnbound ="_201804140420_unbound.tbl";
            ArrayList<PSAI> listUnBound = readPIFromFile(pdblist.get(i),pathUnbound,akhiranUnbound);
        //hitung delta
            ArrayList<Double> deltaASA = countDeltaASA(listUnBound,listBound);
            cs.getComplexlist().get(i).setEpitopeState(deltaASA);
        //tandai epitope atau bukan  
            String epistate= this.printEpitopeState(cs.getComplexlist().get(i));
           // System.out.println(pdblist.get(i)+ " "+ cs.getComplexlist().get(i).getChainIDs()+" "+epistate+"\n");
        }
        
        
        
    }
    
    public ArrayList<Double> countDeltaASA(ArrayList<PSAI> unbound, ArrayList<PSAI> bound ){
    if(unbound.size()==bound.size()){
        ArrayList<Double> delta = new ArrayList();
        for(int i=0;i<bound.size();i++){
            String id= unbound.get(i).getPsaiId();
            String []idc = unbound.get(i).getPsaiId().split("_");
            String res = idc[idc.length-1];
            double unb= calculateRSAByTienMSA(res,unbound.get(i).getASAFromPSAI());
            double b = calculateRSAByTienMSA(res,bound.get(i).getASAFromPSAI());
            delta.add(unb-b);
        }
        return delta;
    }else{
        System.out.println("ukuran data tidak sama");
        return null;
    }
    }
    public ArrayList<PSAI>readPIFromFile(String pdbid, String path,String akhiran){
        ArrayList<PSAI> listPSAI = new ArrayList();
        String fullpath = path+pdbid+"_"+akhiran+".tbl";
        System.out.println(fullpath);
        FileReader fr;
        BufferedReader br;
        String dataline;
        try{
            if(new File(fullpath).exists()){
                fr = new FileReader(new File(fullpath));
            br = new BufferedReader(fr);            
            dataline=br.readLine().trim();
            while(dataline != null){
                //System.out.println(dataline);
                if(!dataline.startsWith("chain id")){
                        dataline = br.readLine().trim();
                }else{
                    dataline = br.readLine().trim();
                   
                    break;
                }                
            }   
            //
           while(dataline!=null){
             
               //String aaId = aa.getPdbid()+"_"+aa.getchainID()+"_"+aa.getresidueID()+"_"+aa.letterToNum();
              
                  PSAI psai = new PSAI();
                  String id = pdbid+"_"+dataline.substring(0,1)+"_"+ dataline.substring(81, 93).trim()+"_"+dataline.substring(94,104).trim();
                  psai.setPsaiId(id);
                  double[] psaiparam = new double[23];
                  int beginIdx = 104;
                  //  System.out.println(dataline.length()+": "+dataline);
                  for(int i=0;i<23;i++){
                        String s=dataline.substring(beginIdx,beginIdx+15).trim();
                        psaiparam[i]=Double.parseDouble(s);
                      //  System.out.println(psaiparam[i]);
                        beginIdx +=16;
                  }
                 //   System.out.println();
                    //this.printArray2D(psaiparam);
                  psai.setPsaiParam(psaiparam);
                  listPSAI.add(psai);
                  dataline = br.readLine().trim();
              }
          
            
            br.close();
              }
            }catch(Exception e){
           // System.out.println(filename);
            System.out.println(e.getMessage());
        }finally{
           // return copyPSAI;
            return listPSAI;
        }
        
    }
    public ArrayList<PSAI> readPIFromFile(String pdbid){
       // System.out.println("\n"+filename);
        ArrayList<PSAI> listPSAI = new ArrayList();
        //ArrayList<PSAI> copyPSAI= new ArrayList();
        //String filename = pdbid+"_201804140420_bound.tbl";
        //String filename = pdbid+"_2019082109_bound.tbl";//1aho_201908210935_bound.tbl
      //  String filename = pdbid+"_201909082102_bound.tbl";//1aho_201908210935_bound.tbl
      //  String filename = pdbid+"_201908061618_bound.tbl";
      //  String filename = pdbid+"_2019082109_bound.tbl";
       // String filename = pdbid+"_201909082102_bound.tbl";
        String filename = pdbid+"_201909082102_bound.tbl";
        //1a2y_201804140420_bound
        //System.out.println(filename);
        FileReader fr;
        BufferedReader br;
        String dataline;
        try{
            if(new File("src/dataset/seppa3/test/psaia/"+filename).exists()){
                   fr= new FileReader(new File("src/dataset/seppa3/test/psaia/"+filename))  ;
//   fr= new FileReader(new File("src/data/struktur_psai/test_seppa3/"+filename))  ;
            // fr = new FileReader(new File("src/dataset/ren2015/tes/binding/"+ filename));
           // fr = new FileReader(new File("src/dataset/seppa3/train/psaia/"+ filename));
           // fr = new FileReader(new File("src/data/struktur_psai/bound/"+filename));
            br = new BufferedReader(fr);            
            dataline=br.readLine().trim();
            while(dataline != null){
               // System.out.println(dataline);
                if(!dataline.startsWith("chain id")){
                        dataline = br.readLine().trim();
                }else{
                    dataline = br.readLine().trim();
                   
                    break;
                }                
            }   
            //
           while(dataline!=null){
             
               //String aaId = aa.getPdbid()+"_"+aa.getchainID()+"_"+aa.getresidueID()+"_"+aa.letterToNum();
              
                  PSAI psai = new PSAI();
                  String id = pdbid+"_"+dataline.substring(0,1)+"_"+ dataline.substring(81, 93).trim()+"_"+this.letterToNum(dataline.substring(94,104).trim());
                  psai.setPsaiId(id);
                  double[] psaiparam = new double[23];
                  int beginIdx = 104;
                  //  System.out.println(dataline.length()+": "+dataline);
                  for(int i=0;i<23;i++){
                        String s=dataline.substring(beginIdx,beginIdx+15).trim();
                        psaiparam[i]=Double.parseDouble(s);
                      //  System.out.println(psaiparam[i]);
                        beginIdx +=16;
                  }
                 //   System.out.println();
                //    this.printArray2D(psaiparam);
                  psai.setPsaiParam(psaiparam);
                  listPSAI.add(psai);
                  dataline = br.readLine().trim();
              }
          
            
            br.close();
              }
            //  copyPSAI = this.copyPIOfChain(chain, listPSAI);
        }catch(Exception e){
           // System.out.println(filename);
            System.out.println(e.getMessage());
        }finally{
           // return copyPSAI;
            return listPSAI;
        }
    }
     private double calculateRSAByTienMSA(String residueName, double area)
  {
      
    if ((residueName.equals("ALA")) || (residueName.equals("A"))) {
      return area*100 / 129.0D;
    }
    if ((residueName.equals("ARG")) || (residueName.equals("R"))) {
      return area*100 / 274.0D;
    }
    if ((residueName.equals("ASN")) || (residueName.equals("N"))) {
      return area*100 / 195.0D;
    }
    if ((residueName.equals("ASP")) || (residueName.equals("D"))) {
      return area*100 / 193.0D;
    }
    if ((residueName.equals("CYS")) || (residueName.equals("C"))) {
      return area*100 / 167.0D;
    }
    if ((residueName.equals("GLN")) || (residueName.equals("Q"))) {
      return area*100 / 223.0D;
    }
    if ((residueName.equals("GLU")) || (residueName.equals("E"))) {
      return area*100 / 225.0D;
    }
    if ((residueName.equals("GLY")) || (residueName.equals("G"))) {
      return area*100 / 104.0D;
    }
    if ((residueName.equals("HIS")) || (residueName.equals("H"))) {
      return area*100 / 224.0D;
    }
    if ((residueName.equals("ILE")) || (residueName.equals("I"))) {
      return area*100 / 197.0D;
    }
    if ((residueName.equals("LEU")) || (residueName.equals("L"))) {
      return area*100 / 201.0D;
    }
    if ((residueName.equals("LYS")) || (residueName.equals("K"))) {
      return area*100 / 236.0D;
    }
    if ((residueName.equals("MET")) || (residueName.equals("M"))) {
      return area*100 / 224.0D;
    }
    if ((residueName.equals("PHE")) || (residueName.equals("F"))) {
      return area*100 / 240.0D;
    }
    if ((residueName.equals("PRO")) || (residueName.equals("P"))) {
      return area*100 / 159.0D;
    }
    if ((residueName.equals("SER")) || (residueName.equals("S"))) {
      return area*100 / 155.0D;
    }
    if ((residueName.equals("THR")) || (residueName.equals("T"))) {
      return area*100 / 172.0D;
    }
    if ((residueName.equals("TRP")) || (residueName.equals("W"))) {
      return area*100 / 285.0D;
    }
    if ((residueName.equals("TYR")) || (residueName.equals("Y"))) {
      return area*100 / 263.0D;
    }
    if ((residueName.equals("VAL")) || (residueName.equals("V"))) {
      return area*100 / 174.0D;
    }
    return -1.0D;
  }
    public  int letterToNum(String letter){
  
      int num=-1;
  if (letter.equals("A")||letter.equals("ALA")) {
      num= 0;
    }
  if(letter.equals("R")||letter.equals("ARG")){
      num= 1;
  } 
  if(letter.equals("N")||letter.equals("ASN")){
      num= 2;
  }
  if(letter.equals("D")||letter.equals("ASP")){
     num= 3;
  }
  if(letter.equals("C")||letter.equals("CYS")){
      num= 4;
  }
  if(letter.equals("Q")||letter.equals("GLN")){
      num= 5;
  }
  if(letter.equals("E")||letter.equals("GLU")){
      num= 6;
  }
  if(letter.equals("G")||letter.equals("GLY")){
      num= 7;
  }
  if(letter.equals("H")||letter.equals("HIS")){
      num= 8;
  }
  if(letter.equals("I")||letter.equals("ILE")){
      num= 9;
  }
  if(letter.equals("L")||letter.equals("LEU")){
      num= 10;
  }
  if(letter.equals("K")||letter.equals("LYS")){
      num= 11;
  }
  if(letter.equals("M")||letter.equals("MET")){
      num= 12;
  }
  if(letter.equals("F")||letter.equals("PHE")){
      num= 13 ;
  }
  if(letter.equals("P")||letter.equals("PRO")){
      num= 14;
  }
  if(letter.equals("S")||letter.equals("SER")){
      num= 15;
  }
  if(letter.equals("T")||letter.equals("THR")){
      num= 16 ;
  }
  if(letter.equals("W")||letter.equals("TRP")){
      num= 17;
  }
  if(letter.equals("Y")||letter.equals("TYR")){
      num= 18;
  }
  if(letter.equals("V")||letter.equals("VAL")){
      num= 19;
  }
  return num;
     
  
  }
}
