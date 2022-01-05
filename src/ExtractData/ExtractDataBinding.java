/*
 * input: 
 */
package ExtractData;

import ComplexStructure.Aminoacid;
import ComplexStructure.Epitope;
import ExtractDataTraining.Benchmark;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author dell
 */
public class ExtractDataBinding {
    private ArrayList<String> pdbid_list;
    private ArrayList<String> chain_list;

    public ArrayList<String> getPdbid_list() {
        return pdbid_list;
    }

    public ArrayList<String> getChain_list() {
        return chain_list;
    }
     public ArrayList<Contact> loadListContact(String list_file){
     //pdbid
     //ag chain
     //ab chain
     ArrayList<Contact> contactList=null;
     try{
         if (new File("src/dataset/"+list_file).exists()) {
             contactList = new ArrayList();
            String dataline = null;
            FileReader fr = null;
            fr = new FileReader(new File("src/dataset/"+list_file));
                BufferedReader br = new BufferedReader(fr);
                dataline = br.readLine();                
                while (dataline != null)
                {
                    System.out.println(dataline);
                    
                    String[] splitdl = dataline.split(":"); 
                    if(splitdl.length==3){
                        Contact c= new Contact();
                        c.pdbid= splitdl[0].trim();
                        c.Ag = splitdl[1].trim();      
                        c.Ab = splitdl[2].trim();
                        contactList.add(c);
                    }
                    dataline = br.readLine();
                }
                br.close();
                
        }
     }catch(Exception e){
     
     }finally{
         return contactList;
     }
     
     
     }
      public void createListDataSetCV(ArrayList<String> pdbtrain,ArrayList<String>chaintrain, ArrayList<String> pdbtestset, ArrayList<String> chaintestset, int part){
        
        int sz = this.pdbid_list.size();
        int set = sz/10;//90/10=9 0;9;18;27
        for(int i=0;i<set;i++){
            int start =i*set;
            int stop = (i+1)*set -1;
           // System.out.println(start+"-->"+stop);
            if(set==part){
            //masukkan ke testset
            for(int j=start;j<=stop;j++){
                pdbtestset.add(this.pdbid_list.get(j));
                chaintestset.add(this.chain_list.get(j));
            }
            }else{
            //masukkan ke trainset
                for(int j=start;j<=stop;j++){
                pdbtrain.add(this.pdbid_list.get(j));
                chaintrain.add(this.chain_list.get(j));
                }
            }
        }
        
                
        
    }
    //load list file
    //ekstract epitop from epi file
     //
    public void loadListFile(String list_file){
    pdbid_list = new ArrayList();
    chain_list = new ArrayList();
    
    //baca file list
    //masukkan dalam list pdbid dan chain
    try{
        if (new File("src/dataset/"+list_file).exists()) {
            String dataline = null;
            FileReader fr = null;
            fr = new FileReader(new File("src/dataset/"+list_file));
                BufferedReader br = new BufferedReader(fr);
                dataline = br.readLine();                
                while (dataline != null)
                {
                   // System.out.println(dataline);
                    String[] splitdl = dataline.split(":"); 
                    if(splitdl.length==2){
                        pdbid_list.add(splitdl[0].trim());
                        chain_list.add(splitdl[1].trim());      
                    }
                    dataline = br.readLine();
                }
                br.close();
                printList(pdbid_list);
                printList(chain_list);
        }
    }catch(IOException e){
            e.printStackTrace();
        }
    
    }
    public void printList(ArrayList<String> list){
        for(String s:list){
            System.out.println(s);
        }
    }
    public ArrayList<String> getResidueIDBinding(ArrayList<Binding> list){
        ArrayList<String> listId= new ArrayList();
        for(Binding b:list){
            listId.add(b.residueId);
        }
        return listId;
    }
    public void printBinding(ArrayList<Binding> list){
        for(Binding b:list){
            System.out.println(b.chain+"\t"+b.residueId+"\t"+b.residue);
        }
    }
    public ArrayList<Epitope> loadEpiFromList(ArrayList<String>pdbidList, ArrayList<String>chainList){
        ArrayList<Epitope> eList = new ArrayList();
    
    for(int i=0;i<pdbidList.size();i++){
        String pdbid = pdbidList.get(i);
        String chain = chainList.get(i);
        ArrayList<Binding> bindlist = this.loadDataBindingFromEpiFile(pdbid, chain);
        System.out.println("jumlah epi "+pdbid+ bindlist.size());
        Epitope e = new Epitope();
        e.setChain(chain);
        e.setPdbid(pdbid);
        ArrayList<String> idxepi = new ArrayList();
        for (Binding b: bindlist){
            idxepi.add(b.residueId);
        }
        e.setIdxEpi(idxepi);
        eList.add(e);
    }
    return eList;
    
        
    }
    public ArrayList<Epitope> loadEpiFromBinding(String filename){
    this.loadListFile(filename);
    ArrayList<Epitope> eList = new ArrayList();
    
    for(int i=0;i<this.pdbid_list.size();i++){
        String pdbid = this.pdbid_list.get(i);
        String chain = this.chain_list.get(i);
        ArrayList<Binding> bindlist = this.loadDataBindingFromEpiFile(pdbid, chain);
        System.out.println("jumlah epi "+pdbid+ bindlist.size());
        Epitope e = new Epitope();
        e.setChain(chain);
        e.setPdbid(pdbid);
        ArrayList<String> idxepi = new ArrayList();
        for (Binding b: bindlist){
            idxepi.add(b.residueId);
        }
        e.setIdxEpi(idxepi);
        eList.add(e);
    }
    return eList;
    
    }
    
    public ArrayList<Epitope> loadBindingFromContact(String filename){
        ArrayList<Epitope> list_epi = new ArrayList();
        ArrayList<Contact> contact=  this.loadListContact(filename);
        
        System.out.println("jumlah contact: "+ contact.size());
        for(int i=0;i<contact.size();i++){
          //  System.out.println(contact.get(i).pdbid +":"+contact.get(i).Ag+":"+contact.get(i).Ab.toString());
            Epitope epi = extractBindingFromContact(contact.get(i));
            System.out.println(epi.getPdbid()+":"+epi.getChain());
            list_epi.add(epi);        
        }
        return list_epi;
    }
    public Epitope extractBindingFromContact(Contact c){
        Epitope e =null;
        ArrayList<String> residuId= new ArrayList();
        ArrayList<String> residu = new ArrayList();
        String dataline = null;
        FileReader fr = null;
        String[] ab = c.Ab.split(",");
       // String name= "src/dataset/contact/"+c.pdbid+"_201908101008_contacts.tbl";
        String name= "src/dataset/seppa3/test/contact/"+c.pdbid+"_201909071958_contacts.tbl";
       // String name= "src/dataset/andersen/contact/"+c.pdbid+".tbl";
       // String name= "src/dataset/seppa3/train/binding/"+c.pdbid+".tbl";
        
        try{
            if (new File(name).exists()) {
                System.out.println(name);
                fr = new FileReader(new File(name));
                BufferedReader br = new BufferedReader(fr);
                dataline = br.readLine();
                boolean state=false;
                Aminoacid aa= new Aminoacid();
                while (dataline != null)
                {                      
                   if(dataline.startsWith(c.Ag.trim())){
                     //create the binding base on chain
                        //   System.out.println(dataline);
                        String sAb = dataline.substring(14,15).trim();
                       // System.out.println(sAb);
                        boolean AbState = false;
                        for(int i=0;i<ab.length;i++){
                            if (sAb.equalsIgnoreCase(ab[i])){
                            AbState = true;
                            }
                        }
                        
                       if(AbState){
                           String id = dataline.substring(3, 7).trim();
                            if(!residuId.contains(id)){
                            residu.add( aa.LetterToletter(dataline.substring(8, 12).trim()));
                            residuId.add(id);
                            System.out.println(dataline);
                       }
                       }
                       
                       
                       
                        
                    }
                   
                dataline = br.readLine();
                }
                br.close();
             e = new Epitope();
             e.setIdxEpi(residuId);
             e.setSequenceAsEpitope(residu);
             e.setPdbid(c.pdbid);
             e.setChain(c.Ag);
            // System.out.println(e.getPdbid()+":"+e.getChain()+":"+e.getIdxEpi());
             return e;
            }
            
        }catch(IOException exc){
            exc.printStackTrace();
            System.out.println(c.pdbid+"not exist");
        }finally{
            return e;
        }
    }
    public ArrayList<Binding> loadDataBindingFromEpilist(String pdbid, String chain){
    //load data dari epilist
    //buat obyek bindListnya
    //
    return null;
    }
    public ArrayList<Binding> loadDataBindingFromEpiFile(String pdbid,String chain){
        ArrayList<Binding> bindList= new ArrayList();
        String dataline = null;
        FileReader fr = null;
        String name= "src/dataset/ren2015/tes/epi/"+pdbid.toUpperCase()+".epi";
        System.out.println(name);
        try{
            if (new File(name).exists()) {
               // System.out.println(pdbid + chain);
                fr = new FileReader(new File(name));
                BufferedReader br = new BufferedReader(fr);
                dataline = br.readLine();
                System.out.println(dataline);
                while (dataline != null)
                {
                     
                    String[] rawDt = dataline.split(" ");
                    if(chain.equalsIgnoreCase(rawDt[0])){
                        Binding bind= new Binding(); 
                        bind.chain= chain;
                        bind.residueId = rawDt[1];
                        bind.residue = rawDt[2];
                        bindList.add(bind);
                    }
                        dataline = br.readLine();
                }
                br.close();
            }
                }catch(IOException e){
            e.printStackTrace();
            System.out.println(pdbid+"not exist");
        }finally{
            return bindList;
        }
    }
     public ArrayList<Binding> loadDataBindingFromFile(String pdbid,String chain){
        
        ArrayList<Binding> bindList= new ArrayList();
        String dataline = null;
        FileReader fr = null;
       // String name= "src/dataset/test_set_binding_residues/"+pdbid+"_201901220118_binding_residues.tbl";
     //   String name= "src/dataset/training_set_binding_residues/"+pdbid+"_201905130838_binding_residues.tbl";//max distance
        //String name= "src/dataset/training_set_binding_residu_45/"+pdbid+"_201905160829_binding_residues.tbl";
        //String name= "src/dataset/training_set_br_4/"+pdbid+"_201905221216_binding_residues.tbl";
        String name= "src/dataset/test_zhang2011/"+pdbid+"_201908031442_binding_residues.tbl";
        try{
            if (new File(name).exists()) {
               // System.out.println(pdbid + chain);
                fr = new FileReader(new File(name));
                BufferedReader br = new BufferedReader(fr);
                dataline = br.readLine();
                
                while (dataline != null)
                {
                    
                    //8 sd 12 : AA
                   if(dataline.startsWith(chain.trim())){
                     //   System.out.println(dataline);
                       Binding bind= new Binding(); 
                        bind.chain= chain;
                        bind.residue= dataline.substring(8, 12).trim();
                        bind.residueId = dataline.substring(3, 7).trim();
                        bindList.add(bind);
                    }
                   
                dataline = br.readLine();
                }
                br.close();
             printBinding(bindList);
            }
            
        }catch(IOException e){
            e.printStackTrace();
            System.out.println(pdbid+"not exist");
        }finally{
            return bindList;
        }
       
    }
}
class Binding{
String chain;
String residueId;
String residue;

}
class Contact{
    String pdbid;
    String Ag;
    String Ab;
}




