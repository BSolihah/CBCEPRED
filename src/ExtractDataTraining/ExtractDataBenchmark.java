/*
 * ekstract data benchmarksequencedata.txt
 */
package ExtractDataTraining;

import ComplexStructure.Epitope;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author OpenGress
 */
public class ExtractDataBenchmark {
    //pdbid,chainid,chainnumber,sequence,epitopeposition
    //simpan setiap item data dalam obyek kelas databenchmark
    public ExtractDataBenchmark(){
        //load from file
    }
    public Epitope getBenchmarkData(String filename){
        String dataline = null;
        FileReader fr = null;
        Epitope epi=null;
        try{
            if (new File("src/datatest/" + filename).exists()) {
                fr = new FileReader(new File("src/datatest/" + filename));
                BufferedReader br = new BufferedReader(fr);
                dataline = br.readLine();
                while (dataline != null)
                {
                // System.out.println(dataline);
                    String[] linestring= dataline.split(" ");
                 //  System.out.println("panjang linestring"+ linestring.length);
                 //  System.out.println( linestring[0]+"\t"+linestring[1]+"\t"+linestring[2]+"\t"+linestring[3]+"\t"+linestring[4]);
                    epi = new Epitope(linestring);
                   //  System.out.println(e.getPdbid()+"\t"+e.getSequence()+"\t"+e.getEpipos());
                    
                    }
                br.close();}
        }catch(IOException e){
            e.printStackTrace();
        }finally{
            return epi;
        }

    }
    public ArrayList<Epitope> loadDatatest(String filename){
        ArrayList<Epitope> listE= new ArrayList();
        String dataline = null;
        FileReader fr = null;
        try{
            if (new File("src/data/" + filename).exists()) {
                fr = new FileReader(new File("src/data/" + filename));
                BufferedReader br = new BufferedReader(fr);
                dataline = br.readLine();
                while (dataline != null)
                {
                // System.out.println(dataline);
                    String[] linestring= dataline.split(" ");
                 //  System.out.println("panjang linestring"+ linestring.length);
                 //  System.out.println( linestring[0]+"\t"+linestring[1]+"\t"+linestring[2]+"\t"+linestring[3]+"\t"+linestring[4]);
                    Epitope e = new Epitope(linestring);
                   //  System.out.println(e.getPdbid()+"\t"+e.getSequence()+"\t"+e.getEpipos());
                    listE.add(e);


                    dataline = br.readLine();
                    }
                br.close();}
        }catch(IOException e){
            e.printStackTrace();
        }finally{
            return listE;
        }

    }
    
    public ArrayList<Epitope> loadBenchmarkData(String filename){
        ArrayList<Epitope> listE= new ArrayList();
        //String filename = "epitopelist.txt";
        //file seqepiexposed berisi sequence, posisi epi dan posisi exposed
        String dataline = null;
        FileReader fr = null;
        try{
            if (new File("src/data/" + filename).exists()) {
                fr = new FileReader(new File("src/data/" + filename));
                BufferedReader br = new BufferedReader(fr);
                dataline = br.readLine();
                while (dataline != null)
                {
                // System.out.println(dataline);
                    String[] linestring= dataline.split(" ");
                 //  System.out.println("panjang linestring"+ linestring.length);
                 //  System.out.println( linestring[0]+"\t"+linestring[1]+"\t"+linestring[2]+"\t"+linestring[3]+"\t"+linestring[4]);
                    Epitope e = new Epitope(linestring);
                   //  System.out.println(e.getPdbid()+"\t"+e.getSequence()+"\t"+e.getEpipos());
                    listE.add(e);


                    dataline = br.readLine();
                    }
                br.close();}
        }catch(IOException e){
            e.printStackTrace();
        }finally{
            return listE;
        }

    }
    public ArrayList<Epitope> loadBenchmarkDataTest(String filename){
        ArrayList<Epitope> listE= new ArrayList();
        //String filename = "epitopelisttestset.txt";
        
        //file seqepiexposed berisi sequence, posisi epi dan posisi exposed
        String dataline = null;
        FileReader fr = null;
        try{
            if (new File("src/data/" + filename).exists()) {
                fr = new FileReader(new File("src/data/" + filename));
                BufferedReader br = new BufferedReader(fr);
                dataline = br.readLine();
                while (dataline != null)
                {
                // System.out.println(dataline);
                    String[] linestring= dataline.split(" ");
                 //  System.out.println("panjang linestring"+ linestring.length);
                 //  System.out.println( linestring[0]+"\t"+linestring[1]+"\t"+linestring[2]+"\t"+linestring[3]+"\t"+linestring[4]);
                    Epitope e = new Epitope(linestring);
                     //System.out.println(e.getPdbid()+"\n"+e.getSequence()+"\n"+e.getEpipos());
                    listE.add(e);


                    dataline = br.readLine();
                    }
                br.close();}
        }catch(IOException e){
            e.printStackTrace();
        }finally{
            return listE;
        }

    }
  
  
    
    
}
