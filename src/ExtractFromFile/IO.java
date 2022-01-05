/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ExtractFromFile;

import ComplexFeature.AAIndex;
import ComplexFeature.AAIndex1;
import ComplexStructure.*;
import ExtractData.CreateComplexStructure;
import bpred.AnalysisOnSurfaceResidue;
import java.io.*;
import java.util.ArrayList;
import java.util.Collection;

/**
 *
 * @author OpenGress
 */
public class IO {
 public void getSegmentFromFile(){
     
 }
 
 
 public void writeAAIndexToFile(AAIndex aai){
     BufferedWriter bw = null;
	FileWriter fw = null;
	try {
            fw = new FileWriter(new File("src/datatraining/I_AAIndex.txt"));
            bw = new BufferedWriter(fw);
            //tulis disini 
          ArrayList<double[]>idx = aai.getAaindex();
            for(int i=0;i<idx.size();i++){
                
                String content = this.getStringFromDA(idx.get(i))+"\n";
                bw.write(content);
            }
		
	} catch (IOException e) {
		e.printStackTrace();
	}
 }
 public String getStringFromDA(double[] arr){
     String data="";
     String val;
     for (int i=0;i<arr.length;i++){
         if(arr[i]== Double.NaN){
             val = "NA";
         }else{
             val = String.valueOf(arr[i]);
         }
         data=data+val+"\t";
     }
     return data;
 }
 public String DAtoString(double[] da){
     String s="";
     for (int i=0;i<da.length;i++){
         s= s+ String.valueOf(da[i])+"\t";
     }
     return s;
 }
 public void writeArrayListStringToFile(String full_path, ArrayList<String> arraylist){
    BufferedWriter bw;
    FileWriter fw;
    try {
            fw = new FileWriter(new File(full_path),true);
            bw = new BufferedWriter(fw);
            for(String s: arraylist)
                bw.write(s);
            
		
	} catch (IOException e) {
		e.printStackTrace();
	}
 }
 public void writeDatasetExposedResidue(){
     BufferedWriter bw;
    FileWriter fw;
    try {
            fw = new FileWriter(new File("src/dataset/exposedResidu.txt"));
            bw = new BufferedWriter(fw);
            
            
		
	} catch (IOException e) {
		e.printStackTrace();
	}
 }
 public void writeLogOddsToFile(ArrayList<String> logodddataset, String filename){
    
    BufferedWriter bw;
    FileWriter fw;
    try {
            fw = new FileWriter(new File("src/datatraining/"+filename+".txt"));
            bw = new BufferedWriter(fw);
            
            for(String logoddstring: logodddataset){
                //System.out.print(aanum[i]+"\t");
                  bw.write(logoddstring);
            }
		
	} catch (IOException e) {
		e.printStackTrace();
	}
 }
public void writeAAIndexOfEpiToFile(){
    String filename = "aaiofEpi.txt";
    BufferedWriter bw;
    FileWriter fw;
    try {
            fw = new FileWriter(new File("src/datatraining/aaiofEpi.txt"));
            bw = new BufferedWriter(fw);
            //tulis disini 
            
            
            AAIndex idxAA = this.extractAAIndexFromFile();
            idxAA.composeAAI();
            
            //ambil daftar aa dari epitope
            AnalysisOnSurfaceResidue osr = new AnalysisOnSurfaceResidue();
            int[] aanum = osr.getAANum();
            
            String stringAAI;
            for(int i=0;i<aanum.length;i++){
                //System.out.print(aanum[i]+"\t");
                stringAAI= this.DAtoString(idxAA.getAAIndexAA(aanum[i]))+"\t"+ 1 +"\n";
                 bw.write(stringAAI);
            }
		
	} catch (IOException e) {
		e.printStackTrace();
	}
    
}
public void writeAAIndexOfNonEpiToFile(){
    
    BufferedWriter bw;
    FileWriter fw;
    try {
            fw = new FileWriter(new File("src/datatraining/aaiofExpResNonEpi.txt"));
            bw = new BufferedWriter(fw);
            //tulis disini 
            
            AAIndex idxAA = this.extractAAIndexFromFile();
            idxAA.composeAAI();
            
            //ambil daftar aa dari epitope
            AnalysisOnSurfaceResidue osr = new AnalysisOnSurfaceResidue();
            osr.doAnalysisOnSurfaceResidu();
            int[] aanumER = osr.getAANumListOfNonEpiExpRes();
            
            String stringAAI;
            for(int i=0;i<aanumER.length;i++){
                //System.out.print(aanum[i]+"\t");
                double[] aai=idxAA.getAAIndexAA(aanumER[i]);
                if(aai!=null){
                    stringAAI= this.DAtoString(aai)+"\t"+ 0 +"\n";
                    bw.write(stringAAI);
                }
                
            }
		
	} catch (IOException e) {
		e.printStackTrace();
	}
}
 public AAIndex extractAAIndexFromFile(){
     String filename = "aaindex1.txt";
        FileReader fr;
        BufferedReader br;
        String dataline;
        AAIndex aaidx= new AAIndex();
        ArrayList<String> H= new ArrayList();
        ArrayList<String> D= new ArrayList();
        try{
            if(new File("src/datatraining/"+ filename).exists()){
            fr = new FileReader(new File("src/datatraining/"+ filename));
            br = new BufferedReader(fr);            
            dataline=br.readLine();
            while(dataline != null){
                if(dataline.startsWith("H")){
                    H.add(dataline.substring(2, dataline.length()).trim());
                }
                if(dataline.startsWith("D")){
                    D.add(dataline.substring(2, dataline.length()).trim());
                }
                if(dataline.startsWith("I")){
                dataline = br.readLine()+br.readLine();
                //dataline = br.readLine();
                dataline=dataline.trim();
                String[] spl = dataline.split(" ");
                double []d = new double[20];
                int j=0;
                for(int i=0;i<spl.length;i++){
                    String s =spl[i].trim();
                    if(s.length()>0){
                        if(s.matches("NA")){
                            d[j]= Double.NaN;
                          //  System.out.println("NA");
                        }else{
                            d[j]=Double.parseDouble(s);
                        }
                        
                        j+=1;
                    }
                 
                }
                //System.out.println(d[19]);
                aaidx.addAAIndexComponent(d);
                
                
                /* for(int i=0;i<spl.length;i++){
                    System.out.println(spl[i]);
                }
                */
              //  aaidx.addAAIndexComponent(idx);
                }
                dataline=br.readLine();
            }
          //  System.out.println(H.get(0)+"sd"+H.get(H.size()-1));
            }
        }catch(Exception e){
           // e.printStackTrace();
        }finally{
            /*
            for(String header: H){
                System.out.println(header);
            }
            for(String desc: D){
                System.out.println(desc);
            }*/
            return aaidx;
        }
 
 }   
    public String[] getSplitStringFromString(String dl){
        String[] splitdl= new String[20];
        String copydl = dl;
        int count=0;
        while(copydl.length()>0 &&count<20){
            int idx = dl.indexOf(" ");
            System.out.println(copydl);
            System.out.println(copydl.length()+";"+idx);
            if(idx<copydl.length()){
                splitdl[count] = copydl.substring(0,idx).trim();
                count+=1;
                
            }
            if(idx>0&&copydl.length()>idx+1){
                    copydl=copydl.substring(idx,copydl.length()).trim();
            }
            
            
            
        }
        return splitdl;
        
    }
    
 public ArrayList<String> listFilesForFolder(final File folder) {
     ArrayList<String> files = new ArrayList();
    for (final File fileEntry : folder.listFiles()) {
        if (!fileEntry.isDirectory()) {
            //listFilesForFolder(fileEntry);
            String fn = fileEntry.getName();
            files.add(fn);
        } 
    }
    return files;
}
 
 public void extractASA(){
    AnalysisOnSurfaceResidue aosr = new AnalysisOnSurfaceResidue ();
    aosr.doAnalysisOnSurfaceResidu();
    
     ArrayList<Chain> c = aosr.getComplexchainlist();
     
     BufferedWriter bw = null;
     FileWriter fw = null;
     String dataline;
     try{
         //fw = new FileWriter(new File("src/dataset/seqepiexposed.txt"),true);
         fw = new FileWriter(new File("src/dataset/datasetASA.txt"),true);
         bw = new BufferedWriter(fw);
        // System.out.println("jumlah chain: "+ c.size());
         for(int i=0;i<c.size();i++){
             Chain cn = c.get(i);
             
            // dataline = cn.getSequenceOfChain()+"\t"+cn.getSequenceOfEpitopeStatus()+"\t"+cn.getExposesResidue()+"\n";
              dataline = cn.getSequenceOfChain()+"\t"+cn.getSequenceOfEpitopeStatus()+"\t"+cn.getExposesResidue()+"\t"+cn.getSecondaryStructure()+"\n";
             //System.out.print(dataline);
             bw.write(dataline);
                        
         }
     }catch(Exception e){
                    System.out.println("can't write");
     }
     
 }
 public void writeExposedorEpiFromContact(String list_file){
      CreateComplexStructure ccs = new CreateComplexStructure();
     // ccs.loadStructureFromContact(list_file);
     ccs.loadStructureFromListFile(list_file);
      ArrayList<Complex> complexList= ccs.getComplexlist();
      BufferedWriter bw = null;
      FileWriter fw = null;
      String dataline;
      try{
        //  fw = new FileWriter(new File("src/dataset/seqepiexposedSeppa3Train.txt"),true);
          fw = new FileWriter(new File("src/dataset/seqepiexposedRen2015.txt"),true);
          bw = new BufferedWriter(fw);
          int totalEpi=0;
          int totalExposed =0;
          for(int i=0;i< complexList.size();i++){
          ArrayList<Aminoacid> listAA= complexList.get(i).getAminoacidVector();
          char chainId = complexList.get(i).getAminoacidVector().get(0).getchainID();
          //System.out.println("chainId"+ chainId);
          Chain c = complexList.get(i).getChainByChianID(chainId);
          dataline = c.getAminoacidSequenceofChain()+" "+ c.getSequenceOfEpitopeStatus()+" "+c.getExposedResidue(0.01)+"\n";
         //pada ren 2015 residu sebagai exposed jika ASA >0
        //  dataline = c.getAminoacidSequenceofChain()+" "+ c.getSequenceOfEpitopeStatus()+" "+c.getAsabasedExposesResidue(0)+"\n";
          //System.out.println(dataline);
          totalEpi += c.getNumEpi();
          totalExposed += c.getNumExposed();
          bw.write(dataline);
          }
          System.out.println("Epi: "+ totalEpi + "Exposed: "+ totalExposed);
           }catch(Exception e){
           
           }
      
      
      
 }
 public void writeExposedOrEpitopeOfSequence(){
     AnalysisOnSurfaceResidue aosr = new AnalysisOnSurfaceResidue ();
     aosr.doAnalysisOnSurfaceResidu();
    
     ArrayList<Chain> c = aosr.getComplexchainlist();
     
     BufferedWriter bw = null;
    FileWriter fw = null;
    String dataline;
     try{
         fw = new FileWriter(new File("src/dataset/seqepiexposed.txt"),true);
        // fw = new FileWriter(new File("src/dataset/seqepiexposedsecstruc.txt"),true);
         bw = new BufferedWriter(fw);
        // System.out.println("jumlah chain: "+ c.size());
         for(int i=0;i<c.size();i++){
             Chain cn = c.get(i);
             
             dataline = cn.getSequenceOfChain()+"\t"+cn.getSequenceOfEpitopeStatus()+"\t"+cn.getExposesResidue()+"\n";
           //   dataline = cn.getSequenceOfChain()+"\t"+cn.getSequenceOfEpitopeStatus()+"\t"+cn.getExposesResidue()+"\t"+cn.getSecondaryStructure()+"\n";
             //System.out.print(dataline);
             bw.write(dataline);
                        
         }
     }catch(Exception e){
                    System.out.println("can't write");
     }
     
 }
 
public void extractASAOfEpiFromFile(){
    File folder = new File("src/datatraining/");
    BufferedWriter bw = null;
    FileWriter fw = null;
    FileReader fr = null;
    String dataline;
    try{
        fw = new FileWriter(new File("src/datatraining/dataasa1/asaepitope.txt"),true);
        bw = new BufferedWriter(fw);
        File[] listOfFiles = folder.listFiles();
        for (int i = 0; i < listOfFiles.length; i++) {
            if (listOfFiles[i].isFile()) {
            
                String filename = listOfFiles[i].getName();
               // System.out.println(filename);
                try{
                    fr = new FileReader(new File("src/datatraining/" + filename));
                    BufferedReader br = new BufferedReader(fr);
                    dataline = br.readLine();
                    while (dataline != null){
                        if(dataline.substring(4, 8).trim().equals("true")){
                    //tulis ke file
                        //System.out.println(dataline);
                        bw.write(dataline+"\n");
                        }
                        dataline = br.readLine();
                    }
                    br.close();
               
          
                }catch(Exception e){
                    System.out.println("file can't read");
                }
                
            } 
        }
         bw.close();
        
    }catch(Exception e){
        System.out.print("can't extract epitope asa from file");
    }
    
}
public void AppendEpitopeAntigenSequenceToFile(String s){
    BufferedWriter bw = null;
    FileWriter fw = null;
    try{
        File file = new File("src/datatraining/epitopofexposedresidu.txt");
        if (!file.exists()) {
            file.createNewFile();
        }
        fw = new FileWriter(file.getAbsoluteFile(), true);
	bw = new BufferedWriter(fw);
	bw.write(s);
			
    }catch(IOException e){
        e.printStackTrace();
    }finally {
		try {
                    if (bw != null)
			bw.close();

                    if (fw != null)
			fw.close();
		} catch (IOException ex) {
			ex.printStackTrace();

		}

	}
}
public void writeAnalysisDistanceOfExposedResidue(){
    AnalysisOnSurfaceResidue aosr = new AnalysisOnSurfaceResidue();
    ArrayList<NearestNeighboor> analnn= aosr.distanceAnalysisOnSurfaceResidu();
    BufferedWriter bw = null;
	FileWriter fw = null;
	try {
            fw = new FileWriter(new File("src/datatraining/" + "distanceAnalysis.txt"));
            bw = new BufferedWriter(fw);
            for(int i=0;i<analnn.size();i++){
                
                String content =analnn.get(i).printNN();
                bw.write(content);
            }
            
        }catch (IOException e) {
		e.printStackTrace();
	} finally {
		try {
                    if (bw != null)
			bw.close();

                    if (fw != null)
			fw.close();
		} catch (IOException ex) {
			ex.printStackTrace();

		}

	}
}
 
    public void writeAnalysisParamASAOfExpossedResidue(ArrayList<Aminoacid> exposedResidue, String filename){
        
        BufferedWriter bw = null;
	FileWriter fw = null;
	try {
            fw = new FileWriter(new File("src/datatraining/" + filename));
            bw = new BufferedWriter(fw);
            //tulis disini 
          
            for(int i=0;i<exposedResidue.size();i++){
                //ASA RSA CN FSE HSE QSE EFSE SFSE
                //String ASA = exposedResidue.get(i).printParamASA();
                String content =String.valueOf(exposedResidue.get(i).letterToNum())+"\t"+ exposedResidue.get(i).getpASA().printParamASA()+"\t"+ exposedResidue.get(i).isAsEpitope()+"\n";
              //  System.out.println(content);
                bw.write(content);
            }
		
	} catch (IOException e) {
		e.printStackTrace();
	} finally {
		try {
                    if (bw != null)
			bw.close();

                    if (fw != null)
			fw.close();
		} catch (IOException ex) {
			ex.printStackTrace();

		}

	}
    }
    public void copyPdbChainFromPdb(String pdbid, String chain, String chain_num){
        BufferedWriter bw = null;
	FileWriter fw = null;
        FileReader fr = null;
        String dataline;
        try{
            // fw = new FileWriter(new File("src/datatraining/pdbchain/" + pdbid +"_"+ chain+chain_num+".pdb"));
             fw = new FileWriter(new File("src/dataset/seppa3/test/nonglyco/chain/" + pdbid +"_"+ chain+chain_num+".pdb"));
             bw = new BufferedWriter(fw);
             String file = "src/dataset/test_seppa3/non_glyco/" + pdbid+".pdb";
             System.out.println(file);
             if (new File("src/dataset/test_seppa3/non_glyco/" + pdbid+".pdb").exists()) {
                fr = new FileReader(new File("src/dataset/test_seppa3/non_glyco/" + pdbid+".pdb"));
                BufferedReader br = new BufferedReader(fr);
                dataline = br.readLine();
                while (dataline != null){
                    
                    if (dataline.startsWith("ATOM")){
                            //tulis dataline ke file
                        String c = dataline.substring(21, 22).trim();
                        if(c.equals(chain)){
                            bw.write(dataline+"\n");
                        }
                        
                    }else {
                            bw.write(dataline+"\n");
                    }
                    
                        dataline= br.readLine();
                    }
                br.close();
                bw.close();
                
                
            }
        }catch(Exception e){
        System.out.println(pdbid +"not found");
        }
        
    }
    public void copyPdbChainFromPdb(String pdbid, String chain){
        BufferedWriter bw = null;
	FileWriter fw = null;
        FileReader fr = null;
        String dataline;
        try{
            // fw = new FileWriter(new File("src/dataset/seppa3/train/chain/" + pdbid +"_"+ chain+".pdb"));
            fw = new FileWriter(new File("src/dataset/seppa3/test/glyco/chain/"+ pdbid +"_"+ chain+".pdb"));
            
             bw = new BufferedWriter(fw);
             if (new File("src/dataset/test_seppa3/non_glyco/"+ pdbid+".pdb").exists()) {
                 System.out.println("File ditemukan"+ pdbid);
                fr = new FileReader(new File("src/dataset/seppa3/test/glyco/" + pdbid+".pdb"));
                BufferedReader br = new BufferedReader(fr);
                dataline = br.readLine();
                while (dataline != null){
                    
                    if (dataline.startsWith("ATOM")){
                            //tulis dataline ke file
                        String c = dataline.substring(21, 22).trim();
                       // System.out.println("chain:"+c);
                        if(c.equalsIgnoreCase(chain)){
                            
                            bw.write(dataline+"\n");
                        }
                        
                    }else {
                           // bw.write(dataline+"\n");
                    }
                    
                        dataline= br.readLine();
                    }
                br.close();
                bw.close();
                
                
            }
        }catch(Exception e){
            System.out.println("File tidak ditemukan"+ pdbid);
        }
        
    }
    public void writeChainToFile(){
        //ambil daftar pdbid dan chain dari file
        //
    }
    
    public void writePdbidChainToFile(){
        ArrayList<Epitope> epilist = this.getEpitopeListFromFile();
        String filename = "daftarpdbid.txt";
        BufferedWriter bw = null;
	FileWriter fw = null;
	try {
            fw = new FileWriter(new File("src/datatraining/" + filename));
            bw = new BufferedWriter(fw);
            //tulis disini 
            for(int i=0;i<epilist.size();i++){
                String content = epilist.get(i).getPdbid()+"\t"+epilist.get(i).getChain()+"\n";
                bw.write(content);
            }
		
	} catch (IOException e) {
		e.printStackTrace();
	} finally {
		try {
                    if (bw != null)
			bw.close();

                    if (fw != null)
			fw.close();
		} catch (IOException ex) {
			ex.printStackTrace();

		}

	}

    }
    public double[][]getBlosumMatrixFromFile(){
        String filename = "blosumforpseudocount.txt";
        String dataline=null;
        FileReader fr = null;
        double[][]blosum = new double[20][20];
        try{
            if(new File("src/datatraining/"+ filename).exists()){
            fr = new FileReader(new File("src/datatraining/"+ filename));
            BufferedReader br = new BufferedReader(fr);
            dataline = br.readLine();
            dataline=br.readLine();
            int line=0;
            while(dataline != null && line<20){
                this.setBlosumAtPositionK(blosum, dataline, line);
                String[] splitdl = dataline.split(";");        
                //System.out.println(idx);
                for(int i=0;i<splitdl.length;i++){
                    blosum[line][i]=Double.parseDouble(splitdl[i]);
                }
                dataline= br.readLine();
                line=line+1;
                
            }
            br.close();
            
            }
        }catch(Exception e){
            e.printStackTrace();
        }finally{
            return blosum;
        }
        
    }
    public void setBlosumAtPositionK(double[][]blosum, String dl,int idx){
        String[] splitdl = dl.split(";");        
        System.out.println(idx);
        for(int i=0;i<splitdl.length;i++){
            blosum[idx][i]=Double.parseDouble(splitdl[i]);
        }
    }
    public void printBlosum(double[][] blosum){
        for(int i=0;i<blosum.length;i++){
            for (int j=0;j<blosum[i].length;j++){
                System.out.print(blosum[i][j]+"\t");
            }
                System.out.println();
        }
    }
    public void getSegment(String seq, String epi){
        ArrayList<String> substring= new ArrayList();
        String[] splitepi= epi.split("0");
             for(int i=0;i<splitepi.length;i++){
                 if(splitepi[i].startsWith("1")){
                     substring.add(splitepi[i]);
                  }
             }
             for(int j=0;j<substring.size();j++){
                 
             }
    }
    public void getEpitopeSegmentFromEpilist(){
        ArrayList<Epitope> epilist = this.getEpitopeListFromFile();
        //baca file banchmark
        String filename = "epitopesegmen.txt";
        String dataline = null;
        FileWriter fw=null;
        BufferedWriter bw=null;
        try{
            fw = new FileWriter(new File("src/datatraining/" + filename));
            bw = new BufferedWriter(fw);
            //tulis disini 
            for(int i=0;i<epilist.size();i++){
                String content = epilist.get(i).getSegmentOfEpi();
                System.out.println(content);
                bw.write(content);
            }
            
        }catch(Exception e){
        System.out.print(e.getMessage());
        }finally {
		try {
                    if (bw != null)
			bw.close();

                    if (fw != null)
			fw.close();
		} catch (IOException ex) {
			ex.printStackTrace();

		}

	}
    }
    public void extractFeature(String filename){
     BufferedWriter bw = null;
     FileReader fr = null;
     String dataline;
     try{
         if (new File("src/datatraining/andersen/" + filename).exists()) {
        fr = new FileReader(new File("src/datatraining/andersen/" + filename));
        BufferedReader br = new BufferedReader(fr);
        dataline = br.readLine();
        while (dataline != null){
            String[] linestring= dataline.split(",");
            System.out.println("jumlah linestring:"+ linestring.length);
        }
            }
     }catch(Exception e){
     
     }
 }
    public ArrayList<Epitope> getEpitopeListFromFile(){
        ArrayList<Epitope> line = new ArrayList();
        //String filename = "benchmarksequencedata.txt";
       // String filename = "bsd_addition.txt";
        String filename = "epitopelist.txt";
        String dataline = null;
        FileReader fr = null;
    try
    {
      if (new File("src/data/" + filename).exists()) {
        fr = new FileReader(new File("src/data/" + filename));
        BufferedReader br = new BufferedReader(fr);
        dataline = br.readLine();
        while (dataline != null)
        {
           // System.out.println(dataline);
             String[] linestring= dataline.split(" ");
             Epitope e = new Epitope(linestring);
//             System.out.print("jumlah epitope"+e.getPdbid()+ ": "+e.getIdxsequenceAsEpitope().size());
             
            
            String pdbid = linestring[0];
            String chain = linestring[1];
            String sequence = linestring[3];
            String asepitope = linestring[4];
            ArrayList<Integer> idxepi = new ArrayList();
            ArrayList<String> seq = new ArrayList();
            for(int i=0;i<asepitope.length();i++){
                String epi = String.valueOf(asepitope.charAt(i));
                if(epi.equalsIgnoreCase("1")){
                String sepi = String.valueOf(sequence.charAt(i));
                idxepi.add(Integer.valueOf(i));
                seq.add(sepi);
                }
                
            }
            e.setIdxsequenceAsEpitope(idxepi);
            e.setSequenceAsEpitope(seq);
            line.add(e);
         //   System.out.println(asepitope);
             
             dataline = br.readLine();
        }
        br.close();
      }
     // System.out.println("jumlah data terbaca: "+ line.size());
      return line;
    }
          
    
    catch (IOException e)
    {
      e.printStackTrace();
      return line;
    }
   
    }
}
