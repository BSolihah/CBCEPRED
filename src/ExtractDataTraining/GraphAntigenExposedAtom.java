/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ExtractDataTraining;

import ComplexFeature.AminoAcidPairByCN;
import ComplexStructure.Atom;
import ComplexStructure.Complex;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author dell
 * digenerate dari atom terekspose
 */
public class GraphAntigenExposedAtom {
    
    
    
    public void generateGraph(Complex c, String path_in_res,String path_in_atom, String filenameres,String filenameatom, double t_distance, double [][]edge_weight, String path_out){        
        
        System.out.println("complex: "+c.getComplexName()+ c.getAminoacidVector().size());
        //baca dari file exposed atom dan exposed residu
        ArrayList<Atom> list_atom_res = loadExposedAtomResiduLevel(path_in_res, filenameres);
        System.out.println("list_atom_res: "+ list_atom_res.get(0).getX()+" "+list_atom_res.get(0).getY()+" "+list_atom_res.get(0).getZ());
        ArrayList<Atom> list_atom_atom = loadExposedAtom(path_in_atom, filenameatom);
        System.out.println("list_atom_atom: "+ list_atom_atom.get(0).getX()+" "+list_atom_atom.get(0).getY()+" "+list_atom_atom.get(0).getZ()+" "+list_atom_atom.get(0).getAtomName() );
        concatResiduInfo(list_atom_res,list_atom_atom);
        //System.out.println(list_atom_res.get(0).getAtomName());
        ArrayList<String> pirlist = generateEdgeList(list_atom_res, t_distance,edge_weight);
        //ArrayList<String> pirlist = generateEdgeCNWList(c,list_atom_res, t_distance,edge_weight);
        //ArrayList<String> pirlist = generateEdgeCNList(c,list_atom_res, t_distance);
        //hitung jarak antar atom terekspose jika memenuhi treshold maka definisikan edge
        
        if(pirlist.size()>0){
            System.out.println("pirlist terbentuk: "+ pirlist.size());
            BufferedWriter bw;
            FileWriter fw;
            try{
              //  fw = new FileWriter(new File("src/dataset/graphexposedaa/"+eoc.getPdbid()+eoc.getChain()+"_"+t_distance+"_"+t_asa+".txt"),true);
              String fn =path_out + filenameatom +"_"+t_distance+".txt";
              System.out.println(fn);
              fw = new FileWriter(new File(fn),true);
                bw = new BufferedWriter(fw);
                for(String s: pirlist){
                   // String s = pir.getPairResidu();
                   
                    bw.write(s+"\n");
                }
                bw.close();
                fw.close();
            }     
            catch(Exception except){
                   System.out.println(except.getMessage());
               }
        }
        
        
        }
    public ArrayList<String> generateEdgeList(ArrayList<Atom> list_atom_res, double t_distance,double [][]edge_weight){
        Set<String> pirlist = new HashSet<String>();
        ArrayList<String> edgelist = new ArrayList();
        for(int i=0;i<list_atom_res.size();i++){
            for(int j=0;j<list_atom_res.size();j++){                
                if(i!=j){
                    Atom a1 = list_atom_res.get(i);
                    Atom a2 = list_atom_res.get(j);
                    if(!a1.getAtomName().equalsIgnoreCase("O")||!a2.getAtomName().equalsIgnoreCase("O")){
                        if (!a1.getResidueID().equalsIgnoreCase(a2.getResidueID())){
                        //hitung jarak
                            double d = distanceBetweenPoint3D(a1,a2);
                           // System.out.println(i+","+j+ ":"+d);
                            if(d <= t_distance){                                
                                //jika memenuhi treshold dan tidak dari residu yang sama maka definisikan edge
                                double w = edge_weight[letterToNum(a1.getResidueName())][letterToNum(a2.getResidueName())];
                                String edge = null;
                                try{
                                    int id1 = Integer.parseInt(a1.getResidueID());
                                    int id2 = Integer.parseInt(a2.getResidueID());

                                    if(id1<=id2){
                                        edge = a1.getResidueID()+"\t"+a2.getResidueID()+"\t"+ String.valueOf(w);
                                    }else{
                                            edge = a2.getResidueID()+"\t"+ a1.getResidueID()+"\t"+String.valueOf(w);
                                    }
                                }catch(Exception e){
                                    System.out.println("except");
                                    if(a1.getResidueID().length()<a2.getResidueID().length()){
                                         edge = a1.getResidueID()+"\t"+a2.getResidueID()+"\t"+ String.valueOf(w);
                                    }else{
                                         edge = a2.getResidueID()+"\t"+a1.getResidueID()+"\t"+ String.valueOf(w);
                                    }                                    
                                }
                               // System.out.println(edge);
                                pirlist.add(edge);
                            }

                        }
                    }
                    }
                        
                }
                
                }
        ArrayList<String> list = new ArrayList<>(pirlist);
        return list;
    }
    public ArrayList<String> generateEdgeCNList(Complex c, ArrayList<Atom> list_atom_res, double t_distance){
        AminoAcidPairByCN aapbycn = new AminoAcidPairByCN(c);
        System.out.println("inisiasi aapbycn");
        Set<String> pirlist = new HashSet<String>();
        ArrayList<String> edgelist = new ArrayList();
        for(int i=0;i<list_atom_res.size();i++){
            for(int j=0;j<list_atom_res.size();j++){                
                if(i!=j){
                    Atom a1 = list_atom_res.get(i);
                    Atom a2 = list_atom_res.get(j);
                    if(!a1.getAtomName().equalsIgnoreCase("O")||!a2.getAtomName().equalsIgnoreCase("O")){
                        if (!a1.getResidueID().equalsIgnoreCase(a2.getResidueID())){
                        //hitung jarak
                            double d = distanceBetweenPoint3D(a1,a2);
                           // System.out.println(i+","+j+ ":"+d);
                            if(d <= t_distance){                                
                                //jika memenuhi treshold dan tidak dari residu yang sama maka definisikan edge
                                //double w = edge_weight[letterToNum(a1.getResidueName())][letterToNum(a2.getResidueName())];
                                String resid_a1 = a1.getResidueID();
                                String resid_a2 = a2.getResidueID();
                                //System.out.println(resid_a1+ "-"+resid_a2);
                                //int cn = aapbycn.calcPairByCN(resid_a1, resid_a2);
                                double w = Double.valueOf(aapbycn.calcPairByCN(resid_a1, resid_a2));
                                //System.out.println("cn: "+cn);
                                String edge = null;
                                try{
                                    int id1 = Integer.parseInt(a1.getResidueID());
                                    int id2 = Integer.parseInt(a2.getResidueID());

                                    if(id1<=id2){
                                       
                                        edge = a1.getResidueID()+"\t"+a2.getResidueID()+"\t"+ String.valueOf(w);
                                    }else{
                                            edge = a2.getResidueID()+"\t"+ a1.getResidueID()+"\t"+String.valueOf(w);
                                    }
                                }catch(Exception e){
                                    System.out.println("except in ordering res id");
                                    if(a1.getResidueID().length()<a2.getResidueID().length()){
                                         edge = a1.getResidueID()+"\t"+a2.getResidueID()+"\t"+ String.valueOf(w);
                                    }else{
                                         edge = a2.getResidueID()+"\t"+a1.getResidueID()+"\t"+ String.valueOf(w);
                                    }                                    
                                }
                               // System.out.println(edge);
                                pirlist.add(edge);
                            }

                        }
                    }
                    }
                        
                }
                
                }
        ArrayList<String> list = new ArrayList<>(pirlist);
        return list;
    }
    public ArrayList<String> generateEdgeCNWList(Complex c, ArrayList<Atom> list_atom_res, double t_distance,double [][]edge_weight){
        AminoAcidPairByCN aapbycn = new AminoAcidPairByCN(c);
        System.out.println("inisiasi aapbycn");
        Set<String> pirlist = new HashSet<String>();
        ArrayList<String> edgelist = new ArrayList();
        for(int i=0;i<list_atom_res.size();i++){
            for(int j=0;j<list_atom_res.size();j++){                
                if(i!=j){
                    Atom a1 = list_atom_res.get(i);
                    Atom a2 = list_atom_res.get(j);
                    if(!a1.getAtomName().equalsIgnoreCase("O")||!a2.getAtomName().equalsIgnoreCase("O")){
                        if (!a1.getResidueID().equalsIgnoreCase(a2.getResidueID())){
                        //hitung jarak
                            double d = distanceBetweenPoint3D(a1,a2);
                           // System.out.println(i+","+j+ ":"+d);
                            if(d <= t_distance){                                
                                //jika memenuhi treshold dan tidak dari residu yang sama maka definisikan edge
                                double waap = edge_weight[letterToNum(a1.getResidueName())][letterToNum(a2.getResidueName())];
                                String resid_a1 = a1.getResidueID();
                                String resid_a2 = a2.getResidueID();
                                //System.out.println(resid_a1+ "-"+resid_a2);
                                //int cn = aapbycn.calcPairByCN(resid_a1, resid_a2);
                                double w = waap*Double.valueOf(aapbycn.calcPairByCN(resid_a1, resid_a2));
                                //System.out.println("cn: "+cn);
                                String edge = null;
                                try{
                                    int id1 = Integer.parseInt(a1.getResidueID());
                                    int id2 = Integer.parseInt(a2.getResidueID());

                                    if(id1<=id2){
                                       
                                        edge = a1.getResidueID()+"\t"+a2.getResidueID()+"\t"+ String.valueOf(w);
                                    }else{
                                            edge = a2.getResidueID()+"\t"+ a1.getResidueID()+"\t"+String.valueOf(w);
                                    }
                                }catch(Exception e){
                                    System.out.println("except in ordering res id");
                                    if(a1.getResidueID().length()<a2.getResidueID().length()){
                                         edge = a1.getResidueID()+"\t"+a2.getResidueID()+"\t"+ String.valueOf(w);
                                    }else{
                                         edge = a2.getResidueID()+"\t"+a1.getResidueID()+"\t"+ String.valueOf(w);
                                    }                                    
                                }
                               // System.out.println(edge);
                                pirlist.add(edge);
                            }

                        }
                    }
                    }
                        
                }
                
                }
        ArrayList<String> list = new ArrayList<>(pirlist);
        return list;
    }
    
    public  int letterToNum(String residuName){
        String letter = residuName;
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
    
    private double distanceBetweenPoint3D(Atom pa, Atom pb)
    {
      double distance = 0.0D;
      double x = pa.getX() - pb.getX();
      double y = pa.getY() - pb.getY();
      double z = pa.getZ() - pb.getZ();
      distance = Math.sqrt(x * x + y * y + z * z);
      return distance;
    }
     public static ArrayList<Atom> loadExposedAtom(String path, String filename){
        ArrayList<Atom> listAtom = new ArrayList();
        try{
            if (new File(path+filename).exists()) {
                FileReader fr = new FileReader(new File(path+filename));
                BufferedReader br = new BufferedReader(fr);
                String dataline = br.readLine(); 
                //System.out.println(dataline);
                while (dataline != null){
                    String[] splitstr = dataline.split(" ");
                    double x = Double.valueOf(splitstr[2]);
                    double y = Double.valueOf(splitstr[3]);
                    double z = Double.valueOf(splitstr[4]);
                    Atom a = new Atom();
                    a.setX(x);
                    a.setY(y);
                    a.setZ(z);
                    a.setResidueName(splitstr[1]);
                    a.setAtomName(splitstr[0]);
                    listAtom.add(a);
                    dataline = br.readLine();
                }
                br.close();

            }
        }
        catch(Exception e){

        }
        return listAtom;
    }
     
     public static ArrayList<Atom> loadExposedAtomResiduLevel(String path, String filename){
        ArrayList<Atom> listAtom = new ArrayList();
        System.out.println("exposed: "+ path+filename);
        try{
            if (new File(path+filename).exists()) {
                FileReader fr = new FileReader(new File(path+filename));
                BufferedReader br = new BufferedReader(fr);
                String dataline = br.readLine(); 
                //System.out.println(dataline);
                while (dataline != null){
                    String[] splitstr = dataline.split(" ");
                    double x = Double.valueOf(splitstr[0]);
                    double y = Double.valueOf(splitstr[1]);
                    double z = Double.valueOf(splitstr[2]);
                    Atom a = new Atom();
                    a.setX(x);
                    a.setY(y);
                    a.setZ(z);
                    a.setResidueName(splitstr[4]);
                    a.setResidueID(splitstr[5]);
                    listAtom.add(a);
                    dataline = br.readLine();
                }
                br.close();

            }
        }
        catch(Exception e){

        }
        return listAtom;
    }
     
     public void concatResiduInfo(ArrayList<Atom> listAtom_res, ArrayList<Atom> listAtom_atom  ){
        // System.out.println("res size:"+ listAtom_res.size());
        // System.out.println("atom size:"+ listAtom_atom.size());
        
         if(listAtom_res.size()==listAtom_atom.size()){
             for (int i=0;i<listAtom_res.size();i++){
                 Atom ar = listAtom_res.get(i);
                 Atom aa = listAtom_atom.get(i);
                 //System.out.println(aa.getAtomName());
                 if(ar.getX()== aa.getX()&& ar.getY()== aa.getY()&& ar.getZ()==aa.getZ()){
                     listAtom_res.get(i).setAtomName(aa.getAtomName());
                        //System.out.println(listAtom.get(i).getAtomName());
                 }
                 
             }
         }
         
         
     }
}
