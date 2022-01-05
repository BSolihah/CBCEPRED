/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package StructureAnalysis;

import ComplexStructure.*;
import bpred.AnalysisOnSurfaceResidue;
import java.util.ArrayList;

/**
 *
 * @author OpenGress
 */
public class StructureBasedAnalysis {
    //hitungContactNumber
    //sphere radius = 13 A
    //input CA residu pusat dan CA komplek
    
    //langkah: menemukan c alpha dalam radius 13A dari calpha
    //merotasikan atom-atom tersebut sehingga align dengan calpha-cbeta atau 
    //calpha-pseudo cbeta dengan sumbu z
    
    //mencari aa tetangga dalam radius d
   
   
   public void calculateParamASA(Aminoacid aa, Complex complex){
       int[]efse = calculateEightSubAreaFSE(aa, complex);                
        int[]sfse = calculateSixSubAreaFSE(aa, complex);
        int[] fse= calculateFSE(aa, complex);
        int[] qse = calculateQSE(aa, complex);
        int[] hse = calculateHSE(qse);
        int cnfromhse = calculateCNFromHSE(hse);
        int cn = ContactNumber(aa, complex);
        ParameterASA pAsa = new ParameterASA();
        pAsa.setCn(cn);
        pAsa.setCnfromHSE(cnfromhse);
        pAsa.setFse(fse);
        pAsa.setHse(hse);
        pAsa.setQse(qse);
        pAsa.setEfse(efse);
        pAsa.setSfse(sfse);
        aa.setpASA(pAsa);
   }
    public int ContactNumber(Aminoacid aa, Complex c){
        int cn=0;
        Atom ca = aa.getCAlpha();
        ArrayList<Atom> listca = c.getCAlphaOfComplex();
        double []dAaToC = new double[listca.size()];
        for(int i=0;i<listca.size();i++){
            dAaToC[i]=distance3DPoint(ca,listca.get(i));
            //print jarak
            if (dAaToC[i]<=13){
                    //System.out.print("i:" +i+ ": "+ dAaToC[i]+"\n");
                cn=cn+1;
            }
            
        }
        return cn;
        
    }
    //analisis jarak residu ke residu tetangganya
    public void distanceAnalysisOnExposedResidue(){
    //ambil residu terexpose
    AnalysisOnSurfaceResidue aosr = new AnalysisOnSurfaceResidue();
    aosr.doAnalysisOnSurfaceResidu();
    ArrayList<Aminoacid> exposedAA= aosr.getExposedAA();
        
    //ambil residu epitope
    //ambil residu nonepitope
    for(int i=0;i<exposedAA.size();i++){
        
    }
        
    }
    //buat fungsi jarak antar residu berdasarkan posisi sembarang atom
    //untuk semua atom dlm residu a hitung jaraknya terhadap residu b
    //jika ditemukan maka residu adalah dekat dengan residu pusat
    public double minDistAtomToAmino(Aminoacid a, Atom x){
        double d= Double.MAX_VALUE;
        ArrayList<Atom> atomA= a.getAtomvector();
        for(Atom atom: atomA){
            double da= this.distance3DPoint(atom, x);
            if(d>da) d=da;
        }
        return d;
    }
    public double minDistancebetweenAA(Aminoacid a, Aminoacid b){
        ArrayList<Atom> atomA= a.getAtomvector();
        double daAtoB = Double.MAX_VALUE;
        for(Atom atom: atomA){
            double datomtoB= this.minDistAtomToAmino(b, atom);
            if(daAtoB>datomtoB) daAtoB=datomtoB;
        }
        return daAtoB;
    }
    //deteksi aa dekat residu pusat
    public ArrayList<Aminoacid> NearAminacid(Aminoacid aa, Complex c,  double d){
        ArrayList<Aminoacid> nearAA= new ArrayList();
        double dist ;
        ArrayList<Aminoacid> aaChain = c.getAminoacidVector();
        for(Aminoacid aac: aaChain){
            dist = this.minDistancebetweenAA(aac, aa);
            if(dist<d)
                nearAA.add(aac);
        }
        return nearAA;
    }
    public ArrayList<String> residuIdStrucWindSmooth(Aminoacid aa, Complex c,  double d){
        ArrayList<Aminoacid> nearaa= NearAminacid(aa,c,d);
        ArrayList<String> idresidu = new ArrayList();
        for(int i=0;i<nearaa.size();i++){
            Aminoacid a= nearaa.get(i);
            idresidu.add(a.getresidueID());
            
        }
        return idresidu;
                
    }
    public ArrayList<Integer> ContactIdOfNeighbor(Aminoacid aa, Complex c){
       // ArrayList<String> listresidcn = new ArrayList();
        ArrayList<Integer> listresidcn = new ArrayList();
        Atom ca = aa.getCAlpha();
        ArrayList<Atom> listca = c.getCAlphaOfComplex();
        double []dAaToC = new double[listca.size()];
        for(int i=0;i<listca.size();i++){
            dAaToC[i]=distance3DPoint(ca,listca.get(i));
            //print jarak
            if (dAaToC[i]<=13){
                    //System.out.print("i:" +i+ ": "+ dAaToC[i]+"\n");
                int atom_no = listca.get(i).getAtomNo();
                //String resid = listca.get(i).getresidueID();
                listresidcn.add(atom_no);
            }
            
        }
        return listresidcn;
    }
    //deteksi residu dalam jarak dA
    public ArrayList<NeighborResidu> NeighboorAminoAcid(Aminoacid aa, Complex c, double d){
        ArrayList<NeighborResidu> neighboorAA= new ArrayList();
        Atom ca = aa.getCAlpha();
        ArrayList<Aminoacid> listaa = c.getAminoacidVector();
        for (int i=0;i<listaa.size();i++){
            double dist = distance3DPoint(ca,listaa.get(i).getCAlpha());
            if (dist<=d){
                NeighborResidu nR = new NeighborResidu();
                nR.setAa(listaa.get(i));
                nR.setDistance(dist);
                nR.setIdxOnComplex(i);
                neighboorAA.add(nR);
            }
        }
        return neighboorAA;
        
    }
    //deteksi residu yang jaraknya dalam radius 3A
     public ArrayList<Atom> listNearNeighboorAtom(Aminoacid aa, Complex c){
        ArrayList<Atom> nearneighboorAtom= new ArrayList();
        Atom ca = aa.getCAlpha();
        ArrayList<Atom> listca = c.getCAlphaOfComplex();
        for (int i=0;i<listca.size();i++){
            double dist = distance3DPoint(ca,listca.get(i));
            if (dist<=4){
                nearneighboorAtom.add(listca.get(i));
            }
        }
        return nearneighboorAtom;
        
    }
    
     //deteksi residu yang jaraknya dalam radius 13A
    public ArrayList<Atom> listNeighboorAtom(Aminoacid aa, Complex c){
        ArrayList<Atom> neighboorAtom= new ArrayList();
        Atom ca = aa.getCAlpha();
        ArrayList<Atom> listca = c.getCAlphaOfComplex();
        for (int i=0;i<listca.size();i++){
            double dist = distance3DPoint(ca,listca.get(i));
            if (dist<=13){
                neighboorAtom.add(listca.get(i));
            }
        }
        return neighboorAtom;
        
    }
    // cbeta_pseudo: merotasikan atom N sebesar -120 terhadap sumbu CA-C
    private AtomCoordinate calculateCbetaPseudo(Aminoacid aa){
    //definisikan vektor dari titik CA ke titik C
        Atom ca = aa.getAtomofAA("CA");
        Atom c = aa.getAtomofAA("C");
        double A = ca.getX();
        double B= ca.getY();
        double C = ca.getZ();
        double u = c.getX()-ca.getX();
        double v= c.getY()-ca.getY();
        double w = c.getZ()-ca.getZ();
        double L = u*u + v*v + w*w;
        double sqrtL = Math.sqrt(L);
        double sintheta = Math.sin(-2.0943951024);
        double costheta = Math.cos(-2.0943951024);
        double satumincostheta = 1- Math.cos(-2.0943951024);
        double m11 = (u*u + (v*v+w*w)*costheta )/L;
        double m12 = (u*v*satumincostheta -w*sqrtL*sintheta)/L;
        double m13 = (u*w*satumincostheta + v*sqrtL*sintheta)/L;
        double m14 = ((A*(v*v+w*w)-u*(B*v+C*w))*satumincostheta + (B*w-C*w)*sqrtL*sintheta)/L;
        double m21 = (u*v*satumincostheta+w*sqrtL*sintheta)/L;
        double m22 = (v*v+(u*u+w*w)*costheta)/L;
        double m23 = (u*w*satumincostheta - u*sqrtL*sintheta)/L;
        double m24 = ((B*(u*u+w*w)-v*(A*u+C*w))*satumincostheta + (C*u-A*w)*sqrtL*sintheta)/L;
        double m31 = (u*w*satumincostheta-v*sqrtL*sintheta)/L;
        double m32= (u*w*satumincostheta + u*sqrtL*sintheta)/L;
        double m33 = (w*w+(u*u+v*v)*costheta)/L;
        double m34 = ((C*(u*u+v*v)-w*(A*u+B*v))*satumincostheta+(A*v-B*u)*sqrtL*sintheta)/L;
                
    //get atom N    
        Atom n= aa.getAtomofAA("N");
        double nx = n.getX();
        double ny = n.getY();
        double nz = n.getZ();
        
        //get cbetapseudo
        AtomCoordinate cbpse = new AtomCoordinate();
        cbpse.setX(nx*m11+ny*m12+nz*m13+m14);
        cbpse.setY(nx*m21+ny*m22+nz*m23+m24);
        cbpse.setZ(nx*m31+ny*m32+nz*m33+m34);
        return cbpse;
        
    }
    private AtomCoordinate getAtomCoordinateofAtom(Atom a){
        AtomCoordinate ac = new AtomCoordinate();
        ac.setX(a.getX());
        ac.setY(a.getY());
        ac.setZ(a.getZ());
        return ac;
    }
    public void printCBetaCoordinateAsOrigin(Aminoacid aa){
        AtomCoordinate cb = this.getAtomCoordinateofAtom(aa.getAtomofAA("CB"));
        if (cb == null){
            AtomCoordinate cbp = this.calculateCbetaPseudo(aa);
            System.out.print(cbp.getX()+", "+cbp.getY()+", "+ cbp.getZ()+"\n");
        } else{
            System.out.print(cb.getX()+", "+cb.getY()+", "+ cb.getZ()+"\n");
        }
    }
    //fungsi mendefinisikan koordinate asal ke koordinate frame ref yang baru
    //input koordinate neighbor atom dan frame ref
    //output representasi baru koordinate neighbor
    public ArrayList<AtomCoordinate> defineCoordOnFrameRef(ArrayList<AtomCoordinate> frameref, ArrayList<Atom> neighboorAtom){
    double m11 = frameref.get(1).getX();
    double m12 = frameref.get(2).getX();
    double m13 = frameref.get(3).getX();
    double m14 = frameref.get(0).getX();
    
    double m21 = frameref.get(1).getY();
    double m22 = frameref.get(2).getY();
    double m23 = frameref.get(3).getY();
    double m24 = frameref.get(0).getY();
    
    double m31 = frameref.get(1).getZ();
    double m32 = frameref.get(2).getZ();
    double m33 = frameref.get(3).getZ();
    double m34 = frameref.get(0).getZ();
    ArrayList<AtomCoordinate> newcoord = new ArrayList();
    for (int i=0;i<neighboorAtom.size();i++){
        AtomCoordinate ac = this.getAtomCoordinateofAtom(neighboorAtom.get(i));
        //System.out.print(ac.getX()+", "+ac.getY()+", "+ ac.getZ()+"\n");
        double newx = m11*ac.getX()+m12*ac.getY()+m13*ac.getZ()+m14;
        double newy = m21*ac.getX()+m22*ac.getY()+m23*ac.getZ()+m24;
        double newz = m31*ac.getX()+m32*ac.getY()+m33*ac.getZ()+m34;
        AtomCoordinate coord = new AtomCoordinate(newx, newy, newz);
      //  System.out.print(coord.getX()+", "+coord.getY()+", "+ coord.getZ()+"\n");
        newcoord.add(coord);
    }
    return newcoord;
    }
    public void printFrameRef(Aminoacid aa){
        ArrayList<AtomCoordinate> frameref = calculateFrameRef(aa);
        for (int i=0; i<frameref.size();i++){
            System.out.print(frameref.get(i).getX()+", "+frameref.get(i).getY()+", "+ frameref.get(i).getZ());
            System.out.print("\n");
        }
    }
    public ArrayList<AtomCoordinate> calculateFrameRef(Aminoacid aa){
        //definisikan e1, e2 dan e3
        //e1: CbN x CbCa / norm (CbNxCbCa)
        AtomCoordinate cb= this.getAtomCoordinateofAtom(aa.getAtomofAA("CB"));
        
        if (cb==null){
            cb= this.calculateCbetaPseudo(aa);
        }
        
        AtomCoordinate n = this.getAtomCoordinateofAtom(aa.getAtomofAA("N"));
        AtomCoordinate ca = this.getAtomCoordinateofAtom(aa.getAtomofAA("CA"));
       // System.out.println("Koordinat Cbeta");
       // System.out.println(cb.getX()+"\t"+cb.getY()+"\t"+cb.getZ());
       // System.out.println("Koordinat N");
       // System.out.println(n.getX()+"\t"+n.getY()+"\t"+n.getZ());
       // System.out.println("Koordinat CA");
       // System.out.println(ca.getX()+"\t"+ca.getY()+"\t"+ca.getZ());
        //vektor e1 = CbN x CbCa
        double CbN_x = n.getX()-cb.getX();
        double CbN_y = n.getY()-cb.getY();
        double CbN_z = n.getZ()-cb.getZ();
        AtomCoordinate CbN = new AtomCoordinate(CbN_x,CbN_y,CbN_z);
       // System.out.println("Koordinat CbN");
       // System.out.println(CbN.getX()+"\t"+CbN.getY()+"\t"+CbN.getZ());
        double CbCa_x = ca.getX()-cb.getX();
        double CbCa_y = ca.getY()-cb.getY();
        double CbCa_z = ca.getZ()-cb.getZ();
        AtomCoordinate CbCa = new AtomCoordinate(CbCa_x,CbCa_y,CbCa_z);
      //  System.out.println("Koordinat CbCa");
       // System.out.println(CbCa.getX()+"\t"+CbCa.getY()+"\t"+CbCa.getZ());
        AtomCoordinate CbNCbCa = crossProduct(CbN,CbCa);
        //System.out.println("crossProduct(CbN,CbCa)");
        //System.out.println(CbNCbCa.getX()+"\t"+CbNCbCa.getY()+"\t"+CbNCbCa.getZ());
        double x2 = Math.pow(CbNCbCa.getX(),2);
        
        double y2=Math.pow(CbNCbCa.getY(),2);
        double z2= Math.pow(CbNCbCa.getZ(),2);
        double normCbNCbCa = Math.sqrt(x2 + y2+z2);
        //System.out.println("normCbNCbCa: "+ normCbNCbCa);
        
        AtomCoordinate e1 = new AtomCoordinate(CbNCbCa);
        e1.setX(e1.getX()/normCbNCbCa);
        e1.setY(e1.getY()/normCbNCbCa);
        e1.setZ(e1.getZ()/normCbNCbCa);
                
        //e2:CbN/norm(CbN)
        double normCbN = Math.sqrt(CbN.getX()*CbN.getX()+CbN.getY()*CbN.getY()+CbN.getZ()*CbN.getZ());
        //System.out.println("normCbN: "+ normCbN);
        AtomCoordinate e2 = new AtomCoordinate(CbN);
        e2.setX(e2.getX()/normCbN);
        e2.setY(CbN.getY()/normCbN);
        e2.setZ(CbN.getZ()/normCbN);
        
        //e3:e1xe2
        AtomCoordinate e3 = crossProduct(e1,e2);
        
        //
        ArrayList<AtomCoordinate> frameRef = new ArrayList();
        frameRef.add(cb);
        frameRef.add(e1);
        frameRef.add(e2);
        frameRef.add(e3);
        return frameRef;
    }
    private AtomCoordinate crossProduct(AtomCoordinate a, AtomCoordinate b){
        //a1= xa; a2=ya; a3=za
        //b1=xb;b2=yb;b3=zb
        //x=a2b3-b2a3 = yazb-ybza
        //y=b1a3-a1b3=xbza-xazb
        //z=a1b2-b1a2=xayb - xbya
        AtomCoordinate atom= new AtomCoordinate();
        atom.setX( a.getY()* b.getZ()-a.getZ()*b.getY());
       atom.setY(a.getZ()*b.getX()- a.getX()*b.getZ());
        atom.setZ(a.getX()*b.getY()-a.getY()*b.getX());
        return atom;
    }
    private double calculateNorm(AtomCoordinate a){
        //hitung norm dengan pytaghoras
        return Math.sqrt(a.getX()*a.getX()+a.getY()*a.getY()+a.getZ()*a.getZ());
    }
    
    
    // qse dari asam amino 
    public int[] calculateQSE(Aminoacid aa, Complex c){
        int [] listquadran = new int [8];
        for(int j=0;j<8;j++){
            listquadran[j]=0;
        }
        ArrayList<Atom> neighbooratom =this.listNeighboorAtom(aa, c);
        //print neighbooratom
        /*
        System.out.println("c alpha tetangga");
        for(Atom a: neighbooratom){
            System.out.println(a.getresidueID()+"_"+a.getresidueName()+"_"+a.getatomName()+":"+a.getX()+"\t"+a.getY()+"\t"+a.getZ());
        }
        * 
        */
        ArrayList<AtomCoordinate> frameref = this.calculateFrameRef(aa);
        //print frameref
        /*
        System.out.println("matrik koordinat qse");
        for(AtomCoordinate ac: frameref){
            System.out.println(ac.getX()+"\t"+ac.getY()+"\t"+ac.getZ());
        }*/
        ArrayList<AtomCoordinate> newneighbooratomcoord= this.defineCoordOnFrameRef(frameref, neighbooratom);
        //print koordinate hasil transformasi
        /*
        System.out.println("Koordinat atom Calpha baru tetangga");
        for(AtomCoordinate ac: newneighbooratomcoord){
            System.out.println(ac.getX()+"\t"+ac.getY()+"\t"+ac.getZ());
        }
        System.out.println("quadrant");
        * 
        */
        for(int i= 0; i<newneighbooratomcoord.size();i++){
            int quadran = detectQuadranOfCoordinateInFrame(frameref.get(0),newneighbooratomcoord.get(i));
          //  System.out.print(quadran +"\t");
            if(quadran == 1){
                listquadran[0]= listquadran[0]+1;
            }else if(quadran == 2){
                listquadran[1]= listquadran[1]+1;
            }else if(quadran == 3){
                listquadran[2]= listquadran[2]+1;
            }else if(quadran == 4){
                listquadran[3]= listquadran[3]+1;
            }else if(quadran == 5){
                listquadran[4]= listquadran[4]+1;
            }else if(quadran == 6){
                listquadran[5]= listquadran[5]+1;
            }else if(quadran == 7){
                listquadran[6]= listquadran[6]+1;
            }else if(quadran == 8){
                listquadran[7]= listquadran[7]+1;
            }
        }
        /*
        //System.out.println();
        System.out.println("list quadran");
        for(int k=0;k<8;k++){
            System.out.print(listquadran[k]+"\t");
        }
        System.out.println();
        * 
        */
        return listquadran;
    }
    //sixSubAreaFSE
    private int detectQuadrantInSixSubAreaFSE(AtomCoordinate origin, AtomCoordinate ac){
    //30 derajat = 0.523599 radian
        int quadran =0;
        //maxdistance = panjang jari-jari kali sudut pembagi area
        double maxdistance = 13*Math.sin(0.523599);
        double uplimit1 = origin.getZ()+ maxdistance;
        double uplimit2 = origin.getZ()+ 2*maxdistance;
        
        double downlimit1 = origin.getZ()-maxdistance;
        double downlimit2 = origin.getZ()-2*maxdistance;
        
         double n = ac.getZ();
         if(n>uplimit2){
            quadran = 1;
        }else if(n>uplimit1&&n<=uplimit2){
            quadran =2;
        }else if(n>origin.getZ()&&n<=uplimit1){
            quadran =3;
        }else if(n< origin.getZ()&& n >= downlimit1){
                quadran =4;
        }else if(n<downlimit1 && n>= downlimit2){
            quadran =5;
        }else if(n<downlimit2){
            quadran =6;
        }
        return quadran;
    }
    
    //eightSubAreaFSE
    private int detectQuadrantInEightSubAreaFSE(AtomCoordinate origin, AtomCoordinate ac){
    //22.5 derajat = 0.3926991 radian (22.5 dejarat)
        int quadran =0;
        //maxdistance = panjang jari-jari kali sudut pembagi area
        double maxdistance = 13*Math.sin(0.3926991);
        double uplimit1 = origin.getZ()+ maxdistance;
        double uplimit2 = origin.getZ()+ 2*maxdistance;
        double uplimit3 = origin.getZ()+ 3*maxdistance;
        double downlimit1 = origin.getZ()-maxdistance;
        double downlimit2 = origin.getZ()-2*maxdistance;
        double downlimit3 = origin.getZ()-3*maxdistance;
         double n = ac.getZ();
         if(n>uplimit3){
            quadran = 1;
        }else if(n>uplimit2&&n<=uplimit3){
            quadran =2;
        }else if(n>uplimit1&&n<=uplimit2){
            quadran = 3;
        }else if(n>origin.getZ()&&n<=uplimit1){
            quadran =4;
        }else if(n< origin.getZ()&& n >= downlimit1){
                quadran =5;
        }else if(n<downlimit1 && n>= downlimit2){
            quadran =6;
        }else if(n<downlimit2 && n >= downlimit3){
            quadran =7;
        }else if(n<downlimit3){
            quadran =8;
        }
        return quadran;
    }
    //fourth sphere exposure
    private int detectQuadranInFSE(AtomCoordinate origin, AtomCoordinate ac){
        int quadran =0;
        double uplimit = origin.getZ()+ 13*Math.sin(0.785398);
        double downlimit = origin.getZ()-13*Math.sin(0.785398);
        double n = ac.getZ();
        if(n>uplimit){
            quadran = 1;
        }else if(n>origin.getZ()&&n<=uplimit){
            quadran =2;
        }else if(n<downlimit){
            quadran = 4;
        }else{
            quadran =3;
        }
        return quadran;
    
    }
    public int calculateCNFromHSE(int[] hse){
        return hse[0]+hse[1];
    }
    public int[]calculateEightSubAreaFSE(Aminoacid aa, Complex c){
        int[] listquadran= new int[8];
        for(int j=0;j<listquadran.length;j++){
            listquadran[j]=0;
        }
        ArrayList<Atom> neighbooratom =this.listNeighboorAtom(aa, c);
        ArrayList<AtomCoordinate> frameref = this.calculateFrameRef(aa);
        if(neighbooratom.size()>0){
            ArrayList<AtomCoordinate> newneighbooratomcoord= this.defineCoordOnFrameRef(frameref, neighbooratom);
            for(int i= 0; i<newneighbooratomcoord.size();i++){
            //int quadran = detectQuadranOfCoordinateInFrame(frameref.get(0),newneighbooratomcoord.get(i));
            //int quadran = this.detectQuadranInFSE(frameref.get(0),newneighbooratomcoord.get(i));
            int quadran = this.detectQuadrantInEightSubAreaFSE(frameref.get(0),newneighbooratomcoord.get(i));
            if(quadran == 1){
                listquadran[0]= listquadran[0]+1;
            }else if(quadran == 2){
                listquadran[1]= listquadran[1]+1;
            }else if(quadran == 3){
                listquadran[2]= listquadran[2]+1;
            }else if(quadran == 4){
                listquadran[3]= listquadran[3]+1;
            }else if(quadran == 5){
                 listquadran[4]= listquadran[4]+1;
            }else if(quadran == 6){
                 listquadran[5]= listquadran[5]+1;
            }else if(quadran == 7){
                 listquadran[6]= listquadran[6]+1;
            }else if(quadran == 8){
                 listquadran[7]= listquadran[7]+1;
            }
        }
        }
        return listquadran;
        
    }
    
    public int[]calculateSixSubAreaFSE(Aminoacid aa, Complex c){
        int[] listquadran= new int[6];
        for(int j=0;j<listquadran.length;j++){
            listquadran[j]=0;
        }
        ArrayList<Atom> neighbooratom =this.listNeighboorAtom(aa, c);
        ArrayList<AtomCoordinate> frameref = this.calculateFrameRef(aa);
        if(neighbooratom.size()>0){
            ArrayList<AtomCoordinate> newneighbooratomcoord= this.defineCoordOnFrameRef(frameref, neighbooratom);
            for(int i= 0; i<newneighbooratomcoord.size();i++){
            //int quadran = detectQuadranOfCoordinateInFrame(frameref.get(0),newneighbooratomcoord.get(i));
            //int quadran = this.detectQuadranInFSE(frameref.get(0),newneighbooratomcoord.get(i));
            int quadran = this.detectQuadrantInSixSubAreaFSE(frameref.get(0),newneighbooratomcoord.get(i));
            if(quadran == 1){
                listquadran[0]= listquadran[0]+1;
            }else if(quadran == 2){
                listquadran[1]= listquadran[1]+1;
            }else if(quadran == 3){
                listquadran[2]= listquadran[2]+1;
            }else if(quadran == 4){
                listquadran[3]= listquadran[3]+1;
            }else if(quadran == 5){
                 listquadran[4]= listquadran[4]+1;
            }else if(quadran == 6){
                 listquadran[5]= listquadran[5]+1;
            }
        }
        }
        return listquadran;
        
    }
    
    public int[]calculateFSE(Aminoacid aa, Complex c){
        int [] listquadran = new int [4];
        for(int j=0;j<4;j++){
            listquadran[j]=0;
        }
        ArrayList<Atom> neighbooratom =this.listNeighboorAtom(aa, c);
        ArrayList<AtomCoordinate> frameref = this.calculateFrameRef(aa);
        if(neighbooratom.size()>0){
            ArrayList<AtomCoordinate> newneighbooratomcoord= this.defineCoordOnFrameRef(frameref, neighbooratom);
        for(int i= 0; i<newneighbooratomcoord.size();i++){
            //int quadran = detectQuadranOfCoordinateInFrame(frameref.get(0),newneighbooratomcoord.get(i));
            int quadran = this.detectQuadranInFSE(frameref.get(0),newneighbooratomcoord.get(i));
            if(quadran == 1){
                listquadran[0]= listquadran[0]+1;
            }else if(quadran == 2){
                listquadran[1]= listquadran[1]+1;
            }else if(quadran == 3){
                listquadran[2]= listquadran[2]+1;
            }else if(quadran == 4){
                listquadran[3]= listquadran[3]+1;
            }
        }
        }
        
        /*
        for(int k=0;k<8;k++){
            System.out.print(listquadran[k]+"\t");
        }
        */
        return listquadran;
    }
    //deteksi kuadran menggunakan QSE
    //input frameref dan titik koordinate
    //output quadran
    private int detectQuadranOfCoordinateInFrame(AtomCoordinate origin,AtomCoordinate ac ){
   /*
    Jika u > Cbx & v >Cbv  & n > Cbn => kuadran 1
    Jika u < Cbx & v >Cbv  & n > Cbn => kuadran 2
    Jika u < Cbx & v <Cbv  & n > Cbn => kuadran 3
    Jika u > Cbx & v <Cbv  & n > Cbn => kuadran 4
    Jika u > Cbx & v >Cbv  & n < Cbn => kuadran 5
    Jika u < Cbx & v >Cbv  & n < Cbn => kuadran 6
    Jika u < Cbx & v <Cbv  & n < Cbn => kuadran 7
    Jika u > Cbx & v <Cbv  & n < Cbn => kuadran 8
    */
        int quadran =0;
        double u = ac.getX();
        double v = ac.getY();
        double n = ac.getZ();
        double Cbx = origin.getX();
        double Cbv = origin.getY();
        double Cbn = origin.getZ();
        if(u > Cbx & v >Cbv  & n > Cbn){
            quadran = 1;
        }else if(u < Cbx & v >Cbv  & n > Cbn){
            quadran =2;
        }else if(u < Cbx & v <Cbv  & n > Cbn){
            quadran =3;
        }else if(u > Cbx & v <Cbv  & n > Cbn){
            quadran = 4;
        }else if(u > Cbx & v >Cbv  & n < Cbn){
            quadran = 5;
        }else if(u < Cbx & v >Cbv  & n < Cbn){
            quadran =6;
        }else if(u < Cbx & v <Cbv  & n < Cbn){
            quadran = 7;
        }else if(u > Cbx & v <Cbv  & n < Cbn){
            quadran = 8;
        }
        return quadran;
        
    }
    public int[] calculateHSE(int[] qse){
        int [] hse = new int[2];
        if(qse.length==8){
            hse[0]= qse[0]+qse[1]+qse[2]+qse[3];
            hse[1]=qse[4]+qse[5]+qse[6]+qse[7];
        }
        return hse;
        
    }
    private double distance3DPoint(Atom pa, Atom pb)
  {
    double distance = 0.0D;
    double x = pa.getX() - pb.getX();
    double y = pa.getY() - pb.getY();
    double z = pa.getZ() - pb.getZ();
    distance = Math.sqrt(x * x + y * y + z * z);
    return distance;
  }
}
