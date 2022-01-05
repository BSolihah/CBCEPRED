/*
 * menggunakan hasil dari ExtractDataBenchmark
 */
package ExtractDataTraining;

import ComplexCalculate.UnboundDSSP;
import ComplexStructure.Complex;
import ComplexStructure.Epitope;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author OpenGress
 */
public class ExtractASA {
    public ExtractASA() throws IOException{
    //ambil data benchmark
        ExtractDataBenchmark edb = new ExtractDataBenchmark();
        ArrayList<Epitope> epilist = edb.loadBenchmarkData("epitopelist.txt");
    //buat obyek UnboundDSSP untuk setiap pdbid
        System.out.println("amino acid vector length");
        for(Epitope e:epilist){
            String pdbid = e.getPdbid();
            String chain = e.getChain();
            String chain_num = e.getNo_chain();
            Complex complex = new Complex(pdbid, chain);
            //complex.getChainByChianID(chain.charAt(0)).setEpitopeStatus(e, chain);
            complex.getChainByChianID(chain.charAt(0)).setAAAsEpitope(e);
            //System.out.println (pdbid+"\t"+complex.getChainByChianID(chain.charAt(0)).getAminoacidVector().size());
           // complex.getChainByChianID(chainid.charAt(0)).getAminoacidVector().size();
          //  System.out.println("panjang sequence"+ e.getSequence().length());
           // complex.setEpitopeStatus(e, chainid.charAt(0));
            UnboundDSSP unbound_dssp = new UnboundDSSP(complex, chain);
            unbound_dssp.writeASAToFile(complex, chain.charAt(0), pdbid);
        }
    //simpan nilai ASAnya untuk semua residu
    }
    
    
    
}
