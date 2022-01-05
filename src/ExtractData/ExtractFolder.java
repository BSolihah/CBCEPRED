/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ExtractData;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author dell
 */
public class ExtractFolder {
    
    public ArrayList<String> listAllFiles(File folder){
         System.out.println("In listAllfiles(File) method");
         ArrayList<String> listfile = new ArrayList();
         File[] fileNames = folder.listFiles();
         for(File file : fileNames){
             // if directory call the same method again
             if(file.isDirectory()){
                 listAllFiles(file);
             }else{
                 listfile.add(file.getName());
                         
             }
         }
         return listfile;
     }
    public void readContent(File file) throws IOException{
         System.out.println("read file " + file.getCanonicalPath() );
         try(BufferedReader br  = new BufferedReader(new FileReader(file))){
               String strLine;
               // Read lines from the file, returns null when end of stream 
               // is reached
               while((strLine = br.readLine()) != null){
                System.out.println("Line is - " + strLine);
               }
         }
     }
     

    
}
