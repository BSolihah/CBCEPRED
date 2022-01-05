/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CBCTop;

import graph.*;

/**
 *
 * @author dell
 */
public class Edge {

        String v1;
        String v2;
        double w;
        public Edge(String s){
            String[] split_s = s.split("\t");
            this.v1 = split_s[0].trim();
            this.v2 = split_s[1].trim();
            this.w = Double.valueOf(split_s[2].trim());
        }
        public String getV1() {
            return v1;
        }
        public String getV2() {
            return v2;
        }
        public double getW() {
            return w;
        }

        @Override
        public String toString() {
            return "Edge{" + "v1=" + v1 + ", v2=" + v2 + ", w=" + w + '}';
        }
        

        
        
        
     
}
