/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package SpatialAnalysis;

/**
 *
 * @author dell
 */
public class Cluster{
        private int idx;
        private int number;
        
        public Cluster(int idx,int number){
            this.idx = idx;
            this.number=number;
        }
        public int getIdx() {
            return idx;
        }

        public void setIdx(int idx) {
            this.idx = idx;
        }

        public int getNumber() {
            return number;
        }

        public void setNumber(int number) {
            this.number = number;
        }

        
        
}
