/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package SpatialAnalysis;

import ComplexStructure.PredictedEpitope;

/**
 *
 * @author dell
 */
public class PairOfPredictedEpitope {
        private PredictedEpitope pe1;
        private PredictedEpitope pe2;
        private double weight;
        public PairOfPredictedEpitope(PredictedEpitope p1, PredictedEpitope p2, double w){
            this.pe1 = p1;
            this.pe2 =p2;
            this.weight=w;
        }

        public PredictedEpitope getPe1() {
            return pe1;
        }

        public PredictedEpitope getPe2() {
            return pe2;
        }

        public double getWeight() {
            return weight;
        }
        
    }