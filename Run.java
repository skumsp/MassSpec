/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package massspec;

import ErrorCorrection.DataSet;
import java.io.IOException;
import java.util.ArrayList;
import org.apache.commons.math3.distribution.NormalDistribution;

/**
 *
 * @author kki8
 */
public class Run {
     public static void main(String[] args) throws IOException 
     {
         char[] nucls = {'A','T','G','C'};
         double[] masses = {329.21,306.17,345.21,305.18};
         int smallCutoff = 0;
         
/*         String sequence = "GATATGATGATGAACTGGTCGCCCACGGCCACC";
         
         Spectrogram spec = new Spectrogram(sequence,nucls,masses,true);
         spec.setSmallCutoff(smallCutoff);
         spec.PrintToFile("test_spec.txt");
         spec.PrintToFileKmers("test_spec_kmers.txt");
         
         double[] frequencies = {10.0};
         Spectrogram[] specs = {spec};
         SpecOperator sop = new SpecOperator();
         
         Spectrogram spec1 = sop.createMixture(specs, frequencies);
         spec1.PrintToFile("test_spec1.txt");*/
         
         String seqFile = "Clone Sequences2.fas";
         DataSet ds = new DataSet(seqFile);
         int nSeqTest = 4;
         Spectrogram[] specs = new Spectrogram[nSeqTest];
         for (int i = 0; i < nSeqTest; i++)
             specs[i] = new Spectrogram(ds.reads.get(i).getNucl(),nucls,masses,smallCutoff);
         
//         double[] frequencies = {10.0,20.0,30.0,40.0};
         double[] frequencies = {2.0,2.0,2.0,94.0};
         
         SpecOperator sop = new SpecOperator();
         Spectrogram mixture = sop.createMixture(specs, frequencies);
         
         double epsilon = 0.005;
         double deviation = 1.0;
         double[] resultEM = sop.expectationMaximization(mixture, specs, deviation, epsilon);
         
         System.out.println();
         
         for (int i = 0; i < nSeqTest; i++)
             System.out.println(100*resultEM[i]);
             
         
     }
    
}
