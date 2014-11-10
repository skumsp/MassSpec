/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package massspec;

import ErrorCorrection.DataSet;
import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author Pavel
 */
public class RunGenotyping {
    public static void main(String[] args) throws IOException 
     {
         char[] nucls = {'A','T','G','C'};
         double[] masses = {329.21,306.17,345.21,305.18};
         int smallCutoff = 0;
         
         int nGenotypes = 2;
         
         ArrayList frequencies = new ArrayList();
/*         double[] freq1= {1.0,2.0,3.0,4.0};
         double[] freq2= {5.0,5.0,5.0,5.0};
         double[] freq3= {20.0,20.0,30.0};*/
         
         double[] freq1= {1.0};
         double[] freq2= {99.0};
         
         frequencies.add(freq1);
         frequencies.add(freq2);
//         frequencies.add(freq3);
          
         
         Spectrogram[] mixturesGenotypes = new Spectrogram[nGenotypes];
         SpecOperator sop = new SpecOperator();
         
         Spectrogram[] cand_real = new Spectrogram[nGenotypes];
         
         for (int i = 0; i < nGenotypes; i++)
         {
             DataSet ds = new DataSet("samples_" + (i+1) + "a.fas");
             double[] fr = (double[]) frequencies.get(i);
             int nSeqTest = fr.length;
             Spectrogram[] specs = new Spectrogram[nSeqTest];
             for (int j = 0; j < nSeqTest; j++)
                specs[j] = new Spectrogram(ds.reads.get(j).getNucl(),nucls,masses,smallCutoff);
             cand_real[i] = specs[0];
             mixturesGenotypes[i] = sop.createMixture(specs, fr);
             
         }
         
         double[] auxFreq = new double[nGenotypes]; 
         for (int i = 0; i < nGenotypes; i++)
             auxFreq[i] = 1.0;
         
         Spectrogram mixture = sop.createMixture(mixturesGenotypes, auxFreq);
         mixture.PrintToFile("mixture_spec.txt");
         
         Spectrogram[] refSpec = new Spectrogram[nGenotypes];
         for (int i = 0; i < nGenotypes; i++)
         {
             DataSet ds = new DataSet("reference_" + (i+1) + "a.fas");
             refSpec[i] = new Spectrogram(ds.reads.get(0).getNucl(),nucls,masses,smallCutoff); 
         }
         
         double epsilon = 0.005;
         double deviation = 1.0;
         double[] resultEM = sop.expectationMaximization(mixture, refSpec, deviation, epsilon);
//         double[] resultEM = sop.expectationMaximization(mixture, cand_real, deviation, epsilon);
         
         System.out.println();
         
         for (int i = 0; i < nGenotypes; i++)
             System.out.println(100*resultEM[i]);
             
         
     }
    
}
