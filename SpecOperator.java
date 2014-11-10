/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package massspec;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.io.*;
import org.apache.commons.math3.distribution.NormalDistribution;


public class SpecOperator {
    Spectrogram createMixture(Spectrogram[] specs, double[] frequencies)
    {
        ArrayList<HashMap> baseCleavage = new ArrayList<HashMap>();
        HashMap nuclMasses = specs[0].nuclMasses;
        int nNucl = specs[0].baseCleavage.size();
        for (int i = 0; i < nNucl; i++)
        {
            HashMap <Double, Double> hm= new HashMap <Double, Double> ();
            for (int j = 0; j < specs.length; j++)
                for (Object o : specs[j].baseCleavage.get(i).entrySet())
                {
                    Map.Entry me = (Map.Entry) o;
                    double mass = (Double) me.getKey();
                    double intencity = (Double) me.getValue();
                    if (hm.containsKey(mass))
                        hm.put(mass, hm.get(mass) + frequencies[j]*intencity);
                    else
                        hm.put(mass, frequencies[j]*intencity);
                }
            baseCleavage.add(hm);
        }
        return new Spectrogram(baseCleavage, nuclMasses);
        
    }
    
    double[] expectationMaximization(Spectrogram spec, Spectrogram[] candidates, double deviation, double epsilon)
    {
        int nPeaksSpec = spec.getNPeaks();
        int nCand = candidates.length;
        double[][] adjMatrixBipart = new double[nCand][nPeaksSpec];
        HashMap peaks = spec.getAllPeaks();
        for (int i = 0; i < nCand; i++)
        {
            HashMap candPeaks = candidates[i].getAllPeaks();
            int j = 0;
            for (Object o : peaks.entrySet())
            {
//                adjMatrixBipart[i][j] = alignPeakPeaks((Map.Entry) o, candPeaks, deviation);
                adjMatrixBipart[i][j] = alignPeakPeaksCutoff((Map.Entry) o, candPeaks, deviation);
                j++;
            }
        }
        
        double[] deg_peaks = new double[nPeaksSpec];
        for (int j = 0; j < nPeaksSpec; j++)
        {
            deg_peaks[j] = 0;
            for (int i = 0; i < nCand; i++)
                deg_peaks[j] +=adjMatrixBipart[i][j]; 
        }
        
        double sumIntens = 0;
        for (Object o : peaks.entrySet())
            sumIntens += (Double) ((Map.Entry) o).getValue();
           
        double[] frequencies = new double[nCand];
        for (int i = 0; i < nCand; i++)
            frequencies[i] = 1.0/nCand;
        double diff = 1.0;
        while (diff > epsilon)
        {
            System.out.println(diff);
            // E-step
            double[][] prob = new double[nCand][nPeaksSpec];
            for (int j = 0; j < nPeaksSpec; j++)
            {
                double sum = 0;
                for (int i = 0; i < nCand; i++)
                    sum += frequencies[i]*adjMatrixBipart[i][j];
/*                if (sum == 0)
                    System.out.println("Sum = 0");*/
                 
                for (int i = 0; i < nCand; i++)
                    if (sum == 0)
                        prob[i][j] = 0;
                    else                     
                        prob[i][j] = (frequencies[i]*adjMatrixBipart[i][j])/sum;
            }
            
            // M-step
            double[] frequencies_old = new double[nCand];
            for (int i = 0; i < nCand; i++)
            {
                frequencies_old[i] = frequencies[i];
                double sum = 0;
                int j = 0;
                for (Object o : peaks.entrySet())
                {
                    sum += prob[i][j]*adjMatrixBipart[i][j]*((Double) ((Map.Entry) o).getValue()); 
                    j++;
                }
                frequencies[i] = sum;
            }
            
            diff = 0;
            for (int i = 0; i < nCand; i++)
                diff+=(frequencies[i] - frequencies_old[i])*(frequencies[i] - frequencies_old[i]);
            diff = Math.sqrt(diff);
            
        }
         
        double sum = 0;
        for (int i = 0; i < nCand; i++)
            sum+=frequencies[i];
         for (int i = 0; i < nCand; i++)
            frequencies[i]/=sum;
        return frequencies;
    }
    
    double alignPeakPeaks(Map.Entry peak, HashMap peaks, double deviation)
    {
        NormalDistribution d = new NormalDistribution(); 
        double mass = (Double) peak.getKey();
        double maxprob = 0;
        for (Object p : peaks.keySet())
        {
            double massCand = (Double) p;
            double z = Math.abs(mass - massCand)/deviation;
            double prob = 1 - d.cumulativeProbability(-z,z);
            if (prob > maxprob)
                maxprob = prob;
        }
        return maxprob;
    }
    
    double alignPeakPeaksCutoff(Map.Entry peak, HashMap peaks, double deviation)
    {
        double mass = (Double) peak.getKey();
        double maxprob = 0;
        for (Object p : peaks.keySet())
        {
            double massCand = (Double) p;
            double z = Math.abs(mass - massCand);
            if (z <= deviation)
                maxprob = 1;
        }
        return maxprob;
    }

}