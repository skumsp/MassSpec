/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package massspec;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;

/**
 *
 * @author Olya
 */
public class Spectrogram {
    public ArrayList <HashMap> baseCleavage;
    public ArrayList <HashMap> kmersMasses;
    public HashMap <Character, Double> nuclMasses;
    public int smallCutoff;
    
    
    Spectrogram(ArrayList <HashMap> localbaseCleavage, HashMap <Character, Double> localnuclMasses)
    {
       baseCleavage=new ArrayList <HashMap> (localbaseCleavage);
       nuclMasses = new HashMap<Character, Double> (localnuclMasses);
    }
    Spectrogram(String nucleotides, char[] nucls, double[] masses, int smallCut)
    {
        smallCutoff = smallCut;
        nuclMasses = new HashMap <Character, Double> ();   //A,C,T,G
        for (int i=0; i<masses.length;i++)
            nuclMasses.put(nucls[i], masses[i]);
        
        baseCleavage=new ArrayList <HashMap> ();

        Set s = nuclMasses.keySet();// set of keys of hash
        for (Object o : s)
        {
            Character nucl = (Character) o;// nucl is one of our bases - A,C,T, or G (order does not matter)
            HashMap <Double, Double> foo= new HashMap <Double, Double> ();

            StringTokenizer st = new StringTokenizer (nucleotides, ""+ nucl);

            while (st.hasMoreTokens())
            {
                String kmer=st.nextToken();
                if (kmer.length()<smallCutoff) 
                    continue;
                double m=calcMassKmer(kmer);
                if(foo.containsKey(m)) 
                    foo.put(m, foo.get(m) + 1);
                else 
                    foo.put(m, 1.0);

            }
            baseCleavage.add(foo);
        }
    }
    Spectrogram(String nucleotides, char[] nucls, double[] masses, int smallCut, boolean storeKmers)
    // <editor-fold defaultstate="collapsed">
    {
        smallCutoff = smallCut;
        nuclMasses = new HashMap <Character, Double> ();   //A,C,T,G
        for (int i=0; i<masses.length;i++)
            nuclMasses.put(nucls[i], masses[i]);
        
        baseCleavage=new ArrayList <HashMap> ();
        kmersMasses=new ArrayList <HashMap> ();

        Set s = nuclMasses.keySet();// set of keys of hash
        for (Object o : s)
        {
            Character nucl = (Character) o;// nucl is one of our bases - A,C,T, or G (order does not matter)
            HashMap <Double, Double> foo= new HashMap <Double, Double> ();
            HashMap<String,Double> koo = new HashMap<String,Double>();

            StringTokenizer st = new StringTokenizer (nucleotides, ""+ nucl);

            while (st.hasMoreTokens())
            {
                String kmer=st.nextToken();
                if (kmer.length()<smallCutoff) 
                    continue;
                double m=calcMassKmer(kmer);
                if(foo.containsKey(m)) 
                    foo.put(m, foo.get(m) + 1);
                else 
                    foo.put(m, 1.0);
                koo.put(kmer, m);

            }
            baseCleavage.add(foo);
            kmersMasses.add(koo);
        }
    }
    // </editor-fold>
    Spectrogram(String addr)
    {
        
    }
    
    void PrintToFile(String nameFile) throws IOException
    {
       FileWriter fw = new FileWriter (nameFile); 
       Set alphabet = nuclMasses.keySet();
       Iterator ir = alphabet.iterator();
       
       for (HashMap hm : baseCleavage)
       {
           
           fw.write(">" + (Character) ir.next() + "\n");
           Set s = hm.entrySet();
           for (Object o : s)
           {
               Map.Entry me = (Map.Entry) o;
               fw.write(me.getKey() + " " + me.getValue() + "\n");                       
           }
       }
       fw.close();
    }
    void PrintToFileKmers(String nameFile) throws IOException
    // <editor-fold>
    {
       FileWriter fw = new FileWriter (nameFile); 
       Set alphabet = nuclMasses.keySet();
       Iterator ir = alphabet.iterator();
       Iterator ir1 = baseCleavage.iterator();
      
       for (HashMap hm : kmersMasses)
       {
           fw.write(">" + (Character) ir.next() + "\n");
           HashMap currcleav = (HashMap) ir1.next();
           Set s = hm.entrySet();
           for (Object o : s)
           {
               Map.Entry me = (Map.Entry) o;
               double m = (Double) me.getValue();
               double intens = (Double) currcleav.get(m);
               fw.write(me.getKey() + " " + me.getValue() + " " + intens +  "\n");                       
           }
       }
       fw.close();
    }
    // </editor-fold>


    double calcMassKmer (String kmer)
    {
        double sum=0;
        for (int i=0; i<kmer.length();i++)
            sum+=nuclMasses.get(kmer.charAt(i));
        return sum;    
    }
    void setSmallCutoff(int c)
    {
        smallCutoff = c;
    }
    int getNPeaks()
    {
        int nPeaks = 0;
        for (HashMap hm : baseCleavage)
            nPeaks += hm.size();
        return nPeaks;
    }
    HashMap<Double,Double> getAllPeaks()
    {
        HashMap<Double,Double> out = new HashMap<Double,Double>();
        for (HashMap hm : baseCleavage)
            for (Object o : hm.entrySet())
            {
                Map.Entry me = (Map.Entry) o;
                double mass = (Double) me.getKey();
                double intencity = (Double) me.getValue();
                if (out.containsKey(mass))
                        out.put(mass, out.get(mass) + intencity);
                else
                        out.put(mass, intencity);
            }
        return out;
    }
    Set getAlphabet()
    {
        return nuclMasses.keySet();
    }
}
