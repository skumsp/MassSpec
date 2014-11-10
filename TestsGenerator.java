/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package massspec;

import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * @author kki8
 */
public class TestsGenerator {
    char [] nucls;
    String seq1;
    String seq2;
    TestsGenerator(char[] n)
    {
        nucls = n;
    }
    String[] Generate(int len, int[] pos) throws IOException
    {
        String[] seq = new String[2];
        seq[0] = "";
        seq[1] = "";
        int j = 0;
        for (int i = 1; i <= len; i++)
        {
            if ((j < pos.length) && (i == pos[j]))
            {
                int k1 = (int) (4*Math.random());
                
                while (nucls[k1] == seq[0].charAt(i-2))
                    k1 = (int) (4*Math.random());
                
                seq[0]+=nucls[k1];
                int k2 = (int) (4*Math.random());
                while ((k1 == k2) || (nucls[k2] == seq[1].charAt(i-2)))
                    k2 = (int) (4*Math.random());
                seq[1]+=nucls[k2];
                j++;
            }
            else
            {
                int k1 = (int) (4*Math.random());
                seq[0]+=nucls[k1];
                seq[1]+=nucls[k1];
            }
        }
        
        seq1 = seq[0];
        seq2 = seq[1];       
        
        return seq;
    }
    void printSequences(String addr) throws IOException
    {
        FileWriter fw = new FileWriter(addr);
        fw.write(">seq0\n");
        fw.write(seq1 + "\n");
        fw.write(">seq1\n");
        fw.write(seq2 + "\n");
        fw.close();
    }
    
}
