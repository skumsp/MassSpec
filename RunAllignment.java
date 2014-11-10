/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package massspec;
import ErrorCorrection.*;
import java.io.IOException;
import java.util.Timer;

/**
 *
 * @author kki8
 */
public class RunAllignment {
    public static void main(String[] args) throws IOException, InterruptedException
     {
         char[] nucls = {'A','T','G','C'};
         double[] masses = {329.21,306.17,345.21,305.18};
         int smallCutoff = 0;
         double epsilon = 0.0001;
         
         int ntests = 20;
         int maxtime = 12000;
         
         long[] times = new long[ntests];
         long totaltime = 0;
         
     
         for (int i = 0; i < ntests; i++)
         {
            TestsGenerator tg = new TestsGenerator(nucls);
            int[] pos = {7,25,48};
            String[] seq =  tg.Generate(60, pos);

            String seqFile = "genTest" + i + ".fas";
            tg.printSequences(seqFile);
            Spectrogram s1 = new Spectrogram(seq[0],nucls,masses,smallCutoff,true);
            Spectrogram s2 = new Spectrogram(seq[1],nucls,masses,smallCutoff,true); 

/*             String seqFile = "genTest8.fas";
            DataSet ds = new DataSet(seqFile,0);
            Spectrogram s1 = new Spectrogram(ds.reads.get(0).getNucl(),nucls,masses,smallCutoff,true);
            Spectrogram s2 = new Spectrogram(ds.reads.get(1).getNucl(),nucls,masses,smallCutoff,true);*/


            s1.PrintToFileKmers(seqFile + "_seq1.seq");
            s2.PrintToFileKmers(seqFile + "_seq2.seq");
            SpecGraph sg = new SpecGraph(s1,s2,epsilon);
            sg.printToFile(seqFile + "_graph.grph");

            String outFile = seqFile + "_flow.lp";
            sg.printMaxFlowProblemCPLEX(outFile);

            String pajekFile = seqFile + "_pajek.net";
            sg.printForPajek(pajekFile);

            String pajekFilePart = seqFile + "_pajek.clu";
            sg.printForPajekPartition(pajekFilePart);

            LPSolver lps = new LPSolver();
            lps.setMaxTime(maxtime);
            
            long startTime = System.currentTimeMillis();
            
            lps.solve(outFile);
            
            
/*            Process p = null;
            int waitTimeout = 3000000;
            Timer timer = new Timer(true);
            InterruptTimerTask interruptTimerTask = 
                new InterruptTimerTask(Thread.currentThread());
            timer.schedule(interruptTimerTask, waitTimeout);
            try {
                p=lps.solve(outFile);
            } catch (InterruptedException e) {
                System.out.println("timeout exeeded");
                p.destroy();
            } finally {
                timer.cancel(); 
            }*/
            
             
            times[i] = System.currentTimeMillis() - startTime;
            totaltime+= times[i];
            System.out.println("Solution time: " + times[i]);
         }


         System.out.println("Average time: " + totaltime/ntests);
         for (int i = 0; i < ntests; i++)
             System.out.print(times[i] + " ");
         System.out.println();
//         System.exit(0); 




     }

}
