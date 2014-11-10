/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package massspec;

import java.io.*;
import java.util.Timer;

/**
 *
 * @author kki8
 */
public class LPSolver {
    int maxtime;
    void setMaxTime(int t)
    {
        maxtime = t;
    }
    public Process solve(String inputFile) throws IOException, InterruptedException
    {
        Runtime run=Runtime.getRuntime(); 
	Process p=null;
        String exPath = "glpk" + File.separator + "glpsol.exe --cpxlp --tmlim " + maxtime + " --output " + inputFile + "_sol.txt " + inputFile;
        p=run.exec(exPath); 
                            
        p.waitFor();
        BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
        System.out.println("Here is the standard output of the command:\n");
        String s;
        while ((s = stdInput.readLine()) != null) 
        {
            System.out.println(s);
        }
        stdInput.close();
        return p;
/*        OutputStream outputStream = p.getOutputStream();
        PrintStream printStream = new PrintStream(outputStream);
        printStream.println();
        printStream.flush();
        printStream.close();*/
    }
    
}
