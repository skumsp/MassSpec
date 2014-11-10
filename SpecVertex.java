/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package massspec;

/**
 *
 * @author kki8
 */
public class SpecVertex {
    double mass;
    double intensity;
    String label;
    int number;
    SpecVertex()
    {
        
    }
    SpecVertex(double m, double i, String c)
    {
        mass = m;
        intensity = i;
        label = c;
    }
    
}
