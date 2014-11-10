/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package massspec;

import java.util.Comparator;

/**
 *
 * @author Pavel
 */
public class VertComparator implements Comparator{

    @Override
    public int compare(Object o1, Object o2) 
    {
        SpecVertex v1 = (SpecVertex) o1;
        SpecVertex v2 = (SpecVertex) o2;
        return v1.number - v2.number;
    }
    
}
