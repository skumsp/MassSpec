/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package massspec;

/**
 *
 * @author kki8
 */
public class SpecEdge {
    SpecVertex begin;
    SpecVertex end;
    double cost;
    double capacity;
    int number;
    String label;
    SpecEdge(SpecVertex v1, SpecVertex v2, double x, double y)
    {
        begin = v1;
        end = v2;
        cost = x;
        capacity = y;
        label = "";
    }
    SpecEdge(SpecVertex v1, SpecVertex v2, double x, double y, String s)
    {
        begin = v1;
        end = v2;
        cost = x;
        capacity = y;
        label = s;
    }
}
