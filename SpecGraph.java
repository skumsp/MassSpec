/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package massspec;
import java.util.*;
import java.io.*;

/**
 *
 * @author kki8
 */
public class SpecGraph {
    ArrayList<SpecVertex> vertices;
    HashMap<SpecVertex,ArrayList> adjLists;
    ArrayList<SpecVertex> vertices1;
    ArrayList<SpecVertex> vertices2;
    ArrayList<SpecVertex> pairvertices1;
    ArrayList<SpecVertex> pairvertices2;
    ArrayList<SpecVertex> auxvertices1; // aux vertices for left fork
    ArrayList<SpecVertex> auxvertices2; // aux vertices for right fork
    ArrayList<SpecVertex> deltvertices1; // aux vertices delta pair left;
    ArrayList<SpecVertex> deltvertices2; // aux vertices delta pair right;
    
    SpecGraph(Spectrogram s1, Spectrogram s2, double epsilon)
    {
        vertices1 = new ArrayList<SpecVertex>();
        vertices2 = new ArrayList<SpecVertex>();
        Set alphabet1 = s1.getAlphabet();
        Set alphabet2 = s2.getAlphabet();
        adjLists = new HashMap<SpecVertex,ArrayList>();
        
        Iterator ir = alphabet1.iterator();
        System.out.println("Creating left part");
        for (HashMap hm : s1.baseCleavage)
        {
            Set s = hm.entrySet();
            String label = "" + (Character) ir.next();
            for (Object o : s)
            {
                Map.Entry me = (Map.Entry) o;
//                label = label + "_" + (Double) me.getKey();
                double m = (Double) me.getKey();
                int intens = ((Double) me.getValue()).intValue();
                for (int i = 0; i < intens; i++)
                    vertices1.add(new SpecVertex((Double) me.getKey(),1,label));
//                vertices1.add(new SpecVertex((Double) me.getKey(),(Double) me.getValue(),label));
            }
        }
        
        System.out.println("Left part size = " + vertices1.size());
        
        System.out.println("Creating right part");
       ir = alphabet2.iterator();
        for (HashMap hm : s2.baseCleavage)
        {
            Set s = hm.entrySet();
            String label = "" + (Character) ir.next();
            for (Object o : s)
            {
                Map.Entry me = (Map.Entry) o;
//                label = label + "_" + (Double) me.getKey();
                double m = (Double) me.getKey();
                int intens = ((Double) me.getValue()).intValue();
                for (int i = 0; i < intens; i++)
                    vertices2.add(new SpecVertex((Double) me.getKey(),1,label));
//                vertices2.add(new SpecVertex((Double) me.getKey(),(Double) me.getValue(),label));
            }
        }
        System.out.println("Right part size = " + vertices2.size());
        
        
        System.out.println("Creating exact match edges");
        for (SpecVertex v1 : vertices1)
            for (SpecVertex v2 : vertices2)
                if ((Math.abs(v1.mass - v2.mass) < epsilon) && (v1.label.equalsIgnoreCase(v2.label)))
                {
                    SpecEdge e = new SpecEdge(v1,v2,0,1);
                    if (adjLists.containsKey(v1))
                    {
                        ArrayList a = adjLists.get(v1);
                        a.add(e);
                        adjLists.put(v1, a);
                    }
                    else
                    {
                        ArrayList<SpecEdge> a = new ArrayList<SpecEdge>();
                        a.add(e);
                        adjLists.put(v1, a);
                    }
                }
        
        
        pairvertices1 = new ArrayList<SpecVertex>();
        pairvertices2 = new ArrayList<SpecVertex>();
        auxvertices1 = new ArrayList<SpecVertex>();
        auxvertices2 = new ArrayList<SpecVertex>();
        deltvertices1 = new ArrayList<SpecVertex>();
        deltvertices2 = new ArrayList<SpecVertex>();
        HashMap nuclmasses = s1.nuclMasses;
        
        // pairvertices for sequence 1
         System.out.println("Creating pairvertices for left side");
        int n = vertices1.size();
        for (int i = 0; i < n; i++)
        {
            for (int j = i+1; j < n; j++)
            {
                System.out.println("Pair (" + i + "," + j + ")");
                if (vertices1.get(i).label.equalsIgnoreCase(vertices1.get(j).label))
                {
//                    if (vertices1.get(i).mass == vertices1.get(j).mass)
//                        continue;
                    
/*                    if ( (Math.abs(vertices1.get(i).mass - 345.21) < epsilon) || (Math.abs(vertices1.get(j).mass - 345.21 ) < epsilon))
                    {                        
                        SpecVertex xxx = vertices1.get(i);
                        SpecVertex yyy = vertices1.get(j);
                        System.out.println();
                    } */
                    
                    for (SpecVertex v : vertices2)
                    {
                        if (v.label.equalsIgnoreCase(vertices1.get(i).label))
                        {
                            for (Object c : alphabet1)
                            {
                                double mmm = vertices1.get(i).mass + vertices1.get(j).mass + (Double) nuclmasses.get((Character) c) - v.mass;
                                if (Math.abs(mmm) < epsilon)
                                {
                                    SpecVertex pairv = new SpecVertex(0,0,"" + vertices1.get(i).mass + "_" + vertices1.get(j).mass);
                                    pairvertices1.add(pairv);
                                    SpecEdge e = new SpecEdge(vertices1.get(i),pairv,0,1);
                                    if (adjLists.containsKey(vertices1.get(i)))
                                    {
                                        ArrayList a = adjLists.get(vertices1.get(i));
                                        a.add(e);
                                        adjLists.put(vertices1.get(i), a);
                                    }
                                    else
                                    {
                                        ArrayList<SpecEdge> a = new ArrayList<SpecEdge>();
                                        a.add(e);
                                        adjLists.put(vertices1.get(i), a);
                                    }
                                    e = new SpecEdge(vertices1.get(j),pairv,0,1);
                                    if (adjLists.containsKey(vertices1.get(j)))
                                    {
                                        ArrayList a = adjLists.get(vertices1.get(j));
                                        a.add(e);
                                        adjLists.put(vertices1.get(j), a);
                                    }
                                    else
                                    {
                                        ArrayList<SpecEdge> a = new ArrayList<SpecEdge>();
                                        a.add(e);
                                        adjLists.put(vertices1.get(j), a);
                                    }
                                    String auxvertnucl1 = "" + (Character) c;
                                    String auxvertnucl2 = "" + v.label.charAt(0);
                                    SpecVertex auxv = new SpecVertex(0,0,"aux_"+auxvertnucl1 + auxvertnucl2 + "_" + vertices1.get(i).mass + "+" + vertices1.get(j).mass + "+m_" + (Character) c + "=" + v.mass);
                                    auxvertices1.add(auxv);
                                    SpecEdge e1 = new SpecEdge(pairv,auxv,0,2);
                                    SpecEdge e2 = new SpecEdge(auxv,v,0,1);
                                    ArrayList a = new ArrayList();
                                    a.add(e1);
                                    adjLists.put(pairv, a);
                                    a = new ArrayList();
                                    a.add(e2);
                                    adjLists.put(auxv, a);
                                                                       
                                }
                            }
                        }
                    }
                }
            }
        }
        
        //pair with empty mass for sequence 1
        System.out.println("Creating pairvertices with empty mass for left side");
        SpecVertex empvert1 = new SpecVertex(0,0,"EmptyVertex1");
        SpecVertex empvert2 = new SpecVertex(0,0,"EmptyVertex2");
        
        for (int i = 0; i < n; i++)
        {
                    System.out.println("vertex " + i);
                    for (SpecVertex v : vertices2)
                    {
                        if (v.label.equalsIgnoreCase(vertices1.get(i).label))
                        {
                            for (Object c : alphabet1)
                            {
                                double mmm = vertices1.get(i).mass + (Double) nuclmasses.get((Character) c) - v.mass;
                                if (Math.abs(mmm) < epsilon)
                                {
                                    SpecVertex pairv = new SpecVertex(0,0,"" + vertices1.get(i).mass + "_emptymass");
                                    pairvertices1.add(pairv);
                                    SpecEdge e = new SpecEdge(vertices1.get(i),pairv,0,1);
                                    if (adjLists.containsKey(vertices1.get(i)))
                                    {
                                        ArrayList a = adjLists.get(vertices1.get(i));
                                        a.add(e);
                                        adjLists.put(vertices1.get(i), a);
                                    }
                                    else
                                    {
                                        ArrayList<SpecEdge> a = new ArrayList<SpecEdge>();
                                        a.add(e);
                                        adjLists.put(vertices1.get(i), a);
                                    }
                                    e = new SpecEdge(empvert1,pairv,0,1);
                                    if (adjLists.containsKey(empvert1))
                                    {
                                        ArrayList a = adjLists.get(empvert1);
                                        a.add(e);
                                        adjLists.put(empvert1, a);
                                    }
                                    else
                                    {
                                        ArrayList<SpecEdge> a = new ArrayList<SpecEdge>();
                                        a.add(e);
                                        adjLists.put(empvert1, a);
                                    }
                                    empvert1.intensity++;
                                    String auxvertnucl1 = "" + (Character) c;
                                    String auxvertnucl2 = "" + v.label.charAt(0);
                                    SpecVertex auxv = new SpecVertex(0,0,"aux_"+auxvertnucl1 + auxvertnucl2 + "_" + vertices1.get(i).mass + "+" + "emptymass" + "+m_" + (Character) c + "=" + v.mass);
                                    auxvertices1.add(auxv);
                                    SpecEdge e1 = new SpecEdge(pairv,auxv,0,2);
                                    SpecEdge e2 = new SpecEdge(auxv,v,0,1);
                                    ArrayList a = new ArrayList();
                                    a.add(e1);
                                    adjLists.put(pairv, a);
                                    a = new ArrayList();
                                    a.add(e2);
                                    adjLists.put(auxv, a);                                                                       
                                }
                            }
                        }
                    }
        }
        
        
        // pairvertices for sequence 2
        System.out.println("Creating pairvertices for right side");
        n = vertices2.size();
        for (int i = 0; i < n; i++)
        {
            for (int j = i+1; j < n; j++)
            {
                System.out.println("Pair (" + i + "," + j + ")");
                if (vertices2.get(i).label.equalsIgnoreCase(vertices2.get(j).label))
                {
//                    if (vertices2.get(i).mass == vertices2.get(j).mass)
//                        continue;
                    
/*                    if ( (Math.abs(vertices2.get(i).mass - 3349.06) < epsilon) || (Math.abs(vertices2.get(j).mass - 3349.06) < epsilon))
                    {                        
                        SpecVertex xxx = vertices2.get(i);
                        SpecVertex yyy = vertices2.get(j);
                        System.out.println();
                    }*/
                    
                    for (SpecVertex v : vertices1)
                    {
                        if (v.label.equalsIgnoreCase(vertices2.get(i).label))
                        {
                            for (Object c : alphabet1)
                            {
                                double mmm = vertices2.get(i).mass + vertices2.get(j).mass + (Double) nuclmasses.get((Character) c) - v.mass;
                                if (Math.abs(mmm) < epsilon)
                                {
                                    SpecVertex pairv = new SpecVertex(0,0,"" + vertices2.get(i).mass + "_" + vertices2.get(j).mass);
                                    pairvertices2.add(pairv);
                                    SpecEdge e = new SpecEdge(pairv,vertices2.get(i),0,1);
                                    ArrayList<SpecEdge> a = new ArrayList<SpecEdge>();
                                    a.add(e);
                                    adjLists.put(pairv, a);
                                    e = new SpecEdge(pairv,vertices2.get(j),0,1);
                                    a = adjLists.get(pairv);
                                    a.add(e);
                                    adjLists.put(pairv, a);
                                    
                                    String auxvertnucl1 = "" + v.label.charAt(0);
                                    String auxvertnucl2 = "" + (Character) c;                                    
                                    SpecVertex auxv = new SpecVertex(0,0,"aux_" + auxvertnucl1 + auxvertnucl2 + "_" + v.mass + "=" + vertices2.get(i).mass + "+" + vertices2.get(j).mass + "+m_" + (Character) c);
                                    auxvertices2.add(auxv);
                                    SpecEdge e1 = new SpecEdge(v,auxv,0,1);
                                    SpecEdge e2 = new SpecEdge(auxv,pairv,0,2);
                                    if (adjLists.containsKey(v))
                                    {
                                        a = adjLists.get(v);
                                        a.add(e1);
                                        adjLists.put(v, a);
                                    }
                                    else
                                    {
                                        a = new ArrayList();
                                        a.add(e1);
                                        adjLists.put(v, a);
                                    }
                                    a = new ArrayList();
                                    a.add(e2);
                                    adjLists.put(auxv, a);
         
                                }
                            }
                        }
                    }
                }
            }
        }
        
        //pair with empty mass for sequence 2
         System.out.println("Creating pairvertices with empty mass for right side");
        n = vertices2.size();
        for (int i = 0; i < n; i++)
        {
                    System.out.println("vertex " + i);
                    for (SpecVertex v : vertices1)
                    {
                        if (v.label.equalsIgnoreCase(vertices2.get(i).label))
                        {
                            for (Object c : alphabet1)
                            {
                                double mmm = vertices2.get(i).mass + (Double) nuclmasses.get((Character) c) - v.mass;
                                if (Math.abs(mmm) < epsilon)
                                {
                                    SpecVertex pairv = new SpecVertex(0,0,"" + vertices2.get(i).mass + "_emptymass");
                                    pairvertices2.add(pairv);
                                    SpecEdge e = new SpecEdge(pairv,vertices2.get(i),0,1);
                                    ArrayList<SpecEdge> a = new ArrayList<SpecEdge>();
                                    a.add(e);
                                    e = new SpecEdge(pairv,empvert2,0,1);
                                    a.add(e);
                                    adjLists.put(pairv, a);
                                    empvert2.intensity++;
                                    
                                    String auxvertnucl1 = "" + v.label.charAt(0);
                                    String auxvertnucl2 = "" + (Character) c;                                    
                                    SpecVertex auxv = new SpecVertex(0,0,"aux_" + auxvertnucl1 + auxvertnucl2 + "_" + v.mass + "=" + vertices2.get(i).mass + "+" + "emptymass" + "+m_" + (Character) c);
                                    auxvertices2.add(auxv);
                                    SpecEdge e1 = new SpecEdge(v,auxv,0,1);
                                    SpecEdge e2 = new SpecEdge(auxv,pairv,0,2);
                                    if (adjLists.containsKey(v))
                                    {
                                        a = adjLists.get(v);
                                        a.add(e1);
                                        adjLists.put(v, a);
                                    }
                                    else
                                    {
                                        a = new ArrayList();
                                        a.add(e1);
                                        adjLists.put(v, a);
                                    }
                                    a = new ArrayList();
                                    a.add(e2);
                                    adjLists.put(auxv, a);
         
                                }
                            }
                        }
                    }
        }
        
        // delta vertices
        System.out.println("Creating delta vertices");
        for (SpecVertex v1 : vertices1)
            for (SpecVertex v2 : vertices2)
            {
                System.out.println("Pair (" + v1.number + "," + v2.number + ")");
                if (v1.label.equalsIgnoreCase(v2.label))
                {
                    for (Object o1 : alphabet1)
                    {
                        boolean found = false;
                        for (Object o2 : alphabet2)
                        {
                            char c1 = (Character) o1;
                            char c2 = (Character) o2;
                            if (c1 == c2)
                                continue;
                            double mm = v1.mass - (Double) nuclmasses.get((Character) o1) - v2.mass + (Double) nuclmasses.get((Character) o2);
                            if (Math.abs(mm) < epsilon)
                            {
                                found = true;
                                SpecVertex deltv1 = new SpecVertex(0,0,"delta_1_" + (Character) o1 + (Character) o2 + "_" + v1.mass + "-m_" + (Character) o1 + "=" +  v2.mass + "-m_" + (Character) o2);
                                deltvertices1.add(deltv1);
                                SpecEdge e = new SpecEdge(v1,deltv1,0,1);
                                if (adjLists.containsKey(v1))
                                {
                                    ArrayList a = adjLists.get(v1);
                                    a.add(e);
                                    adjLists.put(v1, a);
                                }
                                else
                                {
                                    ArrayList<SpecEdge> a = new ArrayList<SpecEdge>();
                                    a.add(e);
                                    adjLists.put(v1, a);
                                }
                                SpecVertex deltv2 = new SpecVertex(0,0,"delta_2_" + (Character) o1 + (Character) o2 + "_" + v1.mass + "-m_" + (Character) o1 + "=" +  v2.mass + "-m_" + (Character) o2);
                                deltvertices2.add(deltv2);
                                e = new SpecEdge(deltv1,deltv2,0,0,"delta");
                                ArrayList<SpecEdge> a = new ArrayList<SpecEdge>();
                                a.add(e);
                                adjLists.put(deltv1, a);
                                e = new SpecEdge(deltv2,v2,0,1);
                                a = new ArrayList<SpecEdge>();
                                a.add(e);
                                adjLists.put(deltv2, a);
                                break;
                            }
                        }
                        if (found)
                            break;
                    }
                }
            }
        
        vertices1.add(empvert1);
        vertices2.add(empvert2);
        
        vertices = new ArrayList<SpecVertex>();
        vertices.addAll(vertices1);
        vertices.addAll(vertices2);
        vertices.addAll(pairvertices1);
        vertices.addAll(pairvertices2);
        vertices.addAll(auxvertices1);
        vertices.addAll(auxvertices2);
        vertices.addAll(deltvertices1);
        vertices.addAll(deltvertices2);
        
        System.out.println("auxvertices1 = " + auxvertices1.size());
        System.out.println("auxvertices2 = " + auxvertices2.size());
        System.out.println("deltvertices1 = " + deltvertices1.size());
        System.out.println("deltvertices2 = " + deltvertices2.size());
        
        int ii = 0;
        int jj = 0;
        for (SpecVertex v1 : auxvertices1)
        {
            ii++;
            jj=0;
            for (SpecVertex v2 : auxvertices2)
            {
                jj++;
                System.out.println("Pair (" + ii + "," + jj + ")/" + "(" + auxvertices1.size() + "," + auxvertices2.size() + ")");
                if (v1.label.substring(4, 6).equalsIgnoreCase(v2.label.substring(4, 6)))
                {
                    SpecEdge e = new SpecEdge(v1,v2,1,1);
                    ArrayList a = adjLists.get(v1);
                    a.add(e);
                    adjLists.put(v1, a);
                }
            }
        }
        int r = deltvertices1.size();
        for (int i = 0; i < r; i++)
            for (int j = i+1; j < r; j++)
            {
                SpecVertex v1 = deltvertices1.get(i);
                SpecVertex v2 = deltvertices1.get(j);
                if (v1.label.substring(8, 10).equalsIgnoreCase(v2.label.substring(8, 10)))
                {
                    SpecVertex w1 = new SpecVertex();
                    SpecVertex w2 = new SpecVertex();
                    for (Object o : adjLists.get(v1))
                        if (deltvertices2.contains(((SpecEdge) o).end))
                        {
                            w1 = ((SpecEdge) o).end;
                            break;
                        }
                    for (Object o : adjLists.get(v2))
                        if (deltvertices2.contains(((SpecEdge) o).end))
                        {
                            w2 = ((SpecEdge) o).end;
                            break;
                        }
                    SpecEdge e = new SpecEdge(v1,w2,0,1);
                    ArrayList a = adjLists.get(v1);
                    a.add(e);
                    adjLists.put(v1, a);
                    e = new SpecEdge(v2,w1,0,1);
                    a = adjLists.get(v2);
                    a.add(e);
                    adjLists.put(v2, a);
                    
                }
            }
        
        SpecVertex source = new SpecVertex(0,0,"source");
        SpecVertex sink = new SpecVertex(0,0,"sink");
        ArrayList a = new ArrayList();
        for (SpecVertex v : vertices1)
            a.add(new SpecEdge(source,v,0,v.intensity));
        adjLists.put(source, a);
        for (SpecVertex v : vertices2)
        {
            a = new ArrayList();
            a.add(new SpecEdge(v,sink,0,v.intensity));
            adjLists.put(v, a);
        }  
        vertices.add(source);
        vertices.add(sink);
        for (SpecVertex v : vertices)
            if (!adjLists.containsKey(v))
                adjLists.put(v, null);
        
        // delete unnecessary auxiliary vertices
        ArrayList<SpecVertex> toRemove = new ArrayList<SpecVertex>();
        for (SpecVertex v : pairvertices1)
        {
            SpecVertex u = ((SpecEdge) adjLists.get(v).get(0)).end;
            ArrayList aaa = adjLists.get(u);
            if (adjLists.get(u).size() == 1)
            {
                adjLists.remove(u);
                adjLists.remove(v);
                toRemove.add(u);
                toRemove.add(v);
                
            }
        }
        
        for (SpecVertex v : auxvertices2)
        {
            int deg = 0;
            for (SpecVertex w : auxvertices1)
            {
                ArrayList adl = adjLists.get(w);
                if (adl == null)
                    continue;
                for (Object o : adl)
                {
                    SpecEdge e = (SpecEdge) o;
                    if (e.end.equals(v))
                        deg++;
                }
            }
            if (deg == 0)
            {
                toRemove.add(v);
                toRemove.add(( (SpecEdge) adjLists.get(v).get(0)).end);
            }
        }
        
        for (SpecVertex v : deltvertices1)
        {
            if (adjLists.get(v).size() == 1)
            {
                SpecVertex u = ((SpecEdge) adjLists.get(v).get(0)).end;
                adjLists.remove(u);
                adjLists.remove(v);
                toRemove.add(u);
                toRemove.add(v);
                
            }
        }
        vertices.removeAll(toRemove);
        pairvertices1.removeAll(toRemove);
        auxvertices1.removeAll(toRemove);
        pairvertices2.removeAll(toRemove);
        auxvertices2.removeAll(toRemove);
        deltvertices1.removeAll(toRemove);
        deltvertices2.removeAll(toRemove);
        for (SpecVertex v : vertices)
        {
            ArrayList<SpecEdge> adj = adjLists.get(v);
            ArrayList<SpecEdge> adj1 = new ArrayList<SpecEdge>();
            if (adj == null)
                continue;
            for (SpecEdge e : adj)
            {
                if (!toRemove.contains(e.end))
                    adj1.add(e);
            }
            if (adj1.size() > 0)
                adjLists.put(v, adj1);
            else
                adjLists.put(v, null);
        }
        for (SpecVertex v : toRemove)
            adjLists.remove(v);
        Set s = adjLists.entrySet();
        int j = 1;
        int iedge = 1;
        for (Object o : s)
        {
            Map.Entry me = (Map.Entry) o;
            ((SpecVertex) me.getKey()).number = j;
            j++;
            if (me.getValue() == null)
                continue;
            Iterator itr = ((ArrayList) me.getValue()).iterator();
            while (itr.hasNext())
            {
                SpecEdge e = (SpecEdge) itr.next();
                e.number = iedge;
                iedge++;
            }
            
        }
    }
    
    void printToFile(String addr) throws IOException
    {
        FileWriter fw = new FileWriter(addr);
                
        Set s = adjLists.entrySet();
        int i = 1;
        for (Object o : s)
        {
            Map.Entry me = (Map.Entry) o;
            fw.write("vertex_" + ((SpecVertex) me.getKey()).number + " : " + ((SpecVertex) me.getKey()).label + "_" + ((SpecVertex) me.getKey()).mass + "\n");
            if (me.getValue() != null)
            {
                Iterator ir = ((ArrayList) me.getValue()).iterator();
                while (ir.hasNext())
                {
                    SpecEdge e = (SpecEdge) ir.next();
                    fw.write("\t\t" + "Vertex_" + e.end.number + " " + e.end.label + "_" + e.end.mass + " " + e.capacity + " " + e.cost + " " + i + "\n");
                    i++;
                }
            }
            else
            {
                fw.write("\t\t empty\n");
            }
        }
        fw.close();
        
    }
    void printMaxFlowProblemMPS(String addr) throws IOException
    {
        FileWriter fw = new FileWriter(addr);
        fw.write("NAME\n");
        fw.write("ROWS\n");
        fw.write(" N OBJ\n");
        for (SpecVertex v : vertices)
            fw.write(" E flow_" + v.label + "\n");
        fw.write(" E flow_source" + "\n");
        fw.write(" E flow_sink" + "\n");
        fw.write("COLUMNS\n");
        fw.write("\t MARK0000 \t 'MARKER' \t'INTORG'\n");
        for (SpecVertex v : vertices1)
        {
            fw.write("\t x_s_" + v.label + "\t" + "flow_source" + "\t" + 1.0 + "\n");
            fw.write("\t x_s_" + v.label + "\t" + "flow_" + v.label + "\t" + 1.0 + "\n");
            if (!adjLists.containsKey(v))
                continue;
            ArrayList a = adjLists.get(v);
            for (Object o:a)
            {
                SpecEdge e = (SpecEdge) o;
                fw.write("\t x_" +e.begin.label + "_" + e.end.label + "\t" + "flow_" + e.begin.label + "\t" + (-1.0) + "\n");
                fw.write("\t x_" +e.begin.label + "_" + e.end.label + "\t" + "flow_" + e.end.label + "\t" + 1.0 + "\n");                
            }
        }
        for (SpecVertex v : vertices2)
        {
            fw.write("\t x_" + v.label + "_t" + "\t" + "flow_sink" + "\t" + 1.0 + "\n");
            fw.write("\t x_" + v.label + "_t" + "\t" + "flow_" + v.label + "\t" + (-1.0) + "\n");
        }
        for (SpecVertex v : pairvertices1)
        {
            if (!adjLists.containsKey(v))
                continue;
            ArrayList a = adjLists.get(v);
            for (Object o:a)
            {
                SpecEdge e = (SpecEdge) o;
                fw.write("\t x_" +e.begin.label + "_" + e.end.label + "\t" + "flow_" + e.begin.label + "\t" + (-1.0) + "\n");
                fw.write("\t x_" +e.begin.label + "_" + e.end.label + "\t" + "flow_" + e.end.label + "\t" + 1.0 + "\n");                
            }
        }
        for (SpecVertex v : pairvertices2)
        {
            if (!adjLists.containsKey(v))
                continue;
            ArrayList a = adjLists.get(v);
            for (Object o:a)
            {
                SpecEdge e = (SpecEdge) o;
                fw.write("\t x_" +e.begin.label + "_" + e.end.label + "\t" + "flow_" + e.begin.label + "\t" + (-1.0) + "\n");
                fw.write("\t x_" +e.begin.label + "_" + e.end.label + "\t" + "flow_" + e.end.label + "\t" + 1.0 + "\n");                
            }
        }
        for (SpecVertex v : auxvertices1)
        {
            if (!adjLists.containsKey(v))
                continue;
            ArrayList a = adjLists.get(v);
            for (Object o:a)
            {
                SpecEdge e = (SpecEdge) o;
                if (e.end.label.startsWith("aux"))
                    fw.write("\t x_" +e.begin.label + "_" + e.end.label + "\t" + "OBJ" + "\t" + e.cost + "\n");
                fw.write("\t x_" +e.begin.label + "_" + e.end.label + "\t" + "flow_" + e.begin.label + "\t" + (-1.0) + "\n");
                fw.write("\t x_" +e.begin.label + "_" + e.end.label + "\t" + "flow_" + e.end.label + "\t" + 1.0 + "\n");                
            }
        }
        for (SpecVertex v : auxvertices2)
        {
            if (!adjLists.containsKey(v))
                continue;
            ArrayList a = adjLists.get(v);
            for (Object o:a)
            {
                SpecEdge e = (SpecEdge) o;
                if (e.end.label.startsWith("aux"))
                    fw.write("\t x_" +e.begin.label + "_" + e.end.label + "\t" + "OBJ" + "\t" + e.cost + "\n");
                fw.write("\t x_" +e.begin.label + "_" + e.end.label + "\t" + "flow_" + e.begin.label + "\t" + (-1.0) + "\n");
                fw.write("\t x_" +e.begin.label + "_" + e.end.label + "\t" + "flow_" + e.end.label + "\t" + 1.0 + "\n");                
            }
        }
        
        fw.write("\t MARK0001 \t 'MARKER' \t 'INTEND'\n");
        fw.write("RHS\n");
        for (SpecVertex v : vertices)
            fw.write("\t flow_" + v.label + "\t" + 0.0 + "\n");
        
        fw.write("BOUNDS\n");
        for (SpecVertex v : vertices)
        {
            if (!adjLists.containsKey(v))
                continue;
            ArrayList a = adjLists.get(v);
            for (Object o : a)
            {
                SpecEdge e = (SpecEdge) o;
                fw.write("\t LO BND \t x_" +e.begin.label + "_" + e.end.label + 0.0 + "\n");
                fw.write("\t UP BND \t x_" +e.begin.label + "_" + e.end.label + e.capacity + "\n");
                
            }
        }
        fw.write("ENDDATA");
        fw.close();
    }
    void printMaxFlowProblemCPLEX(String addr) throws IOException
    {
        FileWriter fw = new FileWriter(addr);
        fw.write("minimize\n");
        fw.write("\t obj: ");
        Set s = adjLists.entrySet();
        int iedge = 1;
        String objfun = "";
        for (Object o : s)
        {
            Map.Entry me = (Map.Entry) o;
            if (me.getValue() == null)
                continue;
            Iterator ir = ((ArrayList) me.getValue()).iterator();
            while (ir.hasNext())
            {
                SpecEdge e = (SpecEdge) ir.next();
                if (e.cost > 0)
                    objfun += " + " + e.cost + " x" + iedge; 
                iedge++;
            }
        }
        objfun = objfun.substring(3);
        fw.write(objfun + "\n");
        fw.write("Subject To\n");
        
        HashMap<SpecVertex,String> flowCons = new HashMap<SpecVertex,String>();
        HashMap<SpecVertex,String> eqPairCons = new HashMap<SpecVertex,String>();
        HashMap<SpecEdge,String> eqDeltaCons = new HashMap<SpecEdge,String>();
        s = adjLists.entrySet();
        iedge = 1;
        for (Object o : s)
        {
            Map.Entry me = (Map.Entry) o;
            if (me.getValue() == null)
                continue;
            SpecVertex v = (SpecVertex) me.getKey();
//            System.out.println(v.number);
            Iterator ir = ((ArrayList) me.getValue()).iterator();
            while (ir.hasNext())
            {
                SpecEdge e = (SpecEdge) ir.next();
                SpecVertex u = e.end;
                if (!(v.label.equalsIgnoreCase("source") && u.label.equalsIgnoreCase("EmptyVertex1")))
                {
                    if (flowCons.containsKey(v))
                    {
                        String constr = flowCons.get(v);
                        constr += " - x" + iedge;
                        flowCons.put(v, constr);                    
                    }
                    else
                    {
                        String constr = " - x" + iedge;
                        flowCons.put(v, constr);
                    }
                }
                if (flowCons.containsKey(u))
                {
                    String constr = flowCons.get(u);
                    constr += " + x" + iedge;
                    flowCons.put(u, constr);                    
                }
                else
                {
                    String constr = "x" + iedge;
                    flowCons.put(u, constr);
                }
                
                if (pairvertices1.contains(u))
                {
                    if (!eqPairCons.containsKey(u))
                    {
                        String constr = "x" + iedge;
                        eqPairCons.put(u, constr);
                    }
                    else
                    {
                        String constr = eqPairCons.get(u);
                        constr+= " - x" + iedge;
                        eqPairCons.put(u, constr);
                    }
                } 
                if (pairvertices2.contains(v))
                {
                    if (!eqPairCons.containsKey(v))
                    {
                        String constr = "x" + iedge;
                        eqPairCons.put(v, constr);
                    }
                    else
                    {
                        String constr = eqPairCons.get(v);
                        constr+= " - x" + iedge;
                        eqPairCons.put(v, constr);
                    }
                }
                
                if (deltvertices1.contains(v) && deltvertices2.contains(u) && e.label.equalsIgnoreCase("delta"))
                {
                    SpecEdge e1 = (SpecEdge) adjLists.get(u).get(0);
                    String constr = "x" + e1.number;
                    for (SpecVertex w : vertices1)
                    {
                        boolean found = false;
                        if (adjLists.get(w) == null)
                            continue;
                        for (Object o1 : adjLists.get(w))
                        {
                            SpecEdge e2 = (SpecEdge) o1;
                            if (e2.end.equals(v))
                            {
                                constr+= " - x" + e2.number;
                                found = true;
                                break;
                            }
                        }
                        if (found)
                            break;
                    }
                    eqDeltaCons.put(e, constr);
                }
                
                iedge++;
            }
        }

        
        int nnonisol = 0;
        for (SpecVertex v : vertices1)
            if ((adjLists.get(v) != null) && (!v.label.equalsIgnoreCase("EmptyVertex1")))
                nnonisol+=v.intensity;
        
        iedge = 1;
        s = adjLists.entrySet();
        
        for (Object o : s)
        {
            Map.Entry me = (Map.Entry) o;
            SpecVertex v = (SpecVertex) me.getKey();
            String constr = flowCons.get(v);
            if (v.label.equalsIgnoreCase("source"))                
                fw.write("\tflowconstr" + v.number + ": " + constr + " = " + (-nnonisol) + "\n");
            else
                if (v.label.equalsIgnoreCase("sink"))                
//                    fw.write("\tflowconstr" + v.number + ": " + constr + " = " + nnonisol + "\n");
                    continue;
                else
                    fw.write("\tflowconstr" + v.number + ": " + constr + " = " + 0 + "\n");
        }
        
        for (SpecVertex v : pairvertices1)
        {
            String constr = eqPairCons.get(v);   
            fw.write("\teqpairconstr" + v.number + ": " + constr + " = " + 0 + "\n");
        }
        for (SpecVertex v : pairvertices2)
        {
            String constr = eqPairCons.get(v);   
            fw.write("\teqpairconstr" + v.number + ": " + constr + " = " + 0 + "\n");
        }
        
        
        s = eqDeltaCons.entrySet();
        for (Object o : s)
        {
            Map.Entry me = (Map.Entry) o;
            SpecEdge e = (SpecEdge) me.getKey();
            String constr = (String) me.getValue();
            fw.write("\teqdeltconstr" + e.number + ": " + constr + " = " + 0 + "\n");
        }
        
        fw.write("Bounds\n");
        iedge = 1;
        s = adjLists.entrySet();
        for (Object o : s)
        {
            Map.Entry me = (Map.Entry) o;
            if (me.getValue() == null)
                continue;
            Iterator ir = ((ArrayList) me.getValue()).iterator();
            while (ir.hasNext())
            {
                SpecEdge e = (SpecEdge) ir.next();
                fw.write("\t0 <= x" + iedge + " <= " + e.capacity + "\n");
                iedge++;
            }
        }
        
        fw.write("Integer\n");
        iedge = 1;
        s = adjLists.entrySet();
        for (Object o : s)
        {
            Map.Entry me = (Map.Entry) o;
            if (me.getValue() == null)
                continue;
            Iterator ir = ((ArrayList) me.getValue()).iterator();
            while (ir.hasNext())
            {
                SpecEdge e = (SpecEdge) ir.next();
                fw.write("\tx" + iedge + "\n");
                iedge++;
            }
        } 
        
        fw.write("End");
        fw.close();
    }
    void printForPajek(String addr) throws IOException
    {
        FileWriter fw = new FileWriter(addr);
        fw.write("*Vertices " + vertices.size() + "\n");
        VertComparator vc = new VertComparator();
        Collections.sort(vertices, vc);
        for (SpecVertex v : vertices)
        {
            fw.write(v.number + " " + v.label + "_" + v.mass + "\n");
        }
        fw.write("*arcs\n");
        for (SpecVertex v : vertices)
        {
            if (v.label.equalsIgnoreCase("EmptyVertex1"))
                System.out.println();
            ArrayList al = adjLists.get(v);
            if (al == null)
                continue;
            for (Object o : al)
            {
                SpecEdge e = (SpecEdge) o;
                if (auxvertices1.contains(e.begin) && auxvertices2.contains(e.end))
                    continue;
                if (e.capacity > 0)
                    fw.write(e.begin.number + " " + e.end.number + " " + e.capacity + "\n");
                else
                    fw.write(e.begin.number + " " + e.end.number + " " + (-1) + "\n");
            }
        }
        fw.write("*Edges\n");
        for (SpecVertex v : auxvertices1)
        {
            ArrayList al = adjLists.get(v);
            if (al == null)
                continue;
            for (Object o : al)
            {
                SpecEdge e = (SpecEdge) o;
                if (auxvertices2.contains(e.end))
                    fw.write(e.begin.number + " " + e.end.number + " " + e.capacity + "\n");
            }
        }
        
        
        fw.close();
    }
    void printForPajekPartition(String addr) throws IOException
    {
        FileWriter fw = new FileWriter(addr);
        fw.write("*Partition Part\n");
        fw.write("*Vertices " + vertices.size() + "\n");
        VertComparator vc = new VertComparator();
        Collections.sort(vertices, vc);
        for (SpecVertex v : vertices)
        {
            int npart = 0;
            if (v.label.equalsIgnoreCase("source"))
                npart = 1;
            if (vertices1.contains(v))
                npart = 2;
            if (pairvertices1.contains(v))
                npart = 3;
            if (auxvertices1.contains(v))
                npart = 4;
            if (auxvertices2.contains(v))
                npart = 5;
            if (pairvertices2.contains(v))
                npart = 6;
            if (vertices2.contains(v))
                npart = 7;
            if (v.label.equalsIgnoreCase("sink"))
                npart = 8;
            if (deltvertices1.contains(v))
                npart = 9;
            if (deltvertices2.contains(v))
                npart = 10;
            fw.write(npart + "\n");
        }
        fw.close();
    }
}
