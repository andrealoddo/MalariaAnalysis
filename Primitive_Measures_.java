import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.filter.*;
import ij.plugin.filter.PlugInFilter;
import ij.measure.*;
import ij.plugin.MeasurementsWriter;
import java.util.*;
import javax.swing.JOptionPane;

public class Primitive_Measures_ implements PlugInFilter {
    ImagePlus imp;
    static String title = "prova";

    protected ResultsTable rt = ResultsTable.getResultsTable();
    boolean doArea;
    boolean doPerimeter;

    boolean doFeret;
    boolean doBreadth;

    boolean doAspRatio;
    boolean doCirc;
    boolean doRoundness;
    boolean doArEquivD;
    boolean doPerEquivD;
    boolean doEquivEllAr;
    boolean doCompactness;
    boolean doShape;

    boolean doIS;
    boolean doMidpoint;
    boolean doDS;

    int [] feret_exstreme_points = {0,0,0,0};
    int [] breadth_exstreme_points = {0,0,0,0};


    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;

        GenericDialog gd = new GenericDialog("Primitive_Measures", IJ.getInstance());
        gd.addStringField("Title: ", title);
        gd.addMessage("Choise measures");

        gd.addCheckbox("SelectAll",false);

        gd.addCheckbox("Area",true);
        gd.addCheckbox("Perimeter",false);


        gd.addCheckbox("Feret",true);
        gd.addCheckbox("Breadth",true);

        gd.addCheckbox("AspRatio",false);
        gd.addCheckbox("Circ", false);
        gd.addCheckbox("Roundness", false);
        gd.addCheckbox("ArEquivD", false);
        gd.addCheckbox("PerEquivD", false);
        gd.addCheckbox("EquivEllAr", false);
        gd.addCheckbox("Compactness", false);
        gd.addCheckbox("Shape", false);

        gd.addCheckbox("IS", false);
        gd.addCheckbox("CG-Midpoint",false);
        gd.addCheckbox("DS-Distance between IS and CG", false);

        //gd.addCheckbox("Show me feret and breadth points.\n Please check Feret and breadth", false);


        gd.showDialog();

        if (gd.wasCanceled())
            return DONE;

        boolean doSelectAll=gd.getNextBoolean ();

        if(doSelectAll){
            doArea = true;
            doPerimeter = true;
            doFeret= true;
            doBreadth= true;
            doAspRatio = true;
            doCirc = true;
            doRoundness = true;
            doArEquivD = true;
            doPerEquivD = true;
            doEquivEllAr = true ;
            doCompactness = true;
            doIS = true;
            doMidpoint = true;
            doDS= true;
        }else{
            doArea = gd.getNextBoolean ();
            doPerimeter = gd.getNextBoolean ();
            doFeret= gd.getNextBoolean ();
            doBreadth= gd.getNextBoolean ();
            doAspRatio=gd.getNextBoolean();
            doCirc=gd.getNextBoolean();
            doRoundness=gd.getNextBoolean();
            doArEquivD=gd.getNextBoolean();
            doPerEquivD=gd.getNextBoolean();
            doEquivEllAr=gd.getNextBoolean();
            doCompactness=gd.getNextBoolean();
            doIS=gd.getNextBoolean();
            doMidpoint = gd.getNextBoolean ();
            doDS = gd.getNextBoolean ();

        }

        return DOES_ALL;

    }

    public void run(ImageProcessor ip) {
        rt.reset();

        rt.incrementCounter();

        int xe = ip.getWidth(); //larghezza e altezza dell'intera immagine
        int ye = ip.getHeight();

        feret_exstreme_points[3]=ye;
        breadth_exstreme_points[2]=xe;

        int [] feret_and_breadth = {breadth(ip, xe, ye), feret(ip, xe, ye)};

        if(doArea){
            int area = area(ip, xe, ye);
            rt.addValue("Area of element", area);
        }

        if(doPerimeter){
            double perimeter =perimeter(ip, xe, ye);
            rt.addValue("Perimeter of element", perimeter);
        }

        if(doMidpoint){//center
            double [] midpoint = midpoint(ip, xe, ye);
            rt.addValue("CG Midpoint X", midpoint[0]);
            rt.addValue("CG Midpoint Y", midpoint[1]);
        }

        if(doFeret){
            rt.addValue("Height/Feret of Element", feret_and_breadth[1]);


        }

        if(doBreadth){
            rt.addValue("Width/Breadth of Element", feret_and_breadth[0]);

        }

        if(doAspRatio){
            double aspRatio = aspRatio(feret_and_breadth);
            rt.addValue("AspRatio", aspRatio);
        }

        if(doCirc){
            double circ = circ((area(ip, xe, ye)), perimeter(ip, xe, ye));
            rt.addValue("Circ", circ);

        }

        if(doRoundness){
            double roundness = roundness(area(ip, xe, ye), feret_and_breadth[1]);
            rt.addValue("Roundness", roundness);

        }

        if(doArEquivD){
            double arEquivD= arEquivD(area(ip, xe, ye));
            rt.addValue("ArEquivD", arEquivD);
        }

        if(doPerEquivD){
            double perEquivD = perEquivD(area(ip, xe, ye));
            rt.addValue("PerEquivD", perEquivD);
        }

        if(doEquivEllAr){
            double equivEllAr = equivEllAr(feret_and_breadth);
            rt.addValue("EquivEllAr", equivEllAr);
        }

        if(doCompactness){
            double compactness = compactness(area(ip, xe, ye), (feret_and_breadth[1]));
            rt.addValue("Compactness", compactness);
        }

        if(doShape){
            double shape = shape (area(ip, xe, ye), perimeter(ip, xe, ye));
            rt.addValue("Shape", shape);
        }


        if(doIS){
            double [] cS = iS (feret_exstreme_points, breadth_exstreme_points);
            rt.addValue("IS-CSx", cS[0]);
            rt.addValue("IS-CSy", cS[1]);

            ip.putPixelValue((int)cS[0],(int)cS[1], 130);
        }

        if(doDS){
            double dS = dS (iS(feret_exstreme_points, breadth_exstreme_points), midpoint(ip, xe, ye));
            rt.addValue("DS", dS);
        }

        for(int i =0; i < 4; i++){
            rt.addValue("feret"+i, feret_exstreme_points[i]);
        }

        for(int i =0; i < 4; i++){
            rt.addValue("breadth"+i, breadth_exstreme_points[i]);
        }

        rt.addResults();
        rt.updateResults();
        rt.show("Results");



    }

    public int area(ImageProcessor ip_for_area, int xe, int ye){
        int area=0;
        for(int x=0; x<xe; x++ ){
            for(int y=0; y<ye; y++) {
                //se il pixel è bianco, fa parte dell'oggetto da analizzare dunque viene conteggiato come "pixel area"
                if (ip_for_area.getPixel(x, y) == 255) {
                    area++;
                }
            }
        }

        return area;
    }

    public double perimeter(ImageProcessor ip_for_perimeter, int xe, int ye){
        int internal_perimeter=0;
        int external_perimeter=0;

        for(int x=0; x<xe; x++ ){
            for(int y=0; y<ye; y++){
                if(ip_for_perimeter.getPixel(x,y)==255){
                    if(ip_for_perimeter.getPixel(x-1,y)==0){
                        internal_perimeter++;
                    } else if (ip_for_perimeter.getPixel(x+1,y)==0) {
                        internal_perimeter++;
                    }else if (ip_for_perimeter.getPixel(x,y-1)==0){
                        internal_perimeter++;
                    } else if (ip_for_perimeter.getPixel(x,y+1)==0){
                        internal_perimeter++;
                    }

                }else {
                    if(ip_for_perimeter.getPixel(x-1,y)==255){
                        external_perimeter++;
                    } else if (ip_for_perimeter.getPixel(x+1,y)==255) {
                        external_perimeter++;
                    }else if (ip_for_perimeter.getPixel(x,y-1)==255){
                        external_perimeter++;
                    } else if (ip_for_perimeter.getPixel(x,y+1)==255){
                        external_perimeter++;
                    }
                }

            }
        }


        return ((double)external_perimeter+(double)internal_perimeter)/2;
    }


    public int feret (ImageProcessor ip_for_feret, int xe, int ye){

        for(int x=0; x<xe; x++ ){
            for(int y=0; y<ye; y++) {
                if (ip_for_feret.getPixel(x, y) == 255) {
                    if (ip_for_feret.getPixel(x - 1, y) == 0) {
                        feret_exstreme_points=feret_exstreme_points(x,y, feret_exstreme_points);
                    } else if (ip_for_feret.getPixel(x + 1, y) == 0) {
                        feret_exstreme_points=feret_exstreme_points(x,y, feret_exstreme_points);
                    } else if (ip_for_feret.getPixel(x, y - 1) == 0) {
                        feret_exstreme_points=feret_exstreme_points(x,y, feret_exstreme_points);
                    } else if (ip_for_feret.getPixel(x, y + 1) == 0) {
                        feret_exstreme_points=feret_exstreme_points(x,y, feret_exstreme_points);
                    }

                }
            }
        }

        if(itIsFeret(((feret_exstreme_points[1]-feret_exstreme_points[3])+1), ((breadth_exstreme_points[0]-breadth_exstreme_points[2])+1)))
        {
            return ((feret_exstreme_points[1]-feret_exstreme_points[3])+1);
        }else{
            int [] tt = feret_exstreme_points;
            feret_exstreme_points=breadth_exstreme_points;
            breadth_exstreme_points=tt;
            return (feret_exstreme_points[1]-feret_exstreme_points[3])+1; //ricalcola
        }


    }

    public int breadth (ImageProcessor ip_for_breadth, int xe, int ye) {

            for (int x = 0; x < xe; x++) {
                for (int y = 0; y < ye; y++) {
                    if (ip_for_breadth.getPixel(x, y) == 255) {
                        if (ip_for_breadth.getPixel(x - 1, y) == 0) {
                            breadth_exstreme_points = breadth_exstreme_points(x, y, breadth_exstreme_points);
                        } else if (ip_for_breadth.getPixel(x + 1, y) == 0) {
                            breadth_exstreme_points = breadth_exstreme_points(x, y, breadth_exstreme_points);
                        } else if (ip_for_breadth.getPixel(x, y - 1) == 0) {
                            breadth_exstreme_points = breadth_exstreme_points(x, y, breadth_exstreme_points);
                        } else if (ip_for_breadth.getPixel(x, y + 1) == 0) {
                            breadth_exstreme_points = breadth_exstreme_points(x, y, breadth_exstreme_points);
                        }

                    }
                }
            }


            if (itIsFeret(((feret_exstreme_points[1] - feret_exstreme_points[3]) + 1), ((breadth_exstreme_points[0] - breadth_exstreme_points[2]) + 1))) {
                int[] tt = feret_exstreme_points;
                feret_exstreme_points = breadth_exstreme_points;
                breadth_exstreme_points = tt;
                return (breadth_exstreme_points[0] - breadth_exstreme_points[2]) + 1;//ricalcola
            } else {
                return (breadth_exstreme_points[0] - breadth_exstreme_points[2]) + 1;
            }

    }


    int[] feret_exstreme_points(int x, int y, int [] feret_exstreme_points){
        if(feret_exstreme_points[1]<y){
            feret_exstreme_points[1]=y;
            feret_exstreme_points[0]=x;
        }

        if(feret_exstreme_points[3]>y){
            feret_exstreme_points[3]=y;
            feret_exstreme_points[2]=x;
        }

        return feret_exstreme_points;
    }

    int[] breadth_exstreme_points(int x, int y, int [] breadth_exstreme_points){
        if(breadth_exstreme_points[0]<x){
            breadth_exstreme_points[1]=y;
            breadth_exstreme_points[0]=x;
        }
        if(breadth_exstreme_points[2]>x){
            breadth_exstreme_points[3]=y;
            breadth_exstreme_points[2]=x;
        }
        return breadth_exstreme_points;
    }

    double aspRatio (int [] feret_and_breadth){
        return feret_and_breadth[1]/feret_and_breadth[0];
    }

    boolean itIsFeret (int feret, int breadth){
        return feret>=breadth;
    }

    double roundness (int area, int feret){
        double pigreco= 3.1415926535;

        return (4 * (double)area )/ pigreco * ((double) feret * feret);

    }

    double circ (int area, double perim){
        double pigreco= 3.1415926535;
        return (4*(double)area)/(perim*perim);
    }

    double arEquivD(int area){
        double pigreco= 3.1415926535;
        return Math.sqrt((4*pigreco)*(double)area);

    }

    double perEquivD(int area){
        double pigreco= 3.1415926535;
        return (double)area/pigreco;
    }

    double compactness(int area, int feret){
        return arEquivD(area)/(double) feret;
    }

    double equivEllAr(int [] feret_and_breadth){
        double pigreco= 3.1415926535;
        return (pigreco*(double)feret_and_breadth[1]*(double)feret_and_breadth[0])/4;
    }

    double shape (int area, double perimeter){
        return (perimeter*perimeter)/(double) area;
    }

    double [] iS (int [] f, int [] b){
        double [] cS = {0,0};

        int xh0=f[0];
        int yh0=f[1];
        int xh1=f[2];
        int yh1=f[3];

        int xw0= b[0];
        int yw0=b[1];
        int xw1=b[2];
        int yw1=b[3];

        double s1= 0.5 *((double) ((xw0-xw1)*(yh0-yw1)-(yw0-yw1)*(xh0-xw1)));
        double s2= 0.5 * ((double) ((xw0-xw1)*(yw1-yh1)-(yw0-yw1)*(xw1-xh1)));

        cS[0]=(xh0+((s1*(xh1-xh0))/(s1+s2)));
        cS[1]=(yh0+((s1*(yh1-yh0))/(s1+s2)));

        rt.addValue("s1", s1);
        rt.addValue("s2", s2);

        return cS;
    }

    public double [] midpoint (ImageProcessor ip_for_midpoint, int xe, int ye){
        int internal_perimeter= 0;
        int external_perimeter= 0;

        int [] sum_internal_outline= {0,0};
        int [] sum_external_outline= {0,0};

        for(int x=0; x<xe; x++ ){
            for(int y=0; y<ye; y++){
                //se il pixel è bianco, fa parte dell'oggetto da analizzare dunque viene conteggiato come "pixel area"
                if(ip_for_midpoint.getPixel(x,y)==255){
                    //perimetro; solo i bordi interni verranno controllati e conteggiati
                    if(ip_for_midpoint.getPixel(x-1,y)==0){
                        internal_perimeter++;
                        sum_internal_outline[0]+=x;
                        sum_internal_outline[1]+=y;
                    } else if (ip_for_midpoint.getPixel(x+1,y)==0) {
                        internal_perimeter++;
                        sum_internal_outline[0]+=x;
                        sum_internal_outline[1]+=y;
                    }else if (ip_for_midpoint.getPixel(x,y-1)==0){
                        internal_perimeter++;
                        sum_internal_outline[0]+=x;
                        sum_internal_outline[1]+=y;
                    } else if (ip_for_midpoint.getPixel(x,y+1)==0){
                        internal_perimeter++;
                        sum_internal_outline[0]+=x;
                        sum_internal_outline[1]+=y;
                    }
                }else {
                    if(ip_for_midpoint.getPixel(x-1,y)==255){
                        sum_external_outline[0]+=x;
                        sum_external_outline[1]+=y;
                        external_perimeter++;
                    } else if (ip_for_midpoint.getPixel(x+1,y)==255) {
                        sum_external_outline[0]+=x;
                        sum_external_outline[1]+=y;
                        external_perimeter++;
                    }else if (ip_for_midpoint.getPixel(x,y-1)==255){
                        sum_external_outline[0]+=x;
                        sum_external_outline[1]+=y;
                        external_perimeter++;
                    } else if (ip_for_midpoint.getPixel(x,y+1)==255){
                        sum_external_outline[0]+=x;
                        sum_external_outline[1]+=y;
                        external_perimeter++;
                    }
                }
            }
        }

        double [] midpoint= {0,0};

        double internal_midpoint_X = (double)sum_internal_outline[0]/internal_perimeter;
        double internal_midpoint_Y = (double)sum_internal_outline[1]/internal_perimeter;

        double external_midpointX = (double)sum_external_outline[0]/external_perimeter;
        double external_midpointY = (double)sum_external_outline[1]/external_perimeter;

        midpoint[0] = (internal_midpoint_X + external_midpointX) /2;
        midpoint[1] = (internal_midpoint_Y + external_midpointY) /2;

        return midpoint;
    }

    double dS(double [] cS, double [] midpoint){
        double xCs= cS[0];
        double yCs= cS[1];
        double xCg= midpoint[0];
        double yCg= midpoint[1];
        
        return  Math.sqrt((Math.pow(xCs-xCg, 2))+(Math.pow(yCs-yCg, 2))); 
    }





}
