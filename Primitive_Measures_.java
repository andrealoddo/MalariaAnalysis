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
    boolean doMidpoint;
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
        gd.addCheckbox("Midpoint",false);

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

        //gd.addCheckbox("Show me feret and breadth points.\n Please check Feret and breadth", false);


        gd.showDialog();

        if (gd.wasCanceled())
            return DONE;

        boolean doSelectAll=gd.getNextBoolean ();

        if(doSelectAll){
            doArea = true;
            doPerimeter = true;
            doMidpoint = true;
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
        }else{
            doArea = gd.getNextBoolean ();
            doPerimeter = gd.getNextBoolean ();
            doMidpoint = gd.getNextBoolean ();

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
            rt.addValue("Midpoint X", midpoint[0]);
            rt.addValue("Midpoint Y", midpoint[1]);
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

        double s1= 1/2*(((double)b[0]-b[2])*((double)f[1]-b[3])-((double)b[1]-b[3])*((double)f[0]-b[2]));
        double s2= 1/2*(((double)b[0]-b[2])*((double)b[3]-f[3])-((double)b[1]-b[3])*((double)b[2]-f[2]));

        //double s1= 1/2*((b[2]-b[0])*(f[3]-b[1])-(b[3]-b[1])*(f[2]-b[0]));
        //double s2= 1/2*((b[2]-b[0])*(b[1]-f[1])-(b[3]-b[1])*(b[0]-f[0]));

        cS[0]=((f[2])+(s1*(f[0]-f[2])/s1+s2));
        cS[1]=((f[3])+(s1*(f[1]-f[3])/s1+s2));

        return cS;
    }



}
