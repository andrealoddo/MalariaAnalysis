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
    boolean doFeret_and_breadth; //feret and breadth
    boolean doAspRatio;
    boolean doCirc;
    boolean doRoundness;
    boolean doArEquivD;
    boolean doPerEquivD;
    boolean doEquivEllAr;
    boolean doCompactness;
    boolean doShape;


    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;

        GenericDialog gd = new GenericDialog("Primitive_Measures", IJ.getInstance());
        gd.addStringField("Title: ", title);
        gd.addMessage("Choise measures");

        gd.addCheckbox("Area",true);
        gd.addCheckbox("Perimeter",false);
        gd.addCheckbox("Midpoint",false);
        gd.addCheckbox("Feret and breadth of element",false);
        gd.addCheckbox("AspRatio",false);
        gd.addCheckbox("Circ", false);
        gd.addCheckbox("Roundness", false);
        gd.addCheckbox("ArEquivD", false);
        gd.addCheckbox("PerEquivD", false);
        gd.addCheckbox("EquivEllAr", false);
        gd.addCheckbox("Compactness", false);
        gd.addCheckbox("Shape", false);


        gd.showDialog();

        if (gd.wasCanceled())
            return DONE;

        doArea = gd.getNextBoolean ();
        doPerimeter = gd.getNextBoolean ();
        doMidpoint = gd.getNextBoolean ();
        doFeret_and_breadth = gd.getNextBoolean ();
        doAspRatio=gd.getNextBoolean();
        doCirc=gd.getNextBoolean();
        doRoundness=gd.getNextBoolean();
        doArEquivD=gd.getNextBoolean();
        doPerEquivD=gd.getNextBoolean();
        doEquivEllAr=gd.getNextBoolean();
        doCompactness=gd.getNextBoolean();


        return DOES_ALL;

    }

    public void run(ImageProcessor ip) {

        rt.reset();

        int xe = ip.getWidth(); //larghezza e altezza dell'intera immagine
        int ye = ip.getHeight();

        rt.incrementCounter();

        if(doArea){
            int area = area(ip, xe, ye);
            rt.addValue("Area of element", area);
        }

        if(doPerimeter){
            double perimeter =perimeter(ip, xe, ye);
            rt.addValue("Perimeter of element", perimeter);
        }

        if(doMidpoint){
            double [] midpoint = midpoint(ip, xe, ye);
            rt.addValue("Midpoint X", midpoint[0]);
            rt.addValue("Midpoint Y", midpoint[1]);
        }

        if(doFeret_and_breadth){
            int [] feret_and_breadth = feret_and_breadth(ip, xe, ye);
            rt.addValue("Width/Breadth of Element", feret_and_breadth[0]);
            rt.addValue("Height/Feret of Element", feret_and_breadth[1]); //sempre alla posizione 1 feret
        }

        if(doAspRatio){
            double aspRatio = aspRatio(feret_and_breadth(ip, xe, ye));
            rt.addValue("AspRatio", aspRatio);
        }

        if(doCirc){
            double circ = circ((area(ip, xe, ye)), perimeter(ip, xe, ye));
            rt.addValue("Circ", circ);

        }

        if(doRoundness){
            int [] feret_and_breadth = feret_and_breadth(ip, xe, ye);
            double roundness = roundness(area(ip, xe, ye), (feret_and_breadth[1]));
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
            double equivEllAr = equivEllAr(feret_and_breadth(ip, xe, ye));
            rt.addValue("EquivEllAr", equivEllAr);
        }

        if(doCompactness){
            int [] feret_and_breadth = feret_and_breadth(ip, xe, ye);
            double compactness = compactness(area(ip, xe, ye), (feret_and_breadth[1]));
            rt.addValue("Compactness", compactness);
        }

        if(doShape){
            double shape = shape (area(ip, xe, ye), perimeter(ip, xe, ye));
            rt.addValue("Shape", shape);
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

    public int [] feret_and_breadth (ImageProcessor ip_for_feret_and_breadth, int xe, int ye){
       int [] supportUp = {0,0}; //per coordinate superiori massime dell'elemento
       int [] supportDown = {0,0}; //per coordinate minime dell'elemento

       supportDown[1]=ye; supportDown[0]=xe;


        for(int x=0; x<xe; x++ ){
            for(int y=0; y<ye; y++) {
                if (ip_for_feret_and_breadth.getPixel(x, y) == 255) {
                    //perimetro; solo i bordi interni verranno controllati e conteggiati
                    if (ip_for_feret_and_breadth.getPixel(x - 1, y) == 0) {
                        supportDown=updateSupportDown(x, y, supportDown); //richiamo funzione di aggiornamento coordinate per calcolo dell'altezza e della larghezza dell'elemento
                        supportUp=updateSupportUp(x, y, supportUp ); //solo nei bordi interni per evitare caos

                    } else if (ip_for_feret_and_breadth.getPixel(x + 1, y) == 0) {
                        supportDown=updateSupportDown(x, y, supportDown); //richiamo funzione di aggiornamento coordinate per calcolo dell'altezza e della larghezza dell'elemento
                        supportUp=updateSupportUp(x, y, supportUp );

                    } else if (ip_for_feret_and_breadth.getPixel(x, y - 1) == 0) {
                        supportDown=updateSupportDown(x, y, supportDown); //richiamo funzione di aggiornamento coordinate per calcolo dell'altezza e della larghezza dell'elemento
                        supportUp=updateSupportUp(x, y, supportUp );
                    } else if (ip_for_feret_and_breadth.getPixel(x, y + 1) == 0) {
                        supportDown=updateSupportDown(x, y, supportDown); //richiamo funzione di aggiornamento coordinate per calcolo dell'altezza e della larghezza dell'elemento
                        supportUp=updateSupportUp(x, y, supportUp );
                    }

                }
            }
        }

        int heightElement= (supportUp[1]-supportDown[1])+1 ; //calcolo altezza e larghezza dell'elemento: il +1 è necessario perchè il primo pixel viene settato a 0 (0,0)
        int widthElement = (supportUp[0]-supportDown[0])+1;

        if(itIsFeret(heightElement, widthElement)){
            int [] feret_and_breadth = {widthElement, heightElement };
            return feret_and_breadth; //0 larghezza 1 altezza
        }else{
            int [] feret_and_breadth = {heightElement, widthElement};
            return feret_and_breadth; //0 larghezza 1 altezza
        }

    }

    int[] updateSupportDown(int x, int y, int [] supportDown){
        if(supportDown[1]>y){supportDown[1]=y; }
        if(supportDown[0]>x){supportDown[0]=x; }

        return supportDown;
    }

    int[] updateSupportUp(int x, int y, int [] supportUp){
        if(supportUp[1]<y){supportUp[1]=y;}
        if(supportUp[0]<x){supportUp[0]=x;}

        return supportUp;
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
}
