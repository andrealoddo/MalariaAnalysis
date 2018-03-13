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
    boolean doHeight_and_width;


    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;

        GenericDialog gd = new GenericDialog("Primitive_Measures", IJ.getInstance());
        gd.addStringField("Title: ", title);
        gd.addMessage("Choise measures");

        gd.addCheckbox("Area",true);
        gd.addCheckbox("Perimeter",false);
        gd.addCheckbox("Midpoint",false);
        gd.addCheckbox("Height and width of element",false);

        gd.showDialog();

        if (gd.wasCanceled())
            return DONE;

        doArea = gd.getNextBoolean ();
        doPerimeter = gd.getNextBoolean ();
        doMidpoint = gd.getNextBoolean ();
        doHeight_and_width = gd.getNextBoolean ();


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

        if(doHeight_and_width){
            int [] height_and_width = height_and_width(ip, xe, ye);
            rt.addValue("Width of Element", height_and_width[0]);
            rt.addValue("Height of Element", height_and_width[1]);
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

    public int [] height_and_width (ImageProcessor ip_for_height_and_width, int xe, int ye){
       int [] supportUp = {0,0}; //per coordinate superiori massime dell'elemento
       int [] supportDown = {0,0}; //per coordinate minime dell'elemento

       supportDown[1]=ye; supportDown[0]=xe;


        for(int x=0; x<xe; x++ ){
            for(int y=0; y<ye; y++) {
                if (ip_for_height_and_width.getPixel(x, y) == 255) {
                    //perimetro; solo i bordi interni verranno controllati e conteggiati
                    if (ip_for_height_and_width.getPixel(x - 1, y) == 0) {
                        supportDown=updateSupportDown(x, y, supportDown); //richiamo funzione di aggiornamento coordinate per calcolo dell'altezza e della larghezza dell'elemento
                        supportUp=updateSupportUp(x, y, supportUp ); //solo nei bordi interni per evitare caos

                    } else if (ip_for_height_and_width.getPixel(x + 1, y) == 0) {
                        supportDown=updateSupportDown(x, y, supportDown); //richiamo funzione di aggiornamento coordinate per calcolo dell'altezza e della larghezza dell'elemento
                        supportUp=updateSupportUp(x, y, supportUp );

                    } else if (ip_for_height_and_width.getPixel(x, y - 1) == 0) {
                        supportDown=updateSupportDown(x, y, supportDown); //richiamo funzione di aggiornamento coordinate per calcolo dell'altezza e della larghezza dell'elemento
                        supportUp=updateSupportUp(x, y, supportUp );
                    } else if (ip_for_height_and_width.getPixel(x, y + 1) == 0) {
                        supportDown=updateSupportDown(x, y, supportDown); //richiamo funzione di aggiornamento coordinate per calcolo dell'altezza e della larghezza dell'elemento
                        supportUp=updateSupportUp(x, y, supportUp );
                    }

                }
            }
        }

        int heightElement= (supportUp[1]-supportDown[1])+1 ; //calcolo altezza e larghezza dell'elemento: il +1 è necessario perchè il primo pixel viene settato a 0 (0,0)
        int widthElement = (supportUp[0]-supportDown[0])+1;
        int [] height_and_width = {widthElement, heightElement};

        return height_and_width; //0 larghezza 1 altezza

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
}
