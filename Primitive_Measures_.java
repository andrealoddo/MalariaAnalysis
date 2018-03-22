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

/*Plugin-di-prova*/
/*Il plugin è suddiviso in tre parti principali:
* una parte di setup e di inizializzazione
* una seconda parte di run in cui vengono richiamate le funzioni
* della terza parte per cui sono volte alle measures*/

public class Primitive_Measures_ implements PlugInFilter {
    ImagePlus imp;//elemento fondamentale per il PluginFilter
    static String title = "prova";

    protected ResultsTable rt = ResultsTable.getResultsTable();

    /*booleani per la checkbox in cui l'utente è invitato a selezionare cosa vuole calcolare .. da rivedere*/
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
 //inizializzazione dei punti estremi del feret e del breadth
    int [] feret_exstreme_points = {0,0,0,0};
    int [] breadth_exstreme_points = {0,0,0,0};

    int boundaryPixelsX =0;
    int boundaryPixelsY =0;

    double pigreco= 3.1415926535;


    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
//creazione della finestra di dialogo per checkbox
        GenericDialog gd = new GenericDialog("Primitive_Measures", IJ.getInstance());
        gd.addStringField("Title: ", title);
        gd.addMessage("Choise measures");

        //valori da immettere nella checkbox

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


        gd.showDialog(); //show

        if (gd.wasCanceled())
            return DONE;

        boolean doSelectAll=gd.getNextBoolean ();  //se viene selezionata la voce "seleziona tutti"

        if(doSelectAll){ // vengono settati i booleani a true
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
        }else{ //altrimenti pescati uno per uno con un certo ordine
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
/*seconda parte: run. Qui vengono richiamate le funzioni necessarie grazie
* al settaggio dei booleani. è necessario un ImagePrccessor.*/
    public void run(ImageProcessor ip) {
        rt.reset(); //reset della finestra per pulirla

        rt.incrementCounter(); //incremeneto riga (in futurò potra servire

        int xe = ip.getWidth(); //larghezza e altezza dell'intera immagine
        int ye = ip.getHeight();

        feret_exstreme_points[3]=ye; //inizializzazione dei punti minori
        breadth_exstreme_points[2]=xe; //inizializzazione dei punti minori

        int [] feret_and_breadth = {breadth(ip, xe, ye), feret(ip, xe, ye)}; //il calcolo del feret e del breath è necessario

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
        
        /*il +1 è necessario perchè il primo pixel dell'immagine viene conteggiato ad 0 dunque se io taglio l'immagine 
        * scopro che vale uno in più rispetto al valore dato senza questo increment
        * perchè la valutazione viene fatta sul conteggio dei pixel e valutando le coordinate*/

        if(doFeret){
            rt.addValue("Height/Feret of Element", feret_and_breadth[1]+1);
        }

        if(doBreadth){
            rt.addValue("Width/Breadth of Element", feret_and_breadth[0]+1);
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
            ip.putPixelValue((int)cS[0],(int)cS[1], 130); //colorazione per capire quale punto è stato visualizzato
        }

        if(doDS){
            double dS = dS (iS(feret_exstreme_points, breadth_exstreme_points), midpoint(ip, xe, ye));
            rt.addValue("DS", dS); //distanza tra cS e cG
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

        /*viene calcolato sia il perimetro interno che quello esterno e applicata una media*/

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


        return ((double)external_perimeter+(double)internal_perimeter)/2; //media tra i due perimetri
    }


    public int feret (ImageProcessor ip_for_feret, int xe, int ye){

        /*viene calcolata con l'utilizzo del confronto della y con un'altra funzione esterna.
        * viene effettuato un controllo a fine ciclo per verificare se è stato trovato il feret o il breadth*/

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

        /*se il feret è minore del breadth allora è stato calcolato il breadth dunque si invertono le coordinate*/

        if(itIsFeret(((feret_exstreme_points[1]-feret_exstreme_points[3])), ((breadth_exstreme_points[0]-breadth_exstreme_points[2]))))
        {
            return ((feret_exstreme_points[1]-feret_exstreme_points[3]));
        }else{
            int [] tt = feret_exstreme_points;
            feret_exstreme_points=breadth_exstreme_points;
            breadth_exstreme_points=tt;
            return (feret_exstreme_points[1]-feret_exstreme_points[3]); //ricalcola
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


            if (itIsFeret(((feret_exstreme_points[1] - feret_exstreme_points[3])), ((breadth_exstreme_points[0] - breadth_exstreme_points[2])))) {
                int[] tt = feret_exstreme_points;
                feret_exstreme_points = breadth_exstreme_points;
                breadth_exstreme_points = tt;
                return (breadth_exstreme_points[0] - breadth_exstreme_points[2]);//ricalcola
            } else {
                return (breadth_exstreme_points[0] - breadth_exstreme_points[2]);
            }

    }


    int[] feret_exstreme_points(int x, int y, int [] feret_exstreme_points){
        /*supporto del calcolo delle coordinate del feret*/
        if(feret_exstreme_points[1]<=y){
            feret_exstreme_points[1]=y;
            feret_exstreme_points[0]=x;
        }

        if(feret_exstreme_points[3]>=y){
            feret_exstreme_points[3]=y;
            feret_exstreme_points[2]=x;
        }

        return feret_exstreme_points;
    }

    int[] breadth_exstreme_points(int x, int y, int [] breadth_exstreme_points){
        /*supporto del calcolo delle coordinate del breadth*/
        if(breadth_exstreme_points[0]<=x){
            breadth_exstreme_points[1]=y;
            breadth_exstreme_points[0]=x;
        }
        if(breadth_exstreme_points[2]>=x){
            breadth_exstreme_points[3]=y;
            breadth_exstreme_points[2]=x;
        }
        return breadth_exstreme_points;
    }

    boolean itIsFeret (int feret, int breadth){
        //funziona booleana di supporto per il feret
        return feret>=breadth;
    }

    double aspRatio (int [] feret_and_breadth){
        //formula presa dal foglip
        return feret_and_breadth[1]/feret_and_breadth[0];
    }


    double roundness (int area, int feret){
     return (4 * (double)area )/ pigreco * ((double) feret * feret);

    }

    double circ (int area, double perim){
        return (4*(double)area)/(perim*perim);
    }

    double arEquivD(int area){
        return Math.sqrt((4*pigreco)*(double)area);

    }

    double perEquivD(int area){
        return (double)area/pigreco;
    }

    double compactness(int area, int feret){
        return arEquivD(area)/(double) feret;
    }

    double equivEllAr(int [] feret_and_breadth){
        return (pigreco*(double)feret_and_breadth[1]*(double)feret_and_breadth[0])/4;
    }

    double shape (int area, double perimeter){
        return (perimeter*perimeter)/(double) area;
    }

    double [] iS (int [] f, int [] b){
        /*calcolo della intersezione tra lenght e width con le coordinate del feret e breadth*/
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

        /*calcolo del centro della massa
        * valutazione solo nei pixel bianchi, pixel dell'oggetto analizzato
        * ricalcolo del pixel interno e del pixel esterno
        * si conteggiano i valori delle y e delle x interni ed esterni con i due vettori contatori sum_internal_outline e sum_external_outline
         * dove 0 è la x e 1 è la y*/

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

        boundaryPixelsX=sum_internal_outline[0];
        boundaryPixelsY=sum_internal_outline[1];

        //si fa la media

        midpoint[0] = (internal_midpoint_X + external_midpointX) /2;
        midpoint[1] = (internal_midpoint_Y + external_midpointY) /2;

        return midpoint;
    }

    double dS(double [] cS, double [] midpoint){
        /*distanza tra i pixel tra il punto di intersezione e il punto del centro*/
        double xCs= cS[0];
        double yCs= cS[1];
        double xCg= midpoint[0];
        double yCg= midpoint[1];

        return  Math.sqrt((Math.pow(xCs-xCg, 2))+(Math.pow(yCs-yCg, 2)));
    }





}
