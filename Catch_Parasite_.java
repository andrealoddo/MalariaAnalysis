import ij.*;
import ij.measure.*;
import ij.plugin.*;
import ij.plugin.filter.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.filter.Analyzer;
import ij.plugin.MeasurementsWriter;


/*PlugIn di prova per poter aggiungere funzionalit? ad Partiche Analyzer.
 * la classe Parassite_Plugin implementa l'interfaccia PlugIn*/
/*la classe estende la classe madre Analyzer Particles che a sua volta si appoggia alla classe Analyzer.*/

/*da settare in SET MEASUREMENTS:
-area
-Centroid
-CenterOfMass
-Perimeter
-Fit Ellipse
-Feret's diameter
-Skewness
-ShapeDescriptor
*/


public class Catch_Parasite_ implements PlugIn {
    /*booleani per la checkbox in cui l'utente ? invitato a selezionare cosa vuole calcolare .. da rivedere*/

    //ConvexHull data
    boolean doConvexArea;
    boolean doConvexPerimeter;

    boolean doAspRatio;
    boolean doCirc;
    boolean doRoundness;
    boolean doArEquivD;
    boolean doPerEquivD;
    boolean doEquivEllAr;
    boolean doCompactness;
    boolean doConvexity;
    boolean doShape;
    boolean doRFactor;
    boolean doArBBox;
    boolean doRectang;
    boolean doHaralickRatio;//non funziona
    boolean doElongation;
    boolean doIS;
    boolean doMidpoint;
    boolean doDS;

    //pigreco necessario per alcuni calcoli

    double pigreco= 3.1415926535;

    public void run(String arg) {
        //Restituisce un riferimento all'immagine attiva o visualizza un messaggio
        // di errore e interrompe il plugin o la macro se non ci sono immagini aperte
        ImagePlus imp = IJ.getImage();
        analyzeParticles(imp); //metodo presente dopo
    }

    public void analyzeParticles(ImagePlus imp) {
        //crea un nuovo oggetto Parassite
        //crea flags, che prende setup
        //genericDialog();
        Parasite pa = new Parasite();
        int flags = pa.setup("", imp); //runna
        if (flags == PlugInFilter.DONE)
          return; //se ? tutto ok runna
        pa.run(imp.getProcessor());
        Analyzer.getResultsTable().show("Results");

    }

    /*Finestra di dialogo:
    la funzione è volta ad aprire una finestra di dialogo per cui l'utente setta tutte le misure che vuole implementare.
    * */
    public boolean genericDialog(){
        GenericDialog gd = new GenericDialog("Parassite Prova", IJ.getInstance());
        gd.addStringField("Title: ", "Catch parasite");
        gd.addMessage("Choise measures");

        gd.addCheckbox("SelectAll",true);

        gd.addCheckbox("Convex Area", false);
        gd.addCheckbox("Convex Perimeter", false);

        gd.addCheckbox("AspRatio",false);
        gd.addCheckbox("Circ", false);
        gd.addCheckbox("Roundness", false);
        gd.addCheckbox("ArEquivD", false);
        gd.addCheckbox("PerEquivD", false);
        gd.addCheckbox("EquivEllAr", false);
        gd.addCheckbox("Compactness", false);
        gd.addCheckbox("Convexity", false);
        gd.addCheckbox("Shape", false);

        gd.addCheckbox("RFactor", false);

        gd.addCheckbox("ArBBox", false);
        gd.addCheckbox("Rectang", false);
        gd.addCheckbox("HaralickRatio", false);
        gd.addCheckbox("Elongation", false);
        gd.showDialog(); //show

        if (gd.wasCanceled())
        return false;

        boolean doSelectAll=gd.getNextBoolean ();  //se viene selezionata la voce "seleziona tutti"

        if(doSelectAll){ // vengono settati i booleani a true
            doConvexArea = true;
            doConvexPerimeter=true;
            doAspRatio = true;
            doCirc = true;
            doRoundness = true;
            doArEquivD = true;
            doPerEquivD = true;
            doEquivEllAr = true ;
            doCompactness = true;
            doConvexity =true;
            doShape = true;
            doRFactor = true;
            doArBBox =true;
            doRectang = true;
            doHaralickRatio=true;
            doElongation=true;

        }else{ //altrimenti pescati uno per uno con un certo ordine
            doConvexArea = gd.getNextBoolean();
            doConvexPerimeter = gd.getNextBoolean();
            doAspRatio=gd.getNextBoolean();
            doCirc=gd.getNextBoolean();
            doRoundness=gd.getNextBoolean();
            doArEquivD=gd.getNextBoolean();
            doPerEquivD=gd.getNextBoolean();
            doEquivEllAr=gd.getNextBoolean();
            doCompactness=gd.getNextBoolean();
            doConvexity=gd.getNextBoolean();
            doShape=gd.getNextBoolean();
            doRFactor=gd.getNextBoolean();
            doArBBox =gd.getNextBoolean();
            doRectang = gd.getNextBoolean();
            doHaralickRatio = gd.getNextBoolean();
            doElongation = gd.getNextBoolean();

        }

        return true;


    }


    //da Analyzer ma rivisitato dopo il confronto con i dati calcolati da MatLab
    final double getArea(Polygon p) {
        if (p==null) return Double.NaN;
        int carea = 0;
        int iminus1;
        /*FUNZIONAMENTO
        * tratta il polygon in maniera tale da scomporlo in più triangoli.*/
        /*ciclo originale
        for (int i=0; i<p.npoints; i++) {
            iminus1 = i-1;
            if (iminus1<0) iminus1=p.npoints-1;
            carea += (p.xpoints[i]+p.xpoints[iminus1])*(p.ypoints[i]-p.ypoints[iminus1]);
        }*/

        for(int i =0; i<p.npoints-1;i++){
            iminus1=i-1;
            if(iminus1<0) iminus1=p.npoints-1;
            carea+=(p.xpoints[i]+p.xpoints[iminus1])*(p.ypoints[i]-p.ypoints[iminus1]);
        }

        //lo triangolizza?
        //return (double) p.npoints;
        return (Math.abs(carea/2.0));
    }

    final double getPerimeter (Polygon p){
        if(p==null) return  Double.NaN;
        /*FUNZIONAMENTO
         * somma le distanze da un punto A ad un punto B
         * */

        double cperimeter = 0.0;
        int iminus1;

        for(int i=0; i < p.npoints-1; i++){
            iminus1=i-1;
            if(iminus1<0) iminus1 = p.npoints-1;

            double x1 =(double) p.xpoints[i];
            double y1 =(double) p.ypoints[i];

            double x2 =(double) p.xpoints[iminus1];
            double y2 =(double) p.ypoints[iminus1];

            cperimeter +=Math.sqrt(((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1)));
        }

        return cperimeter;
    }

    //cuore
    /*? una classe interna dichiarata nella classe Parassite_Plugin*/

    class Parasite extends ParticleAnalyzer {

        int old_measures = Analyzer.getMeasurements();
        /*int necessary_measures =
                ADD_TO_OVERLAY+
                AREA+
                AREA_FRACTION+
                CENTER_OF_MASS+
                CENTROID+
                CIRCULARITY+
                FERET+
                MEAN+
                MEDIAN+
                MIN_MAX+
                MODE+
                PERIMETER+
                SHAPE_DESCRIPTORS+
                STD_DEV+
                STACK_POSITION;*/
        int necessary_measures =
                ALL_STATS
                -KURTOSIS
                -NaN_EMPTY_CELLS
                -INVERT_Y
                -INTEGRATED_DENSITY
                +SHAPE_DESCRIPTORS
                ;


        //*protected static final int NOTHING=0, OUTLINES=1, BARE_OUTLINES=2, ELLIPSES=3, MASKS=4, ROI_MASKS=5,
        //            OVERLAY_OUTLINES=6, OVERLAY_MASKS=7;*//
        @Override
        public boolean showDialog(){
            boolean flag = super.showDialog();
            //setto solo quello che mi serve
            Analyzer.setMeasurements(necessary_measures);
            //Analyzer.setMeasurements(AREA);
            return genericDialog();
        }


        //SAVE RESULT
        @Override
        protected void saveResults(ImageStatistics stats, Roi roi) {
            super.saveResults(stats, roi);
            ImageProcessor ip = imp.getProcessor();
            ip.setRoi(roi);
            float[] meanRadius= rt.getColumn(ResultsTable.MEAN);
            float[] standardDeviationOfRadii= rt.getColumn(ResultsTable.STD_DEV);
            float[] areas = rt.getColumn(ResultsTable.AREA);
            float[] ferets = rt.getColumn(ResultsTable.FERET);
            float [] breadths =  rt.getColumn(ResultsTable.MIN_FERET);
            float[] perims = rt.getColumn(ResultsTable.PERIMETER);
            //double convexArea = roi!=null?getArea(roi.getConvexHull()):stats.pixelCount;
            double convexArea = getArea(roi.getConvexHull());
            double convexPerimeter = getPerimeter(roi.getConvexHull());
            //double peri = roi!=null?getPerimeter(roi.getConvexHull()):stats.pixelCount;


            for (int i = 0; i < areas.length; i++) {
                if(doConvexArea)
                    rt.addValue("*ConvexArea", convexArea);
                if(doConvexPerimeter)
                    rt.addValue("*PerimeterConvexHull", convexPerimeter);
                if(doAspRatio)
                    rt.addValue("*AspRatio", ferets[i]/breadths[i]);
                if (doCirc)
                    rt.addValue("*Circ", 4 * areas[i] / (perims[i] * perims[i]));
                if (doRoundness)
                    rt.addValue("*Roundness", (areas[i] * 4) / ((pigreco) * (ferets[i] * ferets[i])));
                if(doArEquivD)
                    rt.addValue("*ArEquivD", Math.sqrt((4 * pigreco) * areas[i]));
                if(doPerEquivD)
                    rt.addValue("*PerEquivD", areas[i]/pigreco);
                if(doEquivEllAr)
                    rt.addValue("*EquivEllAr", (pigreco*ferets[i]*breadths[i])/4);
                if(doCompactness)
                    rt.addValue("*Compactness", (Math.sqrt((4 * pigreco) * areas[i]))/ferets[i]);
                if(doConvexity)
                    rt.addValue("*Convexity", convexPerimeter/perims[i]);
                if(doShape)
                    rt.addValue("*Shape", (perims[i]*perims[i])/areas[i]);
                if(doRFactor)
                    rt.addValue("*RFactor", convexArea/(ferets[i]*pigreco));
                if(doArBBox)
                    rt.addValue("*ArBBox", ferets[i]*breadths[i]);
                if(doRectang)
                    rt.addValue("*Rectang", areas[i]/(ferets[i]*breadths[i]));
                if(doHaralickRatio)
                   // rt.addValue("*Haralick ratio", meanRadius[i]/standardDeviationOfRadii[i]);
                if(doElongation)
                    rt.addValue("*Elongation", (perims[i]*perims[i])/(4*pigreco*areas[i]));

            }

            Analyzer.setMeasurements(old_measures);
        }
    }
}
