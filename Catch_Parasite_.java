import ij.*;
import ij.measure.*;
import ij.plugin.*;
import ij.plugin.filter.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.filter.Analyzer;
import ij.plugin.MeasurementsWriter;
import ij.plugin.frame.ThresholdAdjuster;
import ij.plugin.ChannelSplitter;
import ij.plugin.Thresholder;
import java.awt.*;

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
    boolean doMINRMAXR;
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
    boolean doModRatio;
    boolean doSphericity;
    boolean doHaralickRatio;//non funziona
    boolean doElongation;
    boolean doNormPeriIndex;


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
        gd.addMessage("The measures area, perim, feret, feret min are necessary for the PlugIn");

        gd.addCheckbox("SelectAll",true);

        gd.addCheckbox("Convex Area", false);
        gd.addCheckbox("Convex Perimeter", false);
        gd.addCheckbox("MinR and MaxR", false);

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
        gd.addCheckbox("ModRatio", false);
        gd.addCheckbox("Sphericity", false);
        gd.addCheckbox("HaralickRatio", false);
        gd.addCheckbox("Elongation", false);
        gd.addCheckbox("NormPeriIndex", false);
        gd.showDialog(); //show

        if (gd.wasCanceled())
        return false;

        boolean doSelectAll=gd.getNextBoolean ();  //se viene selezionata la voce "seleziona tutti"

        if(doSelectAll){ // vengono settati i booleani a true
            doConvexArea = true;
            doConvexPerimeter=true;
            doMINRMAXR=true;
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
            doModRatio= true;
            doSphericity =true;
            doHaralickRatio=true;
            doElongation=true;
            doNormPeriIndex=true;
        }else{ //altrimenti pescati uno per uno con un certo ordine
            doConvexArea = gd.getNextBoolean();
            doConvexPerimeter = gd.getNextBoolean();
            doMINRMAXR=gd.getNextBoolean();
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
            doModRatio = gd.getNextBoolean();
            doSphericity= gd.getNextBoolean();
            doHaralickRatio = gd.getNextBoolean();
            doElongation = gd.getNextBoolean();
            doNormPeriIndex= gd.getNextBoolean();
        }

        return true;


    }

    //da sistemare
    public void setImage(){
        ThresholdAdjuster th = new ThresholdAdjuster();
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

    /*double [] getMinRMaxR(float xCenterOfMass, float yCenterOfMass){
        double [] minR_maxR;

        return minR_maxR;
    }*/



    //cuore
    /*? una classe interna dichiarata nella classe Parassite_Plugin*/

    class Parasite extends ParticleAnalyzer {

        int old_measures = Analyzer.getMeasurements();
        int necessary_measures =
                ALL_STATS +SHAPE_DESCRIPTORS
                -KURTOSIS
                -NaN_EMPTY_CELLS
                -INVERT_Y
                -INTEGRATED_DENSITY
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
            float[] areas = rt.getColumn(ResultsTable.AREA);
            float[] ferets = rt.getColumn(ResultsTable.FERET);
            float[] breadths =  rt.getColumn(ResultsTable.MIN_FERET);
            float[] perims = rt.getColumn(ResultsTable.PERIMETER);
            //double convexArea = roi!=null?getArea(roi.getConvexHull()):stats.pixelCount;
            double convexArea = getArea(roi.getConvexHull());
            double convexPerimeter = getPerimeter(roi.getConvexHull());
            //double peri = roi!=null?getPerimeter(roi.getConvexHull()):stats.pixelCount;

            for (int i = 0; i < areas.length; i++) {
               // double [] minR_maxR = getMinRMaxR(stats.xCenterOfMass, stats.yCenterOfMass);
                if(doConvexArea) //Area of the convex hull polygon
                    rt.addValue("*ConvexArea", convexArea);
                if(doConvexPerimeter) //Perimeter of the convex hull polygon
                    rt.addValue("*PerimeterConvexHull", convexPerimeter);
                //MINR E MAXR PROVA, inteso come raggio del cerchio più piccolo e più grande prendendo come diametro il FERET e il BREADTH
                if(doMINRMAXR){ //
                    rt.addValue("*MinR",breadths[i]/2 ); //DA RIVEDERE Radius of the inscribed circle centred at the middle of mass
                    rt.addValue("*MaxR",ferets[i]/2); //DA RIVEDERE Radius of the enclosing circle centred at the middle of mass
                }
                if(doAspRatio) //Aspect ratio = Feret/Breadth = L/W also called Feret ratio or Eccentricity or Rectangular ratio
                    rt.addValue("*AspRatio", ferets[i]/breadths[i]);
                if (doCirc) //Circularity = 4·?·Area/Perim2  also called Formfactor or Shapefactor
                    rt.addValue("*Circ", 4 * areas[i] / (perims[i] * perims[i]));
                if (doRoundness) //Roundness = 4·Area/(?·Feret2)
                    rt.addValue("*Roundness", (areas[i] * 4) / ((pigreco) * (ferets[i] * ferets[i])));
                if(doArEquivD) //Diameter of a circle with equivalent area,
                    rt.addValue("*ArEquivD", Math.sqrt((4 * pigreco) * areas[i]));
                if(doPerEquivD) //Diameter of a circle with equivalent perimeter,  Area/?
                    rt.addValue("*PerEquivD", areas[i]/pigreco);
                if(doEquivEllAr) //Area of the ellipse with Feret and Breath as major and minor axis,  = (?·Feret·Breadth)/4
                    rt.addValue("*EquivEllAr", (pigreco*ferets[i]*breadths[i])/4);
                if(doCompactness) //Compactness = ?((4/?)·Area)/Feret
                    rt.addValue("*Compactness", (Math.sqrt((4 * pigreco) * areas[i]))/ferets[i]);
                if(doConvexity) //Convexity = Convex_Perim/Perimeter also called rugosity or roughness
                    rt.addValue("*Convexity", convexPerimeter/perims[i]);
                if(doShape) //Shape = Perimeter2/Area also called Thinness ratio
                    rt.addValue("*Shape", (perims[i]*perims[i])/areas[i]);
                if(doRFactor) //RFactor = Convex_Area /(Feret·?)
                    rt.addValue("*RFactor", convexArea/(ferets[i]*pigreco));
                if(doArBBox) //Area of the bounding box along the Feret diameter = Feret·Breadth
                    rt.addValue("*ArBBox", ferets[i]*breadths[i]);
                if(doRectang) //Rectangularity = Area/ArBBox also called Extent
                    rt.addValue("*Rectang", areas[i]/(ferets[i]*breadths[i]));
                if(doModRatio) //Modification ratio = (2·MinR)/Feret
                    rt.addValue("*ModRatio", (breadths[i]/ferets[i])); //2 * MinR / Feret DA RIVEDERE dannorisultati uguali
                if(doSphericity) //Sphericity = MinR/MaxR also called Radius ratio
                    rt.addValue("*Sphericity", (breadths[i]/2)/(ferets[i]/2)); //MinR / MaxR DA RIVEDERE danno risultati ugualu
                if(doElongation) //The inverse of the circularity,  Perim2/(4·?·Area)
                    rt.addValue("*Elongation", (perims[i]*perims[i])/(4*pigreco*areas[i]));
                if(doNormPeriIndex)
                    rt.addValue("*normPeriIndex", (2*Math.sqrt(pigreco*areas[i]))/perims[i]);
                /*Quinn 2014 da Review di Rosado
                normPeriIndex = (2*sqrt(pi*Area))/Perim
                The normalized perimeter index compares the input perimeter to the most compact polygon
                with the same area (equal area circle), meaning you can use it to identify features with irregular
                boundaries.
                The normalized perimeter index uses the equal area circle to normalize the metric.
                */

            }
            /*float[] meanRadius= rt.getColumn(ResultsTable.MEAN);
            float[] standardDeviationOfRadii= rt.getColumn(ResultsTable.STD_DEV);
            if(doHaralickRatio) DA FARE SU PIù OGGE
                // rt.addValue("*Haralick ratio", meanRadius[i]/standardDeviationOfRadii[i]); //da rivedere
                */

            Analyzer.setMeasurements(old_measures); //per risettare misure vecchie
        }
    }
}
