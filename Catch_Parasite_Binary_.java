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
import java.util.ArrayList;
import java.util.Arrays;

/*Catch_Parasite by Tatalessap
 * Catch_Parasite è una classe volta alla misurazione di elementi all'interno di una immagine binaria.
 * Si appoggia al Plugin "ThresholdAdjuster" per poter individuare le zone e gli oggetti da analizzare.
 * Implementa Plugin ma al suo interno vi è una classe Parasite che estende a sua volta il PluginFilter "ParticleAnalyzer"
 * aggiungendo misure non presenti.*/

public class Catch_Parasite_Binary_ implements PlugIn {
    /*Variabili booleani corrispondenti alle misure implementate e aggiunte.*/
    /*Sezione 1*/
    /*ConvexHull*/
    boolean doConvexArea;
    boolean doConvexPerimeter;

    boolean doMINRMAXR;
    boolean doAspRatio;
    boolean doCircularity;
    boolean doRoundness;
    boolean doArEquivD;
    boolean doPerEquivD;
    boolean doEquivEllAr;
    boolean doConcavity;
    boolean doCompactness;
    boolean doConvexity;
    boolean doShape;
    boolean doRFactor;
    boolean doArBBox;
    boolean doRectang;
    boolean doModRatio;
    boolean doSphericity;
    boolean doElongation;
    /*Sezione 2*/
    boolean doNormPeriIndex;
    boolean doHeralickRatio;

    /*piGreco necessario per alcuni calcoli*/
    double pigreco= 3.1415926535;

    /*Metodo run necessaria per i PlugIn
     * Memorizza l'immagine in ingresso
     * richiama il metodo catch_parasite_running dandogli in ingresso l'immagine*/
    public void run(String arg) {
        ImagePlus imp = IJ.getImage();
        catch_parasite_running(imp); //metodo presente dopo
    }

    /*Metodo catch_parasite_running che riceve in ingresso l'immagine
     * creazione di un nuovo oggetto Parasite, implementato come classa interna
     * avvio attraverso setup e controllo
     * avvio attraverso run dell'oggetto Parasite
     * visualizzazione risultati da super - super classe Analyzer*/
    public void catch_parasite_running (ImagePlus imp) {
        Parasite parasite = new Parasite();
        int flags = parasite.setup("", imp);
        /*Controllo*/
        if (flags == PlugInFilter.DONE)
            return;
        parasite.run(imp.getProcessor());
        /*visualizzazione dei risultati*/
        Analyzer.getResultsTable().show("Results");
    }


    /*Inner class Parasite che estende la classe ParticleAnalyzer
     * Al suo interno, attraverso gli Override, si estendono i meotdi già presenti nella classe ParticleAnalyzer*/
    class Parasite extends ParticleAnalyzer {
        //Salvataggio delle misure già settate in SetMeasurements
        int old_measures = Analyzer.getMeasurements();
        //Salvataggio delle nuove misure necessarie per il calcolo delle nuove
        int necessary_measures =
                ALL_STATS+SHAPE_DESCRIPTORS
                        -SKEWNESS
                        -STD_DEV
                        -KURTOSIS
                        -NaN_EMPTY_CELLS
                        -INVERT_Y
                        -INTEGRATED_DENSITY
                ;
        /*showDialog:
        * settaggio della opzione per visualizzare contorni e etichetta numerata
        * settaggio nelle misure necessarie
        * richiamo della genericDialog creata per poter settare le misure da aggiungere*/
        @Override
        public boolean showDialog(){
            boolean flag = super.showDialog();
            super.staticShowChoice = 1; //*protected static final int NOTHING=0, OUTLINES=1, BARE_OUTLINES=2, ELLIPSES=3, MASKS=4, ROI_MASKS=5,
            //                             OVERLAY_OUTLINES=6, OVERLAY_MASKS=7;*//
            //setto solo quello che mi serve
            Analyzer.setMeasurements(necessary_measures);
            if(flag) return genericDialog();
            return false;

        }

        /*saveResults:
        * richiamo della funzione originaria
        * richiamo della funzione aggiuntiva addMeasure a cui vengono passati roi e stats.
        * roi: regione di interesse
        * stats: oggetto di tipo ImageStatistics
        * settaggio delle vecchie misure dopo che si è concluso il calcolo*/
        @Override
        protected void saveResults(ImageStatistics stats, Roi roi) {
            //ThresholdAdjuster th = new ThresholdAdjuster();
            super.saveResults(stats, roi);
            addMeasure(stats, roi, super.rt, super.imp);
            Analyzer.setMeasurements(old_measures); //per risettare misure vecchie
        }


        /*Metodo della creazione della finestra di dialogo, visualizzazione del:
         * - checkbox per selezionare ogni misura
         * - checkbox per selezionare misure singole*/
        public boolean genericDialog() {
            GenericDialog gd = new GenericDialog("Parassite Prova", IJ.getInstance());
            gd.addStringField("Title: ", "Catch parasite");
            gd.addMessage("Choise measures");
            gd.addMessage("The measures area, perim, feret, feret min are necessary for the PlugIn");
            gd.addCheckbox("SelectAll", true); //tutte le misure
            gd.addCheckbox("Convex Area", false);
            gd.addCheckbox("Convex Perimeter", false);
            gd.addCheckbox("MinR and MaxR", false);
            gd.addCheckbox("AspRatio", false);
            gd.addCheckbox("Circ", false);
            gd.addCheckbox("Roundness", false);
            gd.addCheckbox("ArEquivD", false);
            gd.addCheckbox("PerEquivD", false);
            gd.addCheckbox("EquivEllAr", false);
            gd.addCheckbox("Compactness", false);
            gd.addCheckbox("Concavity", false);
            gd.addCheckbox("Convexity", false);
            gd.addCheckbox("Shape", false);
            gd.addCheckbox("RFactor", false);
            gd.addCheckbox("ArBBox", false);
            gd.addCheckbox("Rectang", false);
            gd.addCheckbox("ModRatio", false);
            gd.addCheckbox("Sphericity", false);
            gd.addCheckbox("Elongation", false);
            gd.addCheckbox("NormPeriIndex", false);
            gd.addCheckbox("HeralickRatio", false);
            gd.showDialog(); //show
            if (gd.wasCanceled())
                return false;
            boolean doSelectAll = gd.getNextBoolean();  //se viene selezionata la voce "seleziona tutti"
            if (doSelectAll) { // vengono settati i booleani a true
                doConvexArea = true;
                doConvexPerimeter = true;
                doMINRMAXR = true;
                doAspRatio = true;
                doCircularity = true;
                doRoundness = true;
                doArEquivD = true;
                doPerEquivD = true;
                doEquivEllAr = true;
                doConcavity = true;
                doCompactness = true;
                doConvexity = true;
                doShape = true;
                doRFactor = true;
                doArBBox = true;
                doRectang = true;
                doModRatio = true;
                doSphericity = true;
                doElongation = true;
                doNormPeriIndex = true;
                doHeralickRatio = true;
            } else { //altrimenti pescati uno per uno con un certo ordine
                doConvexArea = gd.getNextBoolean();
                doConvexPerimeter = gd.getNextBoolean();
                doMINRMAXR = gd.getNextBoolean();
                doAspRatio = gd.getNextBoolean();
                doCircularity = gd.getNextBoolean();
                doRoundness = gd.getNextBoolean();
                doArEquivD = gd.getNextBoolean();
                doPerEquivD = gd.getNextBoolean();
                doEquivEllAr = gd.getNextBoolean();
                doCompactness = gd.getNextBoolean();
                doConvexity = gd.getNextBoolean();
                doShape = gd.getNextBoolean();
                doRFactor = gd.getNextBoolean();
                doArBBox = gd.getNextBoolean();
                doRectang = gd.getNextBoolean();
                doModRatio = gd.getNextBoolean();
                doSphericity = gd.getNextBoolean();
                doElongation = gd.getNextBoolean();
                doNormPeriIndex = gd.getNextBoolean();
                doHeralickRatio = gd.getNextBoolean();
            }
            return true;
        }

        /*Riguardante ConvexHull-
         * Rivisitazione del metodo getArea da Analyzer(super-super classe) per cui dati i punti del poligono calcola l'area
         * trattandolo come se fosse composto da danti piccoli triangoli*/
        final double getArea(Polygon p) {
            if (p==null) return Double.NaN;
            int carea = 0;
            int iminus1;

            for(int i =0; i<p.npoints-1;i++){
                iminus1=i-1;
                if(iminus1<0) iminus1=p.npoints-1;
                carea+=(p.xpoints[i]+p.xpoints[iminus1])*(p.ypoints[i]-p.ypoints[iminus1]);
            }
            return (Math.abs(carea/2.0));
        }

        /*Riguardante ConvexHull-
         * Calcolo del perimetro dato i punti del poligono e calcolando le distanze tra i punti*/
        final double getPerimeter (Polygon p){
            if(p==null) return  Double.NaN;

            double cperimeter = 0.0;
            int iminus1;

            for(int i=0; i < p.npoints-1; i++){
                iminus1=i-1;
                if(iminus1<0) iminus1 = p.npoints-1;
                cperimeter += distance(p.xpoints[i], p.ypoints[i], p.xpoints[iminus1], p.ypoints[iminus1]);
            }
            return cperimeter;
        }

        /*Funzione di appoggio per il caolcolo della distanza*/
        double distance(int argx1, int argy1, int argx2, int argy2){
            return Math.sqrt((argx1-argx2)*(argx1-argx2)+(argy1-argy2)*(argy1-argy2));
        }

        /*Riguardante il calcolo del Haralick Ratio
        * Sfruttando la trasformazione in convexHull del roi, si prendono i punti delle coordinate x e y e si lavorano su essi:
        * ovunque la y abbia un valore uguale si considera la x
        * purtroppo essi sono disordinati*/
        double getHeralickRatio(Polygon p){
            /*Essendo disordinati bisogna scorrere gli array delle x e delle y*/
            double sumMean=0.0;
            /*raccoglierà tutti i raggi*/
            double[] radii = new double[p.npoints];
            int xSecond=0;
            int numberOfRadius=0;
            for(int i=0; i<p.npoints-1; i++){
                for(int j=i+1; j<p.npoints-1; j++){
                    if(p.ypoints[j]==p.ypoints[i]){
                        if(p.xpoints[i]<p.xpoints[j]){
                            xSecond=p.xpoints[j];
                        }
                    }
                }
                sumMean+=(xSecond-p.xpoints[i])/2;
                radii[numberOfRadius]=(xSecond-p.xpoints[i])/2;
                numberOfRadius++;
            }

            double media= sumMean/numberOfRadius;
            for(int i=0; i<numberOfRadius; i++){
                sumMean+=Math.abs((radii[i]-media)*(radii[i]-media));
            }
            return media/(Math.sqrt(Math.abs(sumMean/numberOfRadius-1)));
        }

        /*Metodo addMeasure necessario per il calcolo delle misure aggiuntive a cui vengono passati i valori già trattati nella classe estesa (ParticleAnalyzer)*/
        void addMeasure(ImageStatistics stats, Roi roi, ResultsTable rt, ImagePlus imp){
            ImageProcessor ip = imp.getProcessor();
            ip.autoThreshold();
            ip.setRoi(roi); //settaggio della zona di interesse
            float[] areas = rt.getColumn(ResultsTable.AREA); //recupero di alcuni valori necessari per poter calcolare le nuove misure
            float[] ferets = rt.getColumn(ResultsTable.FERET); //pescati dalle colonne perchè non risiedono in nessuna variabile richiamabile
            float[] breadths =  rt.getColumn(ResultsTable.MIN_FERET);
            float[] perims = rt.getColumn(ResultsTable.PERIMETER);
            //Riguardante ConvexHull: trattamento di un oggetto Polygon (Riguardante la ROI data)
            double convexArea = getArea(roi.getConvexHull());
            double convexPerimeter = getPerimeter(roi.getConvexHull());
            double heralickRatio = getHeralickRatio(roi.getConvexHull());
            //ciclo per l'aggiunta dei nuovi valori
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
                if (doCircularity) //Circularity = 4·?·Area/Perim2  also called Formfactor or Shapefactor
                    rt.addValue("*Circularity", 4 * areas[i] / (perims[i] * perims[i]));
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
                if(doConcavity) //Concavity ConvexArea-Area
                    rt.addValue("*Concavity", convexArea-areas[i]);
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
                if(doHeralickRatio){
                    rt.addValue("*HeralickRatio", heralickRatio);
                }
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
        }
    }
}
