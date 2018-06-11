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

import javax.swing.plaf.metal.MetalLookAndFeel;
import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.awt.Polygon;
import java.lang.Object;
import java.util.Arrays;

/*Catch_Parasite by Tatalessap
 * Catch_Parasite � una classe volta alla misurazione di elementi all'interno di una immagine binaria.
 * Si appoggia al Plugin "ThresholdAdjuster" per poter individuare le zone e gli oggetti da analizzare.
 * Implementa Plugin ma al suo interno vi � una classe Parasite che estende a sua volta il PluginFilter "ParticleAnalyzer"
 * aggiungendo misure non presenti.*/

public class Catch_Parasite_Binary_ implements PlugIn {
    /*Metodo run necessaria per i PlugIn
     * Memorizza l'immagine in ingresso
     * richiama il metodo catch_parasite_running dandogli in ingresso l'immagine*/
    public void run(String arg) {
        Parasite parasite = new Parasite();
        ImagePlus imp = IJ.getImage();
        int flags = parasite.setup("", imp);
        /*Controllo*/
        if (flags == PlugInFilter.DONE)
            return;
        parasite.run(imp.getProcessor());
        /*visualizzazione dei risultati*/
        Analyzer.getResultsTable().show("Results");
        //catch_parasite_running(imp); //metodo presente dopo

    }

    /*Inner class Parasite che estende la classe ParticleAnalyzer
     * Al suo interno, attraverso gli Override, si estendono i meotdi gi� presenti nella classe ParticleAnalyzer*/
    class Parasite extends ParticleAnalyzer {
        /*Variabili booleani corrispondenti alle misure implementate e aggiunte.*/
        private boolean[] measuresBW = new boolean[23];

        private boolean
                doConvexArea,
                doConvexPerimeter,
                doMINRMAXR,
                doAspRatio,
                doCircularity,
                doRoundness,
                doArEquivD,
                doPerEquivD,
                doEquivEllAr,
                doCompactness,
                doSolidity,
                doConcavity,
                doConvexity,
                doShape,
                doRFactor,
                doArBBox,
                doRectang,
                doModRatio,
                doSphericity,
                doElongation,
                doNormPeriIndex,
                doHaralickRatio,
                doBendingEnergy;

        private boolean[] measuresGrey = new boolean[7];

        private boolean doMean,
                doSkewness,
                doKurtois,
                doMode,
                doMedian,
                doMax,
                doMin;

        /*piGreco necessario per alcuni calcoli*/
        double pigreco = Math.PI;
        //Salvataggio delle misure gi� settate in SetMeasurements
        int old_measures = Analyzer.getMeasurements();
        //Salvataggio delle nuove misure necessarie per il calcolo delle nuove
        int necessary_measures =
                ALL_STATS
                        + SHAPE_DESCRIPTORS
                        - SKEWNESS
                        - STD_DEV
                        - KURTOSIS
                        - NaN_EMPTY_CELLS
                        - INVERT_Y
                        -MEAN
                        -MEDIAN
                        -MODE;
        //**//


        //**//

        /*showDialog:
         * settaggio della opzione per visualizzare contorni e etichetta numerata
         * settaggio nelle misure necessarie
         * richiamo della genericDialog creata per poter settare le misure da aggiungere*/
        @Override
        public boolean showDialog() {
            super.staticShowChoice = 1;
            boolean flag = super.showDialog();
            //*protected static final int NOTHING=0, OUTLINES=1, BARE_OUTLINES=2, ELLIPSES=3, MASKS=4, ROI_MASKS=5,
            //                             OVERLAY_OUTLINES=6, OVERLAY_MASKS=7;*//
            //setto solo quello che mi serve
            Analyzer.setMeasurements(necessary_measures);
            if (flag) return genericDialog();
            return false;

        }

        /*saveResults:
         * richiamo della funzione originaria
         * richiamo della funzione aggiuntiva addMeasure a cui vengono passati roi e stats.
         * roi: regione di interesse
         * stats: oggetto di tipo ImageStatistics
         * settaggio delle vecchie misure dopo che si � concluso il calcolo*/
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
        protected boolean genericDialog() {
            GenericDialog gd = new GenericDialog("Parassite Prova", IJ.getInstance());
            gd.addStringField("Title: ", "Catch parasite");
            gd.addMessage("Choise measures for B&W");

            gd.addCheckbox("SelectAll", true); //tutte le misure

            gd.addCheckbox("Convex Area", false);
            gd.addToSameRow();
            gd.addCheckbox("Convex Perimeter", false);
            gd.addToSameRow();
            gd.addCheckbox("MinR and MaxR", false);

            gd.addCheckbox("AspRatio", false);
            gd.addToSameRow();
            gd.addCheckbox("Circ", false);
            gd.addToSameRow();
            gd.addCheckbox("Roundness", false);

            gd.addCheckbox("ArEquivD", false);
            gd.addToSameRow();
            gd.addCheckbox("PerEquivD", false);
            gd.addToSameRow();
            gd.addCheckbox("EquivEllAr", false);

            gd.addCheckbox("Compactness", false);
            gd.addToSameRow();
            gd.addCheckbox("Solidity", false);
            gd.addToSameRow();
            gd.addCheckbox("Concavity", false);


            gd.addCheckbox("Convexity", false);
            gd.addToSameRow();
            gd.addCheckbox("Shape", false);
            gd.addToSameRow();
            gd.addCheckbox("RFactor", false);


            gd.addCheckbox("ArBBox", false);
            gd.addToSameRow();
            gd.addCheckbox("Rectang", false);
            gd.addToSameRow();
            gd.addCheckbox("ModRatio", false);


            gd.addCheckbox("Sphericity", false);
            gd.addToSameRow();
            gd.addCheckbox("Elongation", false);
            gd.addToSameRow();
            gd.addCheckbox("NormPeriIndex", false);


            gd.addCheckbox("HaralickRatio", false);
            gd.addToSameRow();
            gd.addCheckbox("Bending Energy", false);

            gd.addMessage("Choise measures for grey");
            gd.addCheckbox("SelectAll", true); //tutte le misure

            gd.addCheckbox("Mean", false);
            gd.addToSameRow();
            gd.addCheckbox("Skewness", false);
            gd.addToSameRow();
            gd.addCheckbox("Kurtois", false);

            gd.addCheckbox("Mode", false);
            gd.addToSameRow();
            gd.addCheckbox("Median", false);
            gd.addToSameRow();
            gd.addCheckbox("Max", false);

            gd.addCheckbox("Min", false);


            gd.showDialog(); //show
            if (gd.wasCanceled())
                return false;
            boolean doSelectAllBW = gd.getNextBoolean();  //se viene selezionata la voce "seleziona tutti"
            if (doSelectAllBW) { // vengono settati i booleani a true
                for (int i = 0; i < measuresBW.length; i++) {
                    measuresBW[i] = true;
                }
            } else { //altrimenti pescati uno per uno con un certo ordine
                for (int i = 0; i < measuresBW.length; i++) {
                    measuresBW[i] = gd.getNextBoolean();
                    ;
                }
            }

            boolean doSelectAllGrey = gd.getNextBoolean();  //se viene selezionata la voce "seleziona tutti"
            if (doSelectAllBW) { // vengono settati i booleani a true
                for (int i = 0; i < measuresGrey.length; i++) {
                    measuresGrey[i] = true;
                }
            } else { //altrimenti pescati uno per uno con un certo ordine
                for (int i = 0; i < measuresGrey.length; i++) {
                    measuresGrey[i] = gd.getNextBoolean();
                    ;
                }
            }




            return true;
        }

        private void setMeasuresExtended() {
            doConvexArea = measuresBW[0];
            doConvexPerimeter = measuresBW[1];
            doMINRMAXR = measuresBW[2];
            doAspRatio = measuresBW[3];
            doCircularity = measuresBW[4];
            doRoundness = measuresBW[5];
            doArEquivD = measuresBW[6];
            doPerEquivD = measuresBW[7];
            doEquivEllAr = measuresBW[8];
            doCompactness = measuresBW[9];
            doSolidity = measuresBW[10];
            doConcavity = measuresBW[11];
            doConvexity = measuresBW[12];
            doShape = measuresBW[13];
            doRFactor = measuresBW[14];
            doArBBox = measuresBW[15];
            doRectang = measuresBW[16];
            doModRatio = measuresBW[17];
            doSphericity = measuresBW[18];
            doElongation = measuresBW[19];
            doNormPeriIndex = measuresBW[20];
            doHaralickRatio = measuresBW[21];
            doBendingEnergy = measuresBW[22];

            doMean = measuresGrey[0];
            doSkewness= measuresGrey[1];
            doKurtois=measuresGrey[2];
            doMode=measuresGrey[3];
            doMedian=measuresGrey[4];
            doMax=measuresGrey[5];
            doMin=measuresGrey[6];
        }

        /*Riguardante ConvexHull-
         * Rivisitazione del metodo getArea da Analyzer(super-super classe) per cui dati i punti del poligono calcola l'area
         * trattandolo come se fosse composto da danti piccoli triangoli*/
        protected double getArea(Polygon p) {
            if (p == null) return Double.NaN;
            int carea = 0;
            int iminus1;

            for (int i = 0; i < p.npoints - 1; i++) {
                iminus1 = i - 1;
                if (iminus1 < 0) iminus1 = p.npoints - 1;
                carea += (p.xpoints[i] + p.xpoints[iminus1]) * (p.ypoints[i] - p.ypoints[iminus1]);
            }
            return (Math.abs(carea / 2.0));
        }

        /*Riguardante ConvexHull-
         * Calcolo del perimetro dato i punti del poligono e calcolando le distanze tra i punti*/
        protected final double getPerimeter(Polygon p) {
            if (p == null) return Double.NaN;

            double cperimeter = 0.0;
            int iminus1;

            for (int i = 0; i < p.npoints - 1; i++) {
                iminus1 = i - 1;
                if (iminus1 < 0) iminus1 = p.npoints - 1;
                cperimeter += distance(p.xpoints[i], p.ypoints[i], p.xpoints[iminus1], p.ypoints[iminus1]);
            }
            return cperimeter;
        }

        /*Funzione di appoggio per il caolcolo della distanza*/
        public double distance(int argx1, int argy1, int argx2, int argy2) {
            return Math.sqrt((argx1 - argx2) * (argx1 - argx2) + (argy1 - argy2) * (argy1 - argy2));
        }

        /*METODI PER I VETTORI*/
        public double[] diff(double[] z) {
            double first = z[0];
            double ultimade = z[z.length - 1];
            double[] result = new double[z.length];

            for (int i = 0; i < z.length; i++) {
                if (i == (z.length - 1)) {
                    result[i] = (ultimade - first);
                    //result[i]=Math.abs(ultimade-first);
                } else {
                    result[i] = (z[i + 1] - z[i]);
                    //result[i]=Math.abs(ultimade-first);
                }
            }
            return result;
        }

        public double[] diffVector(double[] a, double[] b) {
            double[] result = new double[a.length];
            for (int i = 0; i < a.length; i++) {
                result[i] = a[i] - b[i];
                //result[i]=Math.abs(a[i]-b[i]);
            }

            return result;
        }

        public double[] moltVector(double[] a, double[] b) { //ok
            double[] result = new double[a.length];
            for (int i = 0; i < a.length; i++) {
                result[i] = a[i] * b[i];
            }

            return result;
        }

        public double[] divVector(double[] a, double[] b) {
            double[] result = new double[a.length];
            for (int i = 0; i < a.length; i++) {
                result[i] = a[i] * Math.pow(b[i], -1);
                //result[i]=Math.abs(a[i]*Math.pow(b[i], -1));
            }
            return result;
        }

        public double[] sumVector(double[] a, double[] b) {
            double[] result = new double[a.length];
            for (int i = 0; i < a.length; i++) {
                result[i] = (a[i] + b[i]);
                //result[i]=Math.abs(a[i]+b[i]);
            }
            return result;
        }

        public double[] elevationVector(double[] a, double elevation) {
            double[] result = new double[a.length];
            for (int i = 0; i < a.length; i++) {
                result[i] = Math.pow(a[i], elevation);
            }
            return result;
        }

        private double getBendingEnergy(Polygon polygon) {
            //SBAGLIATO DA RICONTROLLARE
            double bendingEnergy = 0;
            double[] k = new double[polygon.npoints];
            int i;

            /*Trasformazione in double dei vettori x e y..
             * */
            double[] x = new double[polygon.npoints];
            for (i = 0; i < x.length; i++) {
                x[i] = (double) polygon.xpoints[i];
            }

            double[] y = new double[polygon.npoints];
            for (i = 0; i < y.length; i++) {
                y[i] = (double) polygon.ypoints[i];
            }

            k = divVector(
                    diffVector
                            (
                                    moltVector(
                                            (diff(x)), (diff(diff(y)))
                                    )
                                    ,
                                    moltVector(
                                            (diff(y)), (diff(diff(x)))
                                    )
                            )
                    ,
                    elevationVector(
                            (sumVector(
                                    elevationVector(diff(x), 2),
                                    elevationVector(diff(y), 2))),
                            1.5
                    )

            );


            for (i = 0; i < k.length; i++) {
                bendingEnergy = bendingEnergy + Math.pow(k[i], 2);
            }

            return bendingEnergy * Math.pow(k.length, -1);
        }
        /*Riguardante il calcolo del Haralick Ratio
         * Sfruttando la trasformazione in convexHull del roi, si prendono i punti delle coordinate x e y e si lavorano su essi:
         * ovunque la y abbia un valore uguale si considera la x
         * purtroppo essi sono disordinati*/
        protected double getHeralickRatio(Polygon p) {
            /*Essendo disordinati bisogna scorrere gli array delle x e delle y*/
            double sumMean = 0.0;
            /*raccoglier� tutti i raggi*/
            double[] radii = new double[p.npoints];
            int xSecond = 0;
            int numberOfRadius = 0;
            for (int i = 0; i < p.npoints - 1; i++) {
                for (int j = i + 1; j < p.npoints - 1; j++) {
                    if (p.ypoints[j] == p.ypoints[i]) {
                        if (p.xpoints[i] < p.xpoints[j]) {
                            xSecond = p.xpoints[j];
                        }
                    }
                }
                sumMean += (xSecond - p.xpoints[i]) / 2;
                radii[numberOfRadius] = (xSecond - p.xpoints[i]) / 2;
                numberOfRadius++;
            }

            double media = sumMean / numberOfRadius;
            for (int i = 0; i < numberOfRadius; i++) {
                sumMean += Math.abs((radii[i] - media) * (radii[i] - media));
            }
            return media / (Math.sqrt(Math.abs(sumMean / numberOfRadius - 1)));
        }

        /*Metodo addMeasure necessario per il calcolo delle misure aggiuntive a cui vengono passati i valori gi� trattati nella classe estesa (ParticleAnalyzer)*/
        protected void addMeasure(ImageStatistics stats, Roi roi, ResultsTable rt, ImagePlus imp) {
            ImageProcessor ip = imp.getProcessor();
            ip.autoThreshold();
            setMeasuresExtended();
            ip.setRoi(roi); //settaggio della zona di interesse

            float[] areas = rt.getColumn(ResultsTable.AREA); //recupero di alcuni valori necessari per poter calcolare le nuove misure
            float[] ferets = rt.getColumn(ResultsTable.FERET); //pescati dalle colonne perch� non risiedono in nessuna variabile richiamabile
            float[] breadths = rt.getColumn(ResultsTable.MIN_FERET);
            float[] perims = rt.getColumn(ResultsTable.PERIMETER);
            //Riguardante ConvexHull: trattamento di un oggetto Polygon (Riguardante la ROI data)
            Polygon polygon = roi.getConvexHull();
            double convexArea = getArea(polygon);
            double convexPerimeter = getPerimeter(polygon);

            //ciclo per l'aggiunta dei nuovi valori
            for (int i = 0; i < areas.length; i++) {
                // double [] minR_maxR = getMinRMaxR(stats.xCenterOfMass, stats.yCenterOfMass);
                if (doConvexArea) {//Area of the convex hull polygon
                    rt.addValue("*ConvexArea", convexArea);
                }
                if (doConvexPerimeter) { //Perimeter of the convex hull polygon
                    rt.addValue("*PerimeterConvexHull", convexPerimeter);
                }
                //MINR E MAXR PROVA, inteso come raggio del cerchio pi� piccolo e pi� grande prendendo come diametro il FERET e il BREADTH
                if (doMINRMAXR) { //
                    rt.addValue("*MinR", breadths[i] / 2); //DA RIVEDERE Radius of the inscribed circle centred at the middle of mass
                    rt.addValue("*MaxR", ferets[i] / 2); //DA RIVEDERE Radius of the enclosing circle centred at the middle of mass
                }
                if (doAspRatio) //Aspect ratio = Feret/Breadth = L/W also called Feret ratio or Eccentricity or Rectangular ratio
                    rt.addValue("*AspRatio", ferets[i] / breadths[i]);
                if (doCircularity) //Circularity = 4�?�Area/Perim2  also called Formfactor or Shapefactor
                    rt.addValue("*Circularity", 4 * pigreco * areas[i] / (perims[i] * perims[i]));
                if (doRoundness) //Roundness = 4�Area/(?�Feret2)
                    rt.addValue("*Roundness", (areas[i] * 4) / ((pigreco) * (ferets[i] * ferets[i])));
                if (doArEquivD) //Diameter of a circle with equivalent area,
                    rt.addValue("*ArEquivD", Math.sqrt((4 * pigreco) * areas[i]));
                if (doPerEquivD) //Diameter of a circle with equivalent perimeter,  Area/?
                    rt.addValue("*PerEquivD", areas[i] / pigreco);
                if (doEquivEllAr) //Area of the ellipse with Feret and Breath as major and minor axis,  = (?�Feret�Breadth)/4
                    rt.addValue("*EquivEllAr", (pigreco * ferets[i] * breadths[i]) / 4);
                if (doCompactness) //Compactness = ?((4/?)�Area)/Feret
                    rt.addValue("*Compactness", (Math.sqrt((4 / pigreco) * areas[i])) / ferets[i]);
                if (doSolidity) //Compactness = ?((4/?)�Area)/Feret
                    rt.addValue("*Solidity", (areas[i] / convexArea));
                if (doConcavity) //Concavity ConvexArea-Area
                    rt.addValue("*Concavity", convexArea - areas[i]);
                if (doConvexity) //Convexity = Convex_Perim/Perimeter also called rugosity or roughness
                    rt.addValue("*Convexity", convexPerimeter / perims[i]);
                if (doShape) //Shape = Perimeter2/Area also called Thinness ratio
                    rt.addValue("*Shape", (perims[i] * perims[i]) / areas[i]);
                if (doRFactor) //RFactor = Convex_Area /(Feret�?)
                    rt.addValue("*RFactor", convexArea / (ferets[i] * pigreco));
                if (doArBBox) //Area of the bounding box along the Feret diameter = Feret�Breadth
                    rt.addValue("*ArBBox", ferets[i] * breadths[i]);
                if (doRectang) //Rectangularity = Area/ArBBox also called Extent
                    rt.addValue("*Rectang", areas[i] / (ferets[i] * breadths[i]));
                if (doModRatio) //Modification ratio = (2�MinR)/Feret
                    rt.addValue("*ModRatio", (breadths[i] / ferets[i])); //2 * MinR / Feret DA RIVEDERE dannorisultati uguali
                if (doSphericity) //Sphericity = MinR/MaxR also called Radius ratio
                    rt.addValue("*Sphericity", (breadths[i] / 2) / (ferets[i] / 2)); //MinR / MaxR DA RIVEDERE danno risultati ugualu
                if (doElongation) //The inverse of the circularity,  Perim2/(4�?�Area)
                    rt.addValue("*Elongation", (perims[i] * perims[i]) / (4 * pigreco * areas[i]));
                if (doNormPeriIndex)
                    rt.addValue("*normPeriIndex", (2 * Math.sqrt(pigreco * areas[i])) / perims[i]);
                if (doHaralickRatio) {
                    double heralickRatio = getHeralickRatio(polygon);
                    rt.addValue("*HaralickRatio", heralickRatio);
                }
                if (doBendingEnergy) {
                    double be = getBendingEnergy(polygon);
                    rt.addValue("*Bending Energy", be);
                    /*
                    Polygon polygon=roi.getConvexHull();

                    double[] x = new double[polygon.npoints];
                    for(int a =0;a<x.length;a++) {
                        x[a] = (double) polygon.xpoints[a];
                    }
                    double [] y = new double[polygon.npoints];
                    for(int c=0;c<y.length;c++) {
                        y[c] = (double) polygon.ypoints[c];
                    }
                    double [] diffX1= diff(x);
                    double [] diffX2= diff(diffX1);
                    double [] diffY1= diff(y);
                    double [] diffY2= diff(diffY1);
                    for(int b=0; b<polygon.npoints; b++){
                        rt.addValue("**X    "+b, x[b]);
                        rt.addValue("**XDIFF1   "+b, diffX1[b]);
                        rt.addValue("**XDIFF2   "+b, diffX2[b]);
                        rt.addValue("**Y    "+b, y[b]);
                        rt.addValue("**YDIFF1   "+b, diffY1[b]);
                        rt.addValue("**YDIFF2   "+b, diffY2[b]);
                       };
                */
                }

                rt.addValue("**", "**");

                if(doMean){
                    rt.addValue("**Mean",stats.mean);
                }
                if(doSkewness){
                    rt.addValue("**Skewness",stats.skewness);
                }
                if(doKurtois){
                    rt.addValue("**Kurtois",stats.kurtosis);
                }

                if(doMode){
                    rt.addValue("**Mode",stats.dmode);
                }

                if (doMedian){
                    rt.addValue("**Median",stats.median);
                }

                if(doMax){
                    rt.addValue("**Max",stats.max);
                }

                if(doMin){
                    rt.addValue("**Min",stats.min );
                }




            }
        }
    }
}
