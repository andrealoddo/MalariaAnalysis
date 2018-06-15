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

/*CatchParasite_Plus2 by Tatalessap
 *
 * Catch_Parasite è una classe volta alla misurazione di elementi all'interno di una immagine (binaria+grey).
 * Si appoggia al Plugin "ThresholdAdjuster" per poter individuare le zone e gli oggetti da analizzare.
 * Implementa Plugin ma al suo interno vi è una classe Parasite che estende a sua volta il PluginFilter "ParticleAnalyzer"
 * aggiungendo misure non presenti in quest'ultimo.*/
public class CatchParasite_Plus2 implements PlugIn {
    
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
    }

    /*Inner class Parasite che estende la classe ParticleAnalyzer
     * Al suo interno, attraverso gli Override, si estendono i meotdi già presenti nella classe ParticleAnalyzer*/
    class Parasite extends ParticleAnalyzer {
        int measure=0;
        private  boolean[] measureImageJ = new boolean[12]; //AREA +AREA_FRACTION +CENTER_OF_MASS +CENTROID
                                                            // +CIRCULARITY +ELLIPSE +FERET +INTEGRATED_DENSITY
                                                            // +MAX_STANDARDS +PERIMETER +SHAPE_DESCRIPTORS +STD_DEV
        /*Vettore per checkBox B&W*/
        private boolean[] measuresBW = new boolean[23];
        /*CheckBox Misure aggiunte B&W*/
        private boolean doConvexArea, doConvexPerimeter, doMINRMAXR,
                doAspRatio, doCircularity, doRoundness,
                doArEquivD, doPerEquivD, doEquivEllAr,
                doCompactness, doSolidity, doConcavity,
                doConvexity, doShape, doRFactor,
                doArBBox, doRectang, doModRatio,
                doSphericity, doElongation, doNormPeriIndex,
                doHaralickRatio, doBendingEnergy;

        /*Vettore per checkBox GREY*/
        private boolean[] measuresGrey = new boolean[10];
        /*CheckBox Misure aggiunte GREY*/
        private boolean doMean, doSkewness, doKurtois,
                doMode, doMedian, doMax,
                doMin, doEntropy, doIntensitySum,
                doSquareIntensitySum;

        /*piGreco necessario per alcuni calcoli*/
        double pigreco = Math.PI;

        //Salvataggio delle misure già settate in SetMeasurements
        int old_measures = Analyzer.getMeasurements();

        /*showDialog:
         * settaggio della opzione per visualizzare contorni e etichetta numerata
         * settaggio nelle misure necessarie
         * richiamo della genericDialog creata per poter settare le misure da aggiungere*/
        @Override
        public boolean showDialog() {
            super.staticShowChoice = 1;
            boolean flag = super.showDialog();
            if (flag) {
                flag= genericDialog();
                setMeasuresExtended();
                Analyzer.setMeasurements(measure);
                return flag;
            }
            return false;
        }

        /*saveResults:
         *richiamo della funzione originale e aggiunta delle nuove misure*/
        protected void saveResults(ImageStatistics stats, Roi roi) {
            super.saveResults(stats, roi);
            /*Settaggio delle misure da aggiungere date dalla checkBox*/
            /*Oggetti e dati necessari per alcune misure*/
            Polygon polygon = roi.getConvexHull();
            double convexArea = getArea(polygon);
            double convexPerimeter = getPerimeter(polygon);
            double [] feret = roi.getFeretValues(); //0 --> feret , 1 --> angle, 2 --> feret min
            double perim= roi.getLength();
            int[] hist = stats.histogram;

            if (doConvexArea) {//Area of the convex hull polygon
                rt.addValue("*ConvexArea", convexArea);
            }
            if (doConvexPerimeter) { //Perimeter of the convex hull polygon
                rt.addValue("*PerimeterConvexHull", convexPerimeter);
            }

            if (doMINRMAXR) { //
                rt.addValue("*MinR", feret[2]/2); //DA RIVEDERE Radius of the inscribed circle centred at the middle of mass
                rt.addValue("*MaxR", feret[0]/2); //DA RIVEDERE Radius of the enclosing circle centred at the middle of mass
            }

            if (doAspRatio) //Aspect ratio = Feret/Breadth = L/W also called Feret ratio or Eccentricity or Rectangular ratio
                rt.addValue("*AspRatio", feret[0]/feret[2]);

            if (doCircularity) //Circularity = 4·?·Area/Perim2  also called Formfactor or Shapefactor
                rt.addValue("*Circularity", (4 * pigreco * stats.area) / (perim*perim));

            if (doRoundness) //Roundness = 4·Area/(?·Feret2)
                rt.addValue("*Roundness", (stats.area*4) / ((pigreco) * (feret[0]*feret[0])));

            if (doArEquivD) //Diameter of a circle with equivalent area,
                rt.addValue("*ArEquivD", Math.sqrt((4 * pigreco) * stats.area));

            if (doPerEquivD) //Diameter of a circle with equivalent perimeter,  Area/?
                rt.addValue("*PerEquivD", stats.area / pigreco);

            if (doEquivEllAr) //Area of the ellipse with Feret and Breath as major and minor axis,  = (?·Feret·Breadth)/4
                rt.addValue("*EquivEllAr", (pigreco * feret[0] * feret[2]) / 4);

            if (doCompactness) //Compactness = ?((4/?)·Area)/Feret
                rt.addValue("*Compactness", (Math.sqrt((4 / pigreco) * stats.area)) / feret[0]);

            if (doSolidity) //Compactness = ?((4/?)·Area)/Feret
                rt.addValue("*Solidity", (stats.area / convexArea));

            if (doConcavity) //Concavity ConvexArea-Area
                rt.addValue("*Concavity", convexArea - stats.area);

            if (doConvexity) //Convexity = Convex_Perim/Perimeter also called rugosity or roughness
                rt.addValue("*Convexity", convexPerimeter / perim);

            if (doShape) //Shape = Perimeter2/Area also called Thinness ratio
                rt.addValue("*Shape", (perim*perim) / stats.area);

            if (doRFactor) //RFactor = Convex_Area /(Feret·?)
                rt.addValue("*RFactor", convexArea / (feret[0] * pigreco));

            if (doArBBox) //Area of the bounding box along the Feret diameter = Feret·Breadth
                rt.addValue("*ArBBox", feret[0] * feret[2]);

            if (doRectang) //Rectangularity = Area/ArBBox also called Extent
                rt.addValue("*Rectang", stats.area / (feret[0] * feret[2]));

            if (doModRatio) //Modification ratio = (2·MinR)/Feret
                rt.addValue("*ModRatio", (feret[2] / feret[0])); //2 * MinR / Feret DA RIVEDERE dannorisultati uguali

            if (doSphericity) //Sphericity = MinR/MaxR also called Radius ratio
                rt.addValue("*Sphericity", (feret[2] / 2) / (feret[0] / 2)); //MinR / MaxR DA RIVEDERE danno risultati ugualu

            if (doElongation) //The inverse of the circularity,  Perim2/(4·?·Area)
                rt.addValue("*Elongation", (perim*perim) / (4 * pigreco * stats.area));

            if (doNormPeriIndex)
                rt.addValue("*normPeriIndex", (2 * Math.sqrt(pigreco * stats.area)) / perim);

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

            /**/
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

            if(doEntropy){
                double e = getEntropy(hist, stats.area);
                rt.addValue("**Entropy", e);
            }

            if(doIntensitySum){
                int intensitySum = getIntensitySum(hist);
                rt.addValue("**Intensity sum", intensitySum);
            }

            if(doSquareIntensitySum){
                int intensitySum = getIntensitySum(hist);
                rt.addValue("**Square intensity sum (SqI sum)", Math.sqrt((double) intensitySum));
            }

            Analyzer.setMeasurements(old_measures); //per risettare misure vecchie
        }

        /*Metodo della creazione della finestra di dialogo, visualizzazione del:
         * - checkbox per selezionare ogni misura
         * - checkbox per selezionare misure singole*/
        private boolean genericDialog() {
            GenericDialog gd = new GenericDialog("Parassite Prova", IJ.getInstance());
            gd.addStringField("Title: ", "Catch parasite");
            gd.addMessage("Choise measures for B&W");

            gd.addCheckbox("SelectAll", false); //tutte le misure
            //AREA +AREA_FRACTION +CENTER_OF_MASS +CENTROID
            // +CIRCULARITY +ELLIPSE +FERET +INTEGRATED_DENSITY
            // +MAX_STANDARDS +PERIMETER +SHAPE_DESCRIPTORS
            // +SLICE +STD_DEV

            //Misure ImageJ
            gd.addCheckbox("Area", false);
            gd.addToSameRow();
            gd.addCheckbox("Center of Mass", false);
            gd.addToSameRow();
            gd.addCheckbox("Integrated Density", false);
            gd.addCheckbox("Area_Fraction", false);
            gd.addToSameRow();
            gd.addCheckbox("Centroid", false);
            gd.addToSameRow();
            gd.addCheckbox("Max_Standards", false);
            gd.addCheckbox("Feret", false);
            gd.addToSameRow();
            gd.addCheckbox("Circularity", false);
            gd.addToSameRow();
            gd.addCheckbox("Std dev", false);
            gd.addCheckbox("Perimeter", false);
            gd.addToSameRow();
            gd.addCheckbox("Ellipse", false);
            gd.addToSameRow();
            gd.addCheckbox("Shape Descriptors", false);
            //

            gd.addCheckbox("Convex Area*", false);
            gd.addToSameRow();
            gd.addCheckbox("AspRatio*", false);
            gd.addToSameRow();
            gd.addCheckbox("ArEquivD*", false);
            gd.addCheckbox("Convex Perimeter*", false);
            gd.addToSameRow();
            gd.addCheckbox("Circ*", false);
            gd.addToSameRow();
            gd.addCheckbox("PerEquivD*", false);
            gd.addCheckbox("MinR and MaxR*", false);
            gd.addToSameRow();
            gd.addCheckbox("Roundness*", false);
            gd.addToSameRow();
            gd.addCheckbox("EquivEllAr*", false);
            gd.addCheckbox("Compactness*", false);
            gd.addToSameRow();
            gd.addCheckbox("Solidity*", false);
            gd.addToSameRow();
            gd.addCheckbox("Shape*", false);
            gd.addCheckbox("Convexity*", false);
            gd.addToSameRow();
            gd.addCheckbox("RFactor*", false);
            gd.addToSameRow();
            gd.addCheckbox("ArBBox*", false);
            gd.addCheckbox("Concavity*", false);
            gd.addToSameRow();
            gd.addCheckbox("Rectang*", false);
            gd.addToSameRow();
            gd.addCheckbox("ModRatio*", false);
            gd.addCheckbox("Sphericity*", false);
            gd.addToSameRow();
            gd.addCheckbox("Elongation*", false);
            gd.addToSameRow();
            gd.addCheckbox("NormPeriIndex*", false);
            gd.addCheckbox("HaralickRatio*", false);
            gd.addToSameRow();
            gd.addCheckbox("Bending Energy*", false);

            gd.addMessage("Choise measures for grey");
            gd.addCheckbox("SelectAll", true); //tutte le misure

            gd.addCheckbox("Mean", false);
            gd.addToSameRow();
            gd.addCheckbox("Skewness", false);
            gd.addToSameRow();
            gd.addCheckbox("Intensity Sum*", false);
            gd.addCheckbox("Mode", false);
            gd.addToSameRow();
            gd.addCheckbox("Kurtosis", false);
            gd.addToSameRow();
            gd.addCheckbox("Square Intensity Sum*", false);
            gd.addCheckbox("Median", false);
            gd.addToSameRow();
            gd.addCheckbox("Entropy*", false);
            gd.addCheckbox("Min", false);
            gd.addCheckbox("Max", false);


            gd.showDialog(); //show
            if (gd.wasCanceled())
                return false;

            boolean doSelectAllBW = gd.getNextBoolean();  //se viene selezionata la voce "seleziona tutti"
            boolean spam;

            if (doSelectAllBW)
            {
                for (int i = 0; i < measureImageJ.length; i++){
                    measureImageJ[i]=true;
                    spam = gd.getNextBoolean();
                }

                for (int i = 0; i < measuresBW.length; i++){
                    measuresBW[i]=true;
                    spam = gd.getNextBoolean();
                }
            } else {

                for (int i = 0; i < measureImageJ.length; i++){
                    measureImageJ[i] = gd.getNextBoolean();
                }

                for (int i = 0; i < measuresBW.length; i++){
                    measuresBW[i] = gd.getNextBoolean();
                }

            }

            boolean doSelectAllGrey = gd.getNextBoolean();  //se viene selezionata la voce "seleziona tutti"

            if (doSelectAllGrey) {
                for (int i = 0; i < measuresGrey.length; i++){
                    measuresGrey[i]=true;
                    spam = gd.getNextBoolean();
                }
            }else{
                for (int i = 0; i < measuresGrey.length; i++){
                    measuresGrey[i]=gd.getNextBoolean();
                }
            }

            return true;
        }


        /**/
        private void setMeasuresExtended() {
            //
            if (measureImageJ[0]==true) {measure+=AREA;}

            if (measureImageJ[1]==true) {measure+=CENTER_OF_MASS;}

            if (measureImageJ[2]==true) {measure+=INTEGRATED_DENSITY;}

            if (measureImageJ[3]==true) {measure+=AREA_FRACTION;}

            if (measureImageJ[4]==true) {measure+=CENTROID;}

            if (measureImageJ[5]==true) {measure+=MAX_STANDARDS-MIN_MAX;}

            if (measureImageJ[6]==true) {measure+=FERET;}

            if (measureImageJ[7]==true) {measure+=CIRCULARITY;}

            if (measureImageJ[8]==true) {measure+=STD_DEV;}

            if (measureImageJ[9]==true) {measure+=PERIMETER;}

            if (measureImageJ[10]==true) {measure+=ELLIPSE;}

            if (measureImageJ[11]==true) {measure+=SHAPE_DESCRIPTORS;}

            //measure = measure -SKEWNESS-KURTOSIS-NaN_EMPTY_CELLS-INVERT_Y-MEAN-MEDIAN-MODE-MIN_MAX;

            //MISURE BW
            doConvexArea = measuresBW[0];
            doAspRatio = measuresBW[1];
            doArEquivD = measuresBW[2];
            doConvexPerimeter = measuresBW[3];
            doCircularity = measuresBW[4];
            doPerEquivD = measuresBW[5];
            doMINRMAXR = measuresBW[6];
            doRoundness= measuresBW[7];
            doEquivEllAr = measuresBW[8];
            doCompactness = measuresBW[9];
            doSolidity = measuresBW[10];
            doShape = measuresBW[11];
            doConvexity = measuresBW[12];
            doRFactor = measuresBW[13];
            doArBBox = measuresBW[14];
            doConcavity = measuresBW[15];
            doRectang = measuresBW[16];
            doModRatio = measuresBW[17];
            doSphericity = measuresBW[18];
            doElongation = measuresBW[19];
            doNormPeriIndex = measuresBW[20];
            doHaralickRatio = measuresBW[21];
            doBendingEnergy = measuresBW[22];

            //Misure Grey
            doMean = measuresGrey [0];
            doSkewness = measuresGrey [1]; if(doSkewness) {measure+=SKEWNESS;}
            doIntensitySum = measuresGrey [2];
            doMode = measuresGrey [3];
            doKurtois = measuresGrey [4]; if(doKurtois) {measure+=KURTOSIS;}
            doSquareIntensitySum = measuresGrey [5];
            doMedian = measuresGrey [6]; if(doMedian) {measure+=MEDIAN;}
            doEntropy = measuresGrey [7];
            doMin = measuresGrey [8];
            doMax = measuresGrey [9];
        }

        /*Riguardante ConvexHull-
         * Rivisitazione del metodo getArea da Analyzer(super-super classe) per cui dati i punti del poligono calcola l'area
         * trattandolo come se fosse composto da danti piccoli triangoli*/
        private double getArea(Polygon p) {
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
        private final double getPerimeter(Polygon p) {
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
        private double getHeralickRatio(Polygon p) {
            /*Essendo disordinati bisogna scorrere gli array delle x e delle y*/
            double sumMean = 0.0;
            /*raccoglierà tutti i raggi*/
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

        private double getEntropy(int[] hist, double area){
            double sum=0.0;
            for(int i=0; i< hist.length; i++){
                if(hist[i]!=0) sum += (double) (hist[i] / area) * log2((double) (hist[i] / area)) ;
            }
            return -sum;
        }
        // log2:  Logarithm base 2
        public double log2(double d) {
            return Math.log(d)/Math.log(2.0);
        }

        private int getIntensitySum(int [] hist){
            int sum=0;
            for(int i=0; i< hist.length; i++){
                sum+=i*hist[i];
            }
            return sum;
        }

        /*

        private double getVariance(int [] hist, int area){
            //media
            int sum = 0;
            for(int i=0; i< hist.length; i++){
                sum+= hist[i];
            }

            double media = (double)sum / (double) area;
            //valore medio del valore che assumono i pixel

            sum=0;

            for(int element: hist){
                sum+=(hist)

            }

        }*/

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
    }
}


