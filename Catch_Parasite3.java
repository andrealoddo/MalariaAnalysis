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

import java.lang.String;
import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.awt.Polygon;
import java.lang.Object;
import java.util.Arrays;

/*Catch_Parasite by Tatalessap
 *
 * Catch_Parasite è una classe volta alla misurazione di elementi all'interno di una immagine (binaria+grey).
 * Si appoggia al Plugin "ThresholdAdjuster" per poter individuare le zone e gli oggetti da analizzare.
 * Implementa Plugin ma al suo interno vi è una classe Parasite che estende a sua volta il PluginFilter "ParticleAnalyzer"
 * aggiungendo misure non presenti in quest'ultimo.*/
public class Catch_Parasite3 implements PlugIn {

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

        private  boolean[] measureImageJ = new boolean[12];
        /*Vettore per checkBox B&W*/
        private boolean[] measuresBW = new boolean[22];
        /*CheckBox Misure aggiunte B&W*/
        private boolean doConvexArea, doConvexPerimeter, doMINRMAXR,
                doAspRatio, doRoundness,
                doArEquivD, doPerEquivD, doEquivEllAr,
                doCompactness, doSolidity, doConcavity,
                doConvexity, doShape, doRFactor,
                doArBBox, doRectang, doModRatio,
                doSphericity, doElongation, doNormPeriIndex,
                doHaralickRatio, doBendingEnergy;

        /*Vettore per checkBox GREY*/
        private boolean[] measuresGrey = new boolean[11];
        /*CheckBox Misure aggiunte GREY*/
        private boolean doMean, doSkewness, doKurtois,
                doMode, doMedian, doMax,
                doMin, doEntropy, doIntensitySum,
                doSquareIntensitySum, doUniformity;

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
            //super.staticShowChoice = 1;
            boolean flag = super.showDialog();
            if (flag) {
                flag= genericDialog();
                setMeasuresExtended();
                Analyzer.setMeasurements(measure);
                return flag;
            }
            return false;
        }

        /*Metodo della creazione della finestra di dialogo, visualizzazione del:
         * - checkbox per selezionare ogni misura
         * - checkbox per selezionare misure singole*/
        private boolean genericDialog() {
            GenericDialog gd = new GenericDialog("Parassite Prova", IJ.getInstance());
            gd.addMessage("Catch parasite");
            gd.addMessage("Choise measures for B&W");

            String[] labels = new String[12];
            boolean[] states = new boolean[12];

            gd.addCheckbox("Select All", true);

            labels[0] = "Area                  ";                  states[0] = false;
            labels[1] = "Standard deviation    ";                  states[1] = false;
            labels[2] = "Centroid              ";                  states[2] = false;
            labels[3] = "Center of mass        ";                  states[3] = false;
            labels[4] = "Perimeter             ";                  states[4] = false;
            labels[5] = "Bounding rectangle    ";                  states[5] = false;
            labels[6] = "Fit ellipse           ";                  states[6] = false;
            labels[7] = "Shape descriptors     ";                  states[7] = false;
            labels[8] = "Feret's diameter      ";                  states[8] = false;
            labels[9] = "Integrated density    ";                  states[9] = false;
            labels[10] = "Area_fraction        ";                  states[10] = false;
            labels[11] = "Stack position       ";                  states[11] = false;
            gd.setInsets(1, 0, 0);
            gd.addCheckboxGroup(5, 3, labels, states);

            gd.addMessage("Choise measures impleted");

            String[] labels1 = new String[11];
            boolean[] states1 = new boolean[11];
            gd.addCheckbox("Select All", true);
            labels1[0] = "Convex Area           ";                  states1[0] = false;
            labels1[1] = "AspRatio*             ";                  states1[1] = false;
            labels1[2] = "ArEquivD*             ";                  states1[2] = false;
            labels1[3] = "Convex Perimeter*     ";                  states1[3] = false;
            labels1[4] = "PerEquivD*            ";                  states1[4] = false;
            labels1[5] = "MinR and MaxR*        ";                  states1[5] = false;
            labels1[6] = "Roundness*            ";                  states1[6] = false;
            labels1[7] = "EquivEllAr*           ";                  states1[7] = false;
            labels1[8] = "Compactness*          ";                  states1[8] = false;
            labels1[9] = "Solidity*             ";                  states1[9] = false;
            labels1[10] = "Shape*               ";                  states1[10] = false;
            gd.setInsets(0, 0, 0);
            gd.addCheckboxGroup(5, 3, labels1, states1);
            String[] labels2 = new String[11];
            boolean[] states2 = new boolean[11];
            labels2[0] = "Convexity*              ";                  states2[0] = false;
            labels2[1] = "RFactor*                ";                  states2[1] = false;
            labels2[2] = "ArBBox*                 ";                  states2[2] = false;
            labels2[3] = "Concavity*              ";                  states2[3] = false;
            labels2[4] = "Rectang*                ";                  states2[4] = false;
            labels2[5] = "ModRatio*               ";                  states2[5] = false;
            labels2[6] = "Sphericity*             ";                  states2[6] = false;
            labels2[7] = "Elongation*             ";                  states2[7] = false;
            labels2[8] = "NormPeriIndex*          ";                  states2[8] = false;
            labels2[9] = "HaralickRatio*          ";                  states2[9] = false;
            labels2[10] = "Bending Energy*        ";                  states2[10] = false;
            gd.setInsets(0, 0, 0);
            gd.addCheckboxGroup(5, 3, labels2, states2);

            gd.addMessage("Choise measures for grey");

            String[] labels3 = new String[11];
            boolean[] states3 = new boolean[11];
            gd.addCheckbox("Select All", true);

            labels3[0] = "Mean                      ";                 states3[0] = false;
            labels3[1] = "Skewness                  ";                 states3[1] = false;
            labels3[2] = "Intensity Sum*            ";                 states3[2] = false;
            labels3[3] = "Mode                      ";                 states3[3] = false;
            labels3[4] = "Kurtosis                  ";                 states3[4] = false;
            labels3[5] = "Entropy*                  ";                 states3[5] = false;
            labels3[6] = "Median                    ";                 states3[6] = false;
            labels3[7] = "SqI sum*                  ";                 states3[7] = false;
            labels3[8] = "Min                       ";                 states3[8] = false;
            labels3[9] = "Uniformity*               ";                 states3[9] = false;
            labels3[10] = "Max                      ";                 states3[10] = false;
            gd.setInsets(0, 0, 0);
            gd.addCheckboxGroup(7, 3, labels3, states3);

            gd.showDialog(); //show
            if (gd.wasCanceled())
                return false;

            boolean spam = false;

            //imageJ
            if(gd.getNextBoolean()){
                for (int i = 0 ; i <measureImageJ.length ; i++){
                    measureImageJ[i]=true;
                    spam = gd.getNextBoolean();
                }
            }else{
                for (int i = 0 ; i <measureImageJ.length ; i++){
                    measureImageJ[i]=gd.getNextBoolean();
                }
            }


            //BW
            if(gd.getNextBoolean()){
                for (int i = 0 ; i < measuresBW.length ; i++){
                    measuresBW[i]=true;
                    spam = gd.getNextBoolean();
                }
            }else{
                for (int i = 0 ; i < measuresBW.length ; i++){
                    measuresBW[i]=gd.getNextBoolean();
                }
            }


            //grey
            if(gd.getNextBoolean()){
                for (int i = 0 ; i < measuresGrey.length ; i++){
                    measuresGrey[i]=true;
                    spam = gd.getNextBoolean();
                }
            }else{
                for (int i = 0 ; i < measuresGrey.length ; i++){
                    measuresGrey[i]=gd.getNextBoolean();
                }
            }

            return true;
        }


        /**/
        private void setMeasuresExtended() {
            /*Misure by ImageJ*/
            if(measureImageJ[0]) {measure+=AREA;}
            if(measureImageJ[1]) {measure+=STD_DEV;}
            if(measureImageJ[2]) {measure+=CENTROID;}
            if(measureImageJ[3]) {measure+=CENTER_OF_MASS;}
            if(measureImageJ[4]) {measure+=PERIMETER;}
            if(measureImageJ[5]) {measure+=RECT;}
            if(measureImageJ[6]) {measure+=ELLIPSE;}
            if(measureImageJ[7]) {measure+=SHAPE_DESCRIPTORS;}
            if(measureImageJ[8]) {measure+=FERET;}
            if(measureImageJ[9]) {measure+=INTEGRATED_DENSITY;}
            if(measureImageJ[10]) {measure+=AREA_FRACTION;}
            if(measureImageJ[11]) {measure+=STACK_POSITION;}

            /*Misure BW*/
            doConvexArea = measuresBW[0];
            doAspRatio = measuresBW[1];
            doArEquivD = measuresBW[2];
            doConvexPerimeter = measuresBW[3];
            doPerEquivD = measuresBW[4];
            doMINRMAXR = measuresBW[5];
            doRoundness= measuresBW[6];
            doEquivEllAr = measuresBW[7];
            doCompactness = measuresBW[8];
            doSolidity = measuresBW[9];
            doShape = measuresBW[10];
            doConvexity = measuresBW[11];
            doRFactor = measuresBW[12];
            doArBBox = measuresBW[13];
            doConcavity = measuresBW[14];
            doRectang = measuresBW[15];
            doModRatio = measuresBW[16];
            doSphericity = measuresBW[17];
            doElongation = measuresBW[18];
            doNormPeriIndex = measuresBW[19];
            doHaralickRatio = measuresBW[20];
            doBendingEnergy = measuresBW[21];

            /*Misure Grey*/
            doMean = measuresGrey [0];
            doSkewness = measuresGrey [1]; if(doSkewness) {measure+=SKEWNESS;}
            doIntensitySum = measuresGrey [2];
            doMode = measuresGrey [3];
            doKurtois = measuresGrey [4]; if(doKurtois) {measure+=KURTOSIS;}
            doEntropy = measuresGrey [7];
            doMedian = measuresGrey [6]; if(doMedian) {measure+=MEDIAN;}
            doSquareIntensitySum = measuresGrey [5];
            doMin = measuresGrey [8];
            doUniformity= measuresGrey[9];
            doMax = measuresGrey [9];
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
                rt.addValue("*ConvexPerimeter", convexPerimeter);
            }

            if (doMINRMAXR) { //
                rt.addValue("*MinR", feret[2]/2); //DA RIVEDERE Radius of the inscribed circle centred at the middle of mass
                rt.addValue("*MaxR", feret[0]/2); //DA RIVEDERE Radius of the enclosing circle centred at the middle of mass
            }

            if (doAspRatio) //Aspect ratio = Feret/Breadth = L/W also called Feret ratio or Eccentricity or Rectangular ratio
                rt.addValue("*AspRatio", feret[0]/feret[2]);

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
            }

            /**/
            if(doMean){
                rt.addValue("**Mean",stats.mean);
            }

            /*if(doSkewness){
                rt.addValue("**Swewness",stats.skewness);
            }

            if(doKurtois){
                rt.addValue("**Kurtosis",stats.kurtosis);
            }*/

            if(doMode){
                rt.addValue("**Mode",stats.dmode);
            }

            /*
            if (doMedian){
                rt.addValue("**Median",stats.median);
            }*/

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

            if(doUniformity){
                rt.addValue("**Uniformity", getUniformity(hist, stats.area));
            }



            Analyzer.setMeasurements(old_measures); //per risettare misure vecchie
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

        private int getIntensitySum(int [] hist){
            int sum=0;
            for(int i=0; i< hist.length; i++){
                sum+=i*hist[i];
            }
            return sum;
        }

        private double getUniformity (int [] hist, double area){
            double uniformity=0.0;

            for (int i=0; i< hist.length; i++){
                uniformity += Math.sqrt((double) hist[i]/area);
            }

            return uniformity;

        }

        /*Funzione di appoggio per il caolcolo della distanza*/
        public double distance(int argx1, int argy1, int argx2, int argy2) {
            return Math.sqrt((argx1 - argx2) * (argx1 - argx2) + (argy1 - argy2) * (argy1 - argy2));
        }

        // log2:  Logarithm base 2
        public double log2(double d) {
            return Math.log(d)/Math.log(2.0);
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
