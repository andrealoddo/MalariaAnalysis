import ij.*;
import ij.measure.*;
import ij.plugin.*;
import ij.plugin.filter.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.filter.Analyzer;


/*PlugIn di prova per poter aggiungere funzionalit� ad Partiche Analyzer.
 * la classe Parassite_Plugin implementa l'interfaccia PlugIn*/


public class Parassite_Plugin implements PlugIn {
    /*booleani per la checkbox in cui l'utente � invitato a selezionare cosa vuole calcolare .. da rivedere*/
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
    boolean doArBBox;
    boolean doRectang;
    boolean doConvexArea;

    boolean doIS;
    boolean doMidpoint;
    boolean doDS;

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
        genericDialog();
        Parassite pa = new Parassite();

        int flags = pa.setup("", imp); //runna
        if (flags == PlugInFilter.DONE)
            return; //se � tutto ok runna
        pa.run(imp.getProcessor());
        Analyzer.getResultsTable().show("Results");

    }

    public void genericDialog(){
        GenericDialog gd = new GenericDialog("Parassite Prova", IJ.getInstance());
        gd.addStringField("Title: ", "PROVA");
        gd.addMessage("Choise measures");

        gd.addCheckbox("SelectAll",false);

        gd.addCheckbox("AspRatio",false);
        gd.addCheckbox("Circ", false);
        gd.addCheckbox("Roundness", false);
        gd.addCheckbox("ArEquivD", false);
        gd.addCheckbox("PerEquivD", false);
        gd.addCheckbox("EquivEllAr", false);
        gd.addCheckbox("Compactness", false);
        gd.addCheckbox("Shape", false);
        gd.addCheckbox("ArBBox", false);
        gd.addCheckbox("Rectang", false);

        gd.addCheckbox("Convex Area", false);


        gd.showDialog(); //show

        //if (gd.wasCanceled())
        //return DONE;

        boolean doSelectAll=gd.getNextBoolean ();  //se viene selezionata la voce "seleziona tutti"

        if(doSelectAll){ // vengono settati i booleani a true
            doAspRatio = true;
            doCirc = true;
            doRoundness = true;
            doArEquivD = true;
            doPerEquivD = true;
            doEquivEllAr = true ;
            doCompactness = true;
            doShape = true;
            doArBBox =true;
            doRectang = true;
            doConvexArea = true;
        }else{ //altrimenti pescati uno per uno con un certo ordine
            doAspRatio=gd.getNextBoolean();
            doCirc=gd.getNextBoolean();
            doRoundness=gd.getNextBoolean();
            doArEquivD=gd.getNextBoolean();
            doPerEquivD=gd.getNextBoolean();
            doEquivEllAr=gd.getNextBoolean();
            doCompactness=gd.getNextBoolean();
            doShape=gd.getNextBoolean();
            doArBBox =gd.getNextBoolean();
            doRectang = gd.getNextBoolean();
            doConvexArea = gd.getNextBoolean();
        }


    }


    /*Mi servirà poi
    int ms = 0;
        ms |= ana.AREA+ana.PERIMETER+ana.ELLIPSE+ana.CENTER_OF_MASS;
        ana.setMeasurements(ms);
        ImageStatistics stats = img.getStatistics(ms);*/
    public void convexHull(){
        float perimeterOfHull;
        float AREAofHull = 0;
        float major, minor;

    }

    //da Analyzer
    final double getArea(Polygon p) {
        if (p==null) return Double.NaN;
        int carea = 0;
        int iminus1;
        for (int i=0; i<p.npoints; i++) {
            iminus1 = i-1;
            if (iminus1<0) iminus1=p.npoints-1;
            carea += (p.xpoints[i]+p.xpoints[iminus1])*(p.ypoints[i]-p.ypoints[iminus1]);
        }
        return (Math.abs(carea/2.0));
    }
    //



    //cuore

    /*� una classe interna dichiarata nella classe Parassite_Plugin*/

    class Parassite extends ParticleAnalyzer {
        @Override
        protected void saveResults(ImageStatistics stats, Roi roi) {
            super.saveResults(stats, roi);
            ImageProcessor ip = imp.getProcessor();
            ip.setRoi(roi);
            float[] areas = rt.getColumn(ResultsTable.AREA);
            float[] ferets = rt.getColumn(ResultsTable.FERET);
            float [] breadths =  rt.getColumn(ResultsTable.MIN_FERET);
            float[] perims = rt.getColumn(ResultsTable.PERIMETER);
            for (int i = 0; i < areas.length; i++) {

                if(doAspRatio)
                    rt.addValue("aspRatio", ferets[i]/breadths[i]);
                if (doCirc)
                    rt.addValue("Circ", 4 * areas[i] / (perims[i] * perims[i]));
                if (doRoundness)
                    rt.addValue("Roundness", (areas[i] * 4) / ((pigreco) * (ferets[i] * ferets[i])));
                if(doArEquivD)
                    rt.addValue("ArEquivD", Math.sqrt((4 * pigreco) * areas[i]));
                if(doPerEquivD)
                    rt.addValue("PerEquivD", areas[i]/pigreco);
                if(doEquivEllAr)
                    rt.addValue("EquivEllAr", (pigreco*ferets[i]*breadths[i])/4);
                if(doCompactness)
                    rt.addValue("Compactness", (Math.sqrt((4 * pigreco) * areas[i]))/ferets[i]);
                if(doShape)
                    rt.addValue("Shape", (perims[i]*perims[i])/areas[i]);
                if(doArBBox)
                    rt.addValue("ArBBox", ferets[i]*breadths[i]);
                if(doRectang)
                    rt.addValue("Rectang", areas[i]/(ferets[i]*breadths[i]));

                if(doConvexArea){
                    double convexArea = roi!=null?getArea(roi.getConvexHull()):stats.pixelCount;
                    rt.addValue("ConvexArea", convexArea);
                }

               // Polygon convexRoi = roi.getConvexHull();



            }


        }

    }
}

