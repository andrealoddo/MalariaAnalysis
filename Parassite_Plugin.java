import ij.*;
import ij.measure.*;
import ij.plugin.*;
import ij.plugin.filter.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;

/*PlugIn di prova per poter aggiungere funzionalità ad Partiche Analyzer.
* la classe Parassite_Plugin implementa l'interfaccia PlugIn*/


public class Parassite_Plugin implements PlugIn {

    public void run(String arg) {
        //Restituisce un riferimento all'immagine attiva o visualizza un messaggio
        // di errore e interrompe il plugin o la macro se non ci sono immagini aperte
        ImagePlus imp = IJ.getImage();
        analyzeParticles(imp); //metodo presente dopo
    }

    public void analyzeParticles(ImagePlus imp) {
        //crea un nuovo oggetto Parassite
        //crea flags, che prende setup
        Parassite pa = new Parassite();
        int flags = pa.setup("", imp);
        if (flags == PlugInFilter.DONE)
            return; //se è tutto ok runna
        pa.run(imp.getProcessor());
        Analyzer.getResultsTable().show("Results");
    }
    //cuore

    class Parassite extends ParticleAnalyzer {
        @Override
        protected void saveResults(ImageStatistics stats, Roi roi) {
            double pigreco = 3.14;
            super.saveResults(stats, roi);
            ImageProcessor ip = imp.getProcessor();
            ip.setRoi(roi);
            float[] areas = rt.getColumn(ResultsTable.AREA);
            float[] ferets = rt.getColumn(ResultsTable.FERET);
            float[] perims = rt.getColumn(ResultsTable.PERIMETER);
            for (int i = 0; i < areas.length; i++) {
                rt.addValue("Roundness", (areas[i] * 4) / ((pigreco) * (ferets[i] * ferets[i])));
                rt.addValue("Circ", 4 * areas[i] / (perims[i] * perims[i]));
                rt.addValue("ArEquivD", Math.sqrt((4 * pigreco) * areas[i]));
            }


        }

    }
}
