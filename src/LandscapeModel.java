
import java.io.FileWriter;
import java.io.IOException;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.StringTokenizer;

/**
 * Created by IntelliJ IDEA.
 * User: BurtsevM
 * Date: 19.04.2007
 * Time: 3:36:38
 */
class LandscapeModel {
    public static int nProblems;
    public static int complexity;
    public static int populationSize;
    public static int nLearningTrials;
    public static int nGenes;
    public static int nAlleles;
    public static int version;
    public static int period;
    public static double mutationProbability;
    public static double[] avrgSumScr;
    public static double[] avrgMaxInnFS;
    public static double[] avrgMeanInnFS;
    public static double[] avrgMaxLnFS;
    public static double[][] avrgTermFSPEl;
    public static double[] nMut = {0.0001};
    public static int[] nG = {4000};
    public static int[] nA = {5};
    public static int[] nLT = {0, 5};
    public static boolean noCompensation = false;
    public static boolean sexualReproduction = false;
    public static int nlt = 10;
    public static int nVer;
    public static int vv;
    static Evolution world;

    public static void main(String[] args) {
        nProblems = 100;
        complexity = 20;
        populationSize = 500;
        nLearningTrials = 5;
        nGenes = 2000;
        nAlleles = 5;
        version = 0;
        period = 5000;
        mutationProbability = 0.005;
        nVer = 10;
        readParam();
        //if ((nVer-version)> 1) {
            avrgSumScr = new double[period];
            avrgMaxInnFS = new double[period];
            avrgMaxLnFS = new double[period];
            avrgMeanInnFS = new double[period];
            avrgTermFSPEl = new double[period][nlt];
            for (int i = 0; i < period; i++) {
                avrgSumScr[i] = 0;
                avrgMaxInnFS[i] = 0;
                avrgMaxLnFS[i] = 0;
                avrgMeanInnFS[i] = 0;
                for (int j = 0; j < nlt; j++) avrgTermFSPEl[i][j]=0;
            }
        //}
        /* int[] nG = {4000};
         int[] nA = {5};
         int[] nLT = {0,5};*/
        for (int m = 0; m < nMut.length; m++)
        for (int lt = 0; lt < nLT.length; lt++)
            for (int g = 0; g < nG.length; g++)
                for (int a = 0; a < nA.length; a++) {
                    nLearningTrials = nLT[lt];
                    nGenes = nG[g];
                    nAlleles = nA[a];
                    for (int v = version; v < nVer; v++) {
                        vv = v;
                        mutationProbability = nMut[m];
                        world = new Evolution(v);
                        world.evolve(period);
                    }
                    if ((nVer-version) > 1) writeAvrgFile();
                }
    }

    static void writeAvrgFile() {
        try {
            FileWriter outF = new FileWriter("pr" + nProblems + "cmp" + complexity + "ln" + nLearningTrials + "mut"
                    + mutationProbability + "g" + LandscapeModel.nGenes + "a" + LandscapeModel.nAlleles + ".aver" + nVer + ".txt");
            outF.write("t\tsumScr\tmaxInnSc\tmaxScr\taverInn\tNPhenElTermFS\r\n");
            for (int i = 0; i < period; i++) {
                outF.write(i + "\t" + (avrgSumScr[i] / nVer) + "\t"+ (avrgMaxInnFS[i] / nVer) +
                        "\t" + (avrgMaxLnFS[i] / nVer) + "\t" + (avrgMeanInnFS[i] / nVer));
                for (int s = 0; s < nlt; s++) outF.write("\t"+avrgTermFSPEl[s]);
                outF.write("\r\n");
            }
            outF.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        for (int i = 0; i < period; i++) {
            avrgSumScr[i] = 0;
            avrgMaxInnFS[i] = 0;
            avrgMaxLnFS[i] = 0;
            avrgMeanInnFS[i] = 0;
            for (int j = 0; j < nlt; j++) avrgTermFSPEl[i][j]=0;
        }
    }
    static void readParam() {
        try {
            BufferedReader r = null;
            r = new BufferedReader(new FileReader("param.txt"));
            StringTokenizer strTknr;
            strTknr = new StringTokenizer(r.readLine(), "\t");
            strTknr.nextToken();
            nProblems = new Integer(strTknr.nextToken()).intValue();
            strTknr = new StringTokenizer(r.readLine(), "\t");
            strTknr.nextToken();
            complexity = new Integer(strTknr.nextToken()).intValue();
            strTknr = new StringTokenizer(r.readLine(), "\t");
            strTknr.nextToken();
            populationSize = new Integer(strTknr.nextToken()).intValue();
            strTknr = new StringTokenizer(r.readLine(), "\t");
            strTknr.nextToken();
            version = new Integer(strTknr.nextToken()).intValue();
            strTknr = new StringTokenizer(r.readLine(), "\t");
            strTknr.nextToken();
            period = new Integer(strTknr.nextToken()).intValue();
            strTknr = new StringTokenizer(r.readLine(), "\t");
            int l = strTknr.countTokens() - 1;
            strTknr.nextToken();
            nMut = new double[l];
            //mutationProbability = new Double(strTknr.nextToken()).doubleValue();
            for (int i = 0; i < l; i++) nMut[i] = new Double(strTknr.nextToken()).doubleValue();
            strTknr = new StringTokenizer(r.readLine(), "\t");
            strTknr.nextToken();
            nVer = new Integer(strTknr.nextToken()).intValue();
            strTknr = new StringTokenizer(r.readLine(), "\t");
            l = strTknr.countTokens() - 1;
            strTknr.nextToken();
            nG = new int[l];
            for (int i = 0; i < l; i++) nG[i] = new Integer(strTknr.nextToken()).intValue();
            strTknr = new StringTokenizer(r.readLine(), "\t");
            l = strTknr.countTokens() - 1;
            strTknr.nextToken();
            nA = new int[l];
            for (int i = 0; i < l; i++) nA[i] = new Integer(strTknr.nextToken()).intValue();
            strTknr = new StringTokenizer(r.readLine(), "\t");
            l = strTknr.countTokens() - 1;
            strTknr.nextToken();
            nLT = new int[l];
            for (int i = 0; i < l; i++) nLT[i] = new Integer(strTknr.nextToken()).intValue();
            strTknr = new StringTokenizer(r.readLine(), "\t");
            strTknr.nextToken();
            noCompensation = new Boolean(strTknr.nextToken()).booleanValue();
            strTknr = new StringTokenizer(r.readLine(), "\t");
            strTknr.nextToken();
            sexualReproduction = new Boolean(strTknr.nextToken()).booleanValue();
            r.close();

        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }
}
