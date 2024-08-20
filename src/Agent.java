
import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: BurtsevM
 * Date: 19.04.2007
 * Time: 3:52:37
 * To change this template use File | Settings | File Templates.
 */
public class Agent {
    boolean[][] phenotype;// vector of coordinates at PhenoLandscape
    boolean[] FS; // presence of functional system in phenotype
    int[] genotype;// vector of coordinates at GeneSpace
    public int nLearningTrials;
    int nPr;
    int cmpl;
    int nG;
    int nA;
    public int curNFS;
    public int innNFS;
    public int curNFSEl;
    public double score;

    Agent(int nPrblms, int cmplxt, int nGenes, int nAlleles, Random rand) {
        nPr = nPrblms;
        cmpl = cmplxt;
        nG = nGenes;
        nA = nAlleles;
        innNFS = 0;
        phenotype = new boolean[nPr][cmpl];
        FS = new boolean[nPr];
        resetPhenotype();
        genotype = new int[nG];
        for (int i = 0; i < nG; i++)
            genotype[i] = rand.nextInt(nA);
    }

    Agent(Agent parent, double mutProb, Random rand) {
        nPr = parent.nPr;
        cmpl = parent.cmpl;
        phenotype = new boolean[nPr][cmpl];
        FS = new boolean[nPr];
        nG = parent.nG;
        nA = parent.nA;
        innNFS = parent.innNFS;
        resetPhenotype();
        genotype = new int[nG];
        for (int i = 0; i < nG; i++)
            if (mutProb > rand.nextDouble())
                genotype[i] = rand.nextInt(nA);
            else
                genotype[i] = parent.genotype[i];
    }

    Agent(Agent parent1, Agent parent2, double mutProb, Random rand) {
        nPr = LandscapeModel.nProblems;
        cmpl = LandscapeModel.complexity;
        phenotype = new boolean[nPr][cmpl];
        FS = new boolean[nPr];
        nG = LandscapeModel.nGenes;
        nA = LandscapeModel.nAlleles;
        if (rand.nextDouble() < 0.5)
            innNFS = parent1.innNFS;
        else
            innNFS = parent2.innNFS;
        resetPhenotype();
        genotype = new int[nG];
        //int recPoint = rand.nextInt(LandscapeModel.nGenes);
        for (int i = 0; i < nG; i++)
            if (mutProb > rand.nextDouble())
                genotype[i] = rand.nextInt(nA);
           /* else { //point recombination
                if (i < recPoint)
                    genotype[i] = parent1.genotype[i];
                else
                    genotype[i] = parent2.genotype[i];
            }*/
           else { //uniform recombination
                if (rand.nextDouble() < 0.5)
                    genotype[i] = parent1.genotype[i];
                else
                    genotype[i] = parent2.genotype[i];
            }
    }

    void resetPhenotype() {
        for (int i = 0; i < nPr; i++) {
            for (int j = 0; j < cmpl; j++)
                phenotype[i][j] = false;
            FS[i] = false;
        }
    }
}
