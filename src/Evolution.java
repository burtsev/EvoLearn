
import java.util.Random;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: BurtsevM
 * Date: 19.04.2007
 * Time: 3:55:23
 * To change this template use File | Settings | File Templates.
 */
public class Evolution {
    //public int mutationRange;
    //public int plasticityRange;
    //public int learningRange;
    //public int lethalTrshld;
    public double mutProb;
    public int popSize;
    public int nLearningTrials;
    int nPrblms;
    int cmplxt;
    GeneSpace gSpace;
    Agent[] pop;
    Agent[] newPop;
    Random rand;
    int nInnateFSEl; // number of innate FS elements
    int sumScore;
    int maxSc;
    int minSc;
    int maxInnateSc;
    int sd;
    int gtime;
    double[] averTermFSpEl;
    double averNFSEl;
    double averInnFS;
    int curFSEdistr[];
    int curInnFSdistr[];
    int curFSLrndDistr[];
    int curGeneDistr[][];
    FileWriter outFSEdistr;
    FileWriter outInnFSdistr;
    FileWriter outFSLrndDistr;
    double lastFSHoles;

    Evolution(int seed) {
        rand = new Random(seed);
        sd = seed;
        nPrblms = LandscapeModel.nProblems;
        cmplxt = LandscapeModel.complexity;
        popSize = LandscapeModel.populationSize;
        nLearningTrials = LandscapeModel.nLearningTrials;
        mutProb = LandscapeModel.mutationProbability;
        gSpace = new GeneSpace(nPrblms, cmplxt, rand);
        pop = new Agent[popSize];
        newPop = new Agent[popSize];
        for (int i = 0; i < popSize; i++) pop[i] = new Agent(nPrblms, cmplxt, gSpace.nGenes, gSpace.nAlleles, rand);
        System.out.println("Population initialized...");
        getPhenFill();
    }

    void evolve(int t) {
        try {
            FileWriter outF = new FileWriter("pr" + nPrblms + "cmp" + cmplxt + "ln" + nLearningTrials + "mut"
                    + mutProb + "g" + LandscapeModel.nGenes + "a" +
                    LandscapeModel.nAlleles + ".v" + LandscapeModel.vv + ".txt");
            outF.write("t\tsumScr\tmaxInnSc\tmaxScr\taFSEl\taInnFS\tNPhenElTermFS\r\n");
            openDistrFiles();
            /*FileWriter holeStat = new FileWriter("pr" + nPrblms + "cmp" + cmplxt + "ln" + nLearningTrials + "mut"
                    + mutProb + "g" + LandscapeModel.nGenes + "a" +
                    LandscapeModel.nAlleles + ".v" + LandscapeModel.vv + ".hs.txt");*/
            gSpace.saveG2PMap("pr" + nPrblms + "cmp" + cmplxt + "ln" + nLearningTrials + "mut"
                    + mutProb + "g" + LandscapeModel.nGenes + "a" +
                    LandscapeModel.nAlleles + ".v" + LandscapeModel.vv + ".g2pm.txt");
            for (int time = 0; time < t; time++) {
                gtime = time;
                setPenotypes();
                learn();
                reproduce(LandscapeModel.sexualReproduction);
                System.out.print("v:" + sd + " lt:" + nLearningTrials + " nG:" + LandscapeModel.nGenes +
                        " nA:" + LandscapeModel.nAlleles +
                        " t:" + time + " sScr:" + sumScore + " maxInnSc:" + maxInnateSc + " maxScr:" + maxSc + "  ");
                for (int s = 0; s < LandscapeModel.nlt; s++) System.out.print(" " + averTermFSpEl[s]);
                System.out.println();
                if (LandscapeModel.nVer > 1) {
                    LandscapeModel.avrgSumScr[time] += sumScore;
                    LandscapeModel.avrgMaxInnFS[time] += maxInnateSc;
                    LandscapeModel.avrgMaxLnFS[time] += maxSc;
                    LandscapeModel.avrgMeanInnFS[time] += averInnFS;
                    for (int s = 0; s < LandscapeModel.nlt; s++)
                        LandscapeModel.avrgTermFSPEl[time][s] += averTermFSpEl[s];
                }
                outF.write(time + "\t" + sumScore + "\t" + maxInnateSc + "\t" + maxSc +
                        "\t" + averNFSEl + "\t" + averInnFS);
                for (int s = 0; s < LandscapeModel.nlt; s++) outF.write("\t" + averTermFSpEl[s]);
                outF.write("\r\n");
                writeDistrs();
                //holeStat.write(time + "\t" + lastFSHoles + "\r\n");
                //if (time == 15000) nLearningTrials = 5;
            }
            outF.close();
            //holeStat.close();
            closeDistrFiles();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    void setPenotypes() {
        maxInnateSc = 0;
        averNFSEl = 0;
        averInnFS = 0;
        lastFSHoles = 0;
        averTermFSpEl = new double[LandscapeModel.nlt];
        curFSEdistr = new int[nPrblms * cmplxt];
        curInnFSdistr = new int[nPrblms];
        curGeneDistr = new int[gSpace.nGenes][gSpace.nAlleles];
        for (int i = 0; i < popSize; i++) {
            pop[i].curNFSEl = 0;
            for (int j = 0; j < gSpace.nGenes; j++) {
                pop[i].phenotype[gSpace.gene[j][pop[i].genotype[j]][0]][gSpace.gene[j][pop[i].genotype[j]][1]] = true;
                curGeneDistr[j][pop[i].genotype[j]]++;
            }
            for (int j = 0; j < nPrblms; j++) {
                pop[i].FS[j] = true;
                for (int k = 0; k < cmplxt; k++)
                    if (!pop[i].phenotype[j][k]) {
                        pop[i].FS[j] = false;
                        pop[i].curNFSEl += j * cmplxt;
                        for (int l = 0; l < cmplxt; l++)
                            if (pop[i].phenotype[j][l]) {
                                pop[i].curNFSEl++;
                                averTermFSpEl[0]++;
                            }
                        for (int s = 1; s < LandscapeModel.nlt; s++)
                            if ((j + s) < nPrblms)
                                for (int l = 0; l < cmplxt; l++)
                                    if (pop[i].phenotype[j + s][l]) averTermFSpEl[s]++;
                        averNFSEl += pop[i].curNFSEl;
                        break;
                    }
                if (!pop[i].FS[j]) {
                    pop[i].curNFS = j;
                    curInnFSdistr[j]++;
                    averInnFS += j;
                    break;
                }
            }
            if (pop[i].innNFS < pop[i].curNFS) pop[i].innNFS = pop[i].curNFS;
            if (pop[i].curNFS > maxInnateSc) maxInnateSc = pop[i].curNFS;
            // start:lastFSHoles
            for (int e = 0; e < cmplxt; e++)
                if (!pop[i].phenotype[nPrblms - 1][e]) lastFSHoles++;
            // end:lastFSHoles
        }
        averNFSEl = averNFSEl / popSize / cmplxt;
        averInnFS = averInnFS / popSize;
        lastFSHoles = lastFSHoles / popSize;
        for (int s = 0; s < LandscapeModel.nlt; s++) averTermFSpEl[s] = averTermFSpEl[s] / popSize;
        getPhenElDistr();
    }

    void learn() {
        sumScore = 0;
        maxSc = 0;
        minSc = nPrblms;
        curFSLrndDistr = new int[nPrblms];
        for (int i = 0; i < popSize; i++) {
            if ((LandscapeModel.noCompensation && (pop[i].curNFS == pop[i].innNFS)) ||
                    (!LandscapeModel.noCompensation)) {
                for (int l = 0; l < nLearningTrials; l++) {
                    if (pop[i].curNFS < (nPrblms - 1)) {
                        nInnateFSEl = 0;
                        for (int j = 0; j < cmplxt; j++)
                            if (pop[i].phenotype[pop[i].curNFS][j]) nInnateFSEl++;
                        if (nInnateFSEl == cmplxt) {
                            pop[i].curNFS++;
                            l--;
                        } else if (Math.pow(0.5, (double) (cmplxt - nInnateFSEl)) > rand.nextDouble())
                            pop[i].curNFS++;
                    }
                }
            }
            sumScore += pop[i].curNFS;
            curFSLrndDistr[pop[i].curNFS]++;
            if (pop[i].curNFS > maxSc) maxSc = pop[i].curNFS;
            if (pop[i].curNFS < minSc) minSc = pop[i].curNFS;
        }
    }

    void reproduce(boolean sex) {
        double selC = 1.1; //selection coefficient for the uniform selection
        double pMaxSc = Math.pow(selC, (maxSc - 1)); // normalization coeff for the uniform selection
        if ((maxSc - minSc) > 0)
            for (int i = 0; i < popSize; i++) pop[i].score = Math.pow(selC, (pop[i].curNFS - 1)) / pMaxSc;
        //(double) pop[i].curNFS / maxSc;   // additive selection (decreasing)
        //(double) (pop[i].curNFS - minSc) / (maxSc - minSc); // / sumScore;
        else
            for (int i = 0; i < popSize; i++) pop[i].score = 1;
        if (sex) sexualRep(); else asexualRep();
    }

    void asexualRep() {
        int n = 0;
        if (sumScore > 0) {
            while (n < popSize) {
                int i = rand.nextInt(popSize);
                if (pop[i].score > rand.nextDouble()) {
                    newPop[n] = new Agent(pop[i], mutProb, rand);
                    n++;
                }
            }
        } else {
            for (int i = 0; i < popSize; i++)
                newPop[i] = new Agent(pop[i], mutProb, rand);
            n = popSize;
        }
        for (int i = 0; i < popSize; i++) pop[i] = newPop[i];
    }

    void sexualRep() {
        int n = 0;
        while (n < popSize) {
            if (sumScore > 0) {
                for (int i = 0; i < popSize; i++)
                    if (pop[i].score > rand.nextDouble()) {
                        newPop[n] = new Agent(pop[i], pop[getParent()], mutProb, rand);
                        n++;
                        if (n == popSize) break;
                    }
            } else {
                for (int i = 0; i < popSize; i++)
                    newPop[i] = new Agent(pop[i], mutProb, rand);
                n = popSize;
            }
        }
        for (int i = 0; i < popSize; i++) pop[i] = newPop[i];
    }

    int getParent() {
        boolean parentFound = false;
        int ind = 0;
        while (!parentFound) {
            ind = rand.nextInt(popSize);
            if (pop[ind].score > rand.nextDouble()) {
                parentFound = true;
                break;
            }
        }
        return ind;
    }

    double getPhenFill() {
        double phenFill = 0;
        double sum;
        int maxEls = 0;
        double averEls = 0;
        int curEls;
        boolean noHole = true;
        setPenotypes();
        for (int i = 0; i < popSize; i++) {
            noHole = true;
            sum = 0;
            curEls = 0;
            for (int j = 0; j < nPrblms; j++)
                for (int k = 0; k < cmplxt; k++) {
                    if (pop[i].phenotype[j][k]) {
                        sum++;
                        if (noHole) curEls++;
                    } else
                        noHole = false;
                }
            if (curEls > maxEls) maxEls = curEls;
            averEls += curEls;
            sum = sum / (nPrblms * cmplxt);
            phenFill += sum;
        }
        phenFill = phenFill / popSize;
        averEls = averEls / popSize;
        System.out.println("filling of phen:" + phenFill + " aFSEls:" + averEls + " maxFSEls:" + maxEls);
        return phenFill;
    }

    void openDistrFiles() throws IOException {
        outFSEdistr = new FileWriter("pr" + nPrblms + "cmp" + cmplxt + "ln" + nLearningTrials + "mut"
                + mutProb + "g" + LandscapeModel.nGenes + "a" +
                LandscapeModel.nAlleles + ".v" + LandscapeModel.vv + ".FSEdistr" + ".txt");
        outInnFSdistr = new FileWriter("pr" + nPrblms + "cmp" + cmplxt + "ln" + nLearningTrials + "mut"
                + mutProb + "g" + LandscapeModel.nGenes + "a" +
                LandscapeModel.nAlleles + ".v" + LandscapeModel.vv + ".InnFSdistr" + ".txt");
        outFSLrndDistr = new FileWriter("pr" + nPrblms + "cmp" + cmplxt + "ln" + nLearningTrials + "mut"
                + mutProb + "g" + LandscapeModel.nGenes + "a" +
                LandscapeModel.nAlleles + ".v" + LandscapeModel.vv + ".FSLrndDistr" + ".txt");
    }

    void closeDistrFiles() throws IOException {
        outFSEdistr.close();
        outInnFSdistr.close();
        outFSLrndDistr.close();
    }

    void writeDistrs() throws IOException {
        for (int i = 0; i < curInnFSdistr.length; i++) {
            outInnFSdistr.write(curInnFSdistr[i] + "\t");           
        }
        outInnFSdistr.write("\r\n");
		for (int i = 0; i < curFSLrndDistr.length; i++) {
           outFSLrndDistr.write(curFSLrndDistr[i] + "\t");           
        }
        outFSLrndDistr.write("\r\n");
        // for (int i = 0; i < 200; i++) outFSEdistr.write(curFSEdistr[i] + "\t");
        //outFSEdistr.write("\r\n");
        // if ((gtime % 5000) == 0) writeGDistr();
    }

    void writeGDistr() throws IOException {
        FileWriter file = new FileWriter("pr" + nPrblms + "cmp" + cmplxt + "ln" + nLearningTrials + "mut"
                + mutProb + "g" + LandscapeModel.nGenes + "a" +
                LandscapeModel.nAlleles + ".v" + LandscapeModel.vv + ".gDistr." + gtime + ".txt");
        for (int i = 0; i < gSpace.nGenes; i++) {
            for (int j = 0; j < gSpace.nAlleles; j++) file.write(curGeneDistr[i][j] + "\t");
            file.write("\r\n");
        }
        file.close();
    }

    void getPhenElDistr() {
        for (int i = 0; i < popSize; i++)
            for (int s = 0; s < nPrblms; s++)
                for (int e = 0; e < cmplxt; e++)
                    if (pop[i].phenotype[s][e]) curFSEdistr[s * cmplxt + e]++;
    }
}
