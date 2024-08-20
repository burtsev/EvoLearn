
import java.util.Random;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: BurtsevM
 * Date: 21.04.2007
 * Time: 5:59:00
 */
public class GeneSpace {
    public int nAlleles;
    public int nGenes;
    public int[][][] gene;
    int cmp;
    int pr;
    boolean holes = true;

    GeneSpace(int nPr, int cmplx, Random rand) {
        nGenes = LandscapeModel.nGenes;
        nAlleles = LandscapeModel.nAlleles;
        cmp = cmplx;
        pr = nPr;
        gene = new int[nGenes][nAlleles][2];
        while (holes) {
            for (int i = 0; i < nGenes; i++)
                for (int j = 0; j < nAlleles; j++) {
                    gene[i][j][0] = rand.nextInt(nPr);
                    gene[i][j][1] = rand.nextInt(cmplx);
                }
            holes = checkHoles();
        }
    }

    public void saveG2PMap(String fileName) throws IOException {
        FileWriter g2pm = new FileWriter(fileName);
        for (int i = 0; i < nGenes; i++) {
            for (int j = 0; j < nAlleles; j++) g2pm.write((gene[i][j][0] * cmp + gene[i][j][1]) + "\t");
            g2pm.write("\r\n");
        }
        g2pm.close();
    }

    boolean checkHoles() {
        int[][] t = new int[pr][cmp];
        boolean hole = false;
        for (int i = 0; i < pr; i++)
            for (int j = 0; j < cmp; j++) t[i][j] = 0;
        for (int i = 0; i < nGenes; i++)
            for (int j = 0; j < nAlleles; j++) t[gene[i][j][0]][gene[i][j][1]]++;
        for (int i = 0; i < pr; i++)
            for (int j = 0; j < cmp; j++) if (t[i][j] == 0) hole = true;
        return hole;
    }
}
