
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.FileWriter;
import java.util.StringTokenizer;

class Cleaner {
    static int step = 200;

    public static void main(String[] args) {
        int time = 0;
        try {

            BufferedReader r = null;
            r = new BufferedReader(new FileReader("pr100cmp20ln5mut2.0E-4g4000a5.aver10.txt"));
            StringTokenizer strTknr;
            String s = new String();
            s = r.readLine();
            s = r.readLine();
            FileWriter outF = new FileWriter("clean.txt");
            while (s != null) {
                strTknr = new StringTokenizer(s, "\t");
                time = new Integer(strTknr.nextToken()).intValue();
                if ((time % step) == 0) outF.write(s + "\r\n");
                s = r.readLine();
            }
            r.close();
            outF.close();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }
}
