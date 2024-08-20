

/**
 * Created by IntelliJ IDEA.
 * User: BurtsevM
 * Date: 21.04.2007
 * Time: 5:53:42
 * To change this template use File | Settings | File Templates.
 */
public class PhenoLandscape {
    public int dimension;
    public int maxScore;
    public int size; // distance form origin to boundary for every direction [0,size]
                     // which gives [-size,size] for total lenghth of  one dimension
 //   public int[] score; // [0;maxScore]
    //public ArrayList score;

    PhenoLandscape(){

    }
    public class Point {   // Point on Landscape
        public int[] coord;
        public int score;
        Point(){
            coord = new int[dimension];
        }
    }
}
