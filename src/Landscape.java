
/**
 * Created by IntelliJ IDEA.
 * User: BurtsevM
 * Date: 19.04.2007
 * Time: 3:40:26
 * To change this template use File | Settings | File Templates.
 */
public class Landscape {
    public int dimension;
    public int size; // distance form origin to boundary for every direction [0,size]
                     // which gives [-size,size] for total lenghth of  one dimension
    public double[][] score; // [-1;1]

}
