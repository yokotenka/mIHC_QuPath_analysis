/*
 * ################################### User Input #######################################
 * In the future, this section will become a GUI so that it is easier for the user to use.
 *
 * Before you run this code, you must:
 *      - successfully import the images into QuPath (where all the images are colour channels)
 *      - rename the colour channels
 *      - run stardist on the region of interest
 *
 * Output:
 *      - Currently the output will be printed out onto the log. The threshold is given as a number in the range 0-255
 *        if you would like the normalise value i.e. between 0-1, simply divide the threshold by 255.
 *      - you can also choose to save the kernel density estimations as a .csv file. This way you can visualise the
 *        density which was fitted using tools like excel, python, matlab etc.
 */


/* Marker List
 * In theory we would have all the markers but we limit to two for testing purposes
 */
markers = [
            "CD45",
            "PanCK"
            ]


/* Baseline markers
 * These are the markers which shouldn't appear in a tumour cell. Since CD45 is present only in immune cells, we use
 * this marker as a baseline marker. Of course we can add in more markers to the list but we haven't really understood
 * the advantage of doing this.
 */
def baseLineMarkers = [
        "CD45"
]

// Tumour marker
def tumourMarker = "PanCK"


// Measurement being used
def measurement = ": Cell: Mean"

// Check tumour distribution
// Saves a .csv of the distribution array if true.
def saveDensities = true
def filePath = "/Users/yokote.k/Documents/densities.csv"


/* #####################################################################################################################
 * ######################################## User should not go beyond this point #######################################
 * ###################################################################################################################*/


/*
 * ########################################### Range class #############################################################
 */

/** Represents a range for the data. For an 8-bit image this will be between 0 and 255
 * @author Kenta Yokote
 */
class Range {

    /** The minimum of the range*/
    private Double min;

    /** The maximum of the range*/
    private Double max;

    /** Constructor for the range class
     * @param min Minimum of the range
     * @param max Maximum of the range
     */
    Range(Double min, Double max){
        if (min.compareTo(max) >= 0 ){
            // Do error handling for when min is greater than equal to max
        }
        this.min = min;
        this.max= max;
    }

    /** Getter for the maximum value of the range
     * @return max
     */
    double getMax() {
        return max;
    }

    /** Getter for the minimum value of the range
     * @return min
     */
    double getMin() {
        return min;
    }

    /** Getter for the difference between the min and max values
     * @return max - min
     */
    double getDifference(){
        return max - min;
    }
}

/*
 * ########################################### Normal Distribution Class ###############################################
 */
/** Represtents the normal distribution. Reason for creating own class is because
 * org.apache.commons.math3.distribution.NormalDistribution did not have a method which returns an array of densities.
 * It was faster making a new class which does this.
 * @author Kenta
 */

import org.apache.commons.math3.util.FastMath;

class NormalDistribution {
    /** The mean of the normal distribution */
    private double mean = 0;
    /** The standard deviation of the normal distribution */
    private double standardDeviation = 1;


    /** Constructor
     * @param mean
     * @param standardDeviation
     */
    NormalDistribution(double mean, double standardDeviation){
        this.mean = mean;
        this.standardDeviation = standardDeviation;
    }

    /** Constructor default is the standard normal distribution
     */
    NormalDistribution(){
    }

    /** Constructor only with the standard deviation
     * @param standardDeviation
     */
    NormalDistribution(double standardDeviation){
        this.standardDeviation = standardDeviation;
    }

    /** Getter for the density for a specified range
     @param range The range of interest
     @param n The number of data points
     */
    Double[] getDensityOfRange(Range range, int n){
        Double[] densityArray = new Double[n];

        Double increment = range.getDifference() / n;
        Double currentPoint = (double) range.getMin();

        for (int i=0; i<n; i++){
            densityArray[i] = density(currentPoint);
            currentPoint += increment;
        }
        return densityArray;
    }

    /** Calculates the density at a specific point
     @param x The point at which the density function is being evaluated at.
     */
    Double density(double x){
        double x0 = (x - mean) / standardDeviation;
        double exponent = -0.5 * x0 * x0;
        return FastMath.exp(exponent) / (standardDeviation * FastMath.sqrt(2 * Math.PI));
    }

    /** Setter for the mean
     * @param mean The new mean for the normal distribution.
     */
    void setMean(double mean) {
        this.mean = mean;
    }

    /** Setter for the standard deviation
     * @param standardDeviation
     */
    void setStandardDeviation(double standardDeviation){
        this.standardDeviation = standardDeviation;
    }
}


/*
 * ######################################### Kernel Density Estimation #################################################
 */
/** A class for Kernel Density Estimation. There are no KDE packages that come pre-installed with QuPath hence was
 * easier to write one which didn't require installing another package.
 *
 * Have not yet implemented solutions for the bounded bias.
 *
 * @author Kenta
 */
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics
class KernelDensityEstimation {
    /** The standard deviation which will be used for individual normal distributions */
    private double bandWidth;

    /** The kernel function name. Defaults to the Gaussian kernel. In reality whatever kernel works the same. */
    private String kernelName ="gaussian";

    /** Weighting for each data point. Not implemented yet. */
    private Double[] weights = null;

    /** The range of the data points */
    private Range range = null;

    /** The estimated kernel density array */
    private Double[] estimate;

    /** The number of data points to be included in the estimated array*/
    private int n = 2048;

    /** The data points which the estimation is performed on */
    private double[] x;


    /** Constructor
     * @param x The data points
     * @param min Minimum in the range.
     * @param max Maximum in the range.
     */
    KernelDensityEstimation(double[] x, double min, double max){
        this.x = x;
        setBandWidthMethod("Default");
        this.range = new Range(min, max);
    }

    /** Estimation
     * @return kernelDensity The kernel density estimation for the given data, x.
     */
    Double[] estimate(){
        Double[] kernelDensity = new Double[n];
        Double[] currentDensity;
        boolean isFirstElement = true;

        // Initialise kernel distribution.
        NormalDistribution kernel = new NormalDistribution(bandWidth);

        // Iterate through every datapoint
        for (Double point : x){

            // Set the mean value for the kernel distribution.
            kernel.setMean(point);

            // Get the distribution of the kernel for the given range
            currentDensity = kernel.getDensityOfRange(range, n);

            // Calculate the kernel density
            if (isFirstElement){
                kernelDensity = currentDensity;
                isFirstElement = false;
            } else{
                for (int i=0; i<kernelDensity.length; i++){
                    kernelDensity[i] = kernelDensity[i] + currentDensity[i];
                }
            }
        }

        // Normalise
        for (int i=0; i<kernelDensity.length; i++) kernelDensity[i] = kernelDensity[i] / x.length;

        return kernelDensity;
    }

    /** Private method which will calculate the bandwidth to be used. Probably should create another class for this.
     * */
    private void setBandWidthMethod(String bandWidthName) throws UnsupportedOperationException {

        if (bandWidthName == "Silverman"){
            throw new UnsupportedOperationException("Silverman not yet supported");
        } else if (bandWidthName == "SJ"){
            throw new UnsupportedOperationException("Sheather & Jones not yet supported");
        } else {
            // Can change this if need be.
            DescriptiveStatistics da = new DescriptiveStatistics(x);
            double xStandardDeviation = da.getStandardDeviation();
            double iqr = da.getPercentile(75) - da.getPercentile(25);
            this.bandWidth = 0.9 * Math.min(xStandardDeviation, iqr/1.34) * Math.pow(x.length, -1/5);
        }
    }

    /** Getter for the bandwidth being used.
     * @return bandWidth
     */
    double getBandWidth(){
        return this.bandWidth;
    }
}

/*
 * ############################################# Utility ###############################################################
 */

/**
 * Function to get the local minima of a given array
 * @param arr The array which the local minima will be found
 * @return ArrayList of all the local minima
 */
def findLocalMinimaIndex(Double[] arr) {

    // Empty vector to store points of
    // local maxima and minima
    ArrayList<Double> mn = new ArrayList<Double>();
    int n = arr.length

    // Iterating over all points to check
    // local maxima and local minima
    for(int i = 1; i < arr.length - 1; i++) {
        // Condition for local minima
        if ((arr[i - 1].compareTo(arr[i])) > 0 &&
                (arr[i].compareTo(arr[i + 1]) < 0)) {
            mn.add(i);
        }
    }

    // Checking whether the last point is
    // local maxima or minima or none
    if (arr[n - 1].compareTo(arr[n - 2]) < 0){
        mn.add(n - 1);
    }
    return mn;
}





/*
 * ############################################ Qupath Scripting #######################################################
 */

// Import
import static qupath.lib.gui.scripting.QPEx.*



// Threshold map
def thresholdMap = new HashMap<String, Double>()


// Collect all the data
def server = getCurrentServer()
def imageData = getCurrentImageData()
def cells = getCellObjects()


def n = 2048
def densityArrayList = new ArrayList<ArrayList<Double>>(server.nChannels())
def csvColumnHeaders = []
csvColumnHeaders.add("Intensity")

// x-axis
def increment = 255 / n
def val = 0;
def x = new double[n];
for (int i=0; i < n; i++){
    x[i] = val
    val += increment
}

densityArrayList.add(new ArrayList<>(Arrays.asList(x)))

/*
 * ############################################ Checks #######################################################
 */
println("Checking parameters")
assert !cells.isEmpty() : "No cells detected"
assert !markers.isEmpty() : "No markers entered"
assert !baseLineMarkers.isEmpty() : "No baseline markers entered"

// ############################################### For regular markers #################################################
println("Calculating for non tumour markers")

for (markerName : markers){
    if (!markerName.equals(tumourMarker)) {
        // Column name
        def columnName = markerName + measurement

        // Change to array
        def markerIntensities = cells.stream()
                .mapToDouble({p -> p.getMeasurementList().getMeasurementValue(columnName)})
                .filter({d -> !Double.isNaN(d)})
                .toArray()

        // Instantiate KernelDensityEstimation class
        def estimation = new KernelDensityEstimation(markerIntensities, 0, 255)

        // Estimate
        def numbers = estimation.estimate() as Double[];

        // Find local minima of the density function (inflection of the distribution function)
        def minima = findLocalMinimaIndex(numbers)

        // Check that we have a minima
        assert (minima.size() > 0) : "Minima not detected for " + markerName

        // Get the threshold
        ArrayList<Double> arrayListNumbers = new ArrayList<>(Arrays.asList(numbers))

        // List of thresholds
        def thresholdList = []

        // Maximum Density
        def maxDensity = Collections.max(arrayListNumbers)
        // The intensity of where the max density occurs
        def maxDensityIndex = arrayListNumbers.indexOf(maxDensity)

        // 25% of the max density
        def upperThreshold = maxDensity * 0.25;

        // Filter the threshold list
        for (int i = 0; i < minima.size(); i++) {
            def index = minima[i]
            // Get threshold which is larger than the intensity with the max density, and less than 25% of the
            // max density.
            if ((numbers[index].compareTo(upperThreshold) <= 0) && (x[index].compareTo(x[maxDensityIndex]) >= 0)) {
                thresholdList.add(index)
            }
        }

        // Make sure the threshold list is not empty
        assert (thresholdList.size()) : "Threshold list empty for " + markerName

        // Get the final threshold
        def threshold = x[thresholdList[0]]
        println("Threshold for "+markerName+": "+threshold)

        // Add to the density array list
        densityArrayList.add(arrayListNumbers)
        csvColumnHeaders.add(markerName)
        thresholdMap.put(markerName, threshold)
    }
}



// ############################################## For tumour marker ####################################################
println("Calculating for the tumour marker")

if (!tumourMarker.isEmpty()){
    def columnName = tumourMarker + measurement

    // Get the cells which are positive for the baseline markers
    def filteredTumourIntensities = cells.stream()
            .mapToDouble({p -> p.getMeasurementList().getMeasurementValue(columnName)})
            .filter({d -> !Double.isNaN(d)})
            .collect()

    for (marker : baseLineMarkers){
        def baselineColumnName = marker + measurement
        def baselineMarkerIntensity = cells.stream()
                .mapToDouble({p -> p.getMeasurementList().getMeasurementValue(baselineColumnName)})
                .filter({d -> !Double.isNaN(d)})
                .collect()

        def newList = []

        def baselineThreshold = thresholdMap.get(marker)
        for (int i=0; i < filteredTumourIntensities.size(); i++){

            if (baselineMarkerIntensity[i] > baselineThreshold){
                newList.add(filteredTumourIntensities[i])
            }
        }
        filteredTumourIntensities = newList
    }

    def stats = new DescriptiveStatistics(filteredTumourIntensities as double[])
    def cutOffForTumour = stats.getPercentile(95)

    columnName = tumourMarker + measurement

    // Change to array
    def markerIntensities = cells.stream()
            .mapToDouble({p -> p.getMeasurementList().getMeasurementValue(columnName)})
            .filter({d -> !Double.isNaN(d)})
            .toArray()

    // Instantiate KernelDensityEstimation class
    def estimation = new KernelDensityEstimation(markerIntensities, 0, 255)

    // Estimate
    def numbers = estimation.estimate() as Double[];

    // Find local minima of the density function (inflection of the distribution function)
    def minima = findLocalMinimaIndex(numbers)

    // Check that we have a minima
    assert (minima.size() > 0)

    // Get the threshold
    ArrayList<Double> arrayListNumbers = new ArrayList<>(Arrays.asList(numbers))

    // threshold for tumour marker
    def threshold = 0

    // Filter the threshold list
    def index = minima[0]

    // Get the first threshold and compare to the cutOffForTumour.
    if ((x[index].compareTo(cutOffForTumour) <= 0)) {
        threshold = x[index]
    } else {
        threshold = cutOffForTumour
    }

    thresholdMap.put(tumourMarker, threshold)
    densityArrayList.add(arrayListNumbers)
    csvColumnHeaders.add(tumourMarker)
}

// Print out the thresholds
println(thresholdMap.toString())

// Save the tumour distribution if requested
if (saveDensities) {
    // Write data to csv
    println("Saving csv for tumour marker")
    def br = new BufferedWriter(new FileWriter(filePath))
    def sb = new StringBuilder()

    // Headers
    sb.append(csvColumnHeaders.join(","))
    sb.append("\n")

    for (int j = 0; j < x.length; j++) {
        for (int k=0; k<densityArrayList.size(); k++){
            sb.append(densityArrayList.get(k).get(j))
            sb.append(",")
        }
        sb.append("\n")
    }
    br.write(sb.toString())
    br.close()
    println("Done")
}