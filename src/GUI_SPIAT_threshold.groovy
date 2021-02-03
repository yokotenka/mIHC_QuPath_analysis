// ############################ Import ######################################

import com.sun.javafx.charts.Legend
import javafx.application.Platform
import javafx.collections.FXCollections
import javafx.collections.ObservableList
import javafx.geometry.Insets
import javafx.geometry.Orientation
import javafx.geometry.Pos
import javafx.scene.Node
import javafx.scene.Scene
import javafx.scene.chart.Chart
import javafx.scene.control.*
import javafx.scene.control.cell.PropertyValueFactory
import javafx.scene.layout.GridPane
import javafx.scene.paint.Color
import javafx.scene.shape.Line
import javafx.stage.Stage
import qupath.lib.gui.QuPathGUI
import qupath.lib.gui.scripting.QPEx
import qupath.lib.objects.classes.PathClassFactory
import qupath.lib.scripting.QP

import static qupath.lib.scripting.QP.*;


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
    private Double min

    /** The maximum of the range*/
    private Double max

    /** Constructor for the range class
     * @param min Minimum of the range
     * @param max Maximum of the range
     */
    Range(Double min, Double max){
        if (min.compareTo(max) >= 0 ){
            // Do error handling for when min is greater than equal to max
        }
        this.min = min
        this.max= max
    }

    /** Getter for the maximum value of the range
     * @return max
     */
    double getMax() {
        return max
    }

    /** Getter for the minimum value of the range
     * @return min
     */
    double getMin() {
        return min
    }

    /** Getter for the difference between the min and max values
     * @return max - min
     */
    double getDifference(){
        return max - min
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

import org.apache.commons.math3.util.FastMath

class NormalDistribution {
    /** The mean of the normal distribution */
    private double mean = 0
    /** The standard deviation of the normal distribution */
    private double standardDeviation = 1


    /** Constructor
     * @param mean
     * @param standardDeviation
     */
    NormalDistribution(double mean, double standardDeviation){
        this.mean = mean
        this.standardDeviation = standardDeviation
    }

    /** Constructor default is the standard normal distribution
     */
    NormalDistribution(){
    }

    /** Constructor only with the standard deviation
     * @param standardDeviation
     */
    NormalDistribution(double standardDeviation){
        this.standardDeviation = standardDeviation
    }

    /** Getter for the density for a specified range
     @param range The range of interest
     @param n The number of data points
     */
    Double[] getDensityOfRange(Range range, int n){
        Double[] densityArray = new Double[n]

        Double increment = range.getDifference() / n
        Double currentPoint = (double) range.getMin()

        for (int i=0; i<n; i++){
            densityArray[i] = density(currentPoint)
            currentPoint += increment
        }
        return densityArray
    }

    /** Calculates the density at a specific point
     @param x The point at which the density function is being evaluated at.
     */
    Double density(double x){
        double x0 = (x - mean) / standardDeviation
        double exponent = -0.5 * x0 * x0
        return FastMath.exp(exponent) / (standardDeviation * FastMath.sqrt(2 * Math.PI))
    }

    /** Setter for the mean
     * @param mean The new mean for the normal distribution.
     */
    void setMean(double mean) {
        this.mean = mean
    }

    /** Setter for the standard deviation
     * @param standardDeviation
     */
    void setStandardDeviation(double standardDeviation){
        this.standardDeviation = standardDeviation
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
    private double bandWidth

    /** The kernel function name. Defaults to the Gaussian kernel. In reality whatever kernel works the same. */
    private String kernelName ="gaussian"

    /** Weighting for each data point. Not implemented yet. */
    private Double[] weights = null

    /** The range of the data points */
    private Range range = null

    /** The estimated kernel density array */
    private Double[] estimate

    /** The number of data points to be included in the estimated array*/
    private int n = 2048

    /** The data points which the estimation is performed on */
    private double[] x


    /** Constructor
     * @param x The data points
     * @param min Minimum in the range.
     * @param max Maximum in the range.
     */
    KernelDensityEstimation(double[] x, double min, double max){
        this.x = x
        setBandWidthMethod("Default")
        this.range = new Range(min, max)
    }

    /** Estimation
     * @return kernelDensity The kernel density estimation for the given data, x.
     */
    Double[] estimate(){
        Double[] kernelDensity = new Double[n]
        Double[] currentDensity
        boolean isFirstElement = true

        // Initialise kernel distribution.
        NormalDistribution kernel = new NormalDistribution(bandWidth)

        // Iterate through every datapoint
        for (Double point : x){

            // Set the mean value for the kernel distribution.
            kernel.setMean(point)

            // Get the distribution of the kernel for the given range
            currentDensity = kernel.getDensityOfRange(range, n)

            // Calculate the kernel density
            if (isFirstElement){
                kernelDensity = currentDensity
                isFirstElement = false
            } else{
                for (int i=0; i<kernelDensity.length; i++){
                    kernelDensity[i] = kernelDensity[i] + currentDensity[i]
                }
            }
        }

        // Normalise
        for (int i=0; i<kernelDensity.length; i++) kernelDensity[i] = kernelDensity[i] / x.length

        return kernelDensity
    }

    /** Private method which will calculate the bandwidth to be used. Probably should create another class for this.
     * */
    private void setBandWidthMethod(String bandWidthName) throws UnsupportedOperationException {

        if (bandWidthName == "Silverman"){
            throw new UnsupportedOperationException("Silverman not yet supported")
        } else if (bandWidthName == "SJ"){
            throw new UnsupportedOperationException("Sheather & Jones not yet supported")
        } else {
            // Can change this if need be.
            DescriptiveStatistics da = new DescriptiveStatistics(x)
            double xStandardDeviation = da.getStandardDeviation()
            double iqr = da.getPercentile(75) - da.getPercentile(25)
            this.bandWidth = 0.9 * Math.min(xStandardDeviation, iqr/1.34) * Math.pow(x.length, -1/5)
        }
    }

    /** Getter for the bandwidth being used.
     * @return bandWidth
     */
    double getBandWidth(){
        return this.bandWidth
    }
}


//class BandWidth{
//    private String method
//
//    /** Calculate the bandWidth from the given method and the data
//     * @param method
//     * @param x
//     */
//    static getBandWidth(String method, Double[] x){
//
//        if (bandWidthName == "Silverman"){
//            throw new UnsupportedOperationException("Silverman not yet supported")
//        } else if (bandWidthName == "SJ"){
//            throw new UnsupportedOperationException("Sheather & Jones not yet supported")
//        } else {
//            // Can change this if need be.
//            DescriptiveStatistics da = new DescriptiveStatistics(x)
//            double xStandardDeviation = da.getStandardDeviation()
//            double iqr = da.getPercentile(75) - da.getPercentile(25)
//            this.bandWidth = 0.9 * Math.min(xStandardDeviation, iqr/1.34) * Math.pow(x.length, -1/5)
//        }
//    }
//}

/*
 * ############################################# Utility ###############################################################
 */

def static findLocalMinimaIndex(Double[] arr) {

    // Empty vector to store points of
    // local maxima and minima
    ArrayList<Double> mn = new ArrayList<Double>()
    int n = arr.length

    // Iterating over all points to check
    // local maxima and local minima
    for(int i = 1; i < arr.length - 1; i++) {
        // Condition for local minima
        if ((arr[i - 1].compareTo(arr[i])) > 0 &&
                (arr[i].compareTo(arr[i + 1]) < 0)) {
            mn.add(i)
        }
    }

    // Checking whether the last point is
    // local maxima or minima or none
    if (arr[n - 1].compareTo(arr[n - 2]) < 0){
        mn.add(n - 1)
    }
    return mn
}


// ############################ Marker Class ######################################
/** Marker class. Has name and information on whether it is selected or not. This is the row class
 * for the first two tables.
 */
class Marker {
    private String name;
    private CheckBox selected;

    /** Constructor
     * @param name
     * @param isSelected
     */
    Marker(String name, boolean isSelected){
        this.name = name
        this.selected = new CheckBox()
        this.selected.setSelected(isSelected)
    }

    /** Gettrer for the name
     * @return
     */
    String getName(){
        return name
    }

    /** Setter for the name
     * @param name
     */
    void setName(String name){
        this.name = name
    }

    /** Getter for the checkbox
     * @return
     */
    CheckBox getSelected(){
        return selected
    }

    /** Setter for the checkbox
     * @param selected
     */
    void setSelected(CheckBox selected){
        this.selected = selected
    }

    /** Checks if selected
     * @return
     */
    boolean isSelected(){
        return selected.isSelected()
    }
}


// ############################ Results Row class ######################################
/** Results row class. For displaying the results in th results table
 *
 */
class ResultsRow {
    private String propertyName
    private Double propertyValue

    ResultsRow(String propertyName, Double propertyValue){
        this.propertyName = propertyName
        this.propertyValue = propertyValue
    }

    String getPropertyName(){
        return this.propertyName
    }

    Double getPropertyValue(){
        return this.propertyValue
    }
}

// ############################## Measurement Classifier Class
/**
 * Helper class to calculate & apply thresholds, resulting in object classifications being set.
 * From Nina Tubau's https://github.com/ninatubau/QuPath_scripts/blob/main/scripts/multichannel_analysis.groovy
 */

// ############################ Get current image data ######################################

// Collect all the data
def server = getCurrentServer()
def imageData = getCurrentImageData()
def cells = getCellObjects()
def path = buildFilePath(PROJECT_BASE_DIR,"classifiers", "object_classifiers")



// ############################ Initialisation for Table creation ######################################
// Observable lists for markers to be thresholded
ObservableList<Marker> thresholdMarkers = FXCollections.observableArrayList()
// Oberservable list for baseline markers
ObservableList<Marker> baselineMarkers = FXCollections.observableArrayList()
// Observable list for marker names
ObservableList<String> markerNames = FXCollections.observableArrayList()

// Initialise
for (int i=0; i<server.nChannels(); i++){
    markerName = server.getChannel(i).getName()
    markerNames.add(markerName)
    thresholdMarkers.add(new Marker(server.getChannel(i).getName() as String, true))
    baselineMarkers.add(new Marker(server.getChannel(i).getName() as String , false))
}

//Settings to control the dialog boxes for the GUI
int col = 0
int row_col_0 = 0


// ############################## The Grid ##########################################
def gridPane = new GridPane()
gridPane.setVgap(5)
gridPane.setPadding(new Insets(10, 10, 10, 10));

// ############################ Table for Thresholded Markers ######################################
def requestLabelThreshold = new Label("Select the markers to be thresholded:")
requestLabelThreshold.setFont(javafx.scene.text.Font.font(15))
gridPane.add(requestLabelThreshold,col, row_col_0++, 1, 1)
requestLabelThreshold.setAlignment(Pos.CENTER)

// Initialise columns
TableColumn<Marker, String> nameColThresh = new TableColumn<>("Name")
nameColThresh.setCellValueFactory(new PropertyValueFactory<>("name"))

TableColumn<Marker, CheckBox> selectedColThresh = new TableColumn<>("Select")
selectedColThresh.setCellValueFactory(new PropertyValueFactory<>("selected"))

// Initialise table
TableView<Marker> tableThresh = new TableView<>()
tableThresh.setItems(thresholdMarkers)
tableThresh.getColumns().addAll(nameColThresh, selectedColThresh)

// add to gridpane
gridPane.add(tableThresh, 0, row_col_0++)


// ############################ Table for Baseline Markers ######################################
def requestLabelBaseline = new Label("Select the baseline markers:")
requestLabelBaseline.setFont(javafx.scene.text.Font.font(15))
gridPane.add(requestLabelBaseline, 0 , row_col_0++, 1, 1)
requestLabelBaseline.setAlignment(Pos.CENTER)

// Initialise columns
TableColumn<Marker, String> nameColBaseline = new TableColumn<>("Name")
nameColBaseline.setCellValueFactory(new PropertyValueFactory<>("name"))

TableColumn<Marker, CheckBox> selectedColBaseline = new TableColumn<>("Select")
selectedColBaseline.setCellValueFactory(new PropertyValueFactory<>("selected"))

// Initialise table
TableView<Marker> tableBaseline = new TableView<>()
tableBaseline.setItems(baselineMarkers)
tableBaseline.getColumns().addAll(nameColBaseline, selectedColBaseline)

// add to gridpane
gridPane.add(tableBaseline, 0, row_col_0++)


// ############################ ComboBox for Tumour Marker ######################################
def requestLabelTumourMarker = new Label("Select the tumour marker:")
requestLabelTumourMarker.setFont(javafx.scene.text.Font.font(15))
gridPane.add(requestLabelTumourMarker, 0, row_col_0++, 3, 1)
requestLabelTumourMarker.setAlignment(Pos.CENTER)

// initliase combobox
ComboBox comboBox = new ComboBox(markerNames)

// add to gridpane
gridPane.add(comboBox, 0, row_col_0++)


// ############################ Button to run script ######################################

Button startButton = new Button()
startButton.setText("Run Thresholding")
gridPane.add(startButton, 0, row_col_0++)

// ################################# Separator ##########################################
Separator separator = new Separator(Orientation.VERTICAL);
gridPane.add(separator, 1, 0, 1, row_col_0++)

// ################################# Results Selection ##########################################
def row_col_2 = 0
def requestLabelResults = new Label("Select marker to display results:")
requestLabelResults.setFont(javafx.scene.text.Font.font(15))
gridPane.add(requestLabelResults, 2, row_col_2++)
requestLabelResults.setAlignment(Pos.CENTER)

ComboBox comboBoxResults = new ComboBox(markerNames)
gridPane.add(comboBoxResults, 3, 0)

// ############################ Line Chart ######################################
import javafx.scene.chart.NumberAxis
import javafx.scene.chart.LineChart
import javafx.scene.layout.Pane
import javafx.scene.chart.XYChart.Series
import javafx.scene.chart.XYChart.Data
final NumberAxis xAxis = new NumberAxis();
final NumberAxis yAxis = new NumberAxis();
xAxis.setAutoRanging(false)
xAxis.setUpperBound(255)
xAxis.setLabel("Intensity")
yAxis.setAutoRanging(false)

//creating the chart
final LineChart<Number,Number> lineChart = new LineChart<Number,Number>(xAxis,yAxis);

lineChart.setTitle("Density Estimation");
lineChart.setCreateSymbols(false)
lineChart.setLegendVisible(false)

Line valueMarker = new Line()

Pane pane = new Pane()
pane.getChildren().addAll(lineChart, valueMarker)

gridPane.add(pane, 2, row_col_2++, 2, 1)


// There is definitely a faster way of doing this. Somehting is wrong
Label thresholdLabel = new Label()
Label percentageLabel = new Label()
gridPane.add(thresholdLabel, 2, row_col_2 + 2)
gridPane.add(percentageLabel, 2, row_col_2 + 3)

// ############################ Display window ######################################
Platform.runLater {
    def stage = new Stage()
    stage.initOwner(QuPathGUI.getInstance().getStage())
    stage.setScene(new Scene( gridPane))
    stage.setTitle("SPIAT Thresholding")
    stage.show()
}

// ############################ Action for when button is pressed ######################################
// Threshold map
def thresholdMap = new HashMap<String, Double>()
def percentageMap = new HashMap<String, Double>()

def n = 2048
def densityArrayList = new ArrayList<ArrayList<Double>>(server.nChannels())

def maxDensityList = []
def maxIntensityList = []
// Measurement being used

// MAKE THIS AN OPTION
def measurement = ": Cell: Mean"


// x-axis
def increment = 255 / n
def val = 0
def x = new double[n]
for (int i=0; i < n; i++){
    x[i] = val
    val += increment
}
// ######################## Get the selected data
def selectedMarkers = []
def selectedBaselineMarkers = []
def invalidInputs = []

// For the colour for each of the markers
HashMap<String, Integer> colourHashMap = new HashMap<>()


/*
 *
 *  EDIT so that when pressing button, it remebers what it has already calculated.
 *
 */
resetDetectionClassifications()
startButton.setOnAction((event) -> {

    selectedMarkers = []
    selectedBaselineMarkers = []
    invalidInputs = []
    lineChart.getData().clear()

    def tumourMarker = comboBox.getValue()

    for (int i = 0; i < server.nChannels(); i++) {
        if (thresholdMarkers.get(i).isSelected()) {
            selectedMarkers.add(server.getChannel(i).getName())
        }
        if (baselineMarkers.get(i).isSelected()) {
            selectedBaselineMarkers.add(server.getChannel(i).getName())
        }
        if (baselineMarkers.get(i).isSelected() && !thresholdMarkers.get(i).isSelected()){
            invalidInputs.add(server.getChannel(i).getName())
        }
    }

    // FIX THIS ERROR HANDLING CANNOT .isEmpty() on a null object
    // ############################ Error handling ######################################
    if (comboBox.getSelectionModel().isEmpty() || !invalidInputs.isEmpty()){
        Alert errorAlert = new Alert(Alert.AlertType.ERROR);
        errorAlert.setHeaderText("Input not valid");
        def errorMsg = ""
        if (!invalidInputs.isEmpty()){

            if (invalidInputs.size() == 1){
                errorMsg += "- "+invalidInputs.toString() + " is selected as a baseline marker but is not selected for thresholding.\n"
            } else{
                errorMsg += "- "+invalidInputs.toString() + " are selected as a baseline markers but are not selected for thresholding.\n"
            }

        }
        if (comboBox.getSelectionModel().isEmpty()){
            errorMsg += "- "+"Tumour not selected"
        }

        errorAlert.setContentText(errorMsg)
        errorAlert.showAndWait()
    }

    /*
    * ############################################ Qupath Scripting #######################################################
    */




// ############################################### For regular markers #################################################
    println("Calculating for non tumour markers")

    def markerIndex =0
    for (markerName in selectedMarkers){
        if (!markerName.equals(tumourMarker)) {
            // Column name
            def columnName = markerName + measurement
            // Change to array
            def markerIntensities = cells.stream()
                    .mapToDouble({ p -> p.getMeasurementList().getMeasurementValue(columnName) })
                    .filter({ d -> !Double.isNaN(d) })
                    .toArray()

            maxIntensityList.add(Collections.max(markerIntensities.toList()))

            // Instantiate KernelDensityEstimation class
            def estimation = new KernelDensityEstimation(markerIntensities, 0, 255)

            // Estimate
            def numbers = estimation.estimate() as Double[]

            // Find local minima of the density function (inflection of the distribution function)
            def minima = findLocalMinimaIndex(numbers)

            // Check that we have a minima
            assert (minima.size() > 0): "Minima not detected for " + markerName

            // Get the threshold
            ArrayList<Double> arrayListNumbers = new ArrayList<>(Arrays.asList(numbers))

            // List of thresholds
            def thresholdList = []

            // Maximum Density
            def maxDensity = Collections.max(arrayListNumbers)
            // The intensity of where the max density occurs
            def maxDensityIndex = arrayListNumbers.indexOf(maxDensity)

            maxDensityList.add(maxDensity)

            // 25% of the max density
            def upperThreshold = maxDensity * 0.25

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
            assert (thresholdList.size()): "Threshold list empty for " + markerName

            def threshold = x[thresholdList[0]]

            println("Threshold for " + markerName + ": " + threshold.toString())

            // Add to the density array list
            densityArrayList.add(arrayListNumbers)
            thresholdMap.put(markerName as String, threshold as Double)

            // For the line chart
            Series series = new Series();
            for (int j=0; j < x.size(); j++){
                series.getData().add(new Data(x[j], numbers[j]))
            }
            series.setName(markerName);
            lineChart.getData().add(series)
            lineChart.getData().get(markerIndex++).getNode().setVisible(false)

            Series vertLine = new Series()
            vertLine.getData().add(new Data(thresholdMap.get(markerName), 0))
            vertLine.getData().add(new Data(thresholdMap.get(markerName), maxDensity * 1.25))
            vertLine.setName("Threshold for "+markerName)
            lineChart.getData().add(vertLine)
            lineChart.getData().get(markerIndex).getNode().setVisible(false)
            lineChart.getData().get(markerIndex++).getNode().setStyle("-fx-stroke: lightslategray;")

            def count = 0
            for (markerIntensity in markerIntensities){
                if (markerIntensity > threshold){
                    count++
                }
            }
            def percentage = (count / markerIntensities.size()) * 100
            percentageMap.put(markerName, percentage as Double)





//            def defaultColor = colourHashMap.containsKey(markerName) ? colourHashMap.get(markerName) : new Integer()

            def positive = cells.findAll {it.getMeasurementList().getMeasurementValue(columnName) > threshold}
            positive.each {
                def currentClass = it.getPathClass()
                def pathClass
                // Create a classifier - only assign the color if this is a single classification
                if (currentClass == null) {
                    pathClass = PathClassFactory.getPathClass(markerName)
                }
                else
                    pathClass = PathClassFactory.getDerivedPathClass(currentClass, markerName, null)
                it.setPathClass(pathClass)
            }
        }
    }
// ############################################## For tumour marker ####################################################
    println("Calculating for the tumour marker")
    if (!tumourMarker.isEmpty()){
        def columnName = tumourMarker + measurement

        // ###### Get the cells which are positive for the baseline markers ######
        def filteredTumourIntensities = cells.stream()
                .mapToDouble({p -> p.getMeasurementList().getMeasurementValue(columnName)})
                .filter({d -> !Double.isNaN(d)})
                .collect()

        for (marker in selectedBaselineMarkers){
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

        // Change to array
        def markerIntensities = cells.stream()
                .mapToDouble({p -> p.getMeasurementList().getMeasurementValue(columnName)})
                .filter({d -> !Double.isNaN(d)})
                .toArray()

        maxIntensityList.add(Collections.max(markerIntensities.toList()))


        // ######## Now find KDE like other markers ###############
        def estimation = new KernelDensityEstimation(markerIntensities, 0, 255)

        // Estimate
        def numbers = estimation.estimate() as Double[]

        // Find local minima of the density function (inflection of the distribution function)
        def minima = findLocalMinimaIndex(numbers)

        // Check that we have a minima
        assert (minima.size() > 0)

        // Get the threshold
        ArrayList<Double> arrayListNumbers = new ArrayList<>(Arrays.asList(numbers))

        def maxDensity = Collections.max(arrayListNumbers)
        maxDensityList.add(maxDensity)
        def maxDensityIndex = arrayListNumbers.indexOf(maxDensity)

        // threshold for tumour marker
        def threshold = 0


        // ######### Find threhsold ############
        for (int i=0; i < minima.size(); i++){
            def index = minima[i]
            if (x[index] < x[maxDensityIndex]){
                continue
            }
            // Get the first threshold and compare to the cutOffForTumour.
            if ((x[index].compareTo(cutOffForTumour) <= 0)) {
                threshold = x[index]
                break
            } else {
                threshold = cutOffForTumour
                break
            }
        }

        println("Threshold for "+tumourMarker+": "+threshold)

        // Put the threshold in the map
        thresholdMap.put(tumourMarker, threshold)
        densityArrayList.add(arrayListNumbers)

        // ####### Add the data into the lineChart
        Series series = new Series();
        for (int j=0; j < x.size(); j++){
            series.getData().add(new Data(x[j], numbers[j]))
        }
        series.setName(tumourMarker);
        lineChart.getData().add(series)
        lineChart.getData().get(markerIndex++).getNode().setVisible(false)

        // Vertical line which will show where the threshold is.
        Series vertLine = new Series()
        vertLine.getData().add(new Data(thresholdMap.get(tumourMarker), 0))
        vertLine.getData().add(new Data(thresholdMap.get(tumourMarker), maxDensity * 1.25))
        vertLine.setName("Threshold for "+tumourMarker)
        lineChart.getData().add(vertLine)
        lineChart.getData().get(markerIndex).getNode().setVisible(false)
        lineChart.getData().get(markerIndex++).getNode().setStyle("-fx-stroke: lightslategray;")

        // ###### Calculate percentage
        def count = 0

        // Can improve this one
        for (markerIntensity in markerIntensities){
            if (markerIntensity > threshold){
                count++
            }
        }
        def percentage = (count / markerIntensities.size()) * 100
        percentageMap.put(tumourMarker, percentage as Double)

        // Set class for each cell
        def positive = cells.findAll {it.getMeasurementList().getMeasurementValue(columnName) > threshold}
        positive.each {
            def currentClass = it.getPathClass()
            def pathClass
            // Create a classifier - only assign the color if this is a single classification
            if (currentClass == null) {
                pathClass = PathClassFactory.getPathClass(tumourMarker)
            }
            else
                pathClass = PathClassFactory.getDerivedPathClass(currentClass, tumourMarker, null)
            it.setPathClass(pathClass)
        }
    }
    fireHierarchyUpdate()
});



// ################### Utility for selecting objects #############################

/**Checks if classificationName is in the pathClass recursively
 * @param pathClass
 * @param classificationName
 * @return
 */
boolean checkForSingleClassification(def pathClass, classificationName) {
    if (pathClass == null)
        return false
    if (pathClass.getName() == classificationName)
        return true
    return checkForSingleClassification(pathClass.getParentClass(), classificationName)
}

/** Checks if all the classification names in the array are in the pathclass
 * @param pathClass
 * @param classificationNames
 * @return
 */
boolean checkForClassifications(def pathClass, String...classificationNames) {
    if (classificationNames.length == 0)
        return false
    for (String name : classificationNames) {
        if (!checkForSingleClassification(pathClass, name))
            return false
    }
    return true
}


// ######################### Action upon combobox selection for results #############################

def currentlySelected = getSelectedObjects()

comboBoxResults.setOnAction {

    // Get the selected marker
    String selectedForResults = comboBoxResults.getValue()

    // Get the axes for the line chart
    NumberAxis lineChartYAxis = (NumberAxis) lineChart.getYAxis()
    NumberAxis lineChartXAxis = (NumberAxis) lineChart.getXAxis()


    for (int k=0; k < selectedMarkers.size() * 2 ; k+=2){
        if (lineChart.getData().get(k).getNode().isVisible()){
            // Change previously visible series to invisible
            lineChart.getData().get(k).getNode().setVisible(false)
            lineChart.getData().get(k+1).getNode().setVisible(false)
        }
        if (lineChart.getData().get(k).getName()==selectedForResults){
            // Change the selected marker series to visible
            lineChart.getData().get(k).getNode().setVisible(true)
            lineChart.getData().get(k+1).getNode().setVisible(true)

            // Change the upper bound of the probability to suit the marker
            lineChartYAxis.setUpperBound(maxDensityList.get(k/2 as int) * 1.25)
            lineChartYAxis.setLowerBound(0)

            // Change the upper bound of the intensity to suit the marker
            lineChartXAxis.setUpperBound(Math.ceil(maxIntensityList.get(k/2 as int) * 1.25))
            lineChartXAxis.setLowerBound(0)

            // title change
            lineChart.setTitle("Density Estimation for " + selectedForResults);

            // Change the displayed threshold and percentage
            thresholdLabel.setText("Threshold for "+selectedForResults+": "+thresholdMap.get(selectedForResults))
            percentageLabel.setText("Percentage for "+selectedForResults+" positive cells in tissue: "+percentageMap.get(selectedForResults) + "%")

            // Select the cells with the markers being positive.
            if (currentlySelected != null) {
                imageData.getHierarchy().selectionModel.deselectObjects(currentlySelected)
            }
            def positive = cells.stream().filter({checkForClassifications(it.getPathClass(), selectedForResults)}).collect()
            imageData.getHierarchy().selectionModel.selectObjects(positive)
            currentlySelected = positive
        }
    }
}



