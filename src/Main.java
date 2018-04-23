import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.lang.Math;
public class Main {

    public static void main(String[] args) {
//        testFFT();
//        testTwoDFFT();
//        testSignalsDataGenerator();
//        doCommonSignalsProblem(); //Problem 1, output in data directory
//        doQuestionTwo(); //Question 2: Sums of signals vs products of signals
//        doQuestionThree(); //Question 3: Affect of phase on PSD

        doQuestionFour(); //Question 4: Filtering
    }

    public static void testFFT(){
        ComplexNumber[] testInputs = generateTestFFTData();
        ComplexNumber[] results = fastFourierTransform(testInputs, 1);
        System.out.println("\nThe fft of the data\n");
        for(int i = 0; i < results.length; i++){
            System.out.format("%f + %fi\n", results[i].getzReal(), results[i].getzImaginary());
        }

        ComplexNumber[] inverseResults = fastFourierTransform(results, -1);
        System.out.println("\nThe inverse of the data\n");
        for(int i = 0; i < results.length; i++){
            System.out.format("%f + %fi\n", inverseResults[i].getzReal(), inverseResults[i].getzImaginary());
        }
    }

    public static void testTwoDFFT(){
        System.out.println("Testing 2D FFT...");
        ComplexNumber[][] data = generateTest2DData();
        ComplexNumber[][] results = twoDFFT(data);
        System.out.println("\nReals\n");
        printTwoDArray(results);
    }

    public static void testSignalsDataGenerator(){
        SignalsDataGenerator testGenerator = new SignalsDataGenerator("data/test_output.txt");
        testGenerator.writeResult("Hello world!");
        testGenerator.flushFile();
        testGenerator.clearFile();
        testGenerator.writeResult("Hello world 2!");
        testGenerator.closeFile();

    }

    public static void doCommonSignalsProblem(){
        //Commented out to prevent constant rewriting of the same values to output files, uncomment for functionality
        doFS();
        doGS();
        doPSDEstimates(); //Does the PSD Estimates for f50 and g50, as per question 1b
    }

    private static void doQuestionTwo(){
        //x(t) = v1(t) + v2(t)
        //y(t) = v1(t) * v2(t)

        //v1(t) = sin(2*pi*13*t)
        //v2(t) = sin(2*pi*31*t)

        double[] v1 = generateV(13);
        double[] v2 = generateV(31);

        ComplexNumber[] xt = generateXT(v1, v2);
        ComplexNumber[] yt = generateYT(v1, v2);

        doPSD(xt, "data/xt_psd.txt");
        doPSD(yt, "data/yt_psd.txt");
    }

    //This method generates v in the time domain, where v1[i] = sin(2*pi*f*t), where t=i/512 (the time interval)
    //f - frequency to be calculated
    private static double[] generateV(int frequency){
        final double INTERVAL = 1.0/512.0;
        double[] results = new double[512];
        for(int i = 0; i < results.length; i++){
            double t = i * INTERVAL;
            double answer = Math.sin(2 * Math.PI * frequency * t);
            results[i] = answer;
        }
        return results;
    }

    //This method takes two functions and adds them together, returning an array of complex numbers
    private static ComplexNumber[] generateXT(double[] v1, double[] v2){
        ComplexNumber[] results = new ComplexNumber[v1.length];
        for(int i = 0; i < results.length; i++){
            ComplexNumber answer = new ComplexNumber(v1[i] + v2[i], 0);
            results[i] = answer;
        }
        return results;
    }

    //This method takes two functions and multiplies them together, returning an array of complex numbers
    private static ComplexNumber[] generateYT(double[] v1, double[] v2){
        ComplexNumber[] results = new ComplexNumber[v1.length];
        for(int i = 0; i < results.length; i++){
            ComplexNumber answer = new ComplexNumber(v1[i] * v2[i], 0);
            results[i] = answer;
        }
        return results;
    }

    public static void doFS(){
        runFS(3, "data/fs3.txt");
        runFS(10, "data/fs10.txt");
        runFS(50, "data/fs50.txt");
    }

    public static void runFS(int sTerm, String filename){
        SignalsDataGenerator writer = new SignalsDataGenerator(filename);
        writer.clearFile();
        //Interval for t(time) to increase by
        final double INTERVAL = 1.0/512.0;
        final int ITERATIONS = 512;

        for(int i = 1; i <= ITERATIONS; i++){
            double t = i * INTERVAL;
            double result = 0;
            for(int k = 1; k <= sTerm; k++){
                result += (Math.sin(2.0 * Math.PI * (2*k - 1) * t))/(2*k-1);
            }
            String toWrite = Double.toString(result);
            writer.writeResult(toWrite);
        }

        writer.closeFile();
    }

    public static void doGS(){
        runGS(3, "data/gs3.txt");
        runGS(10, "data/gs10.txt");
        runGS(50, "data/gs50.txt");
    }

    public static void runGS(int sTerm, String filename){
        SignalsDataGenerator writer = new SignalsDataGenerator(filename);
        writer.clearFile();
        //Interval for t(time) to increase by
        final double INTERVAL = 1.0/512.0;
        final int ITERATIONS = 512;

        for(int i = 1; i <= ITERATIONS; i++){
            double t = i * INTERVAL;
            double result = 0;
            for(int k = 1; k <= sTerm; k++){
                result += (Math.sin(2.0 * Math.PI * (2*k) * t))/(2*k);
            }
            String toWrite = Double.toString(result);
            writer.writeResult(toWrite);
        }

        writer.closeFile();
    }
    //Question 1B: Displays only the positive frequencies
    private static void doPSDEstimates(){
        final String FSPSD_FILENAME = "data/fs50PSD.txt";
        final String GSPSD_FILENAME = "data/gs50PSD.txt";

        final String FS50_FILENAME = "data/fs50.txt";
        final String GS50_FILENAME = "data/gs50.txt";

        final String FL_IN_FILENAME = "data/fl_in.txt";
        final String FL_OUT_FILENAME = "data/fl_out.txt";

        final String GL_IN_FILENAME = "data/gl_in.txt";
        final String GL_OUT_FILENAME = "data/gl_out.txt";

        doPSD(FS50_FILENAME, FSPSD_FILENAME);
        doPSD(GS50_FILENAME, GSPSD_FILENAME);

        doPSD(FL_IN_FILENAME, FL_OUT_FILENAME);
        //TODO: Get feedback to verify that fl_in and gl_in data are correct
        //doPSD(GL_IN_FILENAME, GL_OUT_FILENAME);


    }

    private static void doPSD(String inputDataFilename, String outputDataFilename){
        ComplexNumber[] inputNumbers = parseInputs(inputDataFilename);
        ComplexNumber[] fftResults = fastFourierTransform(inputNumbers, 1);
        double[] psdResults = convertFFTToPSD(fftResults);
        outputPSDResults(psdResults, outputDataFilename);
    }

    //Used if we already have the ComplexNumber[] array from some other source
    private static void doPSD(ComplexNumber[] inputNumbers, String outputDataFilename){
        ComplexNumber[]fftResults = fastFourierTransform(inputNumbers, 1);
        double[] psdResults = convertFFTToPSD(fftResults);
        outputPSDResults(psdResults, outputDataFilename);
    }

    //Assumes a file where each line is a single number, to be converted into a complex number.
    //Assumes inputs are real numbers (non-complex or imaginary)
    private static ComplexNumber[] parseInputs(String filename){
        //TODO: Magic number (512)
        ComplexNumber[] numbers = new ComplexNumber[512];
        int count = 0;
        try {
            File file = new File(filename);
            BufferedReader br = new BufferedReader(new FileReader(file));

            String line;
            while((line = br.readLine()) != null){
                if(!line.isEmpty()) {
                    double value = Double.parseDouble(line);
                    ComplexNumber number = new ComplexNumber(value, 0);
                    numbers[count] = number;
                    count++;
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        if(count != numbers.length) {
            System.out.println("Warning: parseInputs() did not parse the expected amount of numbers. Parsed a total of " + Integer.toString(count) + " numbers");
        }
        return numbers;
    }

    private static double[] convertFFTToPSD(ComplexNumber[] fftValues){
        //To convert FFT to PSD, you take the modulus and square it
        double[] results = new double[fftValues.length];
        for(int i = 0; i < results.length; i++){
            double answer = Math.pow(fftValues[i].getModulus(), 2);
            results[i] = answer;
        }
        return results;
    }

    //Only outputs the positive (first half) of PSD results, due to symmetrical nature
    private static void outputPSDResults(double[] values, String filename){
        SignalsDataGenerator writer = new SignalsDataGenerator(filename);
        writer.clearFile();
        for(int i = 0; i < values.length; i++){
            if(values[i] >= 0){
                String toWrite;
                if(values[i] < .0001 && values[i] > -0.0001) {
                    toWrite = "0";
                } else{
                    double value = values[i];
                    toWrite = Double.toString(values[i]);
                }
                writer.writeResult(toWrite);
            }
        }
        writer.closeFile();
    }

    //Method that applies a Fast Fourier Transform to an Nx1 vector of complex numbers
    //Params: inputs - the Nx1 array
    //        d - the direction code. Either 1 for FFT or -1 for inverse FFT
    public static ComplexNumber[] fastFourierTransform(ComplexNumber[] testInputs, int d){
        int k;
        ComplexNumber t;
        ComplexNumber[] inputs = new ComplexNumber[testInputs.length];
        System.arraycopy(testInputs, 0, inputs, 0, testInputs.length);

        //1) Set theta = (-2*pi*d)/N and r = N/2
        double N = inputs.length;
        double theta = (-2 * Math.PI * d)/N;
        int r = (int) Math.floor(N/2);
        //2)For i = 1 to N-1:
        for(int i = 1; i < N; i++){
            //2a: Set w = cos(i*theta) + jsin(i*theta)
            double real = Math.cos(i*theta);
            double imaginary = Math.sin(i*theta);
            ComplexNumber w = new ComplexNumber(real, imaginary);

            //2b: For k = 0 to N-1 do:
            for(k = 0; k < N; k++){
                //b-1: Set u=1
                ComplexNumber u = new ComplexNumber(1, 0);
                //b-2: For m = 0 to r-1 do:
                for(int m = 0; m < r; m++){
                    t = inputs[k+m].subtract(inputs[(k+m+r)]);
                    inputs[k+m] = inputs[k+m].add(inputs[k+m+r]);
                    inputs[k+m+r] = t.multiply(u);
                    u = w.multiply(u);
                }
                //b-3: set k = k+2r
                k = k + 2*r -1; //-1 due to loop incrementation
            }
            //2c: Set i = 2i and r = r/2
            i = 2*i-1; //-1 due to loop incrementation
            r = r/2;
        }
        //3) For i = 0 to N-1 do:
        for(int i=0; i < N; i++){
            //3a: Set r=i and k = 0
            r = i;
            k = 0;
            //3b: For m = 1 to N-1 do:
            int m = 1;
            while(m < N){
                k = 2 * k + (r % 2);
                r = r / 2;
                m = 2 * m;
            }
            //3c: If k > i do:
            if(k > i){
                t = inputs[i];
                inputs[i] = inputs[k];
                inputs[k] = t;
            }
        }
        //4) If d < 0 then for i = 0 to N-1 do:
        if(d < 0){
            for(int i = 0; i < N; i++){
                inputs[i] = inputs[i].divide(N);
            }
        }
        return inputs;
    }

    public static ComplexNumber[][] twoDFFT(ComplexNumber[][] userInputs){
        //Copy array
        ComplexNumber[][] inputs = new ComplexNumber[userInputs.length][userInputs[0].length];
        for(int i = 0; i < inputs.length; i++) {
            System.arraycopy(userInputs[i], 0, inputs[i], 0, inputs[i].length);
        }

        final int STANDARD_FFT = 1;
        ComplexNumber[][] temp = new ComplexNumber[inputs.length][inputs[0].length];
        ComplexNumber finalResult[][] = new ComplexNumber[inputs.length][inputs[0].length];
        for(int k = 0; k < inputs.length; k++){
            //Do rows
            ComplexNumber[] row = new ComplexNumber[inputs[k].length];
            System.arraycopy(inputs[k],0,row,0,row.length);
            ComplexNumber[] rowResult = fastFourierTransform(row, STANDARD_FFT);
            temp[k] = rowResult;
        }

        System.out.println("First pass:");
        printTwoDArray(temp);
        //Do Columns
        for(int n = 0; n < inputs[0].length; n++){
            ComplexNumber[] column = new ComplexNumber[inputs.length];
            for(int i = 0; i < column.length; i++){
                //add values to form the column
                column[i] = temp[i][n];
            }
            System.out.format("\nColumn: ");
            for(int x = 0; x < column.length; x++) {
                System.out.format("%.2f ", column[x].getzReal());
            }
            ComplexNumber[] columnResult = fastFourierTransform(column, STANDARD_FFT);
            System.out.format("\nColumnFFT: ");
            for(int x = 0; x < columnResult.length; x++) {
                System.out.format("%.2f ", columnResult[x].getzReal());
            }
            //Store as column in final result
            for(int i = 0; i < column.length; i++){
                //add values to form the column
                finalResult[i][n] = columnResult[i];
            }
        }
        return finalResult;
    }

    private static void printTwoDArray(ComplexNumber[][] array) {
        for(int i = 0; i < array.length; i++){
            for(int j = 0; j < array[0].length; j++){
                System.out.format("%.2f ", array[i][j].getzReal());
                if(j % 7 == 0 && j != 0){
                    System.out.format("\n");
                }
            }
        }
    }

    private static ComplexNumber[] generateTestFFTData(){
        double[] reals = {26160.0, 19011.0, 18757, 18405, 17888, 14720, 14285,
        17018, 18014, 17119, 16400, 17497, 17846, 15700, 17636, 17181};
        ComplexNumber[] data = new ComplexNumber[reals.length];
        for(int i = 0; i < reals.length; i++){
            data[i] = new ComplexNumber(reals[i], 0);
        }
        return data;
    }

    private static ComplexNumber[][] generateTest2DData() {
        double[][] reals = new double[][]{
                {0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 70, 80, 90, 0, 0, 0},
                {0, 0, 90, 100, 110, 0, 0, 0},
                {0, 0, 110, 120, 130, 0, 0, 0},
                {0, 0, 130, 140, 150, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0}
        };
        ComplexNumber[][] data = new ComplexNumber[reals.length][reals[0].length];
        for (int i = 0; i < reals.length; i++) {
            for(int j = 0; j < reals[0].length; j++){
                ComplexNumber entry = new ComplexNumber(reals[i][j], 0);
                data[i][j] = entry;
            }
        }
        return data;
    }

    private static void doQuestionThree(){
//        examineTimeVariations(); //Question 3a
        examineHT(); //Question 3b
    }

    private static void examineTimeVariations(){
        ComplexNumber[] start = generateSinglePulse(0, 256);
        ComplexNumber[] mid = generateSinglePulse(127, 256);
        ComplexNumber[] random = generateSinglePulse(38, 256);

        doPSD(start, "data/question3/start.txt");
        doPSD(mid, "data/question3/mid.txt");
        doPSD(random, "data/question3/random.txt");
    }

    private static ComplexNumber[] generateSinglePulse(int pulseLocation, int length){
        ComplexNumber[] results = new ComplexNumber[length];
        for(int i = 0; i < length; i++){
            double value = 0.0;
            if(i == pulseLocation){
                value = 1.0;
            }
            results[i] = new ComplexNumber(value, 0);
        }
        return results;
    }

    private static void examineHT(){
        ComplexNumber[] ht = generateHT(0, "data/question3/ht_raw.txt"); //h(t) = sin(20 * PI * t)
        ComplexNumber[] ht_shifted = generateHT(0.43, "data/question3/ht_shifted.txt");

        //Examine the FFTs for Problem 3.b.i
        ComplexNumber[] ht_fft = fastFourierTransform(ht, 1);
        ComplexNumber[] ht_shifted_fft = fastFourierTransform(ht_shifted, 1);
        outputFFTRealsToFile(ht_fft, "data/question3/ht_fft.txt");
        outputFFTRealsToFile(ht_shifted_fft, "data/question3/ht_shifted_fft.txt");


        doPSD(ht, "data/question3/ht_raw_psd.txt");
        doPSD(ht_shifted, "data/question3/ht_shifted_psd.txt");
    }

    //Generates h(t)
    //c - phase shift in function sin(20 * pi * (t-c))
    private static ComplexNumber[] generateHT(double c, String filename){
        SignalsDataGenerator writer = new SignalsDataGenerator(filename);
        ComplexNumber[] results = new ComplexNumber[512];
        final double INTERVAL = 1.0/512.0;
        for(int i = 0; i < results.length; i++){
            double t = i*INTERVAL;
            ComplexNumber answer = new ComplexNumber(Math.sin(20.0 * Math.PI * (t-c)), 0);
            results[i] = answer;
            writer.writeResult(Double.toString(answer.getzReal()));
        }
        writer.closeFile();
        return results;
    }

    //purpose: This method is used to output a ComplexNumber to a file. For use in problem 3.b.i
    private static void outputFFTRealsToFile(ComplexNumber[] values, String filename){
        SignalsDataGenerator writer = new SignalsDataGenerator(filename);
        writer.clearFile();
        for(int i = 0; i < values.length; i++){
            String real;
            if(values[i].getzReal() > -.0000000001 && values[i].getzReal() < .0000000001) {
                real = "0";
            } else{
                real = Double.toString(values[i].getzReal());
            }
            writer.writeResult(real);
        }
        writer.closeFile();
    }

    private static void doQuestionFour(){
        double[] fs50 = parseFS50();
        double[] lowPassFS = lowPassFilter(fs50, 7);
        double[] highPassFS = highPassFilter(fs50, 43);
        double[] bandPassFS = bandPassFilter(fs50, 4, 7); //inclusive
        double[] notchFS = notchFilter(fs50, 4, 7); //inclusive, SUPPRESSED

        System.out.println("Values:");
        System.out.println("fs50      lowPass      highPass      bandPass      notch      ");
        for(int i = 0; i < fs50.length; i++){
            System.out.format("%6.1f  %10.1f  %10.1f  %10.1f  %10.1f\n",fs50[i], lowPassFS[i], highPassFS[i], bandPassFS[i], notchFS[i]);
        }
    }

    //Takes the data in fs50PSD.txt and converts it into a double[]
    private static double[] parseFS50() {
        final String FS50_PSD_FILENAME = "data/fs50PSD.txt";
        double[] results = new double[256];
        try {
            File file = new File(FS50_PSD_FILENAME);
            BufferedReader br = new BufferedReader(new FileReader(file));
            String line;
            for (int i = 0; i < results.length; i++) {
                if ((line = br.readLine()) != null && !line.isEmpty()){
                    results[i] = Double.parseDouble(line);
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return results;
    }

    //Puts a low-pass filter on a PSD
    //inclusiveEndIndex - an index saying what is the last non-zero value to PASS
    private static double[] lowPassFilter(double[] originalInput, int numValuesToPass){
        //Convert into appropriate indexes for filter() method
        int startingIndex = numValuesToPass;
        int endingIndex = 49;

        return filter(originalInput, startingIndex, endingIndex);
    }

    private static double[] highPassFilter(double[] originalInput, int numValuesToPass){
        int startingIndex = 0;
        int endingIndex = 49 - numValuesToPass;

        return filter(originalInput, startingIndex, endingIndex);
    }

    private static double[] bandPassFilter(double[] originalInput, int startingIndexToPass, int endingIndexToPassInclusive){
        //Band pass must be done in two filter passes, due to limitations in current filter method
        int startingIndex = 0;
        int endingIndex = startingIndexToPass - 1;

        //first pass
        double[] tempResults = filter(originalInput, startingIndex, endingIndex);

        int numRemoved = endingIndex - startingIndex + 1;

        startingIndex = endingIndexToPassInclusive + 1 - numRemoved;
        endingIndex = 49;

        //second pass
        return filter(tempResults, startingIndex, endingIndex);
    }

    private static double[] notchFilter(double[] originalInput, int startingIndexToSuppress, int endIndexToSuppress){
        return filter(originalInput, startingIndexToSuppress, endIndexToSuppress);
    }

    //Filters out non-zero values, based on the parameter indexes. These indexes refer to non-zero entries, all zero entries are ignored
    //For example, the first non-zero entry in the originalInput array will be found at inclusiveStartIndex = 0
    //THE VALUES BETWEEN THE INDEXES STATE WHAT TO FILTER OUT
    private static double[] filter(double[] originalInput, int inclusiveStartIndex, int inclusiveEndIndex){
        int non_zero_index = 0;
        double[] answers = new double[originalInput.length];
        for(int i = 0; i < originalInput.length; i++){
            double currentInput = originalInput[i];
            if(currentInput != 0) {
                if (non_zero_index >= inclusiveStartIndex && non_zero_index <= inclusiveEndIndex) {
                    //Filter out this number
                    answers[i] = 0.0;
                } else{
                    answers[i] = originalInput[i];
                }
                non_zero_index++;
            } else {
                answers[i] = 0.0;
            }
        }
        return answers;
    }

}
