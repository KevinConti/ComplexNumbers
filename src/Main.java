import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.lang.Math;
import java.util.Dictionary;

public class Main {

    public static void main(String[] args) {
//        testFFT();
//        testTwoDFFT();
//        testSignalsDataGenerator();
//        doCommonSignalsProblem(); //Problem 1, output in data directory
//        doQuestionTwo(); //Question 2: Sums of signals vs products of signals
//        doQuestionThree(); //Question 3: Affect of phase on PSD
//        doQuestionFour(); //Question 4: Filtering
//        doQuestionFive(); //Question 5: DTMF Tones
//        doQuestionSix(); //Question 6: Correlation and Convolution
        doQuestionSeven(); //Question 7: Two-Dimensional FFT
    }

    public static void testFFT() {
        ComplexNumber[] testInputs = generateTestFFTData();
        ComplexNumber[] results = fastFourierTransform(testInputs, 1);
        System.out.println("\nThe fft of the data\n");
        for (int i = 0; i < results.length; i++) {
            System.out.format("%f + %fi\n", results[i].getzReal(), results[i].getzImaginary());
        }

        ComplexNumber[] inverseResults = fastFourierTransform(results, -1);
        System.out.println("\nThe inverse of the data\n");
        for (int i = 0; i < results.length; i++) {
            System.out.format("%f + %fi\n", inverseResults[i].getzReal(), inverseResults[i].getzImaginary());
        }
    }

    public static void testTwoDFFT() {
        System.out.println("Testing 2D FFT...");
        ComplexNumber[][] data = generateTest2DData();
        ComplexNumber[][] results = twoDFFT(data, 1);
        System.out.println("\nReals\n");
        printTwoDArray(results);
    }

    public static void testSignalsDataGenerator() {
        SignalsDataGenerator testGenerator = new SignalsDataGenerator("data/test_output.txt");
        testGenerator.writeResult("Hello world!");
        testGenerator.flushFile();
        testGenerator.clearFile();
        testGenerator.writeResult("Hello world 2!");
        testGenerator.closeFile();

    }

    public static void doCommonSignalsProblem() {
        //Commented out to prevent constant rewriting of the same values to output files, uncomment for functionality
        doFS();
        doGS();
        doPSDEstimates(); //Does the PSD Estimates for f50 and g50, as per question 1b
    }

    private static void doQuestionTwo() {
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
    private static double[] generateV(int frequency) {
        final double INTERVAL = 1.0 / 512.0;
        double[] results = new double[512];
        for (int i = 0; i < results.length; i++) {
            double t = i * INTERVAL;
            double answer = Math.sin(2 * Math.PI * frequency * t);
            results[i] = answer;
        }
        return results;
    }

    //This method takes two functions and adds them together, returning an array of complex numbers
    private static ComplexNumber[] generateXT(double[] v1, double[] v2) {
        ComplexNumber[] results = new ComplexNumber[v1.length];
        for (int i = 0; i < results.length; i++) {
            ComplexNumber answer = new ComplexNumber(v1[i] + v2[i], 0);
            results[i] = answer;
        }
        return results;
    }

    //This method takes two functions and multiplies them together, returning an array of complex numbers
    private static ComplexNumber[] generateYT(double[] v1, double[] v2) {
        ComplexNumber[] results = new ComplexNumber[v1.length];
        for (int i = 0; i < results.length; i++) {
            ComplexNumber answer = new ComplexNumber(v1[i] * v2[i], 0);
            results[i] = answer;
        }
        return results;
    }

    public static void doFS() {
        runFS(3, "data/fs3.txt");
        runFS(10, "data/fs10.txt");
        runFS(50, "data/fs50.txt");
    }

    public static void runFS(int sTerm, String filename) {
        SignalsDataGenerator writer = new SignalsDataGenerator(filename);
        writer.clearFile();
        //Interval for t(time) to increase by
        final double INTERVAL = 1.0 / 512.0;
        final int ITERATIONS = 512;

        for (int i = 1; i <= ITERATIONS; i++) {
            double t = i * INTERVAL;
            double result = 0;
            for (int k = 1; k <= sTerm; k++) {
                result += (Math.sin(2.0 * Math.PI * (2 * k - 1) * t)) / (2 * k - 1);
            }
            String toWrite = Double.toString(result);
            writer.writeResult(toWrite);
        }

        writer.closeFile();
    }

    public static void doGS() {
        runGS(3, "data/gs3.txt");
        runGS(10, "data/gs10.txt");
        runGS(50, "data/gs50.txt");
    }

    public static void runGS(int sTerm, String filename) {
        SignalsDataGenerator writer = new SignalsDataGenerator(filename);
        writer.clearFile();
        //Interval for t(time) to increase by
        final double INTERVAL = 1.0 / 512.0;
        final int ITERATIONS = 512;

        for (int i = 1; i <= ITERATIONS; i++) {
            double t = i * INTERVAL;
            double result = 0;
            for (int k = 1; k <= sTerm; k++) {
                result += (Math.sin(2.0 * Math.PI * (2 * k) * t)) / (2 * k);
            }
            String toWrite = Double.toString(result);
            writer.writeResult(toWrite);
        }

        writer.closeFile();
    }

    //Question 1B: Displays only the positive frequencies
    private static void doPSDEstimates() {
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
        doPSD(GL_IN_FILENAME, GL_OUT_FILENAME);


    }

    private static void doPSD(String inputDataFilename, String outputDataFilename) {
        ComplexNumber[] inputNumbers = parseInputs(inputDataFilename);
        ComplexNumber[] fftResults = fastFourierTransform(inputNumbers, 1);
        double[] psdResults = convertFFTToPSD(fftResults);
        outputPSDResults(psdResults, outputDataFilename);
    }

    //Used if we already have the ComplexNumber[] array from some other source
    private static void doPSD(ComplexNumber[] inputNumbers, String outputDataFilename) {
        ComplexNumber[] fftResults = fastFourierTransform(inputNumbers, 1);
        double[] psdResults = convertFFTToPSD(fftResults);
        outputPSDResults(psdResults, outputDataFilename);
    }

    //Assumes a file where each line is a single number, to be converted into a complex number.
    //Assumes inputs are real numbers (non-complex or imaginary)
    private static ComplexNumber[] parseInputs(String filename) {
        ComplexNumber[] numbers = new ComplexNumber[512];
        int count = 0;
        try {
            File file = new File(filename);
            BufferedReader br = new BufferedReader(new FileReader(file));

            String line;
            while ((line = br.readLine()) != null) {
                if (!line.isEmpty()) {
                    double value = Double.parseDouble(line);
                    ComplexNumber number = new ComplexNumber(value, 0);
                    numbers[count] = number;
                    count++;
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        if (count != numbers.length) {
            System.out.println("Warning: parseInputs() did not parse the expected amount of numbers. Parsed a total of " + Integer.toString(count) + " numbers");
        }
        return numbers;
    }

    private static double[] convertFFTToPSD(ComplexNumber[] fftValues) {
        //To convert FFT to PSD, you take the modulus and square it
        double[] results = new double[fftValues.length];
        for (int i = 0; i < results.length; i++) {
            double answer = Math.pow(fftValues[i].getModulus(), 2);
            results[i] = answer;
        }
        return results;
    }

    //Only outputs the positive (first half) of PSD results, due to symmetrical nature
    private static void outputPSDResults(double[] values, String filename) {
        SignalsDataGenerator writer = new SignalsDataGenerator(filename);
        writer.clearFile();
        for (int i = 0; i < values.length; i++) {
            if (values[i] >= 0) {
                String toWrite;
                if (values[i] < .0001 && values[i] > -0.0001) {
                    toWrite = "0";
                } else {
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
    public static ComplexNumber[] fastFourierTransform(ComplexNumber[] testInputs, int d) {
        int k;
        ComplexNumber t;
        ComplexNumber[] inputs = new ComplexNumber[testInputs.length];
        System.arraycopy(testInputs, 0, inputs, 0, testInputs.length);

        //1) Set theta = (-2*pi*d)/N and r = N/2
        double N = inputs.length;
        double theta = (-2 * Math.PI * d) / N;
        int r = (int) Math.floor(N / 2);
        //2)For i = 1 to N-1:
        for (int i = 1; i < N; i++) {
            //2a: Set w = cos(i*theta) + jsin(i*theta)
            double real = Math.cos(i * theta);
            double imaginary = Math.sin(i * theta);
            ComplexNumber w = new ComplexNumber(real, imaginary);

            //2b: For k = 0 to N-1 do:
            for (k = 0; k < N; k++) {
                //b-1: Set u=1
                ComplexNumber u = new ComplexNumber(1, 0);
                //b-2: For m = 0 to r-1 do:
                for (int m = 0; m < r; m++) {
                    t = inputs[k + m].subtract(inputs[(k + m + r)]);
                    inputs[k + m] = inputs[k + m].add(inputs[k + m + r]);
                    inputs[k + m + r] = t.multiply(u);
                    u = w.multiply(u);
                }
                //b-3: set k = k+2r
                k = k + 2 * r - 1; //-1 due to loop incrementation
            }
            //2c: Set i = 2i and r = r/2
            i = 2 * i - 1; //-1 due to loop incrementation
            r = r / 2;
        }
        //3) For i = 0 to N-1 do:
        for (int i = 0; i < N; i++) {
            //3a: Set r=i and k = 0
            r = i;
            k = 0;
            //3b: For m = 1 to N-1 do:
            int m = 1;
            while (m < N) {
                k = 2 * k + (r % 2);
                r = r / 2;
                m = 2 * m;
            }
            //3c: If k > i do:
            if (k > i) {
                t = inputs[i];
                inputs[i] = inputs[k];
                inputs[k] = t;
            }
        }
        //4) If d < 0 then for i = 0 to N-1 do:
        if (d < 0) {
            for (int i = 0; i < N; i++) {
                inputs[i] = inputs[i].divide(N);
            }
        }
        return inputs;
    }

    public static ComplexNumber[][] twoDFFT(ComplexNumber[][] userInputs, int d) {
        //Copy array
        ComplexNumber[][] inputs = new ComplexNumber[userInputs.length][userInputs[0].length];
        for (int i = 0; i < inputs.length; i++) {
            System.arraycopy(userInputs[i], 0, inputs[i], 0, inputs[i].length);
        }

        int STANDARD_FFT = d;
        ComplexNumber[][] temp = new ComplexNumber[inputs.length][inputs[0].length];
        ComplexNumber finalResult[][] = new ComplexNumber[inputs.length][inputs[0].length];
        for (int k = 0; k < inputs.length; k++) {
            //Do rows
            ComplexNumber[] row = new ComplexNumber[inputs[k].length];
            System.arraycopy(inputs[k], 0, row, 0, row.length);
            ComplexNumber[] rowResult = fastFourierTransform(row, STANDARD_FFT);
            temp[k] = rowResult;
        }

//        System.out.println("First pass:");
//        printTwoDArray(temp);
        //Do Columns
        for (int n = 0; n < inputs[0].length; n++) {
            ComplexNumber[] column = new ComplexNumber[inputs.length];
            for (int i = 0; i < column.length; i++) {
                //add values to form the column
                column[i] = temp[i][n];
            }
            ComplexNumber[] columnResult = fastFourierTransform(column, STANDARD_FFT);
            //Store as column in final result
            for (int i = 0; i < column.length; i++) {
                //add values to form the column
                finalResult[i][n] = columnResult[i];
            }
        }
        return finalResult;
    }

    private static void printTwoDArray(ComplexNumber[][] array) {
        for (int i = 0; i < array.length; i++) {
            for (int j = 0; j < array[0].length; j++) {
                System.out.format("%.2f ", array[i][j].getzReal());
                if (j % 7 == 0 && j != 0) {
                    System.out.format("\n");
                }
            }
        }
    }

    private static ComplexNumber[] generateTestFFTData() {
        double[] reals = {26160.0, 19011.0, 18757, 18405, 17888, 14720, 14285,
                17018, 18014, 17119, 16400, 17497, 17846, 15700, 17636, 17181};
        ComplexNumber[] data = new ComplexNumber[reals.length];
        for (int i = 0; i < reals.length; i++) {
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
            for (int j = 0; j < reals[0].length; j++) {
                ComplexNumber entry = new ComplexNumber(reals[i][j], 0);
                data[i][j] = entry;
            }
        }
        return data;
    }

    private static void doQuestionThree() {
//        examineTimeVariations(); //Question 3a
        examineHT(); //Question 3b
    }

    private static void examineTimeVariations() {
        ComplexNumber[] start = generateSinglePulse(0, 256);
        ComplexNumber[] mid = generateSinglePulse(127, 256);
        ComplexNumber[] random = generateSinglePulse(38, 256);

        doPSD(start, "data/question3/start.txt");
        doPSD(mid, "data/question3/mid.txt");
        doPSD(random, "data/question3/random.txt");
    }

    private static ComplexNumber[] generateSinglePulse(int pulseLocation, int length) {
        ComplexNumber[] results = new ComplexNumber[length];
        for (int i = 0; i < length; i++) {
            double value = 0.0;
            if (i == pulseLocation) {
                value = 1.0;
            }
            results[i] = new ComplexNumber(value, 0);
        }
        return results;
    }

    private static void examineHT() {
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
    private static ComplexNumber[] generateHT(double c, String filename) {
        SignalsDataGenerator writer = new SignalsDataGenerator(filename);
        ComplexNumber[] results = new ComplexNumber[512];
        final double INTERVAL = 1.0 / 512.0;
        for (int i = 0; i < results.length; i++) {
            double t = i * INTERVAL;
            ComplexNumber answer = new ComplexNumber(Math.sin(20.0 * Math.PI * (t - c)), 0);
            results[i] = answer;
            writer.writeResult(Double.toString(answer.getzReal()));
        }
        writer.closeFile();
        return results;
    }

    //purpose: This method is used to output a ComplexNumber to a file. For use in problem 3.b.i
    private static void outputFFTRealsToFile(ComplexNumber[] values, String filename) {
        SignalsDataGenerator writer = new SignalsDataGenerator(filename);
        writer.clearFile();
        for (int i = 0; i < values.length; i++) {
            String real;
            if (values[i].getzReal() > -.0000000001 && values[i].getzReal() < .0000000001) {
                real = "0";
            } else {
                real = Double.toString(values[i].getzReal());
            }
            writer.writeResult(real);
        }
        writer.closeFile();
    }

    private static void doQuestionFour() {
        double[] fs50 = parseFS50();
        ComplexNumber[] fs50_complex = new ComplexNumber[fs50.length];
        for (int i = 0; i < fs50_complex.length; i++) {
            fs50_complex[i] = new ComplexNumber(fs50[i], 0);
        }
        ComplexNumber[] fs50_fft = fastFourierTransform(fs50_complex, 1);
        ComplexNumber[] lowPassFS = lowPassFilter(fs50_fft, 7);
        ComplexNumber[] highPassFS = highPassFilter(fs50_fft, 43);
        ComplexNumber[] bandPassFS = bandPassFilter(fs50_fft, 4, 7); //inclusive
        ComplexNumber[] notchFS = notchFilter(fs50_fft, 4, 7); //inclusive, SUPPRESSED

        //Convert signals back to time domain
        double[] fs50_time = convertFilteredSignalToTimeDomain(fs50_fft);
        double[] lowPass_time = convertFilteredSignalToTimeDomain(lowPassFS);
        double[] highPass_time = convertFilteredSignalToTimeDomain(highPassFS);
        double[] bandPass_time = convertFilteredSignalToTimeDomain(bandPassFS);
        double[] notch_time = convertFilteredSignalToTimeDomain(notchFS);

        System.out.println("fs50 after reconversion:");
        for (int i = 0; i < fs50.length / 2; i++) {
            System.out.format("%6.10f\n", fs50_time[i]);
        }

        //output to files for graphing in Excel
        outputQuestion4(lowPass_time, "data/question4/lp_filter.txt");
        outputQuestion4(highPass_time, "data/question4/hp_filter.txt");
        outputQuestion4(bandPass_time, "data/question4/bp_filter.txt");
        outputQuestion4(notch_time, "data/question4/notch_filter.txt");
    }

    //Takes the data in fs50PSD.txt and converts it into a double[]
    private static double[] parseFS50() {
        final String FS50_FILENAME = "data/fs50.txt";
        double[] results = new double[512];
        try {
            File file = new File(FS50_FILENAME);
            BufferedReader br = new BufferedReader(new FileReader(file));
            String line;
            for (int i = 0; i < results.length; i++) {
                if ((line = br.readLine()) != null && !line.isEmpty()) {
                    results[i] = Double.parseDouble(line);
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return results;
    }

    //Puts a low-pass filter
    //int numValuesToPass - passes only this many through
    private static ComplexNumber[] lowPassFilter(ComplexNumber[] originalInput, int numValuesToPass) {
        ComplexNumber[] filter = new ComplexNumber[originalInput.length];

        //Create the filter
        for (int i = 0; i < filter.length; i++) {
            double filterValue;
            ComplexNumber currentValue = originalInput[i];
            if (!realIsZero(currentValue)) {
                //If the value is non-zero:
                if (numValuesToPass > 0) {
                    //if we need to pass more frequencies
                    filterValue = 1.0;
                    numValuesToPass--;
                } else {
                    //We don't need to pass more frequencies
                    filterValue = 0;
                }
            } else {
                //The value is a zero
                if (numValuesToPass > 0) {
                    //We are currently applying filters of 1
                    filterValue = 1;
                } else {
                    //We are applying filters of 0
                    filterValue = 0;
                }
            }
            filter[i] = new ComplexNumber(filterValue, 0);
        }

        ComplexNumber[] results = new ComplexNumber[filter.length];
        for (int i = 0; i < filter.length; i++) {
            results[i] = originalInput[i].multiply(filter[i]);
        }

        return results;
    }

    //Passes the upper 'numValuesToPass' number of values
    private static ComplexNumber[] highPassFilter(ComplexNumber[] originalInput, int numValuesToPass) {
        ComplexNumber[] filter = new ComplexNumber[originalInput.length];
        for (int i = filter.length - 1; i >= 0; i--) {
            double filterValue;
            ComplexNumber currentValue = originalInput[i];
            if (!realIsZero(currentValue)) {
                //If the value is non-zero:
                if (numValuesToPass > 0) {
                    //if we need to pass more frequencies
                    filterValue = 1.0;
                    numValuesToPass--;
                } else {
                    //We don't need to pass more frequencies
                    filterValue = 0;
                }
            } else {
                //The value is a zero
                if (numValuesToPass > 0) {
                    //We are currently applying filters of 1
                    filterValue = 1;
                } else {
                    //We are applying filters of 0
                    filterValue = 0;
                }
            }
            filter[i] = new ComplexNumber(filterValue, 0);
        }
        ComplexNumber[] results = new ComplexNumber[filter.length];
        for (int i = 0; i < filter.length; i++) {
            results[i] = originalInput[i].multiply(filter[i]);
        }

        return results;
    }


    private static ComplexNumber[] bandPassFilter(ComplexNumber[] originalInput, int startingIndexToPass, int endingIndexToPassInclusive) {
        ComplexNumber[] filter = new ComplexNumber[originalInput.length];
        int non_zero_index = 0;
        for (int i = 0; i < filter.length; i++) {
            double filterValue;
            ComplexNumber currentValue = originalInput[i];
            if (!realIsZero(currentValue)) {
                //Value is non-zero
                non_zero_index++;
            }
            filterValue = inRange(startingIndexToPass, endingIndexToPassInclusive, non_zero_index);
            filter[i] = new ComplexNumber(filterValue, 0);
        }
        ComplexNumber[] results = new ComplexNumber[filter.length];
        for (int i = 0; i < filter.length; i++) {
            results[i] = originalInput[i].multiply(filter[i]);
        }

        return results;
    }

    private static double inRange(int startingIndexToPass, int endingIndexToPassInclusive, int non_zero_index) {
        if (non_zero_index >= startingIndexToPass && non_zero_index <= endingIndexToPassInclusive) {
            //pass this value
            return 1.0;
        } else {
            //zero out this value
            return 0.0;
        }
    }

    private static ComplexNumber[] notchFilter(ComplexNumber[] originalInput, int startingIndexToSuppress, int endIndexToSuppress) {
        ComplexNumber[] filter = new ComplexNumber[originalInput.length];
        int non_zero_index = 0;
        for (int i = 0; i < filter.length; i++) {
            double filterValue;
            ComplexNumber currentValue = originalInput[i];
            if (!realIsZero(currentValue)) {
                //Value is non-zero
                non_zero_index++;
            }
            filterValue = inRange(startingIndexToSuppress, endIndexToSuppress, non_zero_index);

            //Flips value so we suppress instead of passing
            if (filterValue == 1.0) {
                filterValue = 0;
            } else {
                filterValue = 1.0;
            }
            filter[i] = new ComplexNumber(filterValue, 0);
        }
        ComplexNumber[] results = new ComplexNumber[filter.length];
        for (int i = 0; i < filter.length; i++) {
            results[i] = originalInput[i].multiply(filter[i]);
        }

        return results;
    }

    private static void outputComplexToFile(ComplexNumber[] values, String filepath) {
        SignalsDataGenerator writer = new SignalsDataGenerator(filepath);
        writer.clearFile();
        String toWrite;
        for (int i = 0; i < values.length; i++) {
            String zReal = Double.toString(values[i].getzReal());
            String zImg = Double.toString(values[i].getzImaginary());
            toWrite = zReal + " + " + zImg + "i";
            writer.writeResult(toWrite);
        }
        writer.closeFile();
    }

    private static boolean realIsZero(ComplexNumber currentValue) {
        if (currentValue.getzReal() > -.000001 && currentValue.getzReal() < .000001) {
            return true;
        } else return false;
    }

    //Takes a signal and applies an inverted FFT, and returns the real numbers of that value
    private static double[] convertFilteredSignalToTimeDomain(ComplexNumber[] inputs) {
        double[] results = new double[inputs.length];

        //Take filtered signal and apply FFT^-1
        ComplexNumber[] invertedSignal = fastFourierTransform(inputs, -1);

        //Return reals
        for (int i = 0; i < results.length; i++) {
            results[i] = invertedSignal[i].getzReal();
        }

        return results;
    }

    private static void outputQuestion4(double[] output, String filepath) {
        SignalsDataGenerator writer = new SignalsDataGenerator(filepath);
        writer.clearFile();
        for (int i = 0; i < output.length; i++) {
            writer.writeResult(Double.toString(output[i]));
        }
        writer.closeFile();
    }

    private static void doQuestionFive() {
        double[] a1 = parseToneDataFile("data/question5/toneDataA1.txt");
        double[] b1 = parseToneDataFile("data/question5/toneDataB1.txt");

        ComplexNumber[] a1_fft = doubleToFFT(a1);
        ComplexNumber[] b1_fft = doubleToFFT(b1);

        outputComplexToFile(a1_fft, "data/question5/a1_fft.txt");
        outputComplexToFile(b1_fft, "data/question5/b1_fft.txt");

        double[] a1_psd = convertFFTToPSD(a1_fft);
        double[] b1_psd = convertFFTToPSD(b1_fft);

        outputPSDResults(a1_psd, "data/question5/a1_psd.txt");
        outputPSDResults(b1_psd, "data/question5/b1_psd.txt");

        String a1_DTMF = mapCorrespondingFrequencies(a1_psd);
        String b1_DTMF = mapCorrespondingFrequencies(b1_psd);

        System.out.println("a1 corresponds to: " + a1_DTMF);
        System.out.println("b1 corresponds to: " + b1_DTMF);
    }

    private static double[] parseToneDataFile(String filepath) {
        double[] results = new double[4096];
        try {
            BufferedReader br = new BufferedReader(new FileReader(filepath));
            for (int i = 0; i < results.length; i++) {
                String line = br.readLine();
                if (line != null) {
                    results[i] = Double.parseDouble(line);
                } else throw new Exception("null line in parseToneDataFile()");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return results;
    }

    //Takes a double input as a real value and creates a complex number, converts it via regular FFT, and returns the result
    private static ComplexNumber[] doubleToFFT(double[] inputs) {
        ComplexNumber[] complexValues = new ComplexNumber[inputs.length];
        for (int i = 0; i < inputs.length; i++) {
            complexValues[i] = new ComplexNumber(inputs[i], 0);
        }
        return fastFourierTransform(complexValues, 1);
    }

    //Takes an psd and returns the corresponding two DTMF tones as Strings
    private static String mapCorrespondingFrequencies(double[] psdInputs) {

        int[] peakIndices = findTwoMaximae(psdInputs); //returns 66 and 153
        System.out.printf("Peak Indices: %d and %d\n", peakIndices[0], peakIndices[1]);
        int[] noisyValues = calculateNoisyValues(peakIndices); //Returns 710 and 1647
        System.out.printf("Noisy Values: %d and %d\n", noisyValues[0], noisyValues[1]);
        int[] matrixIndices = calculateMatrixIndices(noisyValues); //Returns 0,4
        System.out.printf("Matrix Indices: %d and %d\n", matrixIndices[0], matrixIndices[1]);
        String DTMF_Value = calculateDTMF_Values(matrixIndices); //Returns A and 4
        System.out.printf("DTMF Value: %s\n", DTMF_Value);


        return DTMF_Value;
    }

    //Finds the two largest values in the array of doubles, and returns their indices
    private static int[] findTwoMaximae(double[] psdInputs) {
        int[] largestIndices = new int[2];
        double[] largestValues = new double[2];

        //Initialization loop
        for (int i = 0; i < largestIndices.length; i++) {
            largestIndices[i] = -1;
            largestValues[i] = -1.0;
        }

        //Note, this is written as O(n^2). It could probably be done much faster, but as this program only needs to run one time it is fine
        //for now.
        for (int i = 0; i < psdInputs.length / 2; i++) {
            if (psdInputs[i] >= largestValues[0]) {
                largestValues[0] = psdInputs[i];
                largestIndices[0] = i;
            }
        }

        //Next, find the second-largest index, while making sure to avoid getting nearby bins from the first largest
        for (int i = 0; i < psdInputs.length / 2; i++) {
            if (i < largestIndices[0] - 10 || i > largestIndices[0] + 10) {
                if (psdInputs[i] >= largestValues[1]) {
                    largestValues[1] = psdInputs[i];
                    largestIndices[1] = i;
                }
            }
        }

        if (largestIndices[0] > largestIndices[1]) {
            int temp = largestIndices[0];
            largestIndices[0] = largestIndices[1];
            largestIndices[1] = temp;
        }
        return largestIndices;
    }

    //Applies the formula fk = (k*fs)/N with hard-coded values for values known at design-time
    //Returns - The two most dominant frequencies in Hz
    private static int[] calculateNoisyValues(int[] indices) {
        int[] answers = new int[2];

        for (int i = 0; i < answers.length; i++) {
            answers[i] = (indices[i] * 44100) / 4096;
        }
        return answers;
    }

    private static int[] calculateMatrixIndices(int[] noisyValues) {
        int[] indices = new int[noisyValues.length];

        for (int j = 0; j < noisyValues.length; j++) {
            int[] rows = new int[]{697, 770, 852, 941};
            int[] columns = new int[]{1209, 1336, 1477, 1633};
            for (int i = 0; i < rows.length; i++) {
                if (j == 0) {
                    rows[i] = Math.abs(rows[i] - noisyValues[j]);
                } else {
                    columns[i] = Math.abs(columns[i] - noisyValues[j]);
                }
            }

            int currentSmallestValue = 10000000;
            for (int i = 0; i < rows.length; i++) {
                if (j == 0) {
                    if (rows[i] < currentSmallestValue) {
                        currentSmallestValue = rows[i];
                        indices[j] = i;
                    }
                } else {
                    if (columns[i] < currentSmallestValue) {
                        currentSmallestValue = columns[i];
                        indices[j] = i;
                    }
                }
            }
        }
        return indices;
    }

    private static String calculateDTMF_Values(int[] matrix_indices) {
        String[][] DTMF_matrix = new String[][]{
                {"1", "2", "3", "A"},
                {"4", "5", "6", "B"},
                {"7", "8", "9", "C"},
                {"*", "0", "#", "D"}
        };
        return DTMF_matrix[matrix_indices[0]][matrix_indices[1]];
    }

    private static void doQuestionSix() {
        do6A(); //6A: Correlation
        do6B(); //6B: Convolution
    }

    private static void do6A() {
        double[] response_signal = parseResponseSignal("data/question6/response_samples.txt");
        double[] pulse_signal = parsePulseSignal("data/question6/pulse_samples.txt");
        double[] padded_pulse_signal = padPulseSignal(pulse_signal, response_signal.length);
        double[] real_corr = calculateCorrelation(response_signal, padded_pulse_signal);
        int max_corr_index = findMaximumCorrelation(real_corr);
        double[] corrected_pulse_signal = bufferPulseSignal(pulse_signal, max_corr_index, real_corr.length);
        outputRealsToFile(corrected_pulse_signal, "data/question6/corrected_pulse_signal.txt");
        outputRealsToFile(real_corr, "data/question6/correlation.txt");
        print6AResults(max_corr_index);
    }

    private static double[] parseResponseSignal(String filepath) {
        double[] response_signal = parseDoublesFromFile(filepath, 1024);
        return response_signal;
    }

    private static double[] parsePulseSignal(String filepath) {
        return parseDoublesFromFile(filepath, 52);
    }

    private static double[] parseDoublesFromFile(String filepath, int numDoubles) {
        double[] values = new double[numDoubles];
        String line;
        try {
            BufferedReader br = new BufferedReader(new FileReader(filepath));
            for (int i = 0; i < values.length; i++) {
                line = br.readLine();
                values[i] = Double.parseDouble(line);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return values;
    }

    private static double[] padPulseSignal(double[] originalSignal, int finalSize) {
        double[] paddedSignal = new double[finalSize];
        int numZerosNeeded = finalSize - originalSignal.length;
        for (int i = 0; i < finalSize; i++) {
            if (numZerosNeeded > 0) {
                if (i < originalSignal.length) {
                    paddedSignal[i] = originalSignal[i];
                } else {
                    paddedSignal[i] = 0.0;
                    numZerosNeeded--;
                }

            }
        }
        return paddedSignal;
    }

    //fastCorrelation is calculated with the following formula:
    //FFT^-1[FFT(x) * ComplexConjugate(FFT(p))]
    //Imaginary values of this calculation can be ignored (and are)
    private static double[] calculateCorrelation(double[] responseSignal, double[] paddedPulseSignal) {
        double[] corr = new double[responseSignal.length];

        ComplexNumber[] X = fastFourierTransform(ComplexNumber.generateFromDoubles(responseSignal), 1);
        ComplexNumber[] P = fastFourierTransform(ComplexNumber.generateFromDoubles(paddedPulseSignal), 1);
        ComplexNumber[] PConjugate = new ComplexNumber[P.length];
        for (int i = 0; i < PConjugate.length; i++) {
            PConjugate[i] = P[i].conjugate();
        }
        ComplexNumber[] XP = ComplexNumber.massMultiply(X, PConjugate);
        ComplexNumber[] XPInverse = fastFourierTransform(XP, -1);

        for (int i = 0; i < corr.length; i++) {
            corr[i] = XPInverse[i].getzReal();
        }
        return corr;
    }

    private static int findMaximumCorrelation(double[] correlationValues) {
        int indexOfMax = -1;
        double valueOfMax = -99999999.9;

        for (int i = 0; i < correlationValues.length; i++) {
            double currentValue = correlationValues[i];
            if (currentValue > valueOfMax) {
                indexOfMax = i;
                valueOfMax = currentValue;
            }
        }
        return indexOfMax;
    }

    //Creates an array of doubles that is 0 up to maxCorrIndex, the values of nonPaddedPulseSignal, followed by 0s to the maximum
    private static double[] bufferPulseSignal(double[] nonPaddedPulseSignal, int maxCorrIndex, int finalSize) {
        double[] correctedPulseSignal = new double[finalSize];
        int transferIndex = 0; //To keep track of where we are in nonPaddedPulseSignal
        //First, pad the beginning
        for (int i = 0; i < correctedPulseSignal.length; i++) {
            if (i >= maxCorrIndex && i < maxCorrIndex + nonPaddedPulseSignal.length) {
                //add original values
                correctedPulseSignal[i] = nonPaddedPulseSignal[transferIndex];
                transferIndex++;
            } else {
                correctedPulseSignal[i] = 0.0;
            }
        }
        return correctedPulseSignal;
    }

    private static void outputRealsToFile(double[] reals, String filepath) {
        SignalsDataGenerator writer = new SignalsDataGenerator(filepath);
        writer.clearFile();
        for (double real : reals) {
            writer.writeResult(Double.toString(real));
        }
        writer.closeFile();
    }

    private static void print6AResults(int maxCorrIndex) {
        final int VELOCITY_IN_SEA_WATER = 1500;
        final double DT = 1.0 / 50000.0;
        final double TWO_MILLISECONDS = 0.002;

        double answer = (maxCorrIndex / 2.0 * DT + TWO_MILLISECONDS) * VELOCITY_IN_SEA_WATER;

        System.out.printf("Conclusion to 6A:\nThe signal took %d intervals to return.\nWith a %f delay, " +
                "and an interval of %f, this results in:\n" +
                "%.3f meters away\n", maxCorrIndex, TWO_MILLISECONDS, DT, answer);
    }

    private static void do6B(){
        double[] returnedSignal = parseResponseSignal("data/question6/response_samples.txt");
        double[] smoothedSignal = calculateConvolution(returnedSignal);
        outputRealsToFile(smoothedSignal, "data/question6/smoothed_response_signal.txt");
    }

    //Formula: FFT^-1[FFT(pulse)*FFT(input)]
    private static double[] calculateConvolution(double[] returnedSignal){
        double[] pulse = generateNPointPulse(6, returnedSignal.length);
        ComplexNumber[] returnedSignalComplex = ComplexNumber.generateFromDoubles(returnedSignal);
        ComplexNumber[] pulseComplex = ComplexNumber.generateFromDoubles(pulse);
        ComplexNumber[] returned_signal_fft = fastFourierTransform(returnedSignalComplex, 1);
        ComplexNumber[] pulse_fft = fastFourierTransform(pulseComplex, 1);

        ComplexNumber[] multiplied = ComplexNumber.massMultiply(pulse_fft, returned_signal_fft);

        ComplexNumber[] convolutionComplex = fastFourierTransform(multiplied, -1);

        double[] convolution_reals = new double[convolutionComplex.length];
        for(int i = 0; i < convolution_reals.length; i++){
            convolution_reals[i] = convolutionComplex[i].getzReal();
        }
        return convolution_reals;
    }

    private static double[] generateNPointPulse(int numPoints, int pulseSize){
        double[] pulse = new double[pulseSize];

        for(int i = 0; i < pulseSize; i++){
            if(i < numPoints){
                pulse[i] = 1.0/numPoints;
            } else{
                pulse[i] = 0.0;
            }
        }
        return pulse;
    }

    private static void doQuestionSeven(){
        Picture returnSignal = setupReturnPicture();
        Picture pulseSignal = setupPulsePicture();

        returnSignal.save(new File("data/question7/returnSignal.jpg"));
        pulseSignal.save(new File("data/question7/pulseSignal.jpg"));

        double[][] correlation = do2dCorrelation(returnSignal, pulseSignal);
        double[][] scaledCorrelation = scaleCorrelation(correlation);

        output2dRealsToFile(correlation, "data/question7/correlation.txt");
        output2dRealsToFile(scaledCorrelation, "data/question7/scaled_correlation.txt");

        Picture correlationPicture = outputCorrelationPicture(scaledCorrelation);
        correlationPicture.save(new File("data/question7/corrPicture.jpg"));

        double[][] magnitudes = correlationMagnitude(correlationPicture);
        Picture magnitudePicture = outputMagnitudePicture(magnitudes);
        magnitudePicture.save(new File("data/question7/correlation_magnitude.jpg"));
    }

    private static Picture setupReturnPicture(){
        return setupPicture(180, 220, 110, 140, 30, 90);
    }

    private static Picture setupPulsePicture(){
        return setupPicture(0, 0, 30, 120, 15, 90);
    }

    private static Picture setupPicture(int whiteXStart, int whiteYStart, int whiteWidth, int whiteHeight, int blackWidth, int blackHeight){
        Picture picture = new Picture(512, 512);
        setToBlack(picture);
        createBigRectangle(picture, whiteXStart, whiteYStart, whiteWidth, whiteHeight);
        int fullWidth = whiteXStart+whiteWidth;
        int top_c_bar_end = whiteYStart + ((whiteHeight - blackHeight) / 2);
        int bottom_c_bar_start = whiteYStart + whiteHeight - ((whiteHeight-blackHeight)/2);
        blackoutRectanglePortion(picture, fullWidth, top_c_bar_end, bottom_c_bar_start, blackWidth);

        return picture;
    }

    private static void setToBlack(Picture picture){
        for(int i = 0; i < picture.width(); i++){
            for(int j = 0; j < picture.height(); j++){
                picture.set(i, j, Color.BLACK);
            }
        }
    }

    private static void createBigRectangle(Picture picture, int xStart, int yStart, int width, int height){
        for(int i = 0; i < picture.height(); i++){
            for(int j = 0; j < picture.width(); j++){
                if(i >= xStart && i < xStart+width){
                    if(j >= yStart && j < yStart+height){
                        picture.set(i, j, Color.WHITE);
                    }
                }
            }
        }
    }

    private static void blackoutRectanglePortion(Picture picture, int fullWidth, int topCBarEnd, int bottomCBarStart, int blackWidth){
        for(int i = 0; i < picture.width(); i++){
            for(int j = 0; j < picture.height(); j++){
                if(i >= fullWidth - blackWidth && i < fullWidth){
                    if(j >= topCBarEnd && j < bottomCBarStart){
                        picture.set(i, j, Color.BLACK);
                    }
                }
            }
        }
    }

    private static double[][] do2dCorrelation(Picture returnSignal, Picture pulseSignal){
        Color[][] returnColorArray = returnSignal.getColorArray();
        Color[][] pulseColorArray = pulseSignal.getColorArray();

        ComplexNumber[][] returnSignalComplex = colorToComplex(returnColorArray);
        ComplexNumber[][] pulseSignalComplex = colorToComplex(pulseColorArray);

        ComplexNumber[][] return_fft = twoDFFT(returnSignalComplex, 1);
        ComplexNumber[][] pulse_fft = twoDFFT(pulseSignalComplex, 1);
        ComplexNumber[][] pulse_fft_conjugate = new ComplexNumber[pulse_fft.length][pulse_fft[0].length];
        for(int i = 0; i < pulse_fft_conjugate.length; i++){
            for(int j = 0; j < pulse_fft_conjugate[0].length; j++){
                pulse_fft_conjugate[i][j] = pulse_fft[i][j].conjugate();
            }
        }

        ComplexNumber[][] XP = new ComplexNumber[return_fft.length][return_fft[0].length];
        for(int i = 0; i < pulse_fft_conjugate.length; i++){
            for(int j = 0; j < pulse_fft_conjugate[0].length; j++){
                XP[i][j] = return_fft[i][j].multiply(pulse_fft_conjugate[i][j]);
            }
        }

        ComplexNumber[][] corrComplex = twoDFFT(XP, -1);

        double[][] corrDoubles = new double[corrComplex.length][corrComplex[0].length];
        for(int i = 0; i < corrComplex.length; i++){
            for(int j = 0; j < corrComplex[0].length; j++){
                corrDoubles[i][j] = corrComplex[i][j].getzReal();
            }
        }

        return corrDoubles;
    }

    private static ComplexNumber[][] colorToComplex(Color[][] colorArray){
        ComplexNumber[][] results = new ComplexNumber[colorArray.length][colorArray[0].length];
        for(int i = 0; i < results.length; i++){
            for(int j = 0; j < results[0].length; j++){
                results[i][j] = new ComplexNumber(colorArray[i][j].getRed(), 0);
            }
        }
        return results;
    }

    private static double[][] scaleCorrelation(double[][] corr){
        double min = find2dMin(corr);
        min = Math.abs(min);
        double[][] translated_corr = new double[corr.length][corr[0].length];
        for(int i = 0; i < translated_corr.length; i++){
            for(int j = 0; j < translated_corr[0].length; j++){
                translated_corr[i][j] = corr[i][j] + min + 1;
            }
        }

        double[][] log_scaled = new double[translated_corr.length][translated_corr[0].length];

        for(int i = 0; i < log_scaled.length; i++){
            for(int j = 0; j < log_scaled[0].length; j++){
                log_scaled[i][j] = Math.log(translated_corr[i][j]);
            }
        }
        double max = find2dMax(log_scaled);
        min = findNonZeroMin(log_scaled);
        System.out.println("min: "+min);
        System.out.println("max: "+max);

        double[][] translated_to_zero = translateToZero(log_scaled, min);
        max = find2dMax(translated_to_zero);
        min = find2dMin(translated_to_zero);
        System.out.println("min after translation: "+min);
        System.out.println("max after translation: "+max);

        double linearFactor = 255.0 / max;
        System.out.println("Linear Factor: "+ linearFactor);


        double[][] results = new double[translated_to_zero.length][translated_to_zero[0].length];
        for(int i = 0; i < results.length; i++){
            for(int j = 0; j < results[0].length; j++){
                results[i][j] = translated_to_zero[i][j] * linearFactor;
            }
        }

        return results;
    }

    private static double find2dMin(double[][] inputs){
        double lowestValue = 9999999;

        for(int i = 0; i < inputs.length; i++){
            for(int j = 0; j < inputs[0].length; j++){
                if(inputs[i][j] < lowestValue){
                    lowestValue = inputs[i][j];
                }
            }
        }
        return lowestValue;
    }

    private static double find2dMax(double[][] inputs){
        double highestValue = -1;

        for(int i = 0; i < inputs.length; i++){
            for(int j = 0; j < inputs[0].length; j++){
                if(inputs[i][j] > highestValue){
                    highestValue = inputs[i][j];
                }
            }
        }
        return highestValue;
    }

    private static void output2dRealsToFile(double[][] reals, String filepath){
        SignalsDataGenerator writer = new SignalsDataGenerator(filepath);
        writer.clearFile();
        for(int i = 0; i < reals.length; i++){
            for(int j = 0; j < reals[0].length; j++){
                writer.writeResult(Double.toString(reals[i][j]));
            }
        }
        writer.closeFile();
    }

    private static Picture outputCorrelationPicture(double[][] pixelValues){
        Picture picture = new Picture(512, 512);

        for(int i = 0; i < picture.width(); i++){
            for(int j = 0; j < picture.height(); j++){
                int rgb = (int) Math.floor(pixelValues[i][j]);
                if(rgb < 0){
                    rgb = 0;
                }
                Color color = new Color(rgb, rgb, rgb);
                picture.set(j, i, color);
            }
        }

        return picture;
    }

    private static double[][] correlationMagnitude(Picture picture){
        Color[][] inputs = picture.getColorArray();
        ComplexNumber[][] complexes = colorToComplex(inputs);
        return calculateMagnitudes(complexes);
    }

    private static double[][] calculateMagnitudes(ComplexNumber[][] inputs){
        double[][] results = new double[inputs.length][inputs[0].length];
        for(int i = 0; i < results.length; i++){
            for(int j = 0; j < results[0].length; j++){
                results[i][j] = inputs[i][j].getModulus();
            }
        }

        return results;
    }

    private static Picture outputMagnitudePicture(double[][] pixelValues){
        double[][] scaled = scaleCorrelation(pixelValues);
        double[][] scaled_with_red = addRed(scaled);

        Picture picture = new Picture(512, 512);

        for(int i = 0; i < picture.width(); i++){
            for(int j = 0; j < picture.height(); j++){
                int red = (int) Math.floor(scaled_with_red[i][j]);
                int green;
                int blue;
                if(red < 0){
                    red = 0;
                    green = 0;
                    blue = 0;
                } else if(red == 999.0){
                    red = 255;
                    green = 0;
                    blue = 0;
                }
                else{
                    green = red;
                    blue = red;
                }
                Color color = new Color(red, green, blue);
                picture.set(j, i, color);
            }
        }

        return picture;
    }

    private static double[][] addRed(double[][] pixelValues){
        double[][] results = new double[pixelValues.length][pixelValues[0].length];
        double max = find2dMax(pixelValues);
        for(int i = 0; i < pixelValues.length; i++){
            for(int j = 0; j < pixelValues[0].length; j++){
                if(pixelValues[i][j] >= (max * .9)){
                    results[i][j] = 999.0; //Marker for red
                } else {
                    results[i][j] = pixelValues[i][j];
                }
            }
        }
        return results;
    }

    private static double findNonZeroMin(double[][] inputs){
        double lowestValue = 9999999;

        for(int i = 0; i < inputs.length; i++){
            for(int j = 0; j < inputs[0].length; j++){
                if(inputs[i][j] < lowestValue){
                    if(inputs[i][j] < -.00001 || inputs[i][j] > 0.00001) { //check for non-zero
                        lowestValue = inputs[i][j];
                    }
                }
            }
        }
        return lowestValue;
    }

    private static double[][] translateToZero(double[][] inputs, double min){
        double[][] results = new double[inputs.length][inputs[0].length];
        for(int i = 0; i < results.length; i++){
            for(int j = 0; j < results[0].length; j++){
                double currentInput = inputs[i][j];
                if(currentInput > .00001 || currentInput < -.00001){
                    //Is NOT zero
                    results[i][j] = currentInput - min;
                } else{
                    //is zero
                    results[i][j] = 0.0;
                }
            }
        }
        return results;
    }
}