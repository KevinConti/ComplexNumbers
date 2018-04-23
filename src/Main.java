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
        doCommonSignalsProblem(); //Problem 1, output in data directory
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

        doPSD(FS50_FILENAME, FSPSD_FILENAME);
        doPSD(GS50_FILENAME, GSPSD_FILENAME);
    }

    private static void doPSD(String inputDataFilename, String outputDataFilename){
        ComplexNumber[] inputNumbers = parseInputs(inputDataFilename);
        ComplexNumber[] fftResults = fastFourierTransform(inputNumbers, 1);
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
        for(int i = 0; i < values.length / 2; i++){
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
}
