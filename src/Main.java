import java.lang.Math;
public class Main {

    public static void main(String[] args) {
        ComplexNumber[] testInputs = generateTestFFTData();
        ComplexNumber[] results = fastFourierTransform(testInputs, 1);
        System.out.println("\nThe fft of the data\n");
        for(int i = 0; i < results.length; i++){
            System.out.format("%f + %fi\n", results[i].getzReal(), results[i].getzImaginary());
        }
    }

    //Method that applies a Fast Fourier Transform to an Nx1 vector of complex numbers
    //Params: inputs - the Nx1 array
    //        d - the direction code. Either 1 for FFT or -1 for inverse FFT
    public static ComplexNumber[] fastFourierTransform(ComplexNumber[] testInputs, int d){
        int k;
        ComplexNumber t;
        ComplexNumber[] inputs = new ComplexNumber[testInputs.length];
        System.arraycopy(testInputs, 0, inputs, 0, testInputs.length);

        System.out.format("The data as complex\n");
        for(int i = 0; i < inputs.length; i++){
            System.out.format("%.2f + %.2fi\n", inputs[i].getzReal(), inputs[i].getzImaginary());
        }

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
                    System.out.printf("\ni=%d; k=%d; m=%d; r=%d", i,k,m,r);
                }
                //b-3: set k = k+2r
                k = k + 2*r -1; //-1 due to loop incrementation
            }
            //2c: Set i = 2i and r = r/2
            i = 2*i-1; //-1 due to loop incrementation
            r = r/2;
        }

        System.out.printf("\nResults from the first loops\n");
        for(int x = 0; x < inputs.length; x++){
            System.out.format("%f + %fi\n", inputs[x].getzReal(), inputs[x].getzImaginary());
        }
        System.out.format("Rearranging the results\n");
        //3) For i = 0 to N-1 do:
        for(int i=0; i < N; i++){
            //3a: Set r=i and k = 0
            r = i;
            k = 0;
            //3b: For m = 1 to N-1 do:
            int m = 1;
            while(m < N){
                System.out.printf("\ni=%d; k=%d; m=%d; r=%d", i,k,m,r);
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

    private static ComplexNumber[] generateTestFFTData(){
        double[] reals = {26160.0, 19011.0, 18757, 18405, 17888, 14720, 14285,
        17018, 18014, 17119, 16400, 17497, 17846, 15700, 17636, 17181};
        ComplexNumber[] data = new ComplexNumber[reals.length];
        for(int i = 0; i < reals.length; i++){
            data[i] = new ComplexNumber(reals[i], 0);
        }
        return data;
    }
}
