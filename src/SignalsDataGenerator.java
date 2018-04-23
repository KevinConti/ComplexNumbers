import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public class SignalsDataGenerator {
    private String filename;
    private FileWriter fw = null;
    private BufferedWriter bw = null;

    //Error is thrown if filename parameter is not sent correctly
    public SignalsDataGenerator(String filename) {
        this.filename = filename;
        try {
            fw = new FileWriter(filename, true);
            bw = new BufferedWriter(fw);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    //Writes the result with a trailing newline character. Returns true if successfully written
    public boolean writeResult(String toWrite){
        boolean isSuccessful = false;

        try {
            bw.write(toWrite);
            bw.newLine();
            isSuccessful = true;
        } catch (IOException e) {
            e.printStackTrace();
        }

        return isSuccessful;
    }
    public void flushFile(){
        try {
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    //Closes the file. Flushes first for output to appear
    public void closeFile(){
        try {
            bw.close();
            fw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    //Removes all text from the file and resets the bw and fw to append state
    public boolean clearFile(){
        boolean isSuccessful = false;
        try{
            fw = new FileWriter(filename);
            bw = new BufferedWriter(fw);
            bw.write("");
            bw.close();
            fw.close();
            fw = new FileWriter(filename, true);
            bw = new BufferedWriter(fw);
            isSuccessful = true;
        } catch (IOException e) {
            e.printStackTrace();
        }
        return isSuccessful;
    }
}
