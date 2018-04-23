public class ComplexNumber {
    private double zReal;
    private double zImaginary;
    private double modulus;
    private double argument;

    public ComplexNumber(double zReal, double zImaginary) {
        this.zReal = zReal;
        this.zImaginary = zImaginary;
    }

    //Instance Methods
    public ComplexNumber add(ComplexNumber number) {
        return new ComplexNumber(this.getzReal() + number.getzReal(), this.getzImaginary() + number.getzImaginary());
    }

    public ComplexNumber subtract(ComplexNumber number) {
        return new ComplexNumber(this.getzReal() - number.getzReal(), this.getzImaginary() - number.getzImaginary());
    }

    public ComplexNumber multiply(ComplexNumber number) {
        double ac = this.getzReal() * number.getzReal();
        double bd = this.getzImaginary() * number.getzImaginary();
        double ad = this.getzReal() * number.getzImaginary();
        double bc = this.getzImaginary() * number.getzReal();
        return new ComplexNumber((ac - bd), (ad + bc));
    }

    public ComplexNumber multiply(double number, boolean isImaginaryMultiplication){
        ComplexNumber mult = this.multiply(number);
        if(!isImaginaryMultiplication){
            //Then this number is just a scalar, do scalar multiplication
            return mult;
        } else{
            return new ComplexNumber(mult.getzImaginary() * -1, mult.getzReal());
        }
    }

    public ComplexNumber multiply(double scalar){
        return new ComplexNumber(this.getzReal() * scalar, this.getzImaginary() * scalar);
    }

    public ComplexNumber divide(ComplexNumber number) {
        double ac = this.getzReal() * number.getzReal();
        double bd = this.getzImaginary() * number.getzImaginary();
        double ad = this.getzReal() * number.getzImaginary();
        double bc = this.getzImaginary() * number.getzReal();
        double real = (ac + bd) / (Math.pow(number.getzReal(), 2) + Math.pow(number.getzImaginary(), 2));
        double imaginary = (bc - ad) / (Math.pow(number.getzReal(), 2) + Math.pow(number.getzImaginary(), 2));
        return new ComplexNumber(real, imaginary);
    }

    public ComplexNumber divide(double scalar){
        return new ComplexNumber(this.getzReal() / scalar, this.getzImaginary() / scalar);
    }

    public ComplexNumber conjugate(){
        return new ComplexNumber(this.getzReal(), this.getzImaginary() * -1);
    }

    private void updateModulus() {
        this.setModulus(Math.sqrt(Math.pow(this.getzReal(), 2) + Math.pow(this.getzImaginary(), 2)));
    }

    private void updateArgument() {
        double argument = Math.toDegrees(Math.atan(this.getzImaginary() / this.getzReal()));
        if(this.getzImaginary() < 0){
            argument *= -1;
        }
        this.setArgument(argument);
    }

    //GETTERS AND SETTERS

    public double getzReal() {
        return zReal;
    }

    public void setzReal(double zReal) {
        this.zReal = zReal;
    }

    public double getzImaginary() {
        return zImaginary;
    }

    public void setzImaginary(double zImaginary) {
        this.zImaginary = zImaginary;
    }

    public double getModulus() {
        updateModulus();
        return modulus;
    }

    public void setModulus(double modulus) {
        this.modulus = modulus;
    }

    public double getArgument() {
        updateArgument();
        return argument;
    }

    public void setArgument(double argument) {
        this.argument = argument;
    }

    @Override
    public String toString() {
        return "ComplexNumber{" +
                "zReal=" + zReal +
                ", zImaginary=" + zImaginary + "i" +
                ", modulus=" + getModulus() +
                ", argument=" + getArgument() +
                '}';
    }
}
