import org.junit.Test;

import static org.junit.Assert.*;

public class ComplexNumberTest {

    @org.junit.Test
    public void add() {
        ComplexNumber one = new ComplexNumber(5, 2);
        ComplexNumber two = new ComplexNumber(3, -7);

        ComplexNumber myAnswer = one.add(two);
        ComplexNumber actualAnswer = new ComplexNumber(8, -5);

        assertEquals(myAnswer.getzReal(), actualAnswer.getzReal(), 0);
        assertEquals(myAnswer.getzImaginary(), actualAnswer.getzImaginary(), 0);
    }

    @org.junit.Test
    public void subtract() {
        ComplexNumber one = new ComplexNumber(2, -3);
        ComplexNumber two = new ComplexNumber(6, -18);

        ComplexNumber myAnswer = one.subtract(two);
        ComplexNumber actualAnswer = new ComplexNumber(-4, 15);

        assertEquals(myAnswer.getzReal(), actualAnswer.getzReal(), 0);
        assertEquals(myAnswer.getzImaginary(), actualAnswer.getzImaginary(), 0);
    }

    @org.junit.Test
    public void multiply() {
        //Test scalar multiplication
        ComplexNumber one = new ComplexNumber(13, 5);
        double scalar = -4;

        ComplexNumber myAnswer = one.multiply(scalar);
        ComplexNumber actualAnswer = new ComplexNumber(-52, -20);

        assertEquals(myAnswer.getzReal(), actualAnswer.getzReal(), 0);
        assertEquals(myAnswer.getzImaginary(), actualAnswer.getzImaginary(), 0);

        //test imaginary number * complex number
        double imaginary = 2;
        one = new ComplexNumber(3, -8);

        myAnswer = one.multiply(imaginary, true);
        actualAnswer = new ComplexNumber(16, 6);

        assertEquals(myAnswer.getzReal(), actualAnswer.getzReal(), 0);
        assertEquals(myAnswer.getzImaginary(), actualAnswer.getzImaginary(), 0);

        //Test complex number multiplication
        ComplexNumber three = new ComplexNumber(1, 4);
        ComplexNumber four = new ComplexNumber(5, 1);

        myAnswer = three.multiply(four);
        actualAnswer = new ComplexNumber(1, 21);

        assertEquals(myAnswer.getzReal(), actualAnswer.getzReal(), 0);
        assertEquals(myAnswer.getzImaginary(), actualAnswer.getzImaginary(), 0);
    }

    @org.junit.Test
    public void divide() {
        ComplexNumber one = new ComplexNumber(20, -4);
        ComplexNumber two = new ComplexNumber(3, 2);

        ComplexNumber myAnswer = one.divide(two);
        ComplexNumber actualAnswer = new ComplexNumber(4, -4);

        assertEquals(myAnswer.getzReal(), actualAnswer.getzReal(), 0);
        assertEquals(myAnswer.getzImaginary(), actualAnswer.getzImaginary(), 0);
    }

    @org.junit.Test
    public void getModulus() {
        ComplexNumber test = new ComplexNumber(4, 3);
        assertEquals(5, test.getModulus(), .001);
    }

    @org.junit.Test
    public void getArgument() {
        ComplexNumber test = new ComplexNumber(4, 3);
        assertEquals(36.87, test.getArgument(), 0.001);
    }

    @Test
    public void conjugate() {
        //Test positive
        ComplexNumber test = new ComplexNumber(20, 30);
        assertEquals(-30, test.conjugate().getzImaginary(), 0);

        //Test negative
        test = new ComplexNumber(-20, -25);
        assertEquals(25.0, test.conjugate().getzImaginary(), 0);
        assertEquals(-20, test.conjugate().getzReal(), 0);

        //test zero
        test = new ComplexNumber(0, 0);
        assertEquals(0.0, test.conjugate().getzImaginary(), 0);
    }
}