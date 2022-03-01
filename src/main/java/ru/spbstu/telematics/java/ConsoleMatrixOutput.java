package ru.spbstu.telematics.java;

public class ConsoleMatrixOutput implements MatrixOutput{
    @Override
    public void print(double[][][] matrix) {
        for (int i = 0; i < matrix.length; ++i){
            System.out.println("Слой №" + (i+1));
            for (int j = 0; j < matrix[i].length; ++j) {
                for (int k = 0; k < matrix[i][j].length; ++k) {
                    if (k != matrix.length - 1) {
                        System.out.print("\t"+matrix[i][j][k]);
                    } else {
                        System.out.print("\t"+matrix[i][j][k]);
                        System.out.println("");
                    }
                }
            }
            System.out.println("");

        }
        System.out.println("\n==========\n");
    }
}
