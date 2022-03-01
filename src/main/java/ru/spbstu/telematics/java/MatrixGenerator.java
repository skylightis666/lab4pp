package ru.spbstu.telematics.java;

import java.util.Random;

public class MatrixGenerator {
    public static double[][][] getRandomMatrix(int temp, int nX, int nY, int nZ) {
        double[][][] tInit = new double[nX][nY][nZ];
        for (int i = 0; i < nX; i++) {
            for (int j = 0; j < nY; j++) {
                for (int k=0; k<nZ; k++){
                    if (i == 0||i==nX||j == 0||j==nY||k == 0||k==nZ) {
                        tInit[i][j][k] = temp;
                    } else {
                        tInit[i][j][k] = new Random().nextInt(10);
                    }
                }

            }
        }
        return tInit;
    }
}
