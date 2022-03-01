package ru.spbstu.telematics.java;

import org.junit.Assert;
import org.junit.Test;


public class AppTest
{
    @Test
    public void сonservationTest() throws InterruptedException {
        int temp = 22;
        int iterations = 40000;
        int nodes = 200;
        int parallel = 6;
        HeatEquationSolver heatEquationSolver = new HeatEquationSolver(parallel);
        double[][][] tInit = MatrixGenerator.getRandomMatrix(temp, nodes, nodes,nodes);
        double sum = 0;
        for (int i = 0; i < nodes; i++) {
            for (int j = 0; j < nodes; j++) {
                for (int k=0; k < nodes; k++) {
                    sum += tInit[i][j][k];
                }
            }
        }

        heatEquationSolver.calculate(tInit, nodes, nodes, nodes, iterations, temp);
        double sumCalc = 0;
        for (int i = 0; i < nodes; i++) {
            for (int j = 0; j < nodes; j++) {
                for (int k=0; k < nodes; k++) {
                    sumCalc += tInit[i][j][k];
                }
            }
        }
        MatrixOutput matrixOutput = new ConsoleMatrixOutput();
        matrixOutput.print(tInit);
        Assert.assertEquals(sum/sumCalc, 1, 0.03);
    }
}