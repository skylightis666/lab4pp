package ru.spbstu.telematics.java;

import java.time.Duration;
import java.time.Instant;

public class Main   {
    public static void main(String[] args) throws InterruptedException {
        int temp = 22;
        int iterations = 40000;
        int nodes = 10;
        int maxParallelism = 3;


        double[][][] TInit = MatrixGenerator.getRandomMatrix(temp, nodes, nodes, nodes);

        for (int i = 1; i < maxParallelism; i++) {
            HeatEquationSolver heatEquationSolver = new HeatEquationSolver(i);
            ConsoleMatrixOutput consoleMatrixOutput = new ConsoleMatrixOutput();
            System.out.println("Number of nodes:" + nodes);
            System.out.println("Number of threads:" + i);
            Instant start = Instant.now();
            heatEquationSolver.calculate(TInit, nodes, nodes, nodes, iterations, temp);
            Instant finish = Instant.now();
            System.out.printf("Working time = %d milliseconds%n", Duration.between(start, finish).toMillis());
        }
    }
}
