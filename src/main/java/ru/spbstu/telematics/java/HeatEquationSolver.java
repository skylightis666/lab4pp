package ru.spbstu.telematics.java;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class HeatEquationSolver {
    private final ExecutorService executorService;

    public HeatEquationSolver(int paralel) {
        executorService = Executors.newFixedThreadPool(paralel);

    }

    void step(double[][][] t0,
              double[][][] t1,
              int sizeX, int sizeY, int sizeZ,
              double dx,
              double dy,
              double dz,
              double dt,
              int T,
              double h, int k) throws InterruptedException {
        CountDownLatch countDownLatch = new CountDownLatch(sizeX);
        for (int i = 0; i < sizeX; ++i) {
            final int ii = i;
            executorService.submit(() -> {
                keyStep(t0, t1, sizeX, sizeY, sizeZ,
                        dx, dy, dz, dt, ii
                        ,k,T-5, h);
                countDownLatch.countDown();
            });
        }
        countDownLatch.await();
    }

    private void keyStep(double[][][] t0,
                         double[][][] t1,
                         int sizeX, int sizeY, int sizeZ,
                         double dx, double dy, double dz,
                         double dt, int a,
                         double kmat, int Ta, double hs) {

        double mxyz = dt/(kmat/8.418e-5*dx*dy*dz);
        for (int k = 0; k < sizeY; ++k) {
            for (int n = 0; n < sizeZ; ++n) {


                if (k == 0)
                    if (n == 0)
                        if (a == 0)
                            t1[k][n][a] = 8*mxyz*(hs*dx*dz/2*(Ta-t0[k][n][a])
                                    +kmat/4*(dz*(t0[k][n+1][a]+t0[k+1][n][a]-2*t0[k][n][a])
                                    +dx*dy/dz*(t0[k][n][a+1]-t0[k][n][a])))+ t0[k][n][a];
                        else if (a == sizeX-1)
                            t1[k][n][a] = 8*mxyz*(hs*dx*dz/2*(Ta-t0[k][n][a])+hs*dx*dy/4*(Ta-t0[k][n][a])
                                    +kmat/4*(dz*(t0[k][n+1][a]+t0[k+1][n][a]-2*t0[k][n][a])
                                    +dx*dy/dz*(t0[k][n][a-1]-t0[k][n][a])))+ t0[k][n][a];
                        else
                            t1[k][n][a] = 4*mxyz*(hs*dx*dz*(Ta-t0[k][n][a]) + kmat/4*(2*dz*(t0[k][n+1][a]
                                    +t0[k+1][n][a]-2*t0[k][n][a])+dx*dy/dz*(t0[k][n][a-1]
                                    +t0[k][n][a+1]-2*t0[k][n][a])))+ t0[k][n][a];
                    else if (n==sizeZ-1)
                        if (a == 0)
                            t1[k][n][a] = 8*mxyz*(hs*dx*dz/2*(Ta-t0[k][n][a])+kmat/4*(dz*(t0[k][n-1][a]
                                    +t0[k+1][n][a]-2*t0[k][n][a])+dx*dy/dz*(t0[k][n][a+1]
                                    -t0[k][n][a])))+ t0[k][n][a];
                        else if (a == sizeX-1)
                            t1[k][n][a] = 8*mxyz*(hs*dx*dz/2*(Ta-t0[k][n][a])+hs*dx*dy/4*(Ta-t0[k][n][a])
                                    +kmat/4*(dz*(t0[k][n-1][a]+t0[k+1][n][a]-2*t0[k][n][a])
                                    +dx*dy/dz*(t0[k][n][a-1]-t0[k][n][a])))+ t0[k][n][a];
                        else
                            t1[k][n][a] = 4*mxyz*(hs*dx*dz*(Ta-t0[k][n][a]) + kmat/4*(2*dz*(t0[k][n-1][a]
                                    +t0[k+1][n][a]-2*t0[k][n][a])+dx*dy/dz*(t0[k][n][a-1]
                                    +t0[k][n][a+1]-2*t0[k][n][a])))+ t0[k][n][a];
                    else
                    if (a == 0)
                        t1[k][n][a] = 4*mxyz*(hs*dx*dz/2*(Ta-t0[k][n][a]) + kmat/4*(dz*(t0[k][n+1][a]
                                +t0[k][n-1][a]-2*t0[k][n][a])+2*dx*dy/dz*(t0[k][n][a+1]-t0[k][n][a])
                                +2*dz*(t0[k+1][n][a]-t0[k][n][a])))+ t0[k][n][a];
                    else if (a == sizeX-1)
                        t1[k][n][a] = 4*mxyz*(hs*dx*dz/2*(Ta-t0[k][n][a])+hs*dx*dy/2*(Ta-t0[k][n][a])
                                + kmat/4*(2*dz*(t0[k+1][n][a]-t0[k][n][a])+2*dx*dy/dz*(t0[k][n][a-1]-t0[k][n][a])
                                +dz*dx*(t0[k][n-1][a]+t0[k][n+1][a]-2*t0[k][n][a])))+ t0[k][n][a];
                    else
                        t1[k][n][a] = 2*mxyz*(hs*dx*dz*(Ta-t0[k][n][a]) + kmat*(dx*dy/(2*dz)*(t0[k][n][a+1]
                                +t0[k][n][a-1]-2*t0[k][n][a]) + dz/2*(t0[k][n+1][a]+t0[k][n-1][a]-2*t0[k][n][a])
                                +dz*(t0[k+1][n][a]-t0[k][n][a]))) + t0[k][n][a];

                else if (k == sizeY-1)

                    if (n == 0)

                        if (a == 0)
                            t1[k][n][a] = 8*mxyz*(hs*dx*dz/2*(Ta-t0[k][n][a])+kmat/4*(dz*(t0[k][n+1][a]
                                    +t0[k-1][n][a]-2*t0[k][n][a])+dx*dy/dz*(t0[k][n][a+1]
                                    -t0[k][n][a])))+t0[k][n][a];
                        else if (a == sizeX -1)
                            t1[k][n][a] = 8*mxyz*(hs*dx*dz/2*(Ta-t0[k][n][a])+hs*dx*dy/4*(Ta-t0[k][n][a])
                                    +kmat/4*(dz*(t0[k][n+1][a]+t0[k-1][n][a]-2*t0[k][n][a])
                                    +dx*dy/dz*(t0[k][n][a-1]-t0[k][n][a])))+t0[k][n][a];
                        else
                            t1[k][n][a] = 4*mxyz*(hs*dx*dz*(Ta-t0[k][n][a]) + kmat/4*(2*dz*(t0[k][n+1][a]
                                    +t0[k-1][n][a]-2*t0[k][n][a])+dx*dy/dz*(t0[k][n][a-1]
                                    +t0[k][n][a+1]-2*t0[k][n][a])))+ t0[k][n][a];

                    else if (n == sizeZ-1)
                        if (a == 0)
                            t1[k][n][a] = 8*mxyz*(hs*dx*dz/2*(Ta-t0[k][n][a])+kmat/4*(dz*(t0[k][n-1][a]
                                    +t0[k-1][n][a]-2*t0[k][n][a])+dx*dy/dz*(t0[k][n][a+1]
                                    -t0[k][n][a])))+t0[k][n][a];
                        else if (a == sizeX-1)
                            t1[k][n][a] = 8*mxyz*(hs*dx*dz/2*(Ta-t0[k][n][a])+hs*dx*dy/4*(Ta-t0[k][n][a])
                                    +kmat/4*(dz*(t0[k][n-1][a]+t0[k-1][n][a]-2*t0[k][n][a])
                                    +dx*dy/dz*(t0[k][n][a-1]-t0[k][n][a])))+t0[k][n][a];
                        else
                            t1[k][n][a] = 4*mxyz*(hs*dx*dz*(Ta-t0[k][n][a]) + kmat/4*(2*dz*(t0[k][n-1][a]
                                    +t0[k-1][n][a]-2*t0[k][n][a])+dx*dy/dz*(t0[k][n][a-1]
                                    +t0[k][n][a+1]-2*t0[k][n][a])))+ t0[k][n][a];

                    else
                    if (a == 0)
                        t1[k][n][a] = 4*mxyz*(hs*dx*dz/2*(Ta-t0[k][n][a]) + kmat/4*(dz*(t0[k][n+1][a]
                                +t0[k][n-1][a]-2*t0[k][n][a])+2*dx*dy/dz*(t0[k][n][a+1]-t0[k][n][a])
                                +2*dz*(t0[k-1][n][a]-t0[k][n][a])))+ t0[k][n][a];
                    else if (a == sizeX-1)
                        t1[k][n][a] = 4*mxyz*(hs*dx*dz/2*(Ta-t0[k][n][a]) + hs*dx*dy/2*(Ta-t0[k][n][a])
                                + kmat/4*(2*dz*(t0[k-1][n][a]-t0[k][n][a])+2*dx*dy/dz*(t0[k][n][a-1]-t0[k][n][a])
                                +dz*dx*(t0[k][n-1][a]+t0[k][n+1][a]-2*t0[k][n][a])))+ t0[k][n][a];
                    else
                        t1[k][n][a] = 2*mxyz*(hs*dx*dz*(Ta-t0[k][n][a]) + kmat*(dx*dy/(2*dz)*(t0[k][n][a+1]
                                +t0[k][n][a-1]-2*t0[k][n][a]) + dz/2*(t0[k][n+1][a]+t0[k][n-1][a]
                                -2*t0[k][n][a])+dz*(t0[k-1][n][a]-t0[k][n][a]))) + t0[k][n][a];


                else

                if (n == 0)
                    if (a == 0)
                        t1[k][n][a] = 4*mxyz*(hs*dx*dz/2*(Ta-t0[k][n][a]) + kmat/4*(dz*(t0[k+1][n][a]+t0[k-1][n][a]
                                -2*t0[k][n][a])+2*dx*dy/dz*(t0[k][n][a+1]-t0[k][n][a])
                                +2*dz*(t0[k][n+1][a]-t0[k][n][a])))+ t0[k][n][a];
                    else if (a == sizeX-1)
                        t1[k][n][a] = 4*mxyz*(hs*dx*dz/2*(Ta-t0[k][n][a]) + hs*dx*dy/2*(Ta-t0[k][n][a])
                                + kmat/4*(dz*(t0[k+1][n][a]+t0[k-1][n][a]-2*t0[k][n][a])+2*dx*dy/dz*(t0[k][n][a-1]
                                -t0[k][n][a])+2*dz*(t0[k][n+1][a]-t0[k][n][a])))+ t0[k][n][a];
                    else
                        t1[k][n][a] = 2*mxyz*(hs*dx*dz*(Ta-t0[k][n][a]) + kmat*(dx*dy/(2*dz)*(t0[k][n][a+1]
                                +t0[k][n][a-1]-2*t0[k][n][a])+dz/2*(t0[k+1][n][a]+t0[k-1][n][a]-2*t0[k][n][a])
                                +dz*(t0[k][n+1][a]-t0[k][n][a])))+t0[k][n][a];

                else if (n == sizeZ-1)

                    if (a == 0)
                        t1[k][n][a] = 4*mxyz*(hs*dx*dz/2*(Ta-t0[k][n][a]) + kmat/4*(dz*(t0[k+1][n][a]+t0[k-1][n][a]
                                -2*t0[k][n][a])+2*dx*dy/dz*(t0[k][n][a+1]-t0[k][n][a])+2*dz*(t0[k][n-1][a]
                                -t0[k][n][a])))+ t0[k][n][a];
                    else if (a == sizeX-1)
                        t1[k][n][a] = 4*mxyz*(hs*dx*dz/2*(Ta-t0[k][n][a]) + hs*dx*dy/2*(Ta-t0[k][n][a])
                                + kmat/4*(dz*(t0[k+1][n][a]+t0[k-1][n][a]-2*t0[k][n][a])+2*dx*dy/dz*(t0[k][n][a-1]
                                -t0[k][n][a])+2*dz*(t0[k][n-1][a]-t0[k][n][a])))+ t0[k][n][a];
                    else
                        t1[k][n][a] = 2*mxyz*(hs*dx*dz*(Ta-t0[k][n][a]) + kmat*(dx*dy/(2*dz)*(t0[k][n][a+1]+t0[k][n][a-1]
                                -2*t0[k][n][a])+dz/2*(t0[k+1][n][a]+t0[k-1][n][a]-2*t0[k][n][a])
                                +dz*(t0[k][n-1][a]-t0[k][n][a])))+t0[k][n][a];
                else

                if (a == 0)

                    t1[k][n][a] = 2*mxyz*kmat*(dz/2*(t0[k][n+1][a] + t0[k][n-1][a] + t0[k+1][n][a]
                            + t0[k-1][n][a]-4*t0[k][n][a]) + dx*dy/dz *(t0[k][n][a+1]
                            -t0[k][n][a]))+ t0[k][n][a];
                else if (a == sizeX-1)

                    t1[k][n][a] = 2*mxyz*(hs*dx*dy*(Ta-t0[k][n][a])+kmat*(dz*(t0[k][n+1][a]+t0[k][n-1][a]
                            +t0[k+1][n][a]+t0[k-1][n][a]-4*t0[k][n][a])
                            +dx*dy/dz*(t0[k][n][a-1]-t0[k][n][a])))+t0[k][n][a];
                else

                    t1[k][n][a] = mxyz*kmat*(dz*(t0[k][n+1][a] +
                            t0[k][n-1][a] + t0[k+1][n][a] + t0[k-1][n][a]-4*t0[k][n][a])
                            + dx*dy/dz *(t0[k][n][a+1]-2*t0[k][n][a]+t0[k][n][a-1]))+t0[k][n][a];
            }
        }

    }

    void calculate(double[][][] tInit, int nX, int nY, int nZ,int iterations, int temp) throws InterruptedException {
        double dx = 1.0 / nX;
        double dy = 1.0 / nY;
        double dz = 1.0/nZ;
        double dt = 0.005;

        double h=30;
        int K = 247;
        int i1 = 0;
        int i2 = i1 + nX;

        double[][][] T0 = new double[nX][nY][nZ];
        double[][][] T1 = new double[nX][nY][nZ];
        for (int i = i1; i < i2; i++) {
            for (int j = 0; j < nY; j++) {
                for (int k=0; k<nZ; k++) {
                    T0[i - i1][j][k] = tInit[i][j][k];
                    T1[i - i1][j][k] = tInit[i][j][k];
                }
            }
        }

        for (int s = 0; s <= iterations; ++s) {

            step(T0, T1, nX, nY, nZ, dx, dy, dz,dt, temp, h, K);
            double[][][] tmp;
            tmp = T0;
            T0 = T1;
            T1 = tmp;

        }
        for (int i = 0; i < nX; i++) {
            for (int j = 0; j < nX; j++) {
                System.arraycopy(T0[i][j], 0, tInit[i][j], 0, nY);
            }
        }

        executorService.shutdown();
    }
}
