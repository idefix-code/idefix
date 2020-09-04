#include "../idefix.hpp"
#include "hydro.hpp"

// Compute the electrical current on faces
void Hydro::CalcCurrent(DataBlock &data) {
    idfx::pushRegion("Hydro::CalcCurrent");
    IdefixArray4D<real> Vc = data.Vc;
    IdefixArray4D<real> Vs = data.Vs;
    IdefixArray4D<real> J = data.J;

    IdefixArray1D<real> dx1 = data.dx[IDIR];
    IdefixArray1D<real> dx2 = data.dx[JDIR];
    IdefixArray1D<real> dx3 = data.dx[KDIR];

    IdefixArray1D<real> r = data.x[IDIR];
    IdefixArray1D<real> rm = data.xl[IDIR];
    IdefixArray1D<real> th = data.x[JDIR];

    idefix_for("CalcCurrent",
                KOFFSET,data.np_tot[KDIR],
                JOFFSET,data.np_tot[JDIR],
                IOFFSET,data.np_tot[IDIR],
                KOKKOS_LAMBDA (int k, int j, int i) {
                    real Bx1_000 = ZERO_F, Bx1_0m0 = ZERO_F, Bx1_00m = ZERO_F;
                    real Bx2_000 = ZERO_F, Bx2_m00 = ZERO_F, Bx2_00m = ZERO_F;
                    real Bx3_000 = ZERO_F, Bx3_m00 = ZERO_F, Bx3_0m0 = ZERO_F;

                    real d12 = ZERO_F, d13 = ZERO_F, d21 = ZERO_F, d23 = ZERO_F, d31 = ZERO_F, d32 = ZERO_F;
                    real Jx = ZERO_F, Jy = ZERO_F, Jz = ZERO_F;

                    Bx1_000 = Vs(BX1s,k,j,i);
                    #if DIMENSIONS >= 2
                        Bx1_0m0 = Vs(BX1s,k,j-1,i);
                        Bx2_000 = Vs(BX2s,k,j,i);
                        Bx2_m00 = Vs(BX2s,k,j,i-1);
                    #elif COMPONENTS >= 2   // DIMENSIONS == 1 here!
                        Bx2_000 = Vc(BX2,k,j,i);
                        Bx2_m00 = Vc(BX2,k,j,i-1);
                        #if COMPONENTS == 3
                            Bx3_000 = Vc(BX3,k,j,i);
                            Bx3_m00 = Vc(BX3,k,j,i-1);
                        #endif
                    #endif
                    #if DIMENSIONS == 3
                        Bx1_00m = Vs(BX1s,k-1,j,i);
                        Bx2_00m = Vs(BX2s,k-1,j,i);
                        Bx3_000 = Vs(BX3s,k,j,i);
                        Bx3_m00 = Vs(BX3s,k,j,i-1);
                        Bx3_0m0 = Vs(BX3s,k,j-1,i);
                    #elif COMPONENTS == 3  && DIMENSIONS == 2 // DIMENSIONS 2 here
                        Bx3_0m0 = Vc(BX3,k,j-1,i);
                    #endif


                    // define geometrical factors
                    D_EXPAND( d12 = d13 = ONE_F / dx1(i);   ,
                              d21 = d23 = ONE_F / dx2(j);   ,
                              d31 = d32 = ONE_F / dx3(k); )

                    real a13Bx3_000 = Bx3_000;
                    real a13Bx3_m00 = Bx3_m00;
                    real a23Bx3_000 = Bx3_000;
                    real a23Bx3_0m0 = Bx3_0m0;
                    real a12Bx2_000 = Bx2_000;
                    real a12Bx2_m00 = Bx2_m00;


                    #if GEOMETRY == CYLINDRICAL
                        d32 = d31 = ZERO_F;
                        d13 = TWO_F / (FABS(r(i))*r(i) - FABS(r(i-1))*r(i-1));
                        #if COMPONENTS == 3
                            a13Bx3_000 = Bx3_000 * FABS(r(i));
                            a13Bx3_m00 = Bx3_m00 * FABS(r(i-1));
                        #endif
                    #elif GEOMETRY == POLAR
                        d23 /= r(i);
                        d21 /= rm(i);
                        d12 = TWO_F / (FABS(r(i))*r(i) - FABS(r(i-1))*r(i-1));
                        #if COMPONENTS >= 2
                            a12Bx2_000 = Bx2_000 * r(i);
                            a12Bx2_m00 = Bx2_m00 * r(i-1);
                        #endif
                    #elif GEOMETRY == SPHERICAL
                        real s = FABS(SIN(th(j)));
                        real sm = HALF_F*( FABS(SIN(th(j))) + FABS(SIN(th(j-1))));
                        D_EXPAND(d12 /= rm(i);   d13 /= rm(i);     ,
                                d21 /= rm(i);   d23 /= r(i)*sm;   ,
                                d32 /= r(i)*sm; d31 /= rm(i)*s;)

                        #if COMPONENTS >= 2
                            a12Bx2_000 = Bx2_000 * r(i);
                            a12Bx2_m00 = Bx2_m00 * r(i-1);
                        #endif
                        #if COMPONENTS == 3
                            a13Bx3_000 = Bx3_000 * r(i));
                            a13Bx3_m00 = Bx3_m00 * r(i-1);
                            #if DIMENSIONS >= 2
                                a23Bx3_000 = Bx3_000 * FABS(SIN(th(j)));
                                a23Bx3_0m0 = Bx3_0m0 * FABS(SIN(th(j-1)));
                            #endif
                        #endif

                    #endif

                    // Compute actual current
                    
                    #if COMPONENTS == 3
                        #if DIMENSIONS == 3
                            Jx +=  - (Bx2_000 - Bx2_00m)*d32;
                            Jy +=  + (Bx1_000 - Bx1_00m)*d31;
                        #endif
                        #if DIMENSIONS >= 2
                            Jx += + (a23Bx3_000 - a23Bx3_0m0)*d23;
                            Jy += - (a13Bx3_000 - a13Bx3_m00)*d13;
                        #endif
                    #endif

                    #if COMPONENTS >= 2
                    Jz += (a12Bx2_000 - a12Bx2_m00)*d12;
                    #endif
                    #if DIMENSIONS >= 2
                    Jz += - (Bx1_000 - Bx1_0m0)*d21;
                    #endif

                    // Store the beast
                    J(IDIR, k, j, i) = Jx;
                    J(JDIR, k, j, i) = Jy;
                    J(KDIR, k, j, i) = Jz;




                });



    idfx::popRegion();
}