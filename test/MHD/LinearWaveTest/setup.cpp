#include "idefix.hpp"
#include "setup.hpp"

int mode;
real epsilon;



// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Output &output) {

  mode = input.Get<int>("Setup","mode",0);
  epsilon = input.Get<real>("Setup","epsilon",0);

}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);
    IdefixHostArray4D<real> A = IdefixHostArray4D<real>("VectorPotential",3,d.np_tot[KDIR],d.np_tot[JDIR],d.np_tot[IDIR]);

    real x,y,z;

    real gamma = data.hydro->eos->GetGamma();
    // The eigenmode
    real R[8];
    real R0;

    // The mode decomposition is described in
    // T.A. Gardiner, J.M. Stone / Journal of Computational Physics 227 (2008) 4123–4141

    // The background state
    real B[8];
    B[RHO] = 1;
    B[VX1] = 0.0;
    B[VX2] = 0.0;
    B[VX3] = 0.0;
    B[BX1] = 1.0;
    B[BX2] = 3.0/2.0;
    B[BX3] = 0.0;
    B[PRS] = 1/gamma;

    int absMode = std::abs(mode);
    int sign = mode/absMode;

    switch(absMode) {
      case 1:
        // Fast mode
        R0 = 1.0/(2*sqrt(5));
        R[RHO] = 2;
        R[VX1] = sign*4.0;
        R[VX2] = -sign*2.0;
        R[VX3] = 0.0;
        R[BX1] = 0.0;
        R[BX2] = 4.0;
        R[BX3] = 0.0;
        R[PRS] = (gamma-1)*(9.0 - B[BX1]*R[BX1] - B[BX2]*R[BX2] - B[BX3]*R[BX3]);
        break;
      case 2:
        // Alfven mode
        R0 = 1.0;
        R[RHO] = 0;
        R[VX1] = 0;
        R[VX2] = 0;
        R[VX3] = -sign;
        R[BX1] = 0.0;
        R[BX2] = 0.0;
        R[BX3] = 1.0;
        R[PRS] = 0.0;
        break;
      case 3:
        // Slow mode
        R0 = 1.0/(2*sqrt(5));
        R[RHO] = 4;
        R[VX1] = 2*sign;
        R[VX2] = 4*sign;
        R[VX3] = 0.0;
        R[BX1] = 0.0;
        R[BX2] = -2.0;
        R[BX3] = 0.0;
        R[PRS] = (gamma-1)*(3.0 - B[BX1]*R[BX1] - B[BX2]*R[BX2] - B[BX3]*R[BX3]);
        break;
      case 4:
        // entropy mode
        R0 = 0.5;
        B[VX1] = 1; // Add a mean flow to the background
        R[RHO] = 2;
        R[VX1] = 0;
        R[VX2] = 0;
        R[VX3] = 0;
        R[BX1] = 0;
        R[BX2] = 0;
        R[BX3] = 0.0;
        R[PRS] = (gamma-1)*(1.0 - B[BX1]*R[BX1] - B[BX2]*R[BX2] - B[BX3]*R[BX3]-B[VX1]*R[RHO]/2.0);
        break;
      default:
        IDEFIX_ERROR("unknown mode");
    }

    real sinA =  2.0/3.0;
    real sinB =  2.0/sqrt(5.0);
    real cosA = sqrt(1-sinA*sinA);
    real cosB = sqrt(1-sinB*sinB);

    // Projection of the unit vectors 1 and 2 onto our base
    real ex1_x = cosA*cosB;
    real ex1_y = cosA*sinB;
    real ex1_z = sinA;

    real ex2_x = -sinB;
    real ex2_y = cosB;
    real ex2_z = 0;

    real ex3_x = -sinA*cosB;
    real ex3_y = -sinA*sinB;
    real ex3_z = cosA;
/*
    real ex1_x = cosA*cosB;
    real ex1_y = -sinB;
    real ex1_z = -sinA*cosB;

    real ex2_x = cosA*sinB;
    real ex2_y = cosB;
    real ex2_z = sinA*sinB;

    real ex3_x = sinA;
    real ex3_y = 0;
    real ex3_z = cosA;
  */

    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                x=d.x[IDIR](i);
                y=d.x[JDIR](j);
                z=d.x[KDIR](k);

                real x1 = cosA*cosB*x + cosA*sinB*y + sinA*z;

                real phase = epsilon*cos(2.0*M_PI*x1);

                d.Vc(RHO,k,j,i) = B[RHO] + R0*R[RHO]*phase;
                d.Vc(PRS,k,j,i) = B[PRS] + R0*R[PRS]*phase;


                d.Vc(VX1,k,j,i) = B[VX1]*ex1_x + B[VX2]*ex2_x + B[VX3]*ex3_x + R0*(R[VX1]*ex1_x + R[VX2]*ex2_x + R[VX3]*ex3_x)*phase;
                d.Vc(VX2,k,j,i) = B[VX1]*ex1_y + B[VX2]*ex2_y + B[VX3]*ex3_y + R0*(R[VX1]*ex1_y + R[VX2]*ex2_y + R[VX3]*ex3_y)*phase;
                d.Vc(VX3,k,j,i) = B[VX1]*ex1_z + B[VX2]*ex2_z + B[VX3]*ex3_z + R0*(R[VX1]*ex1_z + R[VX2]*ex2_z + R[VX3]*ex3_z)*phase;


                // Now the potential vector
                real A1, A2, A3;

                // Ax1 component
                x=d.x[IDIR](i);
                y=d.xl[JDIR](j);
                z=d.xl[KDIR](k);

                x1 = cosA*cosB*x + cosA*sinB*y + sinA*z;
                real x2 = -sinB*x + cosB*y;

                phase = epsilon*sin(2.0*M_PI*x1)/(2.0*M_PI);

                A1 = 0.0;
                A2 = B[BX3]*x1 + R0*R[BX3]*phase;
                A3 = -(B[BX2]*x1 + R0*R[BX2]*phase) + B[BX1]*x2;

                A(IDIR,k,j,i) = A1*ex1_x + A2*ex2_x + A3*ex3_x;

                // Ax2 component
                x=d.xl[IDIR](i);
                y=d.x[JDIR](j);
                z=d.xl[KDIR](k);

                x1 = cosA*cosB*x + cosA*sinB*y + sinA*z;
                x2 = -sinB*x + cosB*y;

                phase = epsilon*sin(2.0*M_PI*x1)/(2.0*M_PI);

                A1 = 0.0;
                A2 = B[BX3]*x1 + R0*R[BX3]*phase;
                A3 = -(B[BX2]*x1 + R0*R[BX2]*phase) + B[BX1]*x2;

                A(JDIR,k,j,i) = A1*ex1_y + A2*ex2_y + A3*ex3_y;

                // Ax3 component
                x=d.xl[IDIR](i);
                y=d.xl[JDIR](j);
                z=d.x[KDIR](k);

                x1 = cosA*cosB*x + cosA*sinB*y + sinA*z;
                x2 = -sinB*x + cosB*y;

                phase = epsilon*sin(2.0*M_PI*x1)/(2.0*M_PI);

                A1 = 0.0;
                A2 = B[BX3]*x1 + R0*R[BX3]*phase;
                A3 = -(B[BX2]*x1 + R0*R[BX2]*phase) + B[BX1]*x2;

                A(KDIR,k,j,i) = A1*ex1_z + A2*ex2_z + A3*ex3_z;
                /*
                d.Vs(BX1s,k,j,i) = B[BX1]*ex1_x;
                d.Vs(BX2s,k,j,i) = B[BX1]*ex1_y;
                d.Vs(BX3s,k,j,i) = B[BX1]*ex1_z;
                */
              }}}
    // Create the staggered field from the magnetic potential
    d.MakeVsFromAmag(A);


    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void MakeAnalysis(DataBlock & data) {

}
