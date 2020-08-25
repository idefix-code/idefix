#include "idefix.hpp"
#include "setup.hpp"

/*********************************************/
/**
Customized random number generator
Allow one to have consistant random numbers
generators on different architectures.
**/
/*********************************************/
real randm(void) {
    const int a    =    16807;
    const int m =    2147483647;
    static int in0 = 13763 + 2417*idfx::prank;
    int q;
    
    /* find random number  */
    q= (int) fmod((double) a * in0, m);
    in0=q;
    
    return((real) ((double) q/(double)m));
}


// User-defined boundaries
void UserdefBoundary(DataBlock& data, int dir, BoundarySide side, real t) {
    if( (dir==IDIR) && (side == left)) {
        IdefixArray4D<real> Vc = data.Vc;
        IdefixArray4D<real> Vs = data.Vs;
        IdefixArray1D<real> x1 = data.x[IDIR];

        int ighost = data.nghost[IDIR];
        idefix_for("UserDefBoundary",0,data.np_tot[KDIR],0,data.np_tot[JDIR],0,ighost,
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        Vc(RHO,k,j,i) = Vc(RHO,k,j,ighost)/ pow(x1(i)/x1(ighost),1.5);
                        Vc(PRS,k,j,i) = Vc(PRS,k,j,ighost)/ pow(x1(i)/x1(ighost),2.5);;
                        if(Vc(VX1,k,j,ighost)>=ZERO_F) Vc(VX1,k,j,i)=ZERO_F;
			else Vc(VX1,k,j,i) = Vc(VX1,k,j,ighost) /( sqrt(x1(i)/x1(ighost)));
                        Vc(VX2,k,j,i) = Vc(VX2,k,j,ighost) /( sqrt(x1(i)/x1(ighost)));
                        Vc(VX3,k,j,i) = Vc(VX3,k,j,ighost) /( sqrt(x1(i)/x1(ighost)));
                        Vs(BX2s,k,j,i) = Vs(BX2s,k,j,ighost);
                        Vs(BX3s,k,j,i) = ZERO_F; 

                    });
    }
    if( dir==JDIR) {
        IdefixArray4D<real> Vc = data.Vc;
        IdefixArray4D<real> Vs = data.Vs;
        int jghost;
        int jbeg,jend;
        if(side == left) {
            jghost = data.beg[JDIR];
            jbeg = 0;
            jend = data.beg[JDIR];
        }
        else {
            jghost = data.end[JDIR]-1;
            jbeg=data.end[JDIR];
            jend=data.np_tot[JDIR];
        }
        idefix_for("UserDefBoundary",0,data.np_tot[KDIR],jbeg,jend,0,data.np_tot[IDIR],
                    KOKKOS_LAMBDA (int k, int j, int i) {
                        Vc(RHO,k,j,i) = Vc(RHO,k,jghost,i);
                        Vc(PRS,k,j,i) = Vc(PRS,k,jghost,i);
                        Vc(VX1,k,j,i) = Vc(VX1,k,jghost,i);
                        Vc(VX2,k,j,i) = ZERO_F;
                        Vc(VX3,k,j,i) = Vc(VX3,k,jghost,i);
                        Vs(BX1s,k,j,i) = ZERO_F;
                        Vs(BX3s,k,j,i) = ZERO_F;

                    });


    }
    
}

void Potential(DataBlock& data, const real t, IdefixArray1D<real>& x1, IdefixArray1D<real>& x2, IdefixArray1D<real>& x3, IdefixArray3D<real>& phi) {
    
    idefix_for("Potential",0,data.np_tot[KDIR], 0, data.np_tot[JDIR], 0, data.np_tot[IDIR],
        KOKKOS_LAMBDA (int k, int j, int i) {
        phi(k,j,i) = -1.0/x1(i);
    });
    
}

// Default constructor
Setup::Setup() {}

// Initialisation routine. Can be used to allocate
// Arrays or variables which are used later on
Setup::Setup(Input &input, Grid &grid, DataBlock &data, Hydro &hydro) {
    // Set the function for userdefboundary
    hydro.EnrollUserDefBoundary(&UserdefBoundary);
    hydro.EnrollGravPotential(&Potential);
    hydro.SetGamma(1.05);
}

// This routine initialize the flow
// Note that data is on the device.
// One can therefore define locally 
// a datahost and sync it, if needed
void Setup::InitFlow(DataBlock &data) {
    // Create a host copy
    DataBlockHost d(data);
    
    // Make vector potential
    IdefixHostArray4D<real> A = IdefixHostArray4D<real>("Setup_VectorPotential", 3, data.np_tot[KDIR], data.np_tot[JDIR], data.np_tot[IDIR]);
    
    
    real x,y,z;

    real vphi,f,r,th;
    
    real epsilon=0.1;
    real beta=1000;
    
    for(int k = 0; k < d.np_tot[KDIR] ; k++) {
        for(int j = 0; j < d.np_tot[JDIR] ; j++) {
            for(int i = 0; i < d.np_tot[IDIR] ; i++) {
                r=d.x[IDIR](i);
                th=d.x[JDIR](j);
                real R=r*sin(th);
                real z=r*cos(th);
                real Vk=1.0/pow(R,0.5);
                
                real cs2=(epsilon*Vk)*(epsilon*Vk);

                d.Vc(RHO,k,j,i) = 1.0/(R*sqrt(R))*exp(1.0/(cs2)*(1/r-1/R));
                d.Vc(PRS,k,j,i) = cs2*d.Vc(RHO,k,j,i);
                d.Vc(VX1,k,j,i) = 0.0;
                d.Vc(VX2,k,j,i) = 1e-1*(0.5-randm());
                d.Vc(VX3,k,j,i) = Vk*sqrt(R/r-2.5*cs2);
                
                d.Vs(BX1s,k,j,i) = 0.0;
                d.Vs(BX2s,k,j,i) = 0.0;
                d.Vs(BX3s,k,j,i) = 0.0;
                
                real B0 = sqrt(2*cs2/(R*sqrt(R))/sqrt(beta));
                
		d.Vs(BX3s,k,j,i) = B0*cos(R/epsilon)*fmax(1-(z*z)/(4*R*R*epsilon*epsilon),ZERO_F);
                d.Vs(BX3s,k,j,i) *= fmax(tanh(10*(R-1.5)),ZERO_F);
                A(IDIR,k,j,i) = 0.0;
                A(JDIR,k,j,i) = 0.0;
                A(KDIR,k,j,i) = B0*epsilon*cos(R/epsilon)*fmax(1-(z*z)/(4*R*R*epsilon*epsilon),ZERO_F);
            }
        }
    }
    
    // Make the field from the vector potential
    //d.MakeVsFromAmag(A);
    
    // Send it all, if needed
    d.SyncToDevice();
}

// Analyse data to produce an output
void Setup::MakeAnalysis(DataBlock & data, real t) {

}




// Do a specifically designed user step in the middle of the integration
void ComputeUserStep(DataBlock &data, real t, real dt) {

}
