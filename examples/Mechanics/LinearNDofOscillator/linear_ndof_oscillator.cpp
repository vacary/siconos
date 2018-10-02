
#include "SiconosKernel.hpp"
//#define TS_PROJ
//#define TS_COMBINED
using namespace std;

#include <math.h>
class ForcedNdofOscillator : public LagrangianDS
{

  
  
public:
  ForcedNdofOscillator(SP::SiconosVector q0, SP::SiconosVector v0, SP::SiconosMatrix M):LagrangianDS(q0,v0,M)
  {
    int nDof = _ndof;
    SP::SiconosMatrix SparseStiffness(new SimpleMatrix(nDof,nDof,Siconos::SPARSE,3*nDof));
 
    SparseStiffness->setValue(0, 0, 1.0);
    SparseStiffness->setValue(0, 1, -1.0);

    for(unsigned int i = 1; i < nDof-1; i++)
    {
      SparseStiffness->setValue(i,i,2.0);
      SparseStiffness->setValue(i,i-1,-1.0);
      SparseStiffness->setValue(i,i+1,-1.0);
    }
    SparseStiffness->setValue(nDof-1, nDof-2,-1.0);
    SparseStiffness->setValue(nDof-1,nDof-1,1.0);
    setKPtr(SparseStiffness);
    _K->display();
  };

  void computeFExt(double time)
  {
    double A=1.0;
    double omega=10.0;
    _fExt->setValue(_ndof-1, A*cos(omega*time));
  };
    
  void computeFExt(double time, SP::SiconosVector fExt)
  {
    double A=1.0;
    double omega=10.0;
    fExt->setValue(_ndof-1, A*cos(omega*time));
  };
  
  
};

TYPEDEF_SPTR(ForcedNdofOscillator)
int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
#include "UserDefinedParameter.hpp"

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." <<endl<<endl;

    SP::SiconosMatrix SparseMass(new SimpleMatrix(nDof,nDof,Siconos::SPARSE,nDof));

    for(unsigned int i = 0; i < nDof; i++)
    {
      SparseMass->setValue(i,i,1.0);
    }
    
    // -- Initial positions and velocities --
    SP::SiconosVector q0(new SiconosVector(nDof,position_init));
    SP::SiconosVector v0(new SiconosVector(nDof,velocity_init));

    // -- The dynamical system --
    SP::ForcedNdofOscillator bar(new ForcedNdofOscillator(q0,v0,SparseMass));



    // -- Set external forces (weight) --
    //SP::SiconosVector weight(new SiconosVector(nDof,-g*rho*S/l));
    SP::SiconosVector weight(new SiconosVector(nDof,0.0));
    bar->setFExtPtr(weight);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = .5;

    // Interaction bar-floor
    //
    SP::SimpleMatrix H(new SimpleMatrix(1,nDof));
    (*H)(0,0) = 1.0;

    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    SP::Relation relation(new LagrangianLinearTIR(H));

    SP::Interaction inter(new Interaction(nslaw, relation));

    // -------------
    // --- Model ---
    // -------------
    SP::NonSmoothDynamicalSystem impactingBar(new NonSmoothDynamicalSystem(t0, T));

    // add the dynamical system in the non smooth dynamical system
    impactingBar->insertDynamicalSystem(bar);

    // link the interaction and the dynamical system
    impactingBar->link(inter,bar);


    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
#ifdef TS_PROJ
    SP::MoreauJeanDirectProjectionOSI OSI(new MoreauJeanDirectProjectionOSI(theta));
    OSI->setDeactivateYPosThreshold(1e-05);
    OSI->setDeactivateYVelThreshold(0.0);
    OSI->setActivateYPosThreshold(1e-09);
    OSI->setActivateYVelThreshold(100.0);
#else
#ifdef TS_COMBINED
    SP::OneStepIntegrator OSI(new MoreauJeanCombinedProjectionOSI(theta));
#else
    SP::MoreauJeanOSI OSI(new MoreauJeanOSI(theta,0.5));
#endif
#endif
    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0,h));

    // -- (3) one step non smooth problem
    SP::OneStepNSProblem osnspb(new LCP());

    // -- (4) Simulation setup with (1) (2) (3)
#ifdef TS_PROJ
    SP::MLCPProjectOnConstraints position(new MLCPProjectOnConstraints());
    SP::TimeSteppingDirectProjection s(new TimeSteppingDirectProjection(impactingBar, t,OSI, osnspb, position,0));
    s->setProjectionMaxIteration(10);
    s->setConstraintTolUnilateral(1e-10);
    s->setConstraintTol(1e-10);


#else
#ifdef TS_COMBINED
    SP::OneStepNSProblem position(new MLCPProjectOnConstraints(SICONOS_MLCP_ENUM));
    SP::TimeSteppingCombinedProjection s(new TimeSteppingCombinedProjection(impactingBar, t,OSI, osnspb, position,2));
    s->setProjectionMaxIteration(500);
    s->setConstraintTolUnilateral(1e-10);
    s->setConstraintTol(1e-10);
#else
    SP::TimeStepping s(new TimeStepping(impactingBar, t,OSI,osnspb));
#endif
#endif

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---


    int N = floor((T-t0)/h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 11;
    SimpleMatrix dataPlot(N,outputSize);

    SP::SiconosVector q = bar->q();
    SP::SiconosVector v = bar->velocity();
    SP::SiconosVector p = bar->p(1);
    SP::SiconosVector lambda = inter->lambda(1);

    dataPlot(0, 0) = impactingBar->t0();
    dataPlot(0,1) = (*q)(0);
    dataPlot(0,2) = (*v)(0);
    dataPlot(0, 3) = (*p)(0);
    dataPlot(0, 4) = (*lambda)(0);
    dataPlot(0,7) = (*q)(nDof-1);
    dataPlot(0,8) = (*v)(nDof-1);
    dataPlot(0,9) = (*q)((nDof)/2);
    dataPlot(0,10) = (*v)((nDof)/2);


    SP::SiconosVector tmp(new SiconosVector(nDof));


//    std::cout <<"potentialEnergy ="<<potentialEnergy << std::endl;
//     std::cout <<"kineticEnergy ="<<kineticEnergy << std::endl;



    // --- Time loop ---
    cout << "====> Start computation ... " <<endl<<endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();

//    while (s->nextTime() < T)
    while(k < N)
    {
      s->computeOneStep();

//       std::cout << "position"  << std::endl;
//       q->display();
//       std::cout << "velocity"  << std::endl;
//       v->display();

// --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k,1) = (*q)(0);
      dataPlot(k,2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0)/h;
      dataPlot(k, 4) = (*lambda)(0);
      dataPlot(k,7) = (*q)(nDof-1);
      dataPlot(k,8) = (*v)(nDof-1);
      dataPlot(k,9) = (*q)((nDof)/2);
      dataPlot(k,10) = (*v)((nDof)/2);




//      std::cout << "q" << std::endl;
//       q->display();


//       std::cout <<"potentialEnergy ="<<potentialEnergy << std::endl;
//       std::cout <<"kineticEnergy ="<<kineticEnergy << std::endl;



      s->nextStep();
      ++show_progress;
      k++;
    }
    cout<<endl << "End of computation - Number of iterations done: "<<k-1<<endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout<<"====> Output file writing ..."<<endl;
    ioMatrix::write("LinearNDofOscillator.dat", "ascii", dataPlot,"noDim");
  }

  catch(SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch(...)
  {
    cout << "Exception caught in ImpactingBarTS.cpp" << endl;
  }



}
