// User-defined main parameters
unsigned int nDof = 3;// degrees of freedom for the beam
double t0 = 1e-8;                   // initial computation time
double T = 50.;                  // final computation time
double h = 1e-2;                // time step
double position_init = 01e-6;      // initial position
double velocity_init =  -.0;      // initial velocity
double epsilon = 0.0;//1e-1;
double theta = 1/2.0 + epsilon;              // theta for MoreauJeanOSI integrator
//theta = 1.0;
double g = 9.81; // Gravity
g=0.0;

