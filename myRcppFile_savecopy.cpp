// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <Rcpp.h> /*integrating C++ with R */
// #include <iostream>
#include <math.h>
#include <stdlib.h>     /* srand, rand */

#include <random> // for Gaussian distribution 
#include <iostream>
//#include <RcppThread.h>
#include <thread>  // C++11 threads
#include <vector>
#include <chrono> /*framework for working with time and dates.*/
#include <list> /* standard tempelate library container */
#include<vector> /* standard tempelate library container */
#include <fstream>
#include<numeric> /* Offers general numeric operations like std::accumulate (sums up elements in a range) and std::inner_product.*/
#include<string> /* The <string> header provides the std::string class, a powerful and flexible way to work with text.*/
#include <complex> /*Supports complex number arithmetic.*/
// #include "utility.cpp"

using std::vector;
using std::string;
using Rcpp::Rcout; /* Rcpp::Rcout
 is an output stream in C++ code using the Rcpp package that allows messages to be synchronized with R's console output. */


const std::string currentDateTime();
const std::string currentDateTime(double time_offset);
double calculate_stddev(std::vector<double>& v);
double calculate_radial_amplitude(std::vector<std::complex<float>>& x, double& dt, double& fcenter); /* Function declaration with arguments x as vector, and dt and fcenter as scalar */
/*This function is intended to compute the amplitude of radial motion from a complex time-domain signal, using the sampling interval dt and a chosen center frequency fcenter to isolate the desired mode.*/

// Some constants unlikely to be changed
const double pi = 3.14159265359;
const double q = 1.602176634e-19;
const double kb = 1.380649e-23;
const double c = 299792458.0;
double m_e1 = 9.1093837139e-31;  // electron mass removed, but not electron binding energy
double m_e2 = 9.1093837139e-31;
const double hbar = 1.05457182e-34;


// Constants for leapfrog algorithm /* Assuming that these parameter do not change on changing the particles to electrons */
const double rt32 = pow(2.0, 1.0 / 3);
const double c1 = 1.0 / 2 * 1 / (2 - rt32);
const double c2 = 1.0 / 2 * (-rt32 / (2 - rt32) + 1 / (2 - rt32));
const double c3 = c2;
const double c4 = c1;

const double d1 = 1.0 / (2 - rt32);
const double d2 = -rt32 / (2 - rt32);
const double d3 = d1;

// e2low may need rework
const double f0_e2Excitation_undetuned = 957.39680e12; /* Natural optical frequency of e2+ ion */
//const double photon_momentum = 2.116056688185771e-27; /* h*f0_e2Excitation_undetuned/c - This constant is the momentum carried by one e2⁺ cooling photon, used to compute the recoil kick given to the ion during laser absorption and emission.*/

double Fheating_p = 0;

struct e2_particle { /*e2_particle holds the complete state of one e2⁺ ion at a single instant in time, and its values are updated as the simulation advances.*/
double x, y, z;
  double vx, vy, vz;
  double vz_prev; 
  bool isExcited;
};

double gete2TotalEnergy(const e2_particle& e2, double qC2U0_e2) {
  return 0.5 * m_e2 * e2.vz * e2.vz
  + qC2U0_e2 * e2.z * e2.z;
}






class SimulationRun {
public:
  SimulationRun();
  void set(std::string msg) { this->msg = msg; }
  std::string greet() { return msg; }
  double getL() { return this->L; }
  
  std::string msg;
  
  Rcpp::NumericMatrix performLoop();
  // void performLeapfrog_x(unsigned long long int& i, double& x, double& y, double& z, double& vx, double& vy, double& vz, double& qC2U0, double& U0, double& m, double& C4, double& B2,
  //     double& E0_sbdrive_magnetron, double& E0_sbdrive_cyclotron, double& f_minus, double& f_z, double& f_plus);
  // void performLeapfrog_y(double& x, double& y, double& z, double& vx, double& vy, double& vz, double& qC2U0, double& U0, double& m, double& C4, double& B2);
  // void performLeapfrog_z(unsigned long long int &i, double &x, double &y, double& z, double &vx, double &vy, double& vz, double &qC2U0, double &U0, double &m, double &C4, double &B2,
  //     double &E0_sbdrive_magnetron, double &E0_sbdrive_cyclotron, double &f_minus, double &f_z, double &f_plus, bool &bool_calculate_xy, double &Find, double &Fheating);
  
  
  void performLeapfrog_z(const unsigned long long int& i, const double& x, const double& y, double& z, const double& vx, const double& vy, double& vz, const double& qC2U0,
                         const double& U0, const double& m, const double& C4, const double& B2,
                         const double& E0_sbdrive_magnetron, const double& E0_sbdrive_cyclotron, const double& f_minus, const double& f_z, const double& f_plus, const bool& bool_calculate_xy,
                         const double& Find, const double& Fheating);
  
  
  void updateVariables();
  
  // When adding new variable, it must be added 2-3x: here, in the constructor and possibly in rcpp_module
  
  // Simulation run properties
  unsigned long long int number_of_steps;
  double dt;
  unsigned long int save_every_xth;
  bool bool_simulate_e1; /* declaring a boolean variable, if true, proton dynamics will be simulated */
bool bool_simulate_e2; /* declaring a boolean variable, if true, e2 dynamics will be simulated */
//bool bool_simulate_resonator;
//bool bool_calculate_xy_e1; /* declaring a boolean variable, if true, radial motion of proton will be evaluated for simulation */
//bool bool_calculate_xy_e2;
int threadNo; /* nteger label used to distinguish parallel simulation runs, control their start timing, random seeds, and logging output.*/
double t_heating_on = 0; /* declares and initializes a time-control parameter that determines when external heating is turned on during the simulation.*/

// Particle properties
double T_init_e1;
double phase_ini_axial_e1;
double p_absdetuning;
//double rho_ini_minus_e1;
//double rho_ini_plus_e1;
//double phase_ini_minus_e1;
//double phase_ini_plus_e1;


double T_init_e2;
double phase_ini_axial_e2;
double e2_absdetuning;
unsigned int number_of_e2_ions;
bool bool_calculate_single_e2;
//double rho_ini_minus_e2;
//double rho_ini_plus_e2;
//double phase_ini_minus_e2;
//double phase_ini_plus_e2;
bool bool_phase_ini_axial_e2_random;

// Resonator properties and particle-resonator interaction
double R;
double L;
double C;
double f_det; /* resonance frequency of the axial LC detection circuit and serves as the reference frequency to which the proton and e2⁺ axial motions are tuned or detuned.*/
double T_det; /* effective temperature of detector */
double noise_stddev;
//bool bool_ResonatorForceOn;
double phase_ini_det;
double C_common = 5.5e-12; /* common coupling capacitor */
double Q_C_common_ini = 0;
bool bool_simulate_capacitive_coupling = true;

// Trap properties
double Dz_e1 = 3.2e-3;
double delta_TR_p;
double B0;
double B2_p;
double B2_e2;
double Dz_e2 = 3.2e-3;
double delta_TR_e2;
double stability_voltage_e1 = 0;
double stability_voltage_e2 = 0;
double Dz_cap_p = 3.2e-3;
double Dz_cap_e2 = 3.2e-3;
double C2_e1 = 18510;
double C2_e2 = 18510;
double D4_e1 = 1.283e9;
double D4_e2 = 1.283e9;

// External drives
double heatingCoeff;
double E0_sbdrive_magnetron_p;
double E0_sbdrive_cyclotron_p;
double E0_sbdrive_magnetron_e2;
double E0_sbdrive_cyclotron_e2;






private:
  double Fx_sbdrive_cyclotron_p;
  double Fz_sbdrive_cyclotron_p;
  double Fx_sbdrive_magnetron_p;
  double Fz_sbdrive_magnetron_p;
  
  double Fx_sbdrive_cyclotron_e2;
  double Fz_sbdrive_cyclotron_e2;
  double Fx_sbdrive_magnetron_e2;
  double Fz_sbdrive_magnetron_e2;
  
  double Fx_C4_p;
  double Fx_B2_p;
  double Fy_C4_p;
  double Fy_B2_p;
  double Fz_C4_p;
  double Fz_B2_p;
  
  double Fx_C4_e2;
  double Fx_B2_e2;
  double Fy_C4_e2;
  double Fy_B2_e2;
  double Fz_C4_e2;
  double Fz_B2_e2;
  
  
  // Electron 1 initialization
  double f_z_e1;
  double qC2U0_e1;
  double U0_e1;
  double f_plus_e1;
  double f_minus_e1;
  double wplus_e1;
  double wminus_e1;
  
  // Electron 2 initialization
  double f_z_e2;
  double qC2U0_e2;
  double U0_e2;
  double f_plus_e2;
  double f_minus_e2;
  double wplus_e2;
  double wminus_e2;
  
  // No initialization in constructor needed:
  
  double x1, x2, x3;
  double y1, y2, y3;
  double z1, z2, z3;
  
  double vx1, vx2, vx3;
  double vy1, vy2, vy3;
  double vz1, vz2, vz3;
  
  double acc_x1, acc_x2, acc_x3;
  double acc_y1, acc_y2, acc_y3;
  double acc_z1, acc_z2, acc_z3;
  
  
  // helper variables so call by reference is possible
  double Find_e2_coll;
  double Fheating_e2_coll;
  
  // double point20a, point20b, point21a, point21b, point22a, point22b, point23a, point23b;
  // double block20 = 0, block21 = 0, block22 = 0, block23 = 0;
  
};


SimulationRun::SimulationRun() {
  // initialize the variables with the standard values (sometimes any value)
  // If these are not the ones desired, overwrite them in R
  
  Rcout.precision(8);
  
  // Simulation run properties
  number_of_steps = 10;
  dt = 1e-9;
  save_every_xth = 1;
  bool_simulate_e1 = true;
  bool_simulate_e2 = true;
 // bool_simulate_resonator = false;
  //bool_calculate_xy_e1 = false;
//  bool_calculate_xy_e2 = false;
  threadNo = 0;
  
  // Particle properties
  T_init_e1 = 4;
  phase_ini_axial_e1 = 0;
  p_absdetuning = 0;
  //rho_ini_minus_e1 = 1e-10;
  //rho_ini_plus_e1 = 1e-10;
//  phase_ini_minus_e1 = 45 * pi / 180.0;
//  phase_ini_plus_e1 = 0;
  
  T_init_e2 = 4;
  phase_ini_axial_e2 = 90 * pi / 180.0;
  e2_absdetuning = 0;
//  rho_ini_minus_e2 = 1e-10;
 // rho_ini_plus_e2 = 1e-10;
//  phase_ini_minus_e2 = 45 * pi / 180.0;
//  phase_ini_plus_e2 = 0;
  number_of_e2_ions = 1;
  bool_calculate_single_e2 = true;
  bool_phase_ini_axial_e2_random = false;
  
  // Resonator properties and particle-resonator interaction
  R = 236e6;
  L = 3e-3;
  f_det = 626951.500;
  C = 1.0 / (2 * 4 * pi * pi * f_det * f_det) * (sqrt(1.0 / (L * L) - 4 * pi * pi * f_det * f_det / (R * R)) + 1.0 / L);
  //C = 2.14808329e-11; // frequency shift due to R is accounted for here. f_Det = 626951.500
  T_det = 8;
  phase_ini_det = 0;
  noise_stddev = sqrt(2 * kb * T_det / (R * dt));
//  bool_ResonatorForceOn = false;
  
  
  
  // Trap properties
  delta_TR_p = 0;
  B0 = 1.8;
  B2_p = 0;
  B2_e2 = 0;
  
  delta_TR_e2 = 0;
  
  
  
  
  // External drives
  heatingCoeff = 0;
  E0_sbdrive_magnetron_p = 0;
  E0_sbdrive_cyclotron_p = 0;
  E0_sbdrive_magnetron_e2 = 0;
  E0_sbdrive_cyclotron_e2 = 0;
  
  
  // private members which need to be initialized
  Fx_sbdrive_cyclotron_p = 0;
  Fz_sbdrive_cyclotron_p = 0;
  Fx_sbdrive_magnetron_p = 0;
  Fz_sbdrive_magnetron_p = 0;
  
  Fx_sbdrive_cyclotron_e2 = 0;
  Fz_sbdrive_cyclotron_e2 = 0;
  Fx_sbdrive_magnetron_e2 = 0;
  Fz_sbdrive_magnetron_e2 = 0;
  
  Fx_C4_p = 0;
  Fx_B2_p = 0;
  Fy_C4_p = 0;
  Fy_B2_p = 0;
  Fz_C4_p = 0;
  Fz_B2_p = 0;
  
  Fx_C4_e2 = 0;
  Fx_B2_e2 = 0;
  Fy_C4_e2 = 0;
  Fy_B2_e2 = 0;
  Fz_C4_e2 = 0;
  Fz_B2_e2 = 0;
  
  // e1 initialization
  f_z_e1 = f_det + p_absdetuning;
  qC2U0_e1 = 2 * pi * pi * m_e1 * pow(f_z_e1, 2.0);
  U0_e1 = qC2U0_e1 / (q * C2_e1); // for C4
  f_plus_e1 = 1 / (2 * pi) * (q / m_e1 * B0 / 2.0 + sqrt(pow(q / m_e1 * B0 / 2, 2.0) - qC2U0_e1 / m_e1));
  f_minus_e1 = 1 / (2 * pi) * (q / m_e1 * B0 / 2.0 - sqrt(pow(q / m_e1 * B0 / 2, 2.0) - qC2U0_e1 / m_e1));
  wplus_e1 = f_plus_e1 * 2 * pi;
  wminus_e1 = f_minus_e1 * 2 * pi;
  
  // e2 initialization
  f_z_e2 = f_det + e2_absdetuning;
  qC2U0_e2 = 2 * pi * pi * m_e2 * pow(f_z_e2, 2.0);
  U0_e2 = qC2U0_e2 / (q * C2_e1);
  f_plus_e2 = 1 / (2 * pi) * (q / m_e2 * B0 / 2.0 + sqrt(pow(q / m_e2 * B0 / 2, 2.0) - qC2U0_e2 / m_e2));
  f_minus_e2 = 1 / (2 * pi) * (q / m_e2 * B0 / 2.0 - sqrt(pow(q / m_e2 * B0 / 2, 2.0) - qC2U0_e2 / m_e2));
  wplus_e2 = f_plus_e2 * 2 * pi;
  wminus_e2 = f_minus_e2 * 2 * pi;
  
}


// here all Variables which are determined indirectly are determined again
void SimulationRun::updateVariables() {
  
  C = 1.0 / (2 * 4 * pi * pi * f_det * f_det) * (sqrt(1.0 / (L * L) - 4 * pi * pi * f_det * f_det / (R * R)) + 1.0 / L);
  //C = 2.14808329e-11; // frequency shift due to R is accounted for here. f_Det = 626951.500
  noise_stddev = sqrt(2 * kb * T_det / (R * dt));
  
  // e1 initialization
  f_z_e1 = f_det + p_absdetuning;
  qC2U0_e1 = 2 * pi * pi * m_e1 * pow(f_z_e1, 2.0);
  U0_e1 = qC2U0_e1 / (q * C2_e1); // for C4
  f_plus_e1 = 1 / (2 * pi) * (q / m_e1 * B0 / 2.0 + sqrt(pow(q / m_e1 * B0 / 2, 2.0) - qC2U0_e1 / m_e1));
  f_minus_e1 = 1 / (2 * pi) * (q / m_e1 * B0 / 2.0 - sqrt(pow(q / m_e1 * B0 / 2, 2.0) - qC2U0_e1 / m_e1));
  wplus_e1 = f_plus_e1 * 2 * pi;
  wminus_e1 = f_minus_e1 * 2 * pi;
  
  Rcout << "f_z_e1 = " << f_z_e1 << "\n";
  
  // e2 initialization
  f_z_e2 = f_det + e2_absdetuning;
  qC2U0_e2 = 2 * pi * pi * m_e2 * pow(f_z_e2, 2.0);
  U0_e2 = qC2U0_e2 / (q * C2_e1);
  f_plus_e2 = 1 / (2 * pi) * (q / m_e2 * B0 / 2.0 + sqrt(pow(q / m_e2 * B0 / 2, 2.0) - qC2U0_e2 / m_e2));
  f_minus_e2 = 1 / (2 * pi) * (q / m_e2 * B0 / 2.0 - sqrt(pow(q / m_e2 * B0 / 2, 2.0) - qC2U0_e2 / m_e2));
  wplus_e2 = f_plus_e2 * 2 * pi;
  wminus_e2 = f_minus_e2 * 2 * pi;
}



RCPP_MODULE(myModule) {
  using namespace Rcpp;
  
  
  class_<SimulationRun>("SimulationRun")
    
    .default_constructor()
    
    .method("performLoop", &SimulationRun::performLoop)
    
    // Simulation run properties
    .field("number_of_steps", &SimulationRun::number_of_steps)
    .field("dt", &SimulationRun::dt)
    .field("save_every_xth", &SimulationRun::save_every_xth)
    .field("bool_simulate_e1", &SimulationRun::bool_simulate_e1)
    .field("bool_simulate_e2", &SimulationRun::bool_simulate_e2)

    .field("threadNo", &SimulationRun::threadNo)
    .field("t_heating_on", &SimulationRun::t_heating_on)
    
    // Particle  properties
    .field("T_init_e1", &SimulationRun::T_init_e1)
    .field("phase_ini_axial_e1", &SimulationRun::phase_ini_axial_e1)
    .field("p_absdetuning", &SimulationRun::p_absdetuning)
  
    
    
    
    .field("T_init_e2", &SimulationRun::T_init_e2)
    .field("phase_ini_axial_e2", &SimulationRun::phase_ini_axial_e2)
    .field("e2_absdetuning", &SimulationRun::e2_absdetuning)
    .field("number_of_e2_ions", &SimulationRun::number_of_e2_ions)
    .field("bool_calculate_single_e2", &SimulationRun::bool_calculate_single_e2)

    .field("bool_phase_ini_axial_e2_random", &SimulationRun::bool_phase_ini_axial_e2_random)
    
    // Resonator properties  and particle-trap interaction
    .field("R", &SimulationRun::R)
    .field("L", &SimulationRun::L)
    .field("f_det", &SimulationRun::f_det)
    .field("C", &SimulationRun::C)
   // .field("T_det", &SimulationRun::T_det)
    .field("noise_stddev", &SimulationRun::noise_stddev)
 
    .field("C_common", &SimulationRun::C_common)
    .field("Q_C_common_ini", &SimulationRun::Q_C_common_ini)
    .field("bool_simulate_capacitive_coupling", &SimulationRun::bool_simulate_capacitive_coupling)
    
    
    // Trap properties
    .field("Dz_e1", &SimulationRun::Dz_e1)
    .field("Dz_e2", &SimulationRun::Dz_e2)
    .field("delta_TR_p", &SimulationRun::delta_TR_p)
    .field("delta_TR_e2", &SimulationRun::delta_TR_e2)
    .field("B0", &SimulationRun::B0)
    .field("B2_p", &SimulationRun::B2_p)
    .field("B2_e2", &SimulationRun::B2_e2)
    .field("stability_voltage_e1", &SimulationRun::stability_voltage_e1)
    .field("stability_voltage_e2", &SimulationRun::stability_voltage_e2)
    .field("Dz_cap_p", &SimulationRun::Dz_cap_p)
    .field("Dz_cap_e2", &SimulationRun::Dz_cap_e2)
    .field("C2_e1", &SimulationRun::C2_e1)
    .field("C2_e1", &SimulationRun::C2_e1)
    .field("D4_e1", &SimulationRun::D4_e1)
    .field("D4_e2", &SimulationRun::D4_e2)
    
    
    // External drives
    .field("heatingCoeff", &SimulationRun::heatingCoeff)
    .field("E0_sbdrive_magnetron_p", &SimulationRun::E0_sbdrive_magnetron_p)
    .field("E0_sbdrive_cyclotron_p", &SimulationRun::E0_sbdrive_cyclotron_p)
    .field("E0_sbdrive_magnetron_e2", &SimulationRun::E0_sbdrive_magnetron_e2)
    .field("E0_sbdrive_cyclotron_e2", &SimulationRun::E0_sbdrive_cyclotron_e2)
  
    ;
}


Rcpp::NumericMatrix SimulationRun::performLoop() {
  
  // Some initializations
  if (threadNo < 100) {
    std::this_thread::sleep_for(std::chrono::milliseconds(threadNo * 500));
  }
  
  updateVariables();
  
  Rcout << "Thread " << threadNo << " started,     time: " << currentDateTime() << "\n";
  
  /* Initialize random seed */
  std::random_device rd{};
  std::mt19937 rand_gen{ static_cast<unsigned int>(rd() + time(NULL) + threadNo) };
  // std::mt19937 rand_gen{ static_cast<unsigned int>(threadNo) };
  //  Rcout << "Note: Thread not random \n";
  
  std::normal_distribution<double> gauss_distribution(0.0, 1.0);
  std::uniform_real_distribution<double> uniform_perc_dis(0.0, 100.0);
  
  std::mt19937_64 random_generator;
  //cxx::ziggurat_normal_distribution<double> ziggurat_normal{0.0, 1.0};
  //std::cout << ziggurat_normal(random) << '\n';
  
  double x_e1 = 0;
  double y_e1 = 0;
  double z_e1 = 0;
  double vx_e1 = 0;
  double vy_e1 = 0;
  double vz_e1 = 0;
  
  
  Rcout << "T_init_e1 = " << T_init_e1 << "\n";
  
  
  if (bool_simulate_e1) {
    z_e1 = sqrt(kb * T_init_e1 / qC2U0_e1) * sin(phase_ini_axial_e1); // q C2 U0 zmax^2 = kb T
    vz_e1 = sqrt(2 * kb * T_init_e1 / m_e1) * cos(phase_ini_axial_e1); // 1/2 m vmax^2 = kb T
    
    

  }
  Rcout << "z_e1 = " << z_e1 << "\n";
  
  
  
  //double acc_z_e21, acc_z_e22, acc_z_e23;
  //double z_e21, z_e22, z_e23;
  //double e2.vz1, e2.vz2, e2.vz3;
  
  
  
  e2_particle e2;
  
 
  e2.z = 0;
 
  e2.vz = 0;
  e2.vz_prev = 0;
  
  
  if (bool_phase_ini_axial_e2_random) {
    phase_ini_axial_e2 = uniform_perc_dis(rand_gen) / 100.0 * 2 * pi;
  }
  
  e2.z = sqrt(kb * T_init_e2 / qC2U0_e2) * sin(phase_ini_axial_e2);
  e2.vz = sqrt(2 * kb * T_init_e2 / m_e2) * cos(phase_ini_axial_e2);
  
  
  
  // Forces and resonator stuff
  double Fheating_e2 = 0;
  double Fheating_e2_withSine = 0; // saves computation time
  
  double Ures = 0;
  double Iind_e1 = 0, Iind_e2 = 0;
  double I_n = sqrt(2 * kb * T_det / L) * cos(phase_ini_det);
  double J_n = sqrt(2 * kb * T_det / (L * C)) * sin(phase_ini_det); // J is derivative of I 
  double I_np1 = I_n, J_np1 = J_n;
  double Iext_n = 0;
  //double Iext_np1 = 0;
  
  double Find_p = 0, Find_e2 = 0;
  
  double Inoise = 0;
  
  double C4_p = D4_e1 * delta_TR_p;
  double C4_e2 = D4_e2 * delta_TR_e2;
  
  
  
  double x = 0, y = 0, z = 0;
  double vx = 0, vy = 0, vz = 0;
  

  
  double delta_Ekin = 0;
  
  double Qsum_C_common = Q_C_common_ini;
  double F_capacitive_coupling_p = 0, F_capacitive_coupling_e2 = 0;
  double U_C_common = 0;
  
  double F_voltageJitter_p = 0;
  double U_jitter_p = 0;
  
  double F_voltageJitter_e2 = 0;
  double U_jitter_e2 = 0;
  double sqrt_steps_for_60s_inverse = 1.0 / sqrt(60.0 / dt);
  
  double Iind_cap_p = 0;
  double Iind_cap_e2 = 0;
  
  
  
  Rcpp::NumericMatrix ps_data(number_of_steps / save_every_xth, 7);
  
  Rcpp::checkUserInterrupt();
  
  bool bool_estimate_given = false;
  const clock_t begin_time = clock();
  
  for (unsigned long long int i = 0; i < number_of_steps; i++) {
    
   
    if (i % save_every_xth == 0) {
      
    
      
      ps_data(i / save_every_xth, 0) = z_e1;
      ps_data(i / save_every_xth, 1) = vz_e1;
      ps_data(i / save_every_xth, 2) = Ures;
      ps_data(i / save_every_xth, 3) = I_np1;
      ps_data(i / save_every_xth, 4) = e2.z;
      ps_data(i / save_every_xth, 5) = e2.vz;
      ps_data(i / save_every_xth, 6) = gete2TotalEnergy(e2, qC2U0_e2);
      
      
    
      Rcpp::checkUserInterrupt();
      
      if (!bool_estimate_given) {
        if (float(clock() - begin_time) / CLOCKS_PER_SEC > 5) {
          Rcout << std::setprecision(3) << "Thread " << threadNo << ": In the first 5s " << i / 1e9 << "bio steps were made, estimated total time: " << number_of_steps * 1.0 / i * 5 + 5.0 << "s or "
                << number_of_steps * 1.0 / i * 5.0 / 60 + 5.0 / 60 << "min or " << number_of_steps * 1.0 / i * 5.0 / 3600 + 5.0 / 3600 << "h, ";
          Rcout << " estimated time of finish: " << currentDateTime(float(clock() - begin_time) / CLOCKS_PER_SEC * (number_of_steps * 1.0 / i - 1.0)) << "\n";
          bool_estimate_given = true;
        }
      }
      
    }
 
    
    if (stability_voltage_e1 != 0) {
      U_jitter_p += stability_voltage_e1 * U0_e1 * sqrt_steps_for_60s_inverse * gauss_distribution(rand_gen);
      F_voltageJitter_p = 2*q*C2_e1 * U_jitter_p;
    }
    if (stability_voltage_e2 != 0) {
      U_jitter_e2 += stability_voltage_e2 * U0_e2 * sqrt_steps_for_60s_inverse * gauss_distribution(rand_gen);
     F_voltageJitter_e2 = q*C2_e1 * U_jitter_e2;
    }
    
    
   
    Inoise = noise_stddev * gauss_distribution(rand_gen);
   
    if (bool_simulate_e1) {
      Iind_e1 = q / Dz_e1 * vz_e1;
      Iind_cap_p = q / Dz_cap_p * vz_e1;
    }
    
    if (bool_simulate_e2) {
     
        Iind_e2 = q / Dz_e2 * e2.vz;
        Iind_cap_e2 = q / Dz_cap_e2 * e2.vz;
      
    }
    
    
    I_n = I_np1;
    J_n = J_np1;
    Iext_n = Inoise + Iind_e1 + Iind_e2;
   
  
    if (bool_simulate_capacitive_coupling) {
      Qsum_C_common += (Iind_cap_p - Iind_cap_e2) * dt;
      U_C_common = Qsum_C_common / C_common;
      F_capacitive_coupling_p = q / Dz_cap_p * U_C_common;
      F_capacitive_coupling_e2 = -q / Dz_cap_e2 * U_C_common;
    }
    
  
    if (bool_simulate_e1) {
      x = x_e1;
      y = y_e1;
      z = z_e1;
      vx = vx_e1;
      vy = vy_e1;
      vz = vz_e1;
    
      // ------------------------------------------------------------------------------------------------------------------------- Electron 1 z calculation
      // saves computation time
      if (E0_sbdrive_magnetron_p != 0) {
        //Fx_sbdrive_magnetron = q * E0_sbdrive_magnetron * z * sin(2 * pi * (f_minus + f_z) * dt * i);
        Fz_sbdrive_magnetron_p = q * E0_sbdrive_magnetron_p * x * sin(2 * pi * (f_minus_e1 + f_z_e1) * dt * i);
      }
      
      if (E0_sbdrive_cyclotron_p != 0) {
        //Fx_sbdrive_cyclotron = q * E0_sbdrive_cyclotron * z * sin(2 * pi * (f_plus - f_z) * dt * i);
        Fz_sbdrive_cyclotron_p = q * E0_sbdrive_cyclotron_p * x * sin(2 * pi * (f_plus_e1 - f_z_e1) * dt * i);
      }
      
      if (C4_p == 0) {
        Fz_C4_p = 0;
      }
      if (B2_p == 0 ) {
        Fz_B2_p = 0;
      }
      
      
      
      z1 = z + c1 * vz * dt;
      
      if (C4_p != 0 ) {
        Fz_C4_p = (-4) * q * C4_p * U0_e1 * z1 * z1 * z1;
      }
    
     
      acc_z1 = (q * C2_e1 * (U0_e1 + U_jitter_p) * (-2) * z1 + Find_p + Fz_C4_p + Fz_B2_p + Fz_sbdrive_cyclotron_p + Fz_sbdrive_magnetron_p + F_capacitive_coupling_p) / m_e1;
      vz1 = vz + d1 * acc_z1 * dt;
      z2 = z1 + c2 * vz1 * dt;
      
      // here
      
      if (C4_p != 0 ) {
        Fz_C4_p = (-4) * q * C4_p * U0_e1 * z2 * z2 * z2;
      }
     
    
      acc_z2 = (q * C2_e1 * (U0_e1 + U_jitter_p) * (-2) * z2 + Find_p + Fz_C4_p + Fz_B2_p + Fz_sbdrive_cyclotron_p + Fz_sbdrive_magnetron_p + F_capacitive_coupling_p) / m_e1;
      vz2 = vz1 + d2 * acc_z2 * dt;
      z3 = z2 + c3 * vz2 * dt;
      
      if (C4_p != 0 ) {
        Fz_C4_p = (-4) * q * C4_p * U0_e1 * z3 * z3 * z3;
      }
     
     
      acc_z3 = (q * C2_e1 * (U0_e1 + U_jitter_p) * (-2) * z3 + Find_p + Fz_C4_p + Fz_B2_p + Fz_sbdrive_cyclotron_p + Fz_sbdrive_magnetron_p + F_capacitive_coupling_p) / m_e1;
      vz3 = vz2 + d3 * acc_z3 * dt;
      
      z = z3 + c4 * vz3 * dt;
      vz = vz3;
      
      z_e1 = z;
      vz_e1 = vz;
      
      // ------------------------------------------------------------------------------------------------------------------------------------- Proton z calculation end
      
    }
    
    
    
    // ---- 
    if (bool_simulate_e2) {// Note that for e2ryllium its center-of-mass!
      
   
      
      // effective cloud: see wineland 1975 paper
      if (bool_calculate_single_e2) {
        
        if (heatingCoeff > 0 && i * dt > t_heating_on) { // saves computation time,the sin is very costly (experimentally proven)
          Fheating_e2_withSine = heatingCoeff * sin(2 * 2 * pi * f_det * i * dt);
        }
        
      
        x = e2.x;
        y = e2.y;
        z = e2.z;
        vx = e2.vx;
        vy = e2.vy;
        vz = e2.vz;
        // point4d = clock();
     
        
        Fheating_e2 = Fheating_e2_withSine * e2.z;
        
        e2.vz_prev = e2.vz;
      
        
        // ------------------------------------------------------------------------------------------------------------------------- e2 z calculation
        // saves computation time
        if (E0_sbdrive_magnetron_e2 != 0) {
          //Fx_sbdrive_magnetron = q * E0_sbdrive_magnetron * z * sin(2 * pi * (f_minus + f_z) * dt * i);
          Fz_sbdrive_magnetron_e2 = q * E0_sbdrive_magnetron_e2 * x * sin(2 * pi * (f_minus_e2 + f_z_e2) * dt * i);
        }
        
        if (E0_sbdrive_cyclotron_e2 != 0) {
          //Fx_sbdrive_cyclotron = q * E0_sbdrive_cyclotron * z * sin(2 * pi * (f_plus - f_z) * dt * i);
          Fz_sbdrive_cyclotron_e2 = q * E0_sbdrive_cyclotron_e2 * x * sin(2 * pi * (f_plus_e2 - f_z_e2) * dt * i);
        }
        
        if (C4_e2 == 0) {
          Fz_C4_e2 = 0;
        }
        if (B2_e2 == 0 ) {
          Fz_B2_e2 = 0;
        }
        
        z1 = z + c1 * vz * dt;
        
        if (C4_e2 != 0 ) {
          Fz_C4_e2 = (-4) * q * C4_e2 * U0_e2 * z1 * z1 * z1;
        }
      
        acc_z1 = (q * C2_e1 * (U0_e2 + U_jitter_e2) * (-2) * z1 + Find_e2 + Fz_C4_e2 + Fz_B2_e2 + Fz_sbdrive_cyclotron_e2 + Fz_sbdrive_magnetron_e2 + Fheating_e2 + F_capacitive_coupling_e2) / m_e2;
        vz1 = vz + d1 * acc_z1 * dt;
        z2 = z1 + c2 * vz1 * dt;
        
        // here
        
        if (C4_e2 != 0) {
          Fz_C4_e2 = (-4) * q * C4_e2 * U0_e2 * z2 * z2 * z2;
        }
       
      
        acc_z2 = (q * C2_e1 * (U0_e2 + U_jitter_e2) * (-2) * z2 + Find_e2 + Fz_C4_e2 + Fz_B2_e2 + Fz_sbdrive_cyclotron_e2 + Fz_sbdrive_magnetron_e2 + Fheating_e2 + F_capacitive_coupling_e2) / m_e2;
        vz2 = vz1 + d2 * acc_z2 * dt;
        z3 = z2 + c3 * vz2 * dt;
        
        if (C4_e2 != 0 ) {
          Fz_C4_e2 = (-4) * q * C4_e2 * U0_e2 * z3 * z3 * z3;
        }
      
      
        acc_z3 = (q * C2_e1 * (U0_e2 + U_jitter_e2) * (-2) * z3 + Find_e2 + Fz_C4_e2 + Fz_B2_e2 + Fz_sbdrive_cyclotron_e2 + Fz_sbdrive_magnetron_e2 + Fheating_e2 + F_capacitive_coupling_e2) / m_e2;
        vz3 = vz2 + d3 * acc_z3 * dt;
        
        z = z3 + c4 * vz3 * dt;
        vz = vz3;
        
        e2.z = z;
        e2.vz = vz;
        
        // ------------------------------------------------------------------------------------------------------------------------------------- e2 z calculation end
     
      }
     
    }
    
    
    
    if (i > 0 && i % 1000000000 == 0) {
      Rcout << "Thread " << threadNo << " at step " << i / 1000000000 << "e9, time: " << currentDateTime();
      Rcout << "  ,  estimated time of finish: " << currentDateTime(float(clock() - begin_time) / CLOCKS_PER_SEC * (number_of_steps * 1.0 / i - 1.0)) << "\n";
    }
    
    // test
    
    
    
  } // for-loop end
  
  
  
  Rcout << "Thread " << threadNo << " delta_Ekin/kb = " << delta_Ekin / kb << "\n";
  
  
  
  return ps_data;
  
  
  
  
  
  
}



