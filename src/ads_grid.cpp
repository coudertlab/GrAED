#include <local/gridcalc.hpp>
#include <chrono>

int main(int argc, char* argv[]) {
  std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
  std::string structure_file = argv[1];
  std::string forcefield_path = argv[2];
  double temperature = std::stod(argv[3]);
  double cutoff = std::stod(argv[4]);
  double cutoff_sq = cutoff*cutoff;
  double cutoff_6 = (cutoff_sq)*(cutoff_sq)*(cutoff_sq);
  double inv_cutoff_6 = 1.0/cutoff_6;
  double inv_cutoff_12 = inv_cutoff_6*inv_cutoff_6;
  std::string element_guest_str = argv[5];
  double approx_spacing = std::stod(argv[6]);
  double energy_threshold = 40;
  double access_coeff = 0.8;
  if (argv[8]==NULL) {energy_threshold = std::stod(argv[7]);}
  if (argv[9]==NULL) {energy_threshold = std::stod(argv[7]);access_coeff = std::stod(argv[8]);}

  // Error catch
  if ( temperature < 0 ) {throw std::invalid_argument( "Received negative value for the Temperature" );}
  if ( energy_threshold < 0 ) {throw std::invalid_argument( "Received negative value for the Energy Threshold" );}
  if ( access_coeff > 1 || access_coeff < 0 ) {throw std::invalid_argument( "Accessibility Coefficient out of range (Read the purpose of this coeff)" );}

  // key constants
  double const R = 8.31446261815324e-3; // kJ/mol/K
  double const N_A = 6.02214076e23;    // part/mol

  // Input
  double molar_mass = 0;
  double boltzmann_energy_lj = 0;
  double boltzmann_energy_sq_lj = 0;
  double boltzmann_energy_3_lj = 0;
  double boltzmann_energy_4_lj = 0;
  double sum_exp_energy = 0;
  double mean_grid = 0, std_grid = 0;
  gemmi::Grid<double> grid;

  // read_grid_ccp4(grid, grid_file);

  make_energy_grid_ads_descriptors(structure_file,forcefield_path,temperature,cutoff,element_guest_str,approx_spacing,
  energy_threshold,access_coeff, molar_mass, boltzmann_energy_lj, boltzmann_energy_sq_lj, boltzmann_energy_3_lj, boltzmann_energy_4_lj, mean_grid, std_grid, sum_exp_energy, grid, true);

  boltzmann_energy_3_lj = boltzmann_energy_3_lj/sum_exp_energy;
  boltzmann_energy_4_lj = boltzmann_energy_4_lj/sum_exp_energy;
  // save_grid_ccp4(grid, structure_file, approx_spacing, energy_threshold, element_guest_str);

  std::string structure_name = trim(structure_file);
  double Framework_density = molar_mass/(N_A*grid.unit_cell.volume*1e-30); // g/m3
  double energy_mean = boltzmann_energy_lj/sum_exp_energy;
  double enthalpy = energy_mean - R*temperature;  // kJ/mol
  double henry = sum_exp_energy/(grid.data.size()*R*temperature*Framework_density);    // mol/kg/Pa
  // Implement a standard deviation
  double enthalpy_std = std::sqrt(boltzmann_energy_sq_lj/sum_exp_energy - energy_mean*energy_mean);
  double enthalpy_skew = (boltzmann_energy_3_lj-energy_mean*energy_mean*energy_mean)/(enthalpy_std*enthalpy_std*enthalpy_std) - 3*energy_mean/enthalpy_std;
  double enthalpy_kurt = (boltzmann_energy_4_lj-4*energy_mean*boltzmann_energy_3_lj+3*energy_mean*energy_mean*energy_mean*energy_mean)/(enthalpy_std*enthalpy_std*enthalpy_std*enthalpy_std) + 6*energy_mean*energy_mean/(enthalpy_std*enthalpy_std);
  auto t_end = std::chrono::high_resolution_clock::now();
  double elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end-t_start).count();
  // Structure name, Enthalpy (kJ/mol), Henry coeff (mol/kg/Pa), Accessible Surface Area (m2/cm3), Time (s)
  std::cout << structure_name << "," << enthalpy << "," << enthalpy_std << "," << enthalpy_skew << "," << enthalpy_kurt << "," << mean_grid << "," << std_grid << "," << henry << "," << elapsed_time_ms*0.001 << std::endl;
}
