#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

#include <Eigen/Dense>


#include "header.h"








int main(int argc, char* argv[]){
  // int main(){
  std::cout << std::scientific << std::setprecision(25);
  std::clog << std::scientific << std::setprecision(25);

  // --------------------

  const int nptsx = 40;
  const int nptsy = 40;
  const double minOx = -0.2;
  const double maxOx = 0.6;

  const std::string basedir = "/mnt/hdd_barracuda/llnl/reweight/data/32b_v12/";

  const int irow1 = nptsy/2-1;
  std::string mass = "0p4000";
  // double mass_dummy=0.3;
  if(argc>=5) mass = argv[4];

  const int Nt=8;
  const int Ns=32;


  std::vector<double> xpts(nptsx);
  const double deltax = (maxOx-minOx)/nptsx;
  for(int iptx=0; iptx<nptsx; iptx++) xpts[iptx] = minOx + deltax*(iptx+0.5);

  std::vector<double> betas;
  {
    const std::string filename = basedir + "/m"+mass+"avghist_ibx"+std::to_string(0)+"_iby"+std::to_string(0)+"_nojkmeas.bin";
    const HighFive::File f(filename.c_str(), HighFive::File::ReadOnly);
    betas = f.getDataSet("beta").read<std::vector<double>>();
  }
  const int nbeta_meas = betas.size();


  Eigen::ArrayXXd Gamma_list;
  get_Gamma( Gamma_list, nbeta_meas, nptsx, nptsy, irow1, basedir, mass, Ns, Nt );


  int ibeta =700;
  int expn=8;
  if(argc>=2) ibeta = atoi(argv[1]);
  if(argc>=3) expn = atoi(argv[2]);
  const double delta = 1.0*std::pow(10,-expn);
  const Eigen::ArrayXd yy = Gamma_list.block<1,nptsx>( ibeta, 0 );


  {
    // interpolate
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, nptsx);
    gsl_spline_init(spline, xpts.data(), yy.data(), nptsx);

    auto Veff = [&](const double x) { return gsl_spline_eval(spline, x, acc); };
    auto Veff_prime = [&](const double x) { return gsl_spline_eval_deriv(spline, x, acc); };
    auto Veff_prime2 = [&](const double x) { return gsl_spline_eval_deriv2(spline, x, acc); };

    double xA, xB;
    get_local_minima( xA, xB, xpts, deltax, Veff_prime, Veff_prime2);

    auto Veff_sub = [&](const double x) { return Veff(x) - Veff(xB); };

    double tmax=12.;
    if(argc>=4) tmax = atof(argv[3]);
    std::cout << "# tmax = " << tmax << std::endl;

    double dq_init = 0.125*std::pow(10,-6);
    if(argc>=6) dq_init = atof(argv[5]);
    std::cout << "# dq_init = " << dq_init << std::endl;

    int iter_max=1000;
    if(argc>=7) iter_max = atoi(argv[6]);
    std::cout << "# iter_max = " << iter_max << std::endl;

    double tau=1.0e-4;
    if(argc>=8) tau = atof(argv[7]);
    std::cout << "# tau = " << tau << std::endl;

    const std::string dir = "./fit_params_32c_m"+mass+"/";
    std::filesystem::create_directory( dir );

    std::vector<double> ts, qs;
    std::cout << "# debug. search" << std::endl;
    const double c = search_root( ts, qs, dq_init, delta, xA, xB, Veff_prime, tmax, tau, iter_max );
    std::cout << "# dq_init = " << dq_init << std::endl;
    std::cout << "# c = " << c << std::endl;

    {
      char buffer[50];  // maximum expected length of the float
      std::snprintf(buffer, 50, "dq_init: %.15f", dq_init);
      std::string str(buffer);

      const std::string out = dir+"sol_"+std::to_string(ibeta)+".dat";
      write( out, ts, qs, str );
    }

    {
      const std::string out = dir+"c_"+std::to_string(ibeta)+".dat";
      std::ofstream ofs(out);
      ofs << std::scientific << std::setprecision(25) << std::abs(c) << std::endl;
    }
    {
      const std::string out = dir+"dq_init_"+std::to_string(ibeta)+".dat";
      std::ofstream ofs(out);
      ofs << std::scientific << std::setprecision(25) << dq_init << std::endl;
    }


    // bounce solution
    std::vector<double> rs, dIs;
    {
      // integrand
      for(int i=0; i<ts.size()-1; i++) {
        const double r = 0.5*(ts[i]+ts[i+1]);
        const double x = 0.5*(qs[i]+qs[i+1]);
        const double V = Veff_sub(x);
        const double dxdr = ( qs[i+1]-qs[i] )/tau;

        const double dI = r*r * tau * (dxdr*dxdr + V);
        dIs.push_back( dI );
        rs.push_back( r );
      }

      const std::string out = dir+"integrand_"+std::to_string(ibeta)+".dat";
      write( out, rs, dIs, "no 4Pi" );
    }
    {
      double sum = 0.0;
      for(const double elem : dIs) sum += elem;
      sum *= 4.0*M_PI;

      const std::string out = dir+"Scl_"+std::to_string(ibeta)+".dat";
      std::ofstream ofs(out);
      ofs << std::scientific << std::setprecision(25) << sum << std::endl;
    }

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }


  return 0;
}
