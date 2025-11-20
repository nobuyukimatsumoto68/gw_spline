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








std::vector<double> find_all_roots( const std::function<double(const double)>& f,
                                    const std::function<double(const double)>& f_prime,
                                    const double xmin, const double xmax,
                                    const double resolution,
                                    const int iter_max = 100000,
                                    const double tol=1.0e-12 ){
  std::vector<double> roots;
  for(double x0=xmin; x0<xmax-resolution; x0+=resolution){
    const bool is_pos1 = f(x0) > 0;
    const bool is_pos2 = f(x0+resolution) > 0;
    if( is_pos1 != is_pos2 ){
      double root = x0;
      for(int i=0; i<iter_max; i++){
        const double root_old = root;
        root = root - f(root)/f_prime(root);
        if( std::abs(root-root_old)<tol ) break;
      }
      // std::cout << "root = " << root << " f = " << f(root) << std::endl;
      roots.push_back(root);
    }
  }
  return roots;
}








double force(const double t,
             const double p,
             const double q,
             const std::function<double(const double)>& dV
             ){
  // std::cout << "debug. q = " << q << std::endl
  //           << " dV = " << dV(q) << std::endl;
  return -2.0/t * p + 0.5*dV(q);
}

void leapfrog(const double t,
              double& p,
              double& q,
              const double tau,
              const std::function<double(const double)>& dV
              ){
  double p_half = p + 0.5 * tau * force(t+0.5*tau, p, q, dV);
  double p_half_old = p_half;
  for(int i=0; i<10; i++){
    p_half = p + 0.5 * tau * force(t+0.5*tau, p_half, q, dV);
    if( std::abs(p_half-p_half_old)<1.0e-10) break;
    p_half_old = p_half;
  }

  q = q + tau * p_half;
  p = p_half + 0.5 * tau * force(t+0.5*tau, p_half, q, dV);
}


void integrate( std::vector<double>& ts,
                std::vector<double>& qs,
                const double p_init,
                const double q_init,
                const double xA,
                const double xB,
                const std::function<double(const double)>& dV,
                const double tmax=250.,
                const double tau=1.0e-3
                ){
  ts.clear();
  qs.clear();

  const int nsteps = int(tmax/tau);
  double p = p_init;
  double q = q_init;
  for(int i=0; i<nsteps; i++){
    const double t = i*tau;
    leapfrog(t, p, q, tau, dV);
    ts.push_back(t);
    qs.push_back(q);
    if( q<xA-tau || q>xB+tau ) break;
  }
}


double cost( std::vector<double>& ts,
             std::vector<double>& qs,
             const double q_init,
             const double xA,
             const double xB,
             const std::function<double(const double)>& dV,
             const double tmax=250.,
             const double tau=1.0e-3
             ){
  const double p_init = 0.0;
  integrate( ts, qs, p_init, q_init, xA, xB, dV, tmax, tau );
  // std::cout << "# debug. ts[-1] = " << ts[ts.size()-1] << ", qs[-1] = " << qs[qs.size()-1] << std::endl;
  return qs[qs.size()-1] - xB;
}


double search_root( std::vector<double>& ts,
                    std::vector<double>& qs,
                    double& dq_init,
                    double delta,
                    const double xA,
                    const double xB,
                    const std::function<double(const double)>& dV,
                    const double tmax=250.,
                    const double tau=1.0e-3,
                    const int iter_max=1000,
                    const double tol=1.0e-4,
                    const double tol2=1.0e-17
                    ){
  // assert( tau>tol );
  assert( tau*tau<tol );
  double c = cost( ts, qs, xA+dq_init, xA, xB, dV, tmax, tau );
  double c_old = c;

  for(int i=0; i<iter_max; i++){
    // std::cout << "# debug. pt1 " << std::endl;
    if( c*c_old < 0. ) {
      // std::cout << "# debug. pt1 " << std::endl;
      delta *= 0.5;
    }
    c_old = c;
    // else{
    if(c<0.) {
      while(dq_init-delta < 0.) {
        delta *= 0.5;
        // if(delta<tol2) {
        //   std::cout << "# !!! delta too small break." << std::endl;
        //   break;
        // }
      }
      dq_init -= delta;
      // q_init = std::max(q_init, xA);
      // if( std::abs(q_init-xA)<tol2 ) {
      //   std::cout << "# !!! q_init too close to xA." << std::endl;
      //   break;
      // }
    }
    else dq_init += delta;
    // }
    c = cost( ts, qs, xA+dq_init, xA, xB, dV, tmax, tau );
    std::cout << "# debug. c = " << c
              // << "# debug. c_old = " << c_old
              // << "# bool = " << (c*c_old < 0)
              << " dq_init = " << dq_init
              << " delta = " << delta << std::endl;
    if(std::abs(c)<tol && std::abs(ts[ts.size()-1]-tmax)<tol ) {
      std::cout << "# debug. break." << std::endl;
      break;
    }
    if(delta<tol2) {
      std::cout << "# !!! delta too small break." << std::endl;
      break;
    }
    if( std::abs(c-c_old)<tol2 ){
      std::cout << "# !!! no cost difference break" << std::endl;
      break;
    }
  }

  return c;
}


void write(const std::string& filename,
           const std::vector<double>& ts,
           const std::vector<double>& qs,
           const std::string& comment=""){

  std::ofstream ofs(filename);
  ofs << std::scientific << std::setprecision(25);

  assert( ts.size()==qs.size() );

  ofs << "# " << comment << std::endl;

  for(int i=0; i<ts.size(); i++) {
    ofs << std::setw(50) << ts[i] << "\t"
        << std::setw(50) << qs[i] << std::endl;
  }
}














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

    std::vector<double> tmp;
    betas = f.getDataSet("beta").read<std::vector<double>>();
  }
  const int nbeta_meas = betas.size();


  Eigen::ArrayXXd hist_list(nbeta_meas, nptsx);
  const int ipty = irow1;
  for(int iptx=0; iptx<nptsx; iptx++){
    Eigen::ArrayXd yy1, yy2;
    {
      const std::string filename1 = basedir + "/m"+mass+"avghist_ibx"+std::to_string(iptx)+"_iby"+std::to_string(ipty)+"_nojkmeas.bin";
      const HighFive::File f1(filename1.c_str(), HighFive::File::ReadOnly);

      std::vector<double> tmp;
      Eigen::ArrayXd f, fP;
      tmp = f1.getDataSet("f").read<std::vector<double>>();
      f = Eigen::Map<Eigen::ArrayXd>(tmp.data(), tmp.size());
      tmp = f1.getDataSet("fP").read<std::vector<double>>();
      fP = Eigen::Map<Eigen::ArrayXd>(tmp.data(), tmp.size());
      yy1 = (f-fP).exp();
    }

    {
      const std::string filename2 = basedir + "/m"+mass+"avghist_ibx"+std::to_string(iptx)+"_iby"+std::to_string(nptsy-ipty-1)+"_nojkmeas.bin";
      const HighFive::File f2(filename2.c_str(), HighFive::File::ReadOnly);

      std::vector<double> tmp;
      Eigen::ArrayXd f, fP;
      tmp = f2.getDataSet("f").read<std::vector<double>>();
      f = Eigen::Map<Eigen::ArrayXd>(tmp.data(), tmp.size());
      tmp = f2.getDataSet("fP").read<std::vector<double>>();
      fP = Eigen::Map<Eigen::ArrayXd>(tmp.data(), tmp.size());
      yy2 = (f-fP).exp();
    }

    hist_list.block( 0, iptx, nbeta_meas, 1 ) = 0.5*(yy1+yy2);
  }
  const Eigen::ArrayXXd Gamma_list = -std::pow(1.0*Nt/Ns,3) * hist_list.log();


  int ibeta =700; // # in np.arang(450, 504):
  int expn=8;
  double tmax=12.;
  {
    if(argc>=2) ibeta = atoi(argv[1]);
    if(argc>=3) expn = atoi(argv[2]);
  }
  const double delta = 1.0*std::pow(10,-expn);
  // const int ibeta=854;
  const Eigen::ArrayXd yy = Gamma_list.block<1,nptsx>( ibeta, 0 );
  // for(int iptx=0; iptx<nptsx; iptx++){
  //   std::cout << xpts[iptx] << " " << yy[iptx] << std::endl;
  // }

  {
    // interpolate
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, nptsx);
    gsl_spline_init(spline, xpts.data(), yy.data(), nptsx);

    auto Veff = [&](const double x) {
      return gsl_spline_eval(spline, x, acc);
    };

    auto Veff_prime = [&](const double x) {
      return gsl_spline_eval_deriv(spline, x, acc);
    };

    auto Veff_prime2 = [&](const double x) {
      return gsl_spline_eval_deriv2(spline, x, acc);
    };

    std::vector<double> roots;
    double xA, xB;
    {
      const double xmin= *std::min_element(xpts.begin(), xpts.end()); // *std::min(xpts); // -0.1;
      const double xmax= *std::max_element(xpts.begin(), xpts.end()); // *std::min(xpts); // -0.1;
      const double resolution = 0.25*deltax;
      roots = find_all_roots( Veff_prime, Veff_prime2, xmin, xmax, resolution );

      xA = roots[0];
      xB = roots[roots.size()-1];
      assert( Veff_prime2(xA)>0.0 );
      assert( Veff_prime2(xB)>0.0 );
    }

    auto Veff_sub = [&](const double x) {
      return Veff(x) - Veff(xB);
    };


    const int xsize = 1000;
    const double xmin= *std::min_element(xpts.begin(), xpts.end()); // *std::min(xpts); // -0.1;
    const double xmax= *std::max_element(xpts.begin(), xpts.end()); // *std::min(xpts); // -0.1;
    std::vector<double> xs(xsize+1);
    for(int i=0; i<=xsize; i++) xs[i] = (xmax-xmin)/xsize * i + xmin;

    // // // evaluation
    // for(const double x : xs) std::cout << x << " " << Veff_sub(x) << std::endl;


    {
      if(argc>=4) tmax = atof(argv[3]);
      std::cout << "# tmax = " << tmax << std::endl;

      double dq_init = 0.125*std::pow(10,-6);

      if(argc>=6) dq_init = atof(argv[5]);
      // double q_init = xA + dq_init;
      std::cout << "# dq_init = " << dq_init << std::endl;

      int iter_max=1000;
      if(argc>=7) iter_max = atoi(argv[6]);
      std::cout << "# iter_max = " << iter_max << std::endl;

      double tau=1.0e-4;
      if(argc>=8) tau = atof(argv[7]);
      std::cout << "# tau = " << tau << std::endl;


      const std::string dir = "./fit_params_32c_m"+mass+"/";
      std::filesystem::create_directory( dir );
      double c;
      {
        std::vector<double> ts;
        std::vector<double> qs;
        // const double p_init = 0.0;
        // integrate( ts, qs, p_init, q_init, xA, xB, Veff_prime, tmax, tau );
        // for(int i=0; i<ts.size(); i++){
        //   std::cout << ts[i] << " " << qs[i] << std::endl;
        // }

        std::cout << "# debug. search" << std::endl;
        c = search_root( ts, qs, dq_init, delta, xA, xB, Veff_prime, tmax, tau, iter_max );
        std::cout << "# dq_init = " << dq_init << std::endl;
        std::cout << "# c = " << c << std::endl;



        { // if success
          {
            std::vector<double> Vs;
            for(const double x : xs) Vs.push_back( Veff_sub(x) );

            const std::string out = dir+"Veff_"+std::to_string(ibeta)+".dat";
            write( out, xs, Vs, "VB subtracted" );
          }
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
          std::vector<double> rs;
          std::vector<double> dIs;
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

        }
        // if success
      }

    }

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }


  return 0;
}
