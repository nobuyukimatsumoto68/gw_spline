#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdio>
#include <complex>

// #include <hdf5.h>
// #include "highfive/H5File.hpp"
// #include <highfive/H5File.hpp>
// #include <highfive/H5DataSet.hpp>
// #include <highfive/H5DataSpace.hpp>

#include <algorithm>
#include <filesystem>


#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>


using Idx=std::size_t;
using V3 = Eigen::Vector3d;
using M3 = Eigen::Matrix3d;
using Complex = std::complex<double>;

constexpr int NPatches = 10;
constexpr Idx NIcosVertices = 12;

using PatchEdges = std::array<std::array<Idx,4>, NPatches>;
using IcosVertices = std::array<V3, NIcosVertices>;
using Vertices = std::vector<V3>;

inline V3 project( const V3& x ) { return x/x.norm(); }


V3 circumcenter(const V3& r0, const V3& r1, const V3& r2){
  const V3 r10 = r1 - r0;
  const V3 r20 = r2 - r0;

  const V3 tmp1 = r10.squaredNorm() * r20 - r20.squaredNorm() * r10;
  const V3 cross = r10.cross(r20);
  const V3 numer = tmp1.cross(cross);
  const double denom = 2.0*cross.squaredNorm();

  return numer/denom + r0;
}


constexpr int NIh = 120;


struct Rotation {
  std::array<M3, 3> J;

  Rotation(){
    // eps_{1ij}
    J[0] << 0., 0., 0.,
      0., 0., 1.,
      0., -1., 0.;
    // eps_{2ij}
    J[1] << 0., 0., -1.,
      0., 0., 0.,
      1., 0., 0.;
    // eps_{3ij}
    J[2] << 0., 1., 0.,
      -1., 0., 0.,
      0., 0., 0.;
  }

  M3 operator()( const V3& theta ) const {
    M3 expn = M3::Zero();
    for(int i=0; i<3; i++){
      expn -= theta[i]*J[i];
    }
    return expn.exp();
  }

};


struct FullIcosahedralGroup {

  using MultTable = std::vector<std::vector<int>>;
  using Rep = std::vector<Eigen::MatrixXcd>;

  MultTable tab;
  const int is;
  const int it;
  const int iu;

  std::vector<std::string> words;

  Rep reg;

  FullIcosahedralGroup( const std::string& tablefile,
                        const int is_,
                        const int it_,
                        const int iu_
                        )
    : is(is_)
    , it(it_)
    , iu(iu_)
  {
    ReadMultTable( tablefile );
    GetRegRepFromMultTable();

    const int ist = tab[is][it];

    {
      auto check = reg[is]*reg[is];
      double norm = ( check-Eigen::MatrixXcd::Identity(NIh,NIh) ).norm();
      assert( norm<1.0e-12 );
    }
    {
      auto check = reg[it]*reg[it]*reg[it];
      double norm = ( check-Eigen::MatrixXcd::Identity(NIh,NIh) ).norm();
      assert( norm<1.0e-12 );
    }
    {
      auto check = reg[ist]*reg[ist]*reg[ist]*reg[ist]*reg[ist];
      double norm = ( check-Eigen::MatrixXcd::Identity(NIh,NIh) ).norm();
      assert( norm<1.0e-12 );
    }
    {
      auto check = reg[iu]*reg[is]*reg[iu];
      double norm = ( check-reg[is] ).norm();
      assert( norm<1.0e-12 );
    }
    {
      auto check = reg[iu]*reg[it]*reg[iu];
      double norm = ( check-reg[it] ).norm();
      assert( norm<1.0e-12 );
    }

    GetIndependentWords();
  }


  void ReadMultTable( const std::string& filename ){
    tab.clear();

    // std::cout << "# reading simplicial points" << std::endl;
    std::ifstream file(filename);

    std::string str;
    while (std::getline(file, str)){
      std::istringstream iss(str);
      std::vector<int> row;
      for(int i=0; i<NIh; i++){
        int j; iss >> j;
        row.push_back(j);
      }
      tab.push_back( row );
    }
    assert( tab.size()==NIh );
  }


  void GetRegRepFromMultTable(){
    reg.clear();
    reg.resize(NIh);
    for(int g=0; g<NIh; g++){
      reg[g] = Eigen::MatrixXcd::Zero( NIh, NIh );
      for(int i=0; i<NIh; i++){
        reg[g]( tab[g][i], i ) = 1.0;
      }
    }
  }


  void GenerateRep(Rep& rep,
                   const Eigen::MatrixXcd& ms,
                   const Eigen::MatrixXcd& mt
                   ){
    const int nd=ms.cols();

    rep.clear();
    rep.resize(NIh);

    const auto mu = Eigen::MatrixXcd::Identity( nd, nd );

    for(int ig=0; ig<NIh; ig++){
      rep[ig] = Eigen::MatrixXcd::Identity( nd, nd );
      for(const char y : words[ig]) rep[ig] = rep[ig] * replace( y, ms, mt, mu );
    }

    // check
    for(int i=0; i<NIh; i++){
      for(int g=0; g<NIh; g++){
        assert( (rep[ tab[g][i] ] - rep[g] * rep[i]).norm()<1.0e-12 );
      }}
  }


  int replace( const char x ) const {
    if(x=='s') return is;
    else if(x=='t') return it;
    else if(x=='u') return iu;
    else {
      assert(false);
      return 0;
    }
  }

  Eigen::MatrixXcd replace( const char x,
              const Eigen::MatrixXcd& ms,
              const Eigen::MatrixXcd& mt,
              const Eigen::MatrixXcd& mu ) const {
    if(x=='s') return ms;
    else if(x=='t') return mt;
    else if(x=='u') return mu;
    else {
      assert(false);
      return Eigen::MatrixXcd::Zero(ms.rows(),ms.cols());
    }
  }

  M3 replace( const char x,
              const M3& ms,
              const M3& mt,
              const M3& mu ) const {
    if(x=='s') return ms;
    else if(x=='t') return mt;
    else if(x=='u') return mu;
    else {
      assert(false);
      return M3::Zero();
    }
  }

  void GetIndependentWords(){
    std::vector<std::string> words_tmp;
    std::vector<int> elems;

    words_tmp.push_back("");
    elems.push_back(0);

    words_tmp.push_back("u");
    elems.push_back(iu);

    int counter=0;
    bool flag=true;

    while(flag){
      bool flag2=false;

      for( const char x : {'s', 't'} ){
        std::string word = words_tmp[counter] + x;
        std::vector<int> trsl; // (word.size());
        for(const char y : word) trsl.push_back( replace(y) );

        int elem = 0;
        for(int i=0; i<trsl.size(); i++){
          elem = tab[elem][trsl[i]];
        }

        if(std::find(elems.begin(), elems.end(), elem) == elems.end()) {
          words_tmp.push_back(word);
          elems.push_back(elem);
          flag2 = flag2 || true;
        }
      } // end for x

      flag = flag2 || (counter!=words_tmp.size()-1);
      counter++;
    }

    words.clear();
    words.resize(words_tmp.size());

    for(int i=0; i<elems.size(); i++){
      words[elems[i]] = words_tmp[i];
    }
  }

  std::string print() const {
    std::stringstream ss;
    Idx counter=0;
    for(auto word : words){
      ss << counter << ": " << word << std::endl;
      counter++;
    }
    return ss.str();
  }

  M3 rotation( const int i,
               const M3& ms,
               const M3& mt ) const {
    assert( i<NIh );
    const M3 mu = -M3::Identity();

    M3 mat = M3::Identity();
    for(const char y : words[i]) mat = mat*replace(y, ms, mt, mu);

    return mat;
  }

};








struct RefinedIcosahedron {
  const double varphi = 0.5*( 1.0 + std::sqrt(5.0) );
  const double c1 = 0.5/varphi;
  const double s1 = 0.5*std::sqrt(std::sqrt(5.)*varphi);
  const double c2 = -0.5*varphi;
  const double s2 = 0.5*std::sqrt(std::sqrt(5.)/varphi);
  const double c3 = 1.0/std::sqrt(5.);
  const double s3 = 2.0/std::sqrt(5.);

  const Idx L;
  IcosVertices icos;
  PatchEdges patches;
  Vertices vertices;

  // -----------------------------

  RefinedIcosahedron( const Idx L_ )
    : L(L_)
  {
    set_IcosVertices();
    set_PatchEdges();
    FillVerticesSimpleRefinement();
  }

  // -----------------------------

  Idx NVertices() const { return 10*L*L+2; }
  Idx NLinks() const { return 30*L*L; }
  Idx NFaces() const { return 20*L*L; }

  std::string print() const {
    std::stringstream ss;
    Idx counter=0;
    ss << std::scientific << std::setprecision(25);
    for(const V3& x : vertices){
      ss << counter << ": " << x.transpose() << std::endl;
      counter++;
    }
    return ss.str();
  }

  // first 10L^2 is for the regular points
  Idx idx( const int s, const int xl0, const int xl1, const int L ) const {
    assert( xl0 < L );
    assert( xl1 < L );
    return s*L*L + xl1*L + xl0;
  }

  Idx idxS( const int L ) const {
    return 10*L*L;
  }

  Idx idxN( const int L ) const {
    return 10*L*L + 1;
  }

  // -----------------------------

  void set_IcosVertices(){
    icos[0] << 0.,0.,-1.;

    icos[1] << -s3, 0., -c3;
    icos[2] << -s3*c1, -s3*s1, -c3;
    icos[3] << -s3*c2, -s3*s2, -c3;
    icos[4] << -s3*c2, s3*s2, -c3;
    icos[5] << -s3*c1, s3*s1, -c3;

    icos[6] << s3*c2, -s3*s2, c3;
    icos[7] << s3*c1, -s3*s1, c3;
    icos[8] << s3, 0., c3;
    icos[9] << s3*c1, s3*s1, c3;
    icos[10] << s3*c2, s3*s2, c3;

    icos[11] << 0.,0.,1.;
  }

  void set_PatchEdges(){
    patches[0] = std::array<Idx,4>{1, 0, 6, 2};
    patches[1] = std::array<Idx,4>{2, 0, 7, 3};
    patches[2] = std::array<Idx,4>{3, 0, 8, 4};
    patches[3] = std::array<Idx,4>{4, 0, 9, 5};
    patches[4] = std::array<Idx,4>{5, 0, 10,1};

    patches[5] = std::array<Idx,4>{6, 2, 11, 7};
    patches[6] = std::array<Idx,4>{7, 3, 11, 8};
    patches[7] = std::array<Idx,4>{8, 4, 11, 9};
    patches[8] = std::array<Idx,4>{9, 5, 11, 10};
    patches[9] = std::array<Idx,4>{10,1, 11, 6};
  }

  // -----------------------------

  void FillVerticesSimpleRefinement() { // (pA, pB, pD), (pA, pD, pC): right-handed
    vertices.clear();
    vertices.resize( NVertices() );

    for(int s=0; s<NPatches; s++){
      const V3 pA = icos[patches[s][0]];
      const V3 pB = icos[patches[s][1]];
      const V3 pC = icos[patches[s][2]];
      const V3 pD = icos[patches[s][3]];

      //      pD
      //      /|
      //     / |
      //    /  |
      //   /   |
      //  /    |
      // pA----pB
      const V3 vAB = pB-pA;
      const V3 vBD = pD-pB;
      for(Idx i=0; i<L; i++){
        const double ri = 1.0*i/L;
        for(Idx j=0; j<=i; j++){
          const double rj = 1.0*j/L;
          const V3 r = pA + ri*vAB + rj*vBD;
          vertices[idx(s,i,j,L)] = project( r );
        }}

      // pC----pD
      // |     /
      // |    /
      // |   /
      // |  /
      // | /
      // pA
      const V3 vCD = pD-pC;
      const V3 vAC = pC-pA;
      for(Idx j=0; j<L; j++){
        const double rj = 1.0*j/L;
        for(Idx i=0; i<j; i++){
          const double ri = 1.0*i/L;
          const V3 r = pA + ri*vCD + rj*vAC;
          vertices[idx(s,i,j,L)] = project( r );
        }}
    }

    vertices[idxS(L)] = icos[0];
    vertices[idxN(L)] = icos[NIcosVertices-1];
  }

};


struct RefinedIcosahedronDual {
  const RefinedIcosahedron& lattice;

  Vertices vertices;

  RefinedIcosahedronDual(const RefinedIcosahedron& lattice_)
    : lattice(lattice_)
  {
    FillDualPoints();
  }

  Idx NVertices() const { return lattice.NFaces(); }
  // Idx NLinks() const { return 30*L*L; }
  // Idx NFaces() const { return 20*L*L; }


  void FillDualPoints(){
    vertices.clear();

    for(int s=0; s<NPatches; s++){
      // interior
      for(Idx i=0; i<L-1; i++){
        for(Idx j=0; j<L-1; j++){
          V3 r0 = vertices[idx(s,i,j,L)];
          V3 r1; = vertices[idx(s,i,j,L)];
          V3 r2; = vertices[idx(s,i,j,L)];

          V3 r = project( circumcenter(r0, r1, r2) );
          vertices.push_back( r );
        }
        // j=L-1
        vertices[idx(s,i,j,L)];
      }
      // i=L-1
    }

    assert( vertices.size()==NVertices() );
  }

};



template<class Lattice>
struct Orbits {
  const Lattice& lattice;
  const FullIcosahedralGroup& Ih;
  const Rotation& rot;

  M3 ms, mt;
  std::vector<std::vector<Idx>> orbits;

  Orbits(const Lattice& lattice_,
         const FullIcosahedralGroup& Ih_,
         const Rotation& rot_,
         const std::array<int, 2>& link=std::array<int, 2>{0,1},
         const int ipatch=0 )
    : lattice(lattice_)
    , Ih(Ih_)
    , rot(rot_)
  {
    {
      // s
      const V3 p0 = lattice.icos[link[0]];
      const V3 p1 = lattice.icos[link[1]];
      V3 theta = project( p0+p1 );
      theta *= M_PI;
      ms = rot(theta);
      assert( (ms*ms-M3::Identity()).norm()<1.0e-13 );
    }
    {
      // t
      const V3 pa = lattice.icos[lattice.patches[ipatch][0]];
      const V3 pb = lattice.icos[lattice.patches[ipatch][1]];
      const V3 pd = lattice.icos[lattice.patches[ipatch][3]];
      V3 theta = project( pa+pb+pd );
      theta *= 2.0*M_PI/3.0;

      mt = rot(theta);
      assert( (mt*mt*mt-M3::Identity()).norm()<1.0e-13 );
    }
    {
      const M3 mst = ms*mt;
      assert( (mst*mst*mst*mst*mst-M3::Identity()).norm()<1.0e-13 );
    }

    GetOrbits();
  }

  bool is_identical( const V3& x, const V3& y, const double tol=1.0e-12 ) const {
    const double norm = (x-y).norm();
    return norm < tol;
  }

  void GetOrbits(){
    orbits.clear();
    for(const V3& x : lattice.vertices){
      std::vector<Idx> orbit;

      for(int ig=0; ig<NIh; ig++){
        const M3 mg = Ih.rotation(ig, ms, mt);
        const V3 gx = mg*x;

        for(int iy=0; iy<lattice.NVertices(); iy++){
          const V3 y = lattice.vertices[iy];
          if( is_identical(gx, y) ) {
            orbit.push_back(iy);
            break;
          }
        }
      }
      assert( orbit.size()==NIh );
      orbits.push_back(orbit);
    }
  }

  std::string print() const {
    std::stringstream ss;
    Idx counter=0;
    for(const auto& orbit : orbits){
      ss << counter << ": ";
      for(const auto& ix : orbit){
        ss << ix << " ";
      }
      ss << std::endl;
      counter++;
    }
    return ss.str();
  }
};


struct FundamentalOptimizer{
};




// int main(int argc, char* argv[]){
int main(){
  std::cout << std::scientific << std::setprecision(25);
  std::clog << std::scientific << std::setprecision(25);

  // --------------------
  const Idx L = 2;

  RefinedIcosahedron lattice(L);
  std::cout << "# points" << std::endl
            << lattice.print() << std::endl;

  FullIcosahedralGroup Ih( "multtablemathematica.dat",
                           3, 19, 60 );

  Rotation rot;

  FullIcosahedralGroup::Rep rep;
  std::cout << "# Ih elems" << std::endl;
  std::cout << Ih.print() << std::endl;



  Orbits orbits(lattice, Ih, rot);
  std::cout << "# orbits" << std::endl;
  std::cout << orbits.print() << std::endl;

  return 0;
}
