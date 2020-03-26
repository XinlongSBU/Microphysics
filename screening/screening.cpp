#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_REAL.H>
#include <AMReX_Algorithm.H>
#include <network_properties.H>
#include <screen.H>
#include <screen_data.H>

using namespace amrex;
using namespace scrn;

void screening_init() {

}

void screening_finalize() {

}


void add_screening_factor(const int i,
                          const Real z1, const Real a1, const Real z2, const Real a2) {

  BL_ASSERT(i < NSCREEN);

  scn_facs[i].z1 = z1;
  scn_facs[i].z2 = z2;
  scn_facs[i].a1 = a1;
  scn_facs[i].a2 = a2;

  scn_facs[i].zs13 = std::pow(scn_facs[i].z1 + scn_facs[i].z2, 1.0_rt/3.0_rt);
  scn_facs[i].zs13inv = 1.0_rt/scn_facs[i].zs13;
  scn_facs[i].zhat = std::pow(scn_facs[i].z1 + scn_facs[i].z2, 5.0_rt/3.0_rt) -
                     std::pow(scn_facs[i].z1, 5.0_rt/3.0_rt) - std::pow(scn_facs[i].z2, 5.0_rt/3.0_rt);
  scn_facs[i].zhat2 = std::pow(scn_facs[i].z1 + scn_facs[i].z2, 5.0_rt/12.0_rt) -
                     std::pow(scn_facs[i].z1, 5.0_rt/12.0_rt) - std::pow(scn_facs[i].z2, 5.0_rt/12.0_rt);
  scn_facs[i].lzav = (5.0_rt/3.0_rt) * std::log(scn_facs[i].z1 * scn_facs[i].z2 / (scn_facs[i].z1 + scn_facs[i].z2));
  scn_facs[i].aznut = std::pow(scn_facs[i].z1 * scn_facs[i].z1 *
                              scn_facs[i].z2 * scn_facs[i].z2 *
                              scn_facs[i].a1 * scn_facs[i].a2 /
                              (scn_facs[i].a1 + scn_facs[i].a2), 1.0_rt/3.0_rt);
}

