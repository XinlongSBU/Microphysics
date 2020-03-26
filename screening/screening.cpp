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

  auto& s = scn_facs[i];

  s.z1 = z1;
  s.z2 = z2;
  s.a1 = a1;
  s.a2 = a2;

  s.zs13 = std::pow(s.z1 + s.z2, 1.0_rt/3.0_rt);

  s.zs13inv = 1.0_rt/s.zs13;

  s.zhat  =   std::pow(s.z1 + s.z2, 5.0_rt/3.0_rt)
            - std::pow(s.z1, 5.0_rt/3.0_rt)
            - std::pow(s.z2, 5.0_rt/3.0_rt);

  s.zhat2 =   std::pow(s.z1 + s.z2, 5.0_rt/12.0_rt)
            - std::pow(s.z1, 5.0_rt/12.0_rt)
            - std::pow(s.z2, 5.0_rt/12.0_rt);

  s.lzav = (5.0_rt/3.0_rt) * std::log(s.z1 * s.z2 / (s.z1 + s.z2));

  s.aznut = std::pow(s.z1 * s.z1 * s.z2 * s.z2 *
                     s.a1 * s.a2 / (s.a1 + s.a2), 1.0_rt/3.0_rt);
}

