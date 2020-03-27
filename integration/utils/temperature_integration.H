#ifndef _temperature_integration_H_
#define _temperature_integration_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <extern_parameters.H>
#include <network.H>
#include <burn_type.H>
#include <jacobian_sparsity.H>

namespace TempInt
{
    extern AMREX_GPU_MANAGED bool self_heat;
}

/*
  Sets up the temperature equation. This should be called from
  within the actual_rhs routine but is provided here as a convenience
  since most networks will use the same temperature ODE.
*/

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void temperature_rhs(burn_t const& state, RealArray1D<Burn::neqs>& ydot)
{

    Real cv, cp, cvInv, cpInv;

    if (state.self_heat) {

       /*
       Set up the temperature ODE.  For constant pressure, Dp/Dt = 0, we
       evolve :
          dT/dt = (1/c_p) [ Hnuc ]

       For constant volume, div{U} = 0, and we evolve:
          dT/dt = (1/c_v) [ Hnuc ]

       See low Mach paper III, including Eq. A3 for details.
       Note that we no longer include the chemical potential (dE/dX or dH/dX)
       terms because we believe they analytically should vanish.
       */

       if (do_constant_volume_burn) {

          if (!call_eos_in_rhs && dT_crit < 1.0e19_rt) {

             cv = state.cv + (state.T - state.T_old) * state.dcvdt;

          } else {

             cv = state.cv;

          }

          cvInv = 1.0_rt / cv;

          ydot(Burn::itemp) = ydot(Burn::ienuc) * cvInv;

       } else {

          if (!call_eos_in_rhs && dT_crit < 1.0e19_rt) {

             cp = state.cp + (state.T - state.T_old) * state.dcpdt;

          } else {

             cp = state.cp;

          }

          cpInv = 1.0_rt / cp;

          ydot(Burn::itemp) = ydot(Burn::ienuc) * cpInv;

       }
    }
}


/*
  Sets up the temperature entries in the Jacobian. This should be called from
  within the actual_jac routine but is provided here as a convenience
  since most networks will use the same temperature ODE.
*/
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void temperature_jac(burn_t const& state,
                     RealArray2D<Burn::njrows, Burn::njcols>& jac)
{
    real(rt) :: scratch, cspec, cspecInv

    // Temperature Jacobian elements

    if (state.self_heat) {

       if (do_constant_volume_burn) {

          if (!call_eos_in_rhs && dT_crit < 1.0e19_rt) {

             cspec = state.cv + (state.T - state.T_old) * state.dcvdt;

          } else {

             cspec = state.cv;

          }

       } else {

          if (!call_eos_in_rhs && dT_crit < 1.0e19_rt) {

             cspec = state.cp + (state.T - state.T_old) * state.dcpdt;

          } else {

             cspec = state.cp;

          }

       }

       cspecInv = 1.0_rt / cspec;

       // d(itemp)/d(yi)

       for (int k = 0; k < Species::NumSpec; ++k) {
          get_jac_entry(jac, Burn::ienuc, k, scratch);
          scratch = scratch * cspecInv;
          set_jac_entry(jac, Burn::itemp, k, scratch);
       }

       /*
       d(itemp)/d(temp) -- we get this from the equation for d (denuc / dt) / dT
       since dT/dt = 1/c_x denuc/dt in our formalism
       */

       get_jac_entry(jac, Burn::ienuc, Burn::itemp, scratch);
       scratch = scratch * cspecInv;
       set_jac_entry(jac, Burn::itemp, Burn::itemp, scratch);

       // d(itemp)/d(enuc)

       scratch = 0._rt;
       set_jac_entry(jac, Burn::itemp, Burn::ienuc, scratch);

    }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void temperature_rhs_init()
{

    // Provide a default value, then consult the burning_mode.

    self_heat = true;

    if (burning_mode == 0 || burning_mode == 2) {
        self_heat = false;
    } else if (burning_mode == 1 || burning_mode == 3) {
        self_heat = true;
    } else {
        amrex::Error("Error: unknown burning_mode in temperature_rhs_init()");
    }

}

#endif
