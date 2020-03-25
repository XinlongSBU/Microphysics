#ifndef _burn_type_H_
#define _burn_type_H_

#include <AMReX.H>
#include <extern_parameters.H>
#include <network.H>

// A generic structure holding data necessary to do a nuclear burn.

namespace Burn
{
    // Set the number of independent variables -- this should be
    // temperature, enuc + the number of species which participate
    // in the evolution equations.
    constexpr int neqs = 2 + nspec;

    // for dimensioning the Jacobian
#ifdef REACT_SPARSE_JACOBIAN
    constexpr int njrows = NETWORK_SPARSE_JAC_NNZ;
    constexpr int njcols = 1;
#else
    constexpr int njrows = neqs;
    constexpr int njcols = neqs;
#endif

    // Indices of the temperature and energy variables in the work arrays.
    constexpr int itemp = nspec;
    constexpr int ienuc = nspec + 1;

    struct burn_t
    {
        amrex::Real rho;
        amrex::Real T;
        amrex::Real e;
        amrex::Array1D<amrex::Real, 0, Species::NSpec-1> xn;
#if naux > 0
        amrex::Array1D<amrex::Real, 0, naux-1> aux;
#endif

        amrex::Real cv;
        amrex::Real cp;
        amrex::Real y_e;
        amrex::Real eta;
        amrex::Real cs;
        amrex::Real dx;
        amrex::Real abar;
        amrex::Real zbar;

        // Last temperature we evaluated the EOS at
        amrex::Real T_old;

        // Temperature derivatives of specific heat
        amrex::Real dcvdT;
        amrex::Real dcpdT;

        // Whether we are self-heating or not.
        bool self_heat;

        // Zone index information.
        int i, j, k;

        // diagnostics
        int n_rhs, n_jac;

        // Integration time.
        amrex::Real time;

        // Was the burn successful?
        bool success;

        AMREX_GPU_HOST_DEVICE AMREX_INLINE
        void normalize_abundances()
        {
            amrex::Real sum_x = 0.0e0_rt;

            for (int i = 0; i < Species::NSpec; ++i) {
                xn(i) = amrex::max(small_x, amrex::min(1.0e0_rt, xn(i)));
                sum_x += xn(i);
            }

            for (int i = 0; i < Species::NSpec; ++i)
                xn(i) = xn(i) / sum_x;
        }
    };


    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void eos_to_burn(eos_t const& eos_state, burn_t& burn_state)
    {
        // Given an eos type, copy the data relevant to the burn type.

        burn_state.rho  = eos_state.rho
        burn_state.T    = eos_state.T
        burn_state.e    = eos_state.e

        for (int i = 0; i < Species::NSpec; ++i)
            burn_state.xn(i) = eos_state.xn(i);

#if naux > 0
        for (int i = 0; i < naux; ++i)
            burn_state.aux(i) = eos_state.aux(i);
#endif

        burn_state.cv   = eos_state.cv
        burn_state.cp   = eos_state.cp
        burn_state.y_e  = eos_state.y_e
        burn_state.eta  = eos_state.eta
        burn_state.cs   = eos_state.cs
        burn_state.abar = eos_state.abar
        burn_state.zbar = eos_state.zbar
    }


    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void burn_to_eos(burn_t const& burn_state, eos_t& eos_state)
    {
        // Given a burn type, copy the data relevant to the eos type.

        eos_state.rho  = burn_state.rho;
        eos_state.T    = burn_state.T;
        eos_state.e    = burn_state.e;

        for (int i = 0; i < Species::NSpec; ++i)
            eos_state.xn(i) = burn_state.xn(i);

#if naux > 0
        for (int i = 0; i < naux; ++i)
            eos_state.aux(i) = burn_state.aux(i);
#endif

        eos_state.cv   = burn_state.cv;
        eos_state.cp   = burn_state.cp;
        eos_state.y_e  = burn_state.y_e;
        eos_state.eta  = burn_state.eta;
        eos_state.cs   = burn_state.cs;
        eos_state.abar = burn_state.abar;
        eos_state.zbar = burn_state.zbar;
    }
}

#endif
