module actual_network

  use physical_constants, only: ERG_PER_MeV
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  public

  character (len=32), parameter :: network_name = "pynucastro"

  real(rt), parameter :: avo = 6.0221417930e23_rt
  real(rt), parameter :: c_light = 2.99792458e10_rt
  real(rt), parameter :: enuc_conv2 = -avo*c_light*c_light

  real(rt), parameter :: ev2erg  = 1.60217648740e-12_rt
  real(rt), parameter :: mev2erg = ev2erg * 1.0e6_rt
  real(rt), parameter :: mev2gr  = mev2erg / c_light**2

  real(rt), parameter :: mass_neutron  = 1.67492721184e-24_rt
  real(rt), parameter :: mass_proton   = 1.67262163783e-24_rt
  real(rt), parameter :: mass_electron = 9.10938215450e-28_rt

  integer, parameter :: nrates = 5

  ! Evolution and auxiliary
  integer, parameter :: naux  = 0

  ! Number of nuclear species in the network
  integer, parameter :: nspec = 8

  ! For each rate, we need: rate, drate/dT, screening, dscreening/dT
  integer, parameter :: num_rate_groups = 4

  ! Number of reaclib rates
  integer, parameter :: nrat_reaclib = 5
  integer, parameter :: number_reaclib_sets = 6

  ! Number of tabular rates
  integer, parameter :: nrat_tabular = 0

  ! Binding Energies Per Nucleon (MeV)
  real(rt) :: ebind_per_nucleon(nspec)

  ! aion: Nucleon mass number A
  ! aion_inv: 1 / Nucleon mass number A
  ! zion: Nucleon atomic number Z
  ! nion: Nucleon neutron number N
  ! bion: Binding Energies (ergs)

  ! Nuclides
  integer, parameter :: jn   = 1
  integer, parameter :: jp   = 2
  integer, parameter :: jhe4   = 3
  integer, parameter :: jc12   = 4
  integer, parameter :: jo16   = 5
  integer, parameter :: jne20   = 6
  integer, parameter :: jna23   = 7
  integer, parameter :: jmg23   = 8

  ! Reactions
  integer, parameter :: k_c12_c12__he4_ne20   = 1
  integer, parameter :: k_c12_c12__n_mg23   = 2
  integer, parameter :: k_c12_c12__p_na23   = 3
  integer, parameter :: k_he4_c12__o16   = 4
  integer, parameter :: k_n__p__weak__wc12   = 5

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  real(rt), allocatable, save :: aion(:), aion_inv(:), zion(:), bion(:)
  real(rt), allocatable, save :: nion(:), mion(:), wion(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: aion, aion_inv, zion, bion, nion, mion, wion
#endif

  !$acc declare create(aion, aion_inv, zion, bion, nion, mion, wion)

#ifdef REACT_SPARSE_JACOBIAN
  ! Shape of Jacobian in Compressed Sparse Row format
  integer, parameter   :: NETWORK_SPARSE_JAC_NNZ = 45
  integer, allocatable :: csr_jac_col_index(:), csr_jac_row_count(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: csr_jac_col_index, csr_jac_row_count
#endif
#endif

contains

  subroutine actual_network_init()

    implicit none

    integer :: i

    ! Allocate ion info arrays
    allocate(aion(nspec))
    allocate(aion_inv(nspec))
    allocate(zion(nspec))
    allocate(bion(nspec))
    allocate(nion(nspec))
    allocate(mion(nspec))
    allocate(wion(nspec))

    spec_names(jn)   = "neutron"
    spec_names(jp)   = "hydrogen-1"
    spec_names(jhe4)   = "helium-4"
    spec_names(jc12)   = "carbon-12"
    spec_names(jo16)   = "oxygen-16"
    spec_names(jne20)   = "neon-20"
    spec_names(jna23)   = "sodium-23"
    spec_names(jmg23)   = "magnesium-23"

    short_spec_names(jn)   = "n"
    short_spec_names(jp)   = "h1"
    short_spec_names(jhe4)   = "he4"
    short_spec_names(jc12)   = "c12"
    short_spec_names(jo16)   = "o16"
    short_spec_names(jne20)   = "ne20"
    short_spec_names(jna23)   = "na23"
    short_spec_names(jmg23)   = "mg23"

    ebind_per_nucleon(jn)   = 0.00000000000000e+00_rt
    ebind_per_nucleon(jp)   = 0.00000000000000e+00_rt
    ebind_per_nucleon(jhe4)   = 7.07391500000000e+00_rt
    ebind_per_nucleon(jc12)   = 7.68014400000000e+00_rt
    ebind_per_nucleon(jo16)   = 7.97620600000000e+00_rt
    ebind_per_nucleon(jne20)   = 8.03224000000000e+00_rt
    ebind_per_nucleon(jna23)   = 8.11149300000000e+00_rt
    ebind_per_nucleon(jmg23)   = 7.90111500000000e+00_rt

    aion(jn)   = 1.00000000000000e+00_rt
    aion(jp)   = 1.00000000000000e+00_rt
    aion(jhe4)   = 4.00000000000000e+00_rt
    aion(jc12)   = 1.20000000000000e+01_rt
    aion(jo16)   = 1.60000000000000e+01_rt
    aion(jne20)   = 2.00000000000000e+01_rt
    aion(jna23)   = 2.30000000000000e+01_rt
    aion(jmg23)   = 2.30000000000000e+01_rt

    aion_inv(jn)   = 1.0_rt/1.00000000000000e+00_rt
    aion_inv(jp)   = 1.0_rt/1.00000000000000e+00_rt
    aion_inv(jhe4)   = 1.0_rt/4.00000000000000e+00_rt
    aion_inv(jc12)   = 1.0_rt/1.20000000000000e+01_rt
    aion_inv(jo16)   = 1.0_rt/1.60000000000000e+01_rt
    aion_inv(jne20)   = 1.0_rt/2.00000000000000e+01_rt
    aion_inv(jna23)   = 1.0_rt/2.30000000000000e+01_rt
    aion_inv(jmg23)   = 1.0_rt/2.30000000000000e+01_rt

    zion(jn)   = 0.00000000000000e+00_rt
    zion(jp)   = 1.00000000000000e+00_rt
    zion(jhe4)   = 2.00000000000000e+00_rt
    zion(jc12)   = 6.00000000000000e+00_rt
    zion(jo16)   = 8.00000000000000e+00_rt
    zion(jne20)   = 1.00000000000000e+01_rt
    zion(jna23)   = 1.10000000000000e+01_rt
    zion(jmg23)   = 1.20000000000000e+01_rt

    nion(jn)   = 1.00000000000000e+00_rt
    nion(jp)   = 0.00000000000000e+00_rt
    nion(jhe4)   = 2.00000000000000e+00_rt
    nion(jc12)   = 6.00000000000000e+00_rt
    nion(jo16)   = 8.00000000000000e+00_rt
    nion(jne20)   = 1.00000000000000e+01_rt
    nion(jna23)   = 1.20000000000000e+01_rt
    nion(jmg23)   = 1.10000000000000e+01_rt

    do i = 1, nspec
       bion(i) = ebind_per_nucleon(i) * aion(i) * ERG_PER_MeV
    end do

    ! Set the mass
    mion(:) = nion(:) * mass_neutron + zion(:) * (mass_proton + mass_electron) &
         - bion(:)/(c_light**2)

    ! Molar mass
    wion(:) = avo * mion(:)

    ! Common approximation
    !wion(:) = aion(:)

    !$acc update device(aion, aion_inv, zion, bion, nion, mion, wion)

#ifdef REACT_SPARSE_JACOBIAN
    ! Set CSR format metadata for Jacobian
    allocate(csr_jac_col_index(NETWORK_SPARSE_JAC_NNZ))
    allocate(csr_jac_row_count(nspec + 3)) ! neq + 1

    csr_jac_col_index = [ &
      1, &
      4, &
      9, &
      1, &
      2, &
      4, &
      9, &
      3, &
      4, &
      9, &
      3, &
      4, &
      9, &
      3, &
      4, &
      5, &
      9, &
      4, &
      6, &
      9, &
      4, &
      7, &
      9, &
      4, &
      8, &
      9, &
      1, &
      2, &
      3, &
      4, &
      5, &
      6, &
      7, &
      8, &
      9, &
      1, &
      2, &
      3, &
      4, &
      5, &
      6, &
      7, &
      8, &
      9, &
      10  ]

    csr_jac_row_count = [ &
      1, &
      4, &
      8, &
      11, &
      14, &
      18, &
      21, &
      24, &
      27, &
      36, &
      46  ]
#endif

  end subroutine actual_network_init


  subroutine actual_network_finalize()
    ! Deallocate storage arrays

    if (allocated(aion)) then
       deallocate(aion)
    endif

    if (allocated(aion_inv)) then
       deallocate(aion_inv)
    endif

    if (allocated(zion)) then
       deallocate(zion)
    endif

    if (allocated(bion)) then
       deallocate(bion)
    endif

    if (allocated(nion)) then
       deallocate(nion)
    endif

    if (allocated(mion)) then
       deallocate(mion)
    endif

    if (allocated(wion)) then
       deallocate(wion)
    endif

#ifdef REACT_SPARSE_JACOBIAN
    if (allocated(csr_jac_col_index)) then
       deallocate(csr_jac_col_index)
    endif

    if (allocated(csr_jac_row_count)) then
       deallocate(csr_jac_row_count)
    endif
#endif

  end subroutine actual_network_finalize


  subroutine ener_gener_rate(dydt, enuc)
    ! Computes the instantaneous energy generation rate

    !$acc routine seq

    implicit none

    real(rt) :: dydt(nspec), enuc

    !$gpu

    ! This is basically e = m c**2

    enuc = sum(dydt(:) * mion(:)) * enuc_conv2

  end subroutine ener_gener_rate

end module actual_network
