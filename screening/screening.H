#ifndef _screening_H_
#define _screening_H_

#include <screen.H>

void screening_init();

void screening_finalize();

void add_screening_factor(const int i,
                          const amrex::Real z1, const amrex::Real a1,
                          const amrex::Real z2, const amrex::Real a2);

#endif
