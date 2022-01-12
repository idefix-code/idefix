// ***********************************************************************************
// Idefix MHD astrophysical code
// Copyright(C) 2020-2022 Geoffroy R. J. Lesur <geoffroy.lesur@univ-grenoble-alpes.fr>
// and other code contributors
// Licensed under CeCILL 2.1 License, see COPYING for more information
// ***********************************************************************************

// These shortcuts are largely inspired from K-Athena
// https://gitlab.com/pgrete/kathena
// by P. Grete.

#ifndef LOOP_HPP_
#define LOOP_HPP_

#include <string>
#include "idefix.hpp"
#include "global.hpp"

#define KOKKOS_VECTOR_LENGTH  8


#ifdef INNER_TTR_LOOP
  #define TPINNERLOOP Kokkos::TeamThreadRange
#elif defined INNER_TVR_LOOP
  #define TPINNERLOOP Kokkos::ThreadVectorRange
#else
  #define TPINNERLOOP Kokkos::TeamThreadRange
#endif

typedef Kokkos::TeamPolicy<>               team_policy;
typedef Kokkos::TeamPolicy<>::member_type  member_type;

// 1D default loop pattern
template <typename Function>
inline void idefix_for(const std::string & NAME,
                       const int & IB, const int & IE,
                       Function function) {
  idefix_for(idfx::defaultLoopPattern,NAME,IB,IE,function);
}

// 2D default loop pattern
template <typename Function>
inline void idefix_for(const std::string & NAME,
                       const int & JB, const int & JE,
                       const int & IB, const int & IE,
                       Function function) {
  idefix_for(idfx::defaultLoopPattern,NAME,JB,JE,IB,IE,function);
}

// 3D default loop pattern
template <typename Function>
inline void idefix_for(const std::string & NAME,
                const int & KB, const int & KE,
                const int & JB, const int & JE,
                const int & IB, const int & IE,
                Function function) {
  idefix_for(idfx::defaultLoopPattern,NAME,KB,KE,JB,JE,IB,IE,function);
}

// 4D default loop pattern
template <typename Function>
inline void idefix_for(const std::string & NAME,
                const int & NB, const int & NE,
                const int & KB, const int & KE,
                const int & JB, const int & JE,
                const int & IB, const int & IE,
                Function function) {
  idefix_for(idfx::defaultLoopPattern,NAME,NB,NE,KB,KE,JB,JE,IB,IE,function);
}


// 1D loop
template <typename Function>
inline void idefix_for(LoopPattern lp, const std::string & NAME,
                       const int & IB, const int & IE,
                       Function function) {
  const int NI = IE - IB;
  Kokkos::parallel_for(NAME, NI,
    KOKKOS_LAMBDA (const int& IDX) {
      int i = IDX;
      i += IB;
      function(i);
  });
}


// 2D loop
template <typename Function>
inline void idefix_for(LoopPattern lp, const std::string & NAME,
                       const int & JB, const int & JE,
                       const int & IB, const int & IE,
                       Function function) {
  // Kokkos 1D Range
  if (lp == LoopPattern::RANGE) {
    const int NJ = JE - JB;
    const int NI = IE - IB;
    const int NJNI = NJ * NI;
    Kokkos::parallel_for(NAME, NJNI,
      KOKKOS_LAMBDA (const int& IDX) {
        int j = IDX  / NI;
        int i = IDX - j*NI;
        j += JB;
        i += IB;
        function(j,i);
    });

    // MDRange loops
  } else if (lp == LoopPattern::MDRANGE) {
    Kokkos::parallel_for(NAME,
      Kokkos::MDRangePolicy<Kokkos::Rank<2, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>
        ({JB,IB},{JE,IE}), function);

    // TeamPolicies with single inner loops
  } else if (lp == LoopPattern::TPX || lp == LoopPattern::TPTTRTVR ) {
    const int NJ = JE - JB;
    Kokkos::parallel_for(NAME, team_policy (NJ, Kokkos::AUTO,KOKKOS_VECTOR_LENGTH),
      KOKKOS_LAMBDA (member_type team_member) {
        const int j = team_member.league_rank() + JB;
        Kokkos::parallel_for(TPINNERLOOP<>(team_member,IB,IE),
                             [&] (const int i) {
                               function(j,i);
                            });
    });

    // SIMD FOR loops
  } else if (lp == LoopPattern::SIMDFOR) {
    for (auto j = JB; j < JE; j++)
#pragma omp simd
      for (auto i = IB; i < IE; i++)
        function(j,i);
  } else {
    throw std::runtime_error("Unknown/undefined LoopPattern used.");
  }
}


// 3D loop
template <typename Function>
inline void idefix_for(LoopPattern lp, const std::string & NAME,
                       const int & KB, const int & KE,
                       const int & JB, const int & JE,
                       const int & IB, const int & IE,
                       Function function) {
  // Kokkos 1D Range
  if (lp == LoopPattern::RANGE) {
    const int NK = KE - KB;
    const int NJ = JE - JB;
    const int NI = IE - IB;
    const int NKNJNI = NK*NJ*NI;
    const int NJNI = NJ * NI;
    Kokkos::parallel_for(NAME,NKNJNI,
      KOKKOS_LAMBDA (const int& IDX) {
        int k = IDX / NJNI;
        int j = (IDX - k*NJNI) / NI;
        int i = IDX - k*NJNI - j*NI;
        k += KB;
        j += JB;
        i += IB;
        function(k,j,i);
    });

  // MDRange loops
  } else if (lp == LoopPattern::MDRANGE) {
    Kokkos::parallel_for(NAME,
      Kokkos::MDRangePolicy<Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>
        ({KB,JB,IB},{KE,JE,IE}), function);

  // TeamPolicy with single inner loops
  } else if (lp == LoopPattern::TPX) {
    const int NK = KE - KB;
    const int NJ = JE - JB;
    const int NKNJ = NK * NJ;
    Kokkos::parallel_for(NAME,
      team_policy (NKNJ, Kokkos::AUTO,KOKKOS_VECTOR_LENGTH),
      KOKKOS_LAMBDA (member_type team_member) {
        const int k = team_member.league_rank() / NJ + KB;
        const int j = team_member.league_rank() % NJ + JB;
        Kokkos::parallel_for(TPINNERLOOP<>(team_member,IB,IE),
          [&] (const int i) {
            function(k,j,i);
        });
      });

  // TeamPolicy with nested TeamThreadRange and ThreadVectorRange
  } else if (lp == LoopPattern::TPTTRTVR) {
    const int NK = KE - KB;
    Kokkos::parallel_for(NAME,
      team_policy (NK, Kokkos::AUTO,KOKKOS_VECTOR_LENGTH),
      KOKKOS_LAMBDA (member_type team_member) {
        const int k = team_member.league_rank() + KB;
        Kokkos::parallel_for(
          Kokkos::TeamThreadRange<>(team_member,JB,JE),
          [&] (const int j) {
            Kokkos::parallel_for(
              Kokkos::ThreadVectorRange<>(team_member,IB,IE),
              [&] (const int i) {
                function(k,j,i);
              });
          });
      });

  // SIMD FOR loops
  } else if (lp == LoopPattern::SIMDFOR) {
    for (auto k = KB; k < KE; k++)
      for (auto j = JB; j < JE; j++)
#pragma omp simd
        for (auto i = IB; i < IE; i++)
          function(k,j,i);
  } else {
    throw std::runtime_error("Unknown/undefined LoopPattern used.");
  }
}

// 4D loop
template <typename Function>
inline void idefix_for(LoopPattern lp, const std::string & NAME,
                       const int NB, const int NE,
                       const int KB, const int KE,
                       const int JB, const int JE,
                       const int IB, const int IE,
                       Function function) {
  // Kokkos 1D Range
  if (lp == LoopPattern::RANGE) {
    const int NN = (NE) - (NB);
    const int NK = (KE) - (KB);
    const int NJ = (JE) - (JB);
    const int NI = (IE) - (IB);
    const int NNNKNJNI = NN*NK*NJ*NI;
    const int NKNJNI = NK*NJ*NI;
    const int NJNI = NJ * NI;
    Kokkos::parallel_for(NAME,NNNKNJNI,
      KOKKOS_LAMBDA (const int& IDX) {
        int n = IDX / NKNJNI;
        int k = (IDX - n*NKNJNI) / NJNI;
        int j = (IDX - n*NKNJNI - k*NJNI) / NI;
        int i = IDX - n*NKNJNI - k*NJNI - j*NI;
        n += (NB);
        k += (KB);
        j += (JB);
        i += (IB);
        function(n,k,j,i);
    });

  // MDRange loops
  } else if (lp == LoopPattern::MDRANGE) {
    Kokkos::parallel_for(NAME,
      Kokkos::MDRangePolicy<Kokkos::Rank<4,Kokkos::Iterate::Right, Kokkos::Iterate::Right>>
        ({NB,KB,JB,IB},{NE,KE,JE,IE}), function);

  // TeamPolicy loops
  } else if (lp == LoopPattern::TPX) {
    const int NN = NE - NB;
    const int NK = KE - KB;
    const int NJ = JE - JB;
    const int NKNJ = NK * NJ;
    const int NNNKNJ = NN * NK * NJ;
    Kokkos::parallel_for(NAME,
      team_policy (NNNKNJ, Kokkos::AUTO,KOKKOS_VECTOR_LENGTH),
      KOKKOS_LAMBDA (member_type team_member) {
        int n = team_member.league_rank() / NKNJ;
        int k = (team_member.league_rank() - n*NKNJ) / NJ;
        int j = team_member.league_rank() - n*NKNJ - k*NJ + JB;
        n += NB;
        k += KB;
        Kokkos::parallel_for(
          TPINNERLOOP<>(team_member,IB,IE),
          [&] (const int i) {
            function(n,k,j,i);
          });
      });

  // TeamPolicy with nested TeamThreadRange and ThreadVectorRange
  } else if (lp == LoopPattern::TPTTRTVR) {
    const int NN = NE - NB;
    const int NK = KE - KB;
    const int NNNK = NN * NK;
    Kokkos::parallel_for(NAME,
      team_policy (NNNK, Kokkos::AUTO,KOKKOS_VECTOR_LENGTH),
      KOKKOS_LAMBDA (member_type team_member) {
        int n = team_member.league_rank() / NK + NB;
        int k = team_member.league_rank() % NK + KB;
        Kokkos::parallel_for(
          Kokkos::TeamThreadRange<>(team_member,JB,JE),
          [&] (const int j) {
            Kokkos::parallel_for(
              Kokkos::ThreadVectorRange<>(team_member,IB,IE),
              [&] (const int i) {
                function(n,k,j,i);
              });
          });
      });

  // SIMD FOR loops
  } else if (lp == LoopPattern::SIMDFOR) {
    for (auto n = NB; n < NE; n++)
      for (auto k = KB; k < KE; k++)
        for (auto j = JB; j < JE; j++)
#pragma omp simd
          for (auto i = IB; i < IE; i++)
            function(n,k,j,i);
  } else {
    throw std::runtime_error("Unknown/undefined LoopPattern used.");
  }
}

#endif // LOOP_HPP_
