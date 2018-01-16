! Calculates temperature and free electron density using Newton-Raphson solver
!
!     Equilibrium ionization fractions, optically thin media, based on:
!     Katz, Weinberg & Hernquist, 1996: Astrophysical Journal Supplement v.105, p.19
!
! Units are CGS, **BUT** 6 fractions: ne, nh0, nhp, nhe0, nhep, nhepp
!       are in units of nh (hydrogen number density)
!

module eos_module

  use amrex_fort_module, only : rt => amrex_real
  use iso_c_binding, only: c_double

  implicit none

  ! Routines:
  public  :: nyx_eos_given_RT, nyx_eos_given_RT_vec, nyx_eos_T_given_Re, nyx_eos_T_given_Re_vec, eos_init_small_pres
  public  :: nyx_eos_nh0_and_nhep, iterate_ne, iterate_ne_vec
  private :: ion_n

  real(rt), public :: xacc ! EOS Newton-Raphson convergence tolerance
  real(c_double), public :: vode_rtol, vode_atol_scaled ! VODE integration tolerances

  contains

      subroutine fort_setup_eos_params (xacc_in, vode_rtol_in, vode_atol_scaled_in) &
                                       bind(C, name='fort_setup_eos_params')
        use amrex_fort_module, only : rt => amrex_real
        implicit none
        real(rt), intent(in) :: xacc_in, vode_rtol_in, vode_atol_scaled_in

        xacc = xacc_in
        vode_rtol = vode_rtol_in
        vode_atol_scaled = vode_atol_scaled_in

      end subroutine fort_setup_eos_params

     ! ****************************************************************************

      subroutine eos_init_small_pres(R, T, Ne, P, a)

        use amrex_fort_module, only : rt => amrex_real
        use atomic_rates_module, ONLY: YHELIUM
        use fundamental_constants_module, only: mp_over_kb

        implicit none

        real(rt), intent(  out) :: P
        real(rt), intent(in   ) :: R, T, Ne
        real(rt), intent(in   ) :: a

        real(rt) :: mu

        mu = (1.0d0+4.0d0*YHELIUM) / (1.0d0+YHELIUM+Ne)
        P  = R*T / (mp_over_kB * mu)

      end subroutine eos_init_small_pres

     ! ****************************************************************************

      subroutine nyx_eos_soundspeed(c, R, e)

        use meth_params_module, only: gamma_const, gamma_minus_1

        implicit none

        real(rt), intent(in   ) :: R, e
        real(rt), intent(  out) :: c

        ! sound speed: c^2 = gamma*P/rho
        c = sqrt(gamma_const * gamma_minus_1 *e)

      end subroutine nyx_eos_soundspeed

     ! ****************************************************************************

      subroutine nyx_eos_S_given_Re(S, R, T, Ne, a)

        use bl_constants_module, only: M_PI
        use atomic_rates_module, ONLY: YHELIUM
        use fundamental_constants_module, only: mp_over_kb
        use fundamental_constants_module, only: k_B, hbar, m_proton
        implicit none

        real(rt),          intent(  out) :: S
        real(rt),          intent(in   ) :: R, T, Ne, a

        real(rt) :: mu, dens, t1, t2, t3

        mu = (1.0d0+4.0d0*YHELIUM) / (1.0d0+YHELIUM+Ne)
        dens = R/(a*a*a)

        ! Entropy (per gram) of an ideal monoatomic gas (Sactur-Tetrode equation)
        ! NOTE: this expression is only valid for gamma = 5/3.
        t1 = (mu*m_proton);            t1 = t1*t1*sqrt(t1)
        t2 = (k_B*T);                  t2 = t2*sqrt(t2)
        t3 = (2.0d0*M_PI*hbar*hbar);   t3 = t3*sqrt(t3)

        S = (1.d0 / (mu*mp_over_kB)) * (2.5d0 + log(t1/dens*t2/t3))

      end subroutine nyx_eos_S_given_Re

     ! ****************************************************************************

      subroutine nyx_eos_given_RT(e, P, R, T, Ne, a)

        use atomic_rates_module, ONLY: YHELIUM
        use fundamental_constants_module, only: mp_over_kb
        use meth_params_module, only: gamma_minus_1
        implicit none

        double precision,          intent(  out) :: e, P
        double precision,          intent(in   ) :: R, T, Ne
        double precision,          intent(in   ) :: a

        double precision :: mu

        mu = (1.0d0+4.0d0*YHELIUM) / (1.0d0+YHELIUM+Ne)
        e  = T / (gamma_minus_1 * mp_over_kB * mu)

        P  = gamma_minus_1 * R * e

      end subroutine nyx_eos_given_RT

     ! ****************************************************************************

      subroutine nyx_eos_given_RT_vec(e, P, R, T, Ne, a, veclen)

        use atomic_rates_module, ONLY: YHELIUM
        use fundamental_constants_module, only: mp_over_kb
        use meth_params_module, only: gamma_minus_1
        implicit none

        integer, intent(in) :: veclen
        real(rt), dimension(veclen), intent(  out) :: e, P
        real(rt), dimension(veclen), intent(in   ) :: R, T, Ne
        real(rt),          intent(in   ) :: a

        real(rt), dimension(veclen) :: mu
        integer :: i

        do i = 1, veclen
          mu(i) = (1.0d0+4.0d0*YHELIUM) / (1.0d0+YHELIUM+Ne(i))
          e(i)  = T(i) / (gamma_minus_1 * mp_over_kB * mu(i))
  
          P(i)  = gamma_minus_1 * R(i) * e(i)
        end do

      end subroutine nyx_eos_given_RT_vec

     ! ****************************************************************************

      subroutine nyx_eos_T_given_Re(JH, JHe, T, Ne, R_in, e_in, a, species)

      use atomic_rates_module, ONLY: XHYDROGEN, MPROTON
      use fundamental_constants_module, only: density_to_cgs, e_to_cgs

      ! In/out variables
      integer,    intent(in)    :: JH, JHe
      real(rt),   intent(inout) :: T, Ne
      real(rt),   intent(in   ) :: R_in, e_in
      real(rt),   intent(in   ) :: a
      real(rt), optional, intent(out) :: species(5)

      double precision :: nh, nh0, nhep, nhp, nhe0, nhepp
      double precision :: z, rho, U
      integer          :: NR

      ! This converts from code units to CGS
      rho = R_in * density_to_cgs / a**3
        U = e_in * e_to_cgs
      nh  = rho*XHYDROGEN/MPROTON

      z   = 1.d0/a - 1.d0

      call iterate_ne(JH, Jhe, z, U, T, NR, nh, ne, nh0, nhp, nhe0, nhep, nhepp)

      if (present(species)) then
         species(1) = nh0
         species(2) = nhp
         species(3) = nhe0
         species(4) = nhep
         species(5) = nhepp
      endif

      end subroutine nyx_eos_T_given_Re

     ! ****************************************************************************

      subroutine nyx_eos_T_given_Re_vec(T, Ne, R_in, e_in, a, veclen)

      use amrex_fort_module, only : rt => amrex_real
      use atomic_rates_module, ONLY: XHYDROGEN, MPROTON
      use fundamental_constants_module, only: density_to_cgs, e_to_cgs

      ! In/out variables
      integer, intent(in) :: veclen
      real(rt), dimension(veclen), intent(inout) :: T, Ne
      real(rt), dimension(veclen), intent(in   ) :: R_in, e_in
      real(rt),                    intent(in   ) :: a

      real(rt), dimension(veclen) :: nh, nh0, nhep, nhp, nhe0, nhepp, rho, U
      real(rt) :: z

      ! This converts from code units to CGS
      rho = R_in * density_to_cgs / a**3
        U = e_in * e_to_cgs
      nh  = rho*XHYDROGEN/MPROTON

      z   = 1.d0/a - 1.d0

      call iterate_ne_vec(z, U, T, nh, ne, nh0, nhp, nhe0, nhep, nhepp, veclen)

      end subroutine nyx_eos_T_given_Re_vec

     ! ****************************************************************************

      subroutine nyx_eos_nh0_and_nhep(JH, JHe, z, rho, e, nh0, nhep)
      ! This is for skewers analysis code, input is in CGS

      use atomic_rates_module, only: XHYDROGEN, MPROTON

      ! In/out variables
      integer, intent(in) :: JH, Jhe
      real(rt),           intent(in   ) :: z, rho, e
      real(rt),           intent(  out) :: nh0, nhep

      real(rt) :: nh, nhp, nhe0, nhepp, T, ne
      integer  :: NR

      nh  = rho*XHYDROGEN/MPROTON
      ne  = 1.0d0 ! Guess

      call iterate_ne(JH, JHe, z, e, T, NR, nh, ne, nh0, nhp, nhe0, nhep, nhepp)

      nh0  = nh*nh0
      nhep = nh*nhep

      end subroutine nyx_eos_nh0_and_nhep

     ! ****************************************************************************

      subroutine iterate_ne_vec(z, U, t, nh, ne, nh0, nhp, nhe0, nhep, nhepp, veclen)

      use atomic_rates_module, ONLY: this_z, YHELIUM, BOLTZMANN, MPROTON, TCOOLMAX_R
      use meth_params_module, only: gamma_minus_1
      use amrex_error_module, only: amrex_abort

      integer :: i

      integer, intent(in) :: veclen
      real(rt), intent (in   ) :: z
      real(rt), dimension(veclen), intent(in) :: U, nh
      real(rt), dimension(veclen), intent (inout) :: ne
      real(rt), dimension(veclen), intent (  out) :: t, nh0, nhp, nhe0, nhep, nhepp

      real(rt), parameter :: xacc = 1.0d-6

      integer, dimension(veclen)  :: JH, JHe
      real(rt), dimension(veclen) :: f, df, eps, mu
      real(rt), dimension(veclen) :: nhp_plus, nhep_plus, nhepp_plus
      real(rt), dimension(veclen) :: dnhp_dne, dnhep_dne, dnhepp_dne, dne
      real(rt), dimension(veclen):: U_in, t_in, nh_in, ne_in
      real(rt), dimension(veclen) :: nhp_out, nhep_out, nhepp_out
      integer :: vec_count, orig_idx(veclen)
      integer :: ii
      character(len=128) :: errmsg

      ! Check if we have interpolated to this z
      if (abs(z-this_z) .gt. xacc*z) then
          write(errmsg, *) "iterate_ne_vec(): Wrong redshift! z = ", z, " but this_z = ", this_z
          call amrex_abort(errmsg)
      end if

      ii = 0
      ne(1:veclen) = 1.0d0 ! 0 is a bad guess

      do  ! Newton-Raphson solver
         ii = ii + 1

         ! Ion number densities
         do i = 1, veclen
           mu(i) = (1.0d0+4.0d0*YHELIUM) / (1.0d0+YHELIUM+ne(i))
           t(i)  = gamma_minus_1*MPROTON/BOLTZMANN * U(i) * mu(i)
         end do
         vec_count = 0
         do i = 1, veclen
           if (t(i) .ge. TCOOLMAX_R) then ! Fully ionized plasma
             nhp(i)   = 1.0d0
             nhep(i)  = 0.0d0
             nhepp(i) = YHELIUM
           else
             vec_count = vec_count + 1
             U_in(vec_count) = U(i)
             t_in(vec_count) = t(i)
             nh_in(vec_count) = nh(i)
             ne_in(vec_count) = ne(i)
             orig_idx(vec_count) = i
           endif
         end do

         call ion_n_vec(JH(1:vec_count), &
                    JHe(1:vec_count), &
                    U_in(1:vec_count), &
                    nh_in(1:vec_count), &
                    ne_in(1:vec_count), &
                    nhp_out(1:vec_count), &
                    nhep_out(1:vec_count), &
                    nhepp_out(1:vec_count), &
                    t_in(1:vec_count), &
                    vec_count)
         nhp(orig_idx(1:vec_count)) = nhp_out(1:vec_count)
         nhep(orig_idx(1:vec_count)) = nhep_out(1:vec_count)
         nhepp(orig_idx(1:vec_count)) = nhepp_out(1:vec_count)

         ! Forward difference derivatives
         do i = 1, veclen
           if (ne(i) .gt. 0.0d0) then
              eps(i) = xacc*ne(i)
           else
              eps(i) = 1.0d-24
           endif
         end do
         do i = 1, veclen
           mu(i) = (1.0d0+4.0d0*YHELIUM) / (1.0d0+YHELIUM+ne(i)+eps(i))
           t(i)  = gamma_minus_1*MPROTON/BOLTZMANN * U(i) * mu(i)
         end do
         vec_count = 0
         do i = 1, veclen
           if (t(i) .ge. TCOOLMAX_R) then ! Fully ionized plasma
             nhp_plus(i)   = 1.0d0
             nhep_plus(i)  = 0.0d0
             nhepp_plus(i) = YHELIUM
           else
             vec_count = vec_count + 1
             U_in(vec_count) = U(i)
             t_in(vec_count) = t(i)
             nh_in(vec_count) = nh(i)
             ne_in(vec_count) = ne(i)+eps(i)
             orig_idx(vec_count) = i
           endif
         end do

         call ion_n_vec(JH(1:vec_count), &
                    JHe(1:vec_count), &
                    U_in(1:vec_count), &
                    nh_in(1:vec_count), &
                    ne_in(1:vec_count), &
                    nhp_out(1:vec_count), &
                    nhep_out(1:vec_count), &
                    nhepp_out(1:vec_count), &
                    t_in(1:vec_count), &
                    vec_count)
         nhp_plus(orig_idx(1:vec_count)) = nhp_out(1:vec_count)
         nhep_plus(orig_idx(1:vec_count)) = nhep_out(1:vec_count)
         nhepp_plus(orig_idx(1:vec_count)) = nhepp_out(1:vec_count)

         do i = 1, veclen
           dnhp_dne(i)   = (nhp_plus(i)   - nhp(i))   / eps(i)
           dnhep_dne(i)  = (nhep_plus(i)  - nhep(i))  / eps(i)
           dnhepp_dne(i) = (nhepp_plus(i) - nhepp(i)) / eps(i)
         end do

         do i = 1, veclen
           f(i)   = ne(i) - nhp(i) - nhep(i) - 2.0d0*nhepp(i)
           df(i)  = 1.0d0 - dnhp_dne(i) - dnhep_dne(i) - 2.0d0*dnhepp_dne(i)
           dne(i) = f(i)/df(i)
         end do

         do i = 1, veclen
           ne(i) = max((ne(i)-dne(i)), 0.0d0)
         end do

         if (maxval(abs(dne(1:veclen))) < xacc) exit

         if (ii .gt. 15) &
            STOP 'iterate_ne_vec(): No convergence in Newton-Raphson!'

      enddo

      ! Get rates for the final ne
      do i = 1, veclen
        mu(i) = (1.0d0+4.0d0*YHELIUM) / (1.0d0+YHELIUM+ne(i))
        t(i)  = gamma_minus_1*MPROTON/BOLTZMANN * U(i) * mu(i)
      end do
      vec_count = 0
      do i = 1, veclen
        if (t(i) .ge. TCOOLMAX_R) then ! Fully ionized plasma
          nhp(i)   = 1.0d0
          nhep(i)  = 0.0d0
          nhepp(i) = YHELIUM
        else
          vec_count = vec_count + 1
          U_in(vec_count) = U(i)
          t_in(vec_count) = t(i)
          nh_in(vec_count) = nh(i)
          ne_in(vec_count) = ne(i)
          orig_idx(vec_count) = i
        endif
      end do
      call ion_n_vec(JH(1:vec_count), &
                 JHe(1:vec_count), &
                 U_in(1:vec_count), &
                 nh_in(1:vec_count), &
                 ne_in(1:vec_count), &
                 nhp_out(1:vec_count), &
                 nhep_out(1:vec_count), &
                 nhepp_out(1:vec_count), &
                 t_in(1:vec_count), &
                 vec_count)
      nhp(orig_idx(1:vec_count)) = nhp_out(1:vec_count)
      nhep(orig_idx(1:vec_count)) = nhep_out(1:vec_count)
      nhepp(orig_idx(1:vec_count)) = nhepp_out(1:vec_count)

      ! Neutral fractions:
      do i = 1, veclen
        nh0(i)   = 1.0d0 - nhp(i)
        nhe0(i)  = YHELIUM - (nhep(i) + nhepp(i))
      end do
      end subroutine iterate_ne_vec

     ! ****************************************************************************

      subroutine ion_n_vec(JH, JHe, U, nh, ne, nhp, nhep, nhepp, t, vec_count)

      use amrex_fort_module, only : rt => amrex_real
      use meth_params_module, only: gamma_minus_1
      use atomic_rates_module, ONLY: YHELIUM, MPROTON, BOLTZMANN, &
                                     TCOOLMIN, TCOOLMAX, NCOOLTAB, deltaT, &
                                     AlphaHp, AlphaHep, AlphaHepp, Alphad, &
                                     GammaeH0, GammaeHe0, GammaeHep, &
                                     ggh0, gghe0, gghep

      integer, intent(in) :: vec_count
      integer, dimension(vec_count), intent(in) :: JH, JHe
      real(rt), intent(in   ) :: U(vec_count), nh(vec_count), ne(vec_count)
      real(rt), intent(  out) :: nhp(vec_count), nhep(vec_count), nhepp(vec_count), t(vec_count)
      real(rt) :: ahp(vec_count), ahep(vec_count), ahepp(vec_count), ad(vec_count), geh0(vec_count), gehe0(vec_count), gehep(vec_count)
      real(rt) :: ggh0ne(vec_count), gghe0ne(vec_count), gghepne(vec_count)
      real(rt) :: mu(vec_count), tmp(vec_count), logT(vec_count), flo(vec_count), fhi(vec_count)
      real(rt), parameter :: smallest_val=tiny(1.0d0)
      integer :: j(vec_count), i

      mu(:) = (1.0d0+4.0d0*YHELIUM) / (1.0d0+YHELIUM+ne(:))
      t(:)  = gamma_minus_1*MPROTON/BOLTZMANN * U(:) * mu(:)

      logT(1:vec_count) = dlog10(t(1:vec_count))

      ! Temperature floor
      do i = 1, vec_count
        if (logT(i) .le. TCOOLMIN) logT(i) = TCOOLMIN + 0.5d0*deltaT
      end do

      ! Interpolate rates
      do i = 1, vec_count
        tmp(i) = (logT(i)-TCOOLMIN)/deltaT
        j(i) = int(tmp(i))
        fhi(i) = tmp(i) - j(i)
        flo(i) = 1.0d0 - fhi(i)
        j(i) = j(i) + 1 ! F90 arrays start with 1
      end do

      do i = 1, vec_count
        ahp(i)   = flo(i)*AlphaHp  (j(i)) + fhi(i)*AlphaHp  (j(i)+1)
        ahep(i)  = flo(i)*AlphaHep (j(i)) + fhi(i)*AlphaHep (j(i)+1)
        ahepp(i) = flo(i)*AlphaHepp(j(i)) + fhi(i)*AlphaHepp(j(i)+1)
        ad(i)    = flo(i)*Alphad   (j(i)) + fhi(i)*Alphad   (j(i)+1)
        geh0(i)  = flo(i)*GammaeH0 (j(i)) + fhi(i)*GammaeH0 (j(i)+1)
        gehe0(i) = flo(i)*GammaeHe0(j(i)) + fhi(i)*GammaeHe0(j(i)+1)
        gehep(i) = flo(i)*GammaeHep(j(i)) + fhi(i)*GammaeHep(j(i)+1)
      end do

      do i = 1, vec_count
        if (ne(i) .gt. 0.0d0) then
           ggh0ne(i)   = JH(i)  * ggh0  / (ne(i)*nh(i))
           gghe0ne(i)  = JH(i)  * gghe0 / (ne(i)*nh(i))
           gghepne(i)  = JHe(i) * gghep / (ne(i)*nh(i))
        else
           ggh0ne(i)   = 0.0d0
           gghe0ne(i)  = 0.0d0
           gghepne(i)  = 0.0d0
        endif
      end do

      ! H+
      do i = 1, vec_count
        nhp(i) = 1.0d0 - ahp(i)/(ahp(i) + geh0(i) + ggh0ne(i))
      end do

      ! He+
      do i = 1, vec_count
        if ((gehe0(i) + gghe0ne(i)) .gt. smallest_val) then
  
           nhep(i)  = YHELIUM/(1.0d0 + (ahep(i)  + ad(i)     )/(gehe0(i) + gghe0ne(i)) &
                                  + (gehep(i) + gghepne(i))/ahepp(i))
        else
           nhep(i)  = 0.0d0
        endif
      end do

      ! He++
      do i = 1, vec_count
        if (nhep(i) .gt. 0.0d0) then
           nhepp(i) = nhep(i)*(gehep(i) + gghepne(i))/ahepp(i)
        else
           nhepp(i) = 0.0d0
        endif
      end do

      end subroutine ion_n_vec

     ! ****************************************************************************

      subroutine iterate_ne(JH, JHe, z, U, t, NR, nh, ne, nh0, nhp, nhe0, nhep, nhepp)

      use amrex_error_module, only: amrex_abort
      use atomic_rates_module, only: this_z, YHELIUM

      integer :: i

      integer, intent(in) :: JH, JHe
      integer, intent(out) :: NR
      real(rt), intent (in   ) :: z, U, nh
      real(rt), intent (inout) :: ne
      real(rt), intent (  out) :: t, nh0, nhp, nhe0, nhep, nhepp

      real(rt) :: f, df, eps
      real(rt) :: nhp_plus, nhep_plus, nhepp_plus
      real(rt) :: dnhp_dne, dnhep_dne, dnhepp_dne, dne
      character(len=128) :: errmsg

      ! Check if we have interpolated to this z
      if (abs(z-this_z) .gt. xacc*z) then
          write(errmsg, *) "iterate_ne(): Wrong redshift! z = ", z, " but this_z = ", this_z
          call amrex_abort(errmsg)
      end if

      i = 0
      NR = 0
      ne = 1.0d0 ! 0 is a bad guess
      do  ! Newton-Raphson solver
         i = i + 1

         ! Ion number densities
         call ion_n(JH, JHe, U, nh, ne, nhp, nhep, nhepp, t)

         ! Forward difference derivatives
         if (ne .gt. 0.0d0) then
            eps = xacc*ne
         else
            eps = 1.0d-24
         endif
         call ion_n(JH, JHe, U, nh, (ne+eps), nhp_plus, nhep_plus, nhepp_plus, t)

         NR  = NR + 2

         dnhp_dne   = (nhp_plus   - nhp)   / eps
         dnhep_dne  = (nhep_plus  - nhep)  / eps
         dnhepp_dne = (nhepp_plus - nhepp) / eps

         f   = ne - nhp - nhep - 2.0d0*nhepp
         df  = 1.0d0 - dnhp_dne - dnhep_dne - 2.0d0*dnhepp_dne
         dne = f/df

         ne = max((ne-dne), 0.0d0)

         if (abs(dne) < xacc) exit

         if (i .gt. 10) then
            !$OMP CRITICAL
            print*, "ITERATION: ", i, " NUMBERS: ", z, t, ne, nhp, nhep, nhepp, df
            if (i .gt. 12) &
               STOP 'iterate_ne(): No convergence in Newton-Raphson!'
            !$OMP END CRITICAL
         endif

      enddo

      ! Get rates for the final ne
      call ion_n(JH, JHe, U, nh, ne, nhp, nhep, nhepp, t)
      NR  = NR + 1

      ! Neutral fractions:
      nh0   = 1.0d0 - nhp
      nhe0  = YHELIUM - (nhep + nhepp)
      end subroutine iterate_ne

     ! ****************************************************************************

      subroutine ion_n(JH, JHe, U, nh, ne, nhp, nhep, nhepp, t)

      use meth_params_module,  only: gamma_minus_1
      use atomic_rates_module, only: YHELIUM, MPROTON, BOLTZMANN, &
                                     TCOOLMIN, TCOOLMAX, NCOOLTAB, deltaT, &
                                     AlphaHp, AlphaHep, AlphaHepp, Alphad, &
                                     GammaeH0, GammaeHe0, GammaeHep, &
                                     ggh0, gghe0, gghep

      integer, intent(in) :: JH, JHe
      real(rt), intent(in   ) :: U, nh, ne
      real(rt), intent(  out) :: nhp, nhep, nhepp, t
      real(rt) :: ahp, ahep, ahepp, ad, geh0, gehe0, gehep
      real(rt) :: ggh0ne, gghe0ne, gghepne
      real(rt) :: mu, tmp, logT, flo, fhi
      real(rt), parameter :: smallest_val=tiny(1.0d0)
      integer :: j


      mu = (1.0d0+4.0d0*YHELIUM) / (1.0d0+YHELIUM+ne)
      t  = gamma_minus_1*MPROTON/BOLTZMANN * U * mu

      logT = dlog10(t)
      if (logT .ge. TCOOLMAX) then ! Fully ionized plasma
         nhp   = 1.0d0
         nhep  = 0.0d0
         nhepp = YHELIUM
         return
      endif

      ! Temperature floor
      if (logT .le. TCOOLMIN) logT = TCOOLMIN + 0.5d0*deltaT

      ! Interpolate rates
      tmp = (logT-TCOOLMIN)/deltaT
      j = int(tmp)
      fhi = tmp - j
      flo = 1.0d0 - fhi
      j = j + 1 ! F90 arrays start with 1

      ahp   = flo*AlphaHp  (j) + fhi*AlphaHp  (j+1)
      ahep  = flo*AlphaHep (j) + fhi*AlphaHep (j+1)
      ahepp = flo*AlphaHepp(j) + fhi*AlphaHepp(j+1)
      ad    = flo*Alphad   (j) + fhi*Alphad   (j+1)
      geh0  = flo*GammaeH0 (j) + fhi*GammaeH0 (j+1)
      gehe0 = flo*GammaeHe0(j) + fhi*GammaeHe0(j+1)
      gehep = flo*GammaeHep(j) + fhi*GammaeHep(j+1)

      if (ne .gt. 0.0d0) then
         ggh0ne   = JH  * ggh0  / (ne*nh)
         gghe0ne  = JH  * gghe0 / (ne*nh)
         gghepne  = JHe * gghep / (ne*nh)
      else
         ggh0ne   = 0.0d0
         gghe0ne  = 0.0d0
         gghepne  = 0.0d0
      endif

      ! H+
      nhp = 1.0d0 - ahp/(ahp + geh0 + ggh0ne)

      ! He+
      if ((gehe0 + gghe0ne) .gt. smallest_val) then

         nhep  = YHELIUM/(1.0d0 + (ahep  + ad     )/(gehe0 + gghe0ne) &
                                + (gehep + gghepne)/ahepp)
      else
         nhep  = 0.0d0
      endif

      ! He++
      if (nhep .gt. 0.0d0) then
         nhepp = nhep*(gehep + gghepne)/ahepp
      else
         nhepp = 0.0d0
      endif

      end subroutine ion_n


end module eos_module
