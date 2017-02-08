
      ! *************************************************************************************

      subroutine compute_temp(lo,hi, &
                              state   ,s_l1,s_l2,s_l3, s_h1,s_h2,s_h3, &
                              diag_eos,d_l1,d_l2,d_l3, d_h1,d_h2,d_h3, &
                              comoving_a, print_fortran_warnings)

      use eos_module
      use atomic_rates_module, only: this_z, interp_to_this_z
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEINT, UEDEN, &
                                     TEMP_COMP, NE_COMP, small_temp, heat_cool_type
      use  eos_params_module

      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer         , intent(in   ) :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      integer         , intent(in   ) :: print_fortran_warnings
      double precision, intent(inout) ::    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NVAR)
      double precision, intent(inout) :: diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,2)
      double precision, intent(in   ) :: comoving_a

      integer          :: i,j,k
      double precision :: rhoInv,eint
      double precision :: ke(hi(1)-lo(1)+1),dummy_pres(hi(1)-lo(1)+1)
      double precision :: z
      double precision :: eos_inputs_pos_ueint(hi(1)-lo(1)+1,4)
      double precision :: eos_inputs_neg_ueint(hi(1)-lo(1)+1,4)
      integer :: orig_indices(hi(1)-lo(1)+1,3)
      integer :: pos_eos_count, neg_eos_count
      double precision :: small_temp_vec(hi(1)-lo(1)+1)

      z = 1.d0/comoving_a - 1.d0

      if (heat_cool_type.gt.0) then
          if (z .ne. this_z) &
             call interp_to_this_z(z)
      end if

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               if (state(i,j,k,URHO) <= 0.d0) then
                  print *,'   '
                  print *,'>>> Error: compute_temp ',i,j,k
                  print *,'>>> ... negative density ',state(i,j,k,URHO)
                  print *,'    '
                  call bl_error("Error:: compute_temp_3d.f90 &
                                &:: compute_temp")
               end if
            enddo
         enddo
      enddo

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)

            pos_eos_count = 0
            neg_eos_count = 0

            do i = lo(1),hi(1)
               rhoInv = 1.d0 / state(i,j,k,URHO)

               if (state(i,j,k,UEINT) > 0.d0) then

                   pos_eos_count = pos_eos_count + 1

                   eos_inputs_pos_ueint(pos_eos_count,1) = diag_eos(i,j,k,TEMP_COMP)
                   eos_inputs_pos_ueint(pos_eos_count,2) = diag_eos(i,j,k,NE_COMP)
                   eos_inputs_pos_ueint(pos_eos_count,3) = state(i,j,k,URHO)
                   eos_inputs_pos_ueint(pos_eos_count,4) = state(i,j,k,UEINT)*rhoInv

                   orig_indices(pos_eos_count,1) = i
                   orig_indices(pos_eos_count,2) = j
                   orig_indices(pos_eos_count,3) = k

               else

                   neg_eos_count = neg_eos_count + 1

                   eos_inputs_neg_ueint(neg_eos_count,1) = diag_eos(i,j,k,TEMP_COMP) ! DON'T NEED THIS; GET RID OF IT
                   eos_inputs_neg_ueint(neg_eos_count,2) = diag_eos(i,j,k,NE_COMP)
                   eos_inputs_neg_ueint(neg_eos_count,3) = state(i,j,k,URHO)
                   eos_inputs_neg_ueint(neg_eos_count,4) = state(i,j,k,UEINT)

                   orig_indices(neg_eos_count,1) = i
                   orig_indices(neg_eos_count,2) = j
                   orig_indices(neg_eos_count,3) = k

               end if
             end do

             ! For cells with positive E_int
             call nyx_eos_T_given_Re_vec(eos_inputs_pos_ueint(1:pos_eos_count,1), &
                                         eos_inputs_pos_ueint(1:pos_eos_count,2), &
                                         eos_inputs_pos_ueint(1:pos_eos_count,3), &
                                         eos_inputs_pos_ueint(1:pos_eos_count,4), &
                                         comoving_a, &
                                         pos_eos_count)
             diag_eos(orig_indices(1:pos_eos_count,1),j,k,TEMP_COMP) = eos_inputs_pos_ueint(1:pos_eos_count,1)
             diag_eos(orig_indices(1:pos_eos_count,1),j,k,NE_COMP)   = eos_inputs_pos_ueint(1:pos_eos_count,2)

             ! For cells with negative E_int
             call nyx_eos_given_RT_vec(eos_inputs_neg_ueint(1:neg_eos_count,4), &
                                   dummy_pres(1:neg_eos_count), &
                                   eos_inputs_neg_ueint(1:neg_eos_count,3), &
                                   small_temp_vec(1:neg_eos_count), &
                                   eos_inputs_neg_ueint(1:neg_eos_count,2), &
                                   comoving_a, &
                                   neg_eos_count)

             ke(1:neg_eos_count) = 0.5d0 * (state(orig_indices(1:neg_eos_count,1),j,k,UMX)*state(orig_indices(1:neg_eos_count,1),j,k,UMX) + &
                                   state(orig_indices(1:neg_eos_count,1),j,k,UMY)*state(orig_indices(1:neg_eos_count,1),j,k,UMY) + &
                                   state(orig_indices(1:neg_eos_count,1),j,k,UMZ)*state(orig_indices(1:neg_eos_count,1),j,k,UMZ)) * rhoInv

             diag_eos(orig_indices(1:neg_eos_count,1),j,k,TEMP_COMP) = small_temp_vec(1:neg_eos_count)
             state(orig_indices(1:neg_eos_count,1),j,k,UEINT) = eos_inputs_neg_ueint(1:neg_eos_count,3) * eos_inputs_neg_ueint(1:neg_eos_count,4)
             state(orig_indices(1:neg_eos_count,1),j,k,UEDEN) = eos_inputs_neg_ueint(1:neg_eos_count,4) + ke(1:neg_eos_count)

         enddo
      enddo

      end subroutine compute_temp

      subroutine compute_rho_temp(lo,hi,dx, &
                                     state,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                                  diag_eos,d_l1,d_l2,d_l3,d_h1,d_h2,d_h3, &
                                  rho_ave,rho_T_sum, &
                                  T_sum,T_meanrho_sum,rho_sum,vol_sum,vol_mn_sum)

      use meth_params_module, only : NVAR, URHO, TEMP_COMP

      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer         , intent(in   ) :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      double precision, intent(in   ) :: dx(3)
      double precision, intent(in   ) :: rho_ave
      double precision, intent(in   ) ::    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NVAR)
      double precision, intent(in   ) :: diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,2)
      double precision, intent(inout) :: rho_T_sum, rho_sum, T_sum, T_meanrho_sum
      double precision, intent(inout) :: vol_sum, vol_mn_sum

      integer          :: i,j,k
      double precision :: rho_hi, rho_lo, vol

      vol = dx(1)*dx(2)*dx(3)
      rho_hi = 1.1d0*rho_ave
      rho_lo = 0.9d0*rho_ave
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
                   T_sum =     T_sum + vol*diag_eos(i,j,k,TEMP_COMP)
               rho_T_sum = rho_T_sum + state(i,j,k,URHO)*diag_eos(i,j,k,TEMP_COMP)
                 rho_sum =   rho_sum + state(i,j,k,URHO)
                 if ( (state(i,j,k,URHO) .lt. rho_hi) .and. &
                      (state(i,j,k,URHO) .gt. rho_lo) .and. &
                      (diag_eos(i,j,k,TEMP_COMP) .le. 1.0e5) ) then
                         T_meanrho_sum = T_meanrho_sum + vol*dlog10(diag_eos(i,j,k,TEMP_COMP))
                         vol_mn_sum = vol_mn_sum + vol
                 endif
                 vol_sum = vol_sum + vol
            enddo
         enddo
      enddo

      end subroutine compute_rho_temp

      subroutine compute_max_temp_loc(lo,hi, &
                                      state   ,s_l1,s_l2,s_l3, s_h1,s_h2,s_h3, &
                                      diag_eos,d_l1,d_l2,d_l3, d_h1,d_h2,d_h3, &
                                      max_temp, den_maxt, imax, jmax, kmax)

      use meth_params_module, only : TEMP_COMP, NVAR, URHO

      implicit none
      integer         , intent(in   ) :: lo(3),hi(3)
      integer         , intent(in   ) :: s_l1,s_l2,s_l3,s_h1,s_h2,s_h3
      integer         , intent(in   ) :: d_l1,d_l2,d_l3,d_h1,d_h2,d_h3
      double precision, intent(inout) ::    state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NVAR)
      double precision, intent(inout) :: diag_eos(d_l1:d_h1,d_l2:d_h2,d_l3:d_h3,2)
      double precision, intent(in   ) :: max_temp
      double precision, intent(  out) :: den_maxt
      integer         , intent(inout) :: imax,jmax,kmax

      integer                         :: i,j,k
      double precision                :: one_minus_eps

      one_minus_eps = 1.d0 - 1.d-12

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               if (diag_eos(i,j,k,TEMP_COMP) .ge. one_minus_eps*max_temp) then
                  imax = i
                  jmax = j
                  kmax = k
                  den_maxt = state(i,j,k,URHO)
               end if
            enddo
         enddo
      enddo

      end subroutine compute_max_temp_loc
