! Title: Solving the income fluctuation problem with vectorization/parallelization
! Author: Alessandro Di Nola
! References:
!  (1) https://julia.quantecon.org/dynamic_programming/ifp.html
!  (2) Dynamic Programming on the GPU via JAX
!   QE note https://notes.quantecon.org/submission/622ed4daf57192000f918c61
!===============================================================================!
module mod_numerical
   implicit none

   private
   public :: linspace
   integer, parameter, private :: wp=kind(0.0d0)

contains

   function linspace(my_start, my_stop, n)
      ! Purpose: replicate Matlab function <linspace>
      implicit none
      ! Inputs:
      integer :: n
      real(kind=wp) :: my_start, my_stop
      ! Function result:
      real(kind=wp) :: linspace(n)
      ! Locals:
      integer :: i
      real(kind=wp) :: step, grid(n)

      if (n == 1) then
         grid(n) = (my_start + my_stop)/2.d0
      elseif (n >= 2) then
         grid(1) = my_start
         if (n > 2) then
            step = (my_stop - my_start)/(real(n - 1, 8))
            do i = 2, n - 1
               grid(i) = grid(i - 1) + step
            end do
         end if
         grid(n) = my_stop
      end if
      linspace = grid

   end function linspace

end module mod_numerical
!===============================================================================!

module mod_globals
   implicit none

   integer, parameter, private :: wp=kind(0.0d0)
   real(kind=wp), parameter :: large_negative = -1.0d10
   character(len=*), parameter :: savedir = "output/"
   real(kind=wp) :: R, beta, gamma, a_min, a_max, rho, sig_z
   integer :: n_a, n_z, n_ap, par_fortran

   real(kind=wp), allocatable :: a_grid(:), ap_grid(:), z_grid(:), z_tran(:, :)
   real(kind=wp), allocatable :: a_grid2(:, :), ap_grid2(:, :) !needed for vectorized code

contains

   elemental function f_util(c) result(util)
      ! Purpose: utility function
      ! Assumption: c is positive, this is enforced
      ! elsewhere in the code
      real(kind=wp), intent(in) :: c
      real(kind=wp) :: util

      if (abs(gamma - 1.0d0) < 1.0d-6) then
         util = log(c)
      else
         util = c**(1.0d0 - gamma)/(1.0d0 - gamma)
      end if

   end function f_util

end module mod_globals
!===============================================================================!

module mod_vfi
   use mod_globals
   use omp_lib
   implicit none

   private
   public :: bellman_op, bellman_op_vec, bellman_op_vec2
   integer, parameter, private :: wp=kind(0.0d0)

contains

   subroutine bellman_op(v_new, pol_ap, v)
      ! Purpose: Bellman operator, it does one step of the VFI
      ! This version is loop-based and can be parallelized with OpenMP
      ! For vectorized versions, see subroutine bellman_op_vec, bellman_op_vec2
      implicit none
      ! Inputs:
      real(kind=wp), intent(in) :: v(:, :)
      ! Outputs:
      real(kind=wp), intent(out) :: v_new(:, :), pol_ap(:, :)
      ! Locals:
      real(kind=wp), allocatable :: EV(:, :)
      integer :: a_c, z_c, ap_c, ap_ind
      real(kind=wp) :: a_val, z_val, aprime_val, cons, v_max, v_temp

      ! Compute expected value function EV(a',z) = sum over z' of v(a',z')*z_tran(z,z')
      EV = matmul(v, transpose(z_tran))

      ! Step through the state space
      !$omp parallel if (par_fortran==1) default(shared) private(z_c,a_c,a_val,z_val,v_max,ap_ind,ap_c,aprime_val,cons,v_temp)
      !$omp do collapse(2)
      do z_c = 1, n_z
         do a_c = 1, n_a
            ! Current states
            a_val = a_grid(a_c)
            z_val = z_grid(z_c)
            ! Initialize v_new
            v_max = large_negative
            ap_ind = 0
            ! Choose a' optimally by stepping through all possible values
            do ap_c = 1, n_a
               aprime_val = ap_grid(ap_c)
               cons = R*a_val + z_val - aprime_val
               if (cons > 0.0d0) then
                  v_temp = f_util(cons) + beta*EV(ap_c, z_c)
                  !v_temp = f_util(cons) + beta*sum(v(ap_c,:)*z_tran(z_c,:))
                  if (v_temp > v_max) then
                     v_max = v_temp
                     ap_ind = ap_c
                  end if
               end if
            end do !end a'
            v_new(a_c, z_c) = v_max
            pol_ap(a_c, z_c) = ap_grid(ap_ind)
         end do
      end do
      !$omp enddo
      !$omp end parallel

   end subroutine bellman_op

   subroutine bellman_op_vec(v_new, pol_ap, v)
      ! Purpose: Bellman operator, it does one step of the VFI
      ! This version is partially vectorized (eliminated loop over a')
      ! For a loop-based version, see subroutine bellman_op
      ! For a fully vectorized version, see subroutine bellman_op_vec2
      implicit none
      ! Inputs:
      real(kind=wp), intent(in) :: v(:, :)
      ! Outputs:
      real(kind=wp), intent(out) :: v_new(:, :), pol_ap(:, :)
      ! Locals:
      real(kind=wp), allocatable :: EV(:, :), cons(:), v_temp(:)
      integer :: a_c, z_c, ap_ind, ap_c
      real(kind=wp) :: a_val, z_val

      allocate (cons(n_a), v_temp(n_a))

      ! Compute expected value function EV(a',z) = sum over z' of v(a',z')*z_tran(z,z')
      EV = matmul(v, transpose(z_tran))

      ! Step through the state space
      !$omp parallel if (par_fortran==1) default(shared) private(z_c,a_c,a_val,z_val,ap_ind,ap_c,cons,v_temp)
      !$omp do collapse(2)
      do z_c = 1, n_z
         do a_c = 1, n_a
            ! Current state
            a_val = a_grid(a_c)
            z_val = z_grid(z_c)
            ! Compute consumption
            cons = R*a_val + z_val - ap_grid ! (n_ap,1)
            ! NOTE: where and merge are slower than forall
            ! NOTE: forall and do concurrent are equivalent with ifort
            ! but do concurrent is very slow with ifx!
            !where (cons>0.0d0)
            !    util = f_util(cons)  ! (n_ap,1)
            !elsewhere
            !    util = large_negative
            !end where
            !util = merge(f_util(cons),large_negative,cons>0.0d0)
            ! v_temp = large_negative
            ! do concurrent (ap_c=1:n_a, cons(ap_c)>0.0d0)
            !    v_temp(ap_c) = f_util(cons(ap_c))+beta*EV(ap_c,z_c)
            ! enddo
            v_temp = large_negative
            forall (ap_c=1:n_a, cons(ap_c) > 0.0d0)
               v_temp(ap_c) = f_util(cons(ap_c)) + beta*EV(ap_c, z_c)
            end forall
            ap_ind = maxloc(v_temp, dim=1)
            v_new(a_c, z_c) = v_temp(ap_ind)
            pol_ap(a_c, z_c) = a_grid(ap_ind)
         end do
      end do
      !$omp enddo
      !$omp end parallel

   end subroutine bellman_op_vec

   subroutine bellman_op_vec2(v_new, pol_ap, v)
      ! Purpose: Bellman operator, it does one step of the VFI
      ! This version is partially vectorized (eliminated loops over a and a')
      ! For a loop-based version, see subroutine bellman_op
      implicit none
      ! Inputs:
      real(kind=wp), intent(in) :: v(:, :)
      ! Outputs:
      real(kind=wp), intent(out) :: v_new(:, :), pol_ap(:, :)
      ! Locals:
      real(kind=wp), allocatable :: EV(:, :), cons(:, :), v_temp(:, :)
      integer, allocatable :: ap_ind(:)
      integer :: z_c, ap_c, a_c
      real(kind=wp) :: z_val

      allocate (cons(n_ap, n_a), v_temp(n_ap, n_a))

      ! Compute expected value function EV(a',z) = sum over z' of v(a',z')*z_tran(z,z')
      EV = matmul(v, transpose(z_tran))

      ! Step through the state space
      do z_c = 1, n_z
         ! Current state
         z_val = z_grid(z_c)
         ! Compute consumption
         cons = R*a_grid2 + z_val - ap_grid2 ! (n_ap,n_a)
         v_temp = large_negative
         forall (a_c=1:n_a, ap_c=1:n_ap, cons(ap_c, a_c) > 0.0d0)
            v_temp(ap_c, a_c) = f_util(cons(ap_c, a_c)) + beta*EV(ap_c, z_c)
         end forall
         ! v_temp is a 2-dim array(a',a), we need to find the max along a'
         ap_ind = maxloc(v_temp, dim=1)
         v_new(:, z_c) = maxval(v_temp, dim=1)
         pol_ap(:, z_c) = a_grid(ap_ind)
      end do

   end subroutine bellman_op_vec2

end module mod_vfi
!===============================================================================!

program main
   use mod_globals
   use mod_vfi
   use mod_utilities, only: disp, writescalar, writedim
   use mod_numerical, only: linspace
   use omp_lib
   use, intrinsic :: iso_fortran_env, only : compiler_version
   use, intrinsic :: iso_fortran_env, only : compiler_options

   implicit none

   ! Declarations
   integer, parameter :: wp=kind(0.0d0)
   integer :: istat, it, verbose, maxit, z_c, a_c, method
   real(kind=wp) :: t1, t2, err, tol
   !integer, allocatable :: pol_ap_ind(:,:)
   real(kind=wp), allocatable :: V(:, :), V_new(:, :), pol_ap(:, :), pol_c(:, :)
   namelist /args/ method, verbose, par_fortran

   write (*, *) "Starting program..."

   ! ---------------------- Set numerical parameters ----------------------!
   method = 1       ! 0: loop-based, 1: vectorized, 2: vectorized2
   verbose = 1       ! If 1, print iteration info
   par_fortran = 0       ! If 1, parallelize with OpenMP
   tol = 1.0d-6  ! Tolerance for VFI
   maxit = 1000    ! Maximum number of iterations

   ! ---------------------- Set economic parameter values ----------------------!
   R = 1.03d0 ! Gross interest rate, R = 1+r
   beta = 0.96d0 ! Discount factor
   gamma = 1.0d0  ! Coeff. rel. risk aversion (if 1, log utility)

   ! ---------------------- Set grids and probabilities ------------------------!
   a_min = 0.0d0 ! Minimum value of asset grid
   a_max = 4.0d0 ! Maximum value of asset grid
   n_a = 1000  ! Number of grid points for a (state variable)
   ! ---------------------------------------------------------------------------!
   call command_line()
   write (*, nml=args)
   ! ---------------------------------------------------------------------------!
   n_ap = n_a   ! Number of grid points for a' (choice variable)
   a_grid = linspace(a_min, a_max, n_a)
   ap_grid = a_grid
   !call disp("a_grid = ",a_grid(1:100))
   ! Expand a_grid to 2D
   allocate (a_grid2(n_a, n_a), ap_grid2(n_a, n_a), stat=istat)
   if (istat /= 0) error stop "Allocation of a_grid2 and ap_grid2 failed!"

   do a_c = 1, n_a
      a_grid2(a_c, :) = a_grid !varies along rows
      ap_grid2(:, a_c) = ap_grid !varies along columns
   end do

   !call disp("a_grid2  = ",a_grid2)
   !call disp("ap_grid2 = ",ap_grid2)
   !pause

   n_z = 2 ! Number of grid points for z (shock)
   allocate (z_grid(n_z), z_tran(n_z, n_z), stat=istat)
   if (istat /= 0) error stop "Allocation of z_grid and z_tran failed!"

   z_grid = [0.5d0, 1.0d0]
   ! Transition matrix for z (rows sum to 1)
   z_tran(1, :) = [0.6d0, 0.4d0]
   z_tran(2, :) = [0.05d0, 0.95d0]

   !call disp("z_grid = ",z_grid)
   !call disp("z_tran = ",z_tran)

   ! ---------------------- Initialize value function ----------------------!
   allocate (V(n_a, n_z), V_new(n_a, n_z), pol_ap(n_a, n_z), pol_c(n_a, n_z), stat=istat)
   if (istat /= 0) error stop "Allocation of V and pol_ap failed!"

   V = 0.0d0 ! Initial guess for value function
   it = 1
   err = tol + 1.0d0

   write (*, *) "Starting VFI..."
   write (*, *) "Method = ", method
   t1 = omp_get_wtime()

   do while (err > tol .and. it < maxit)

      if (method == 0) then
         call bellman_op(V_new, pol_ap, V)
      elseif (method == 1) then
         call bellman_op_vec(V_new, pol_ap, V)
      elseif (method == 2) then
         call bellman_op_vec2(V_new, pol_ap, V)
      else
         error stop "Invalid method!"
      end if
      err = maxval(abs(V_new - V))

      if (verbose == 1) then
         write (*, *) "Iteration ", it, " error = ", err
      end if

      ! Update
      V = V_new
      it = it + 1
   end do

   t2 = omp_get_wtime()

   write (*, *) "============================================================="
   write (*,'(1x,*(g0,1x))')                                        &
                'VFI runs for',      real(t2 - t1),     'seconds.', &
                'compiled by',       compiler_version(),            &
                'using the options', compiler_options()
   write (*, *) "============================================================="
   !pause

   ! Compute policy function for consumption
   do z_c = 1, n_z
      pol_c(:, z_c) = R*a_grid + z_grid(z_c) - pol_ap(:, z_c)
   end do

   write (*, *) "Writing results to file..."
   call sub_write()

contains

   subroutine sub_write()
      ! NOTE: The subroutines writescalar and writedim are defined in mod_utilities
      ! If you don't have mod_utilities, just comment out the calls to these subroutines
      call writescalar(n_a, trim(savedir)//'n_a.txt')
      call writescalar(n_z, trim(savedir)//'n_z.txt')

      ! Write grids and probs to files
      call writedim(z_grid, trim(savedir)//'z_grid.txt')
      call writedim(a_grid, trim(savedir)//'a_grid.txt')
      call writedim(V, trim(savedir)//'V.txt')
      call writedim(pol_ap, trim(savedir)//'pol_ap.txt')
      call writedim(pol_c, trim(savedir)//'pol_c.txt')
      !call writedim(mu,trim(savedir)//'mu.txt')

   end subroutine sub_write

   subroutine command_line()
      implicit none
      character(len=255)           :: message ! use for I/O error messages
      character(len=:), allocatable :: string  ! stores command line argument
      integer                      :: ios
      integer                      :: sum

      string = get_namelist()  ! return command line arguments as NAMELIST input
      read (string, nml=args, iostat=ios, iomsg=message) ! internal read of namelist
      if (ios .ne. 0) then
         write (*, '("ERROR:",i0,1x,a)') ios, trim(message)
         write (*, *) 'COMMAND OPTIONS ARE'
         write (*, nml=args)
         stop 1
      end if
   end subroutine command_line

   function get_namelist() result(string)
      character(len=:), allocatable :: string         ! stores command line argument
      integer :: command_line_length
      call get_command(length=command_line_length)   ! get length needed to hold command
      allocate (character(len=command_line_length) :: string)
      call get_command(string)
      ! trim off command name and get command line arguments
      string = adjustl(string)//' '                    ! assuming command verb does not have spaces in it
      string = string(index(string, ' '):)
      string = "&args "//string//" /"                   ! add namelist prefix and terminator
   end function get_namelist
end program main
