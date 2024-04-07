# timing1

Skip to main content
Fortran Discourse
Performance of vectorized code in ifort and ifx

    ​

Performance of vectorized code in ifort and ifx
Apr 7
11m
aledinola
1h

I created a small Fortran example of discrete dynamic programming that occurs quite often in macroeconomics. It is an optimal savings problem for a consumer who received idiosyncratic shocks to their income.

I was interested in comparing loop-based Fortran code vs vectorized code, both in ifort and in ifx. Just to avoid misunderstandings, the definition of “vectorizing” that I am using is the following: I vectorize if I replace this code

do i=1,n
yvec(i) = exp(xvec(i))
enddo

with the following code:

yvec = exp(xvec)

With this in mind, I coded the dynamic programming example in three different subroutines:

    bellman_op, which is the benchmark based entirely on loops (three do loops)
    bellman_op_vec, which is partially vectorized (I eliminated the innermost loop over a’)
    bellman_op_vec2, which is even more vectorized (I eliminated the two innermost loops)

bellman_op is here:

subroutine bellman_op(v_new,pol_ap, v)
    ! Purpose: Bellman operator, it does one step of the VFI
    ! This version is loop-based and can be parallelized with OpenMP
    ! For vectorized versions, see subroutine bellman_op_vec, bellman_op_vec2
    implicit none
    ! Inputs:
    real(8), intent(in) :: v(:,:)
    ! Outputs:
    real(8), intent(out) :: v_new(:,:), pol_ap(:,:)
    ! Locals:
    real(8), allocatable :: EV(:,:)
    integer :: a_c, z_c, ap_c, ap_ind
    real(8) :: a_val, z_val, aprime_val, cons, v_max, v_temp

    ! Compute expected value function EV(a',z) = sum over z' of v(a',z')*z_tran(z,z')
    EV = matmul(v,transpose(z_tran))

    ! Step through the state space
    !$omp parallel if (par_fortran==1) default(shared) private(z_c,a_c,a_val,z_val,&
    !$ v_max,ap_ind,ap_c,aprime_val,cons,v_temp)
    !$omp do collapse(2)
    do z_c=1,n_z
        do a_c=1,n_a
            ! Current states
            a_val = a_grid(a_c)
            z_val = z_grid(z_c)
            ! Initialize v_new
            v_max = large_negative
            ap_ind = 0
            ! Choose a' optimally by stepping through all possible values
            do ap_c=1,n_a
                aprime_val = ap_grid(ap_c)
                cons = R*a_val + z_val - aprime_val
                if (cons>0.0d0) then
                    v_temp = f_util(cons) + beta*EV(ap_c,z_c)
                    !v_temp = f_util(cons) + beta*sum(v(ap_c,:)*z_tran(z_c,:))
                    if (v_temp>v_max) then
                        v_max = v_temp
                        ap_ind = ap_c
                    end if
                endif
            enddo !end a'
            v_new(a_c,z_c)  = v_max
            pol_ap(a_c,z_c) = ap_grid(ap_ind)
        enddo
    enddo
    !$omp enddo
    !$omp end parallel
    
    end subroutine bellman_op

bellman_op_vec:

subroutine bellman_op_vec(v_new,pol_ap, v)
     ! Purpose: Bellman operator, it does one step of the VFI
     ! This version is partially vectorized (eliminated loop over a')
     ! For a loop-based version, see subroutine bellman_op
     ! For a fully vectorized version, see subroutine bellman_op_vec2
     implicit none
     ! Inputs:
     real(8), intent(in) :: v(:,:)
     ! Outputs:
     real(8), intent(out) :: v_new(:,:), pol_ap(:,:)
     ! Locals:
     real(8), allocatable :: EV(:,:), cons(:), v_temp(:)
     integer :: a_c,z_c,ap_ind,ap_c
     real(8) :: a_val,z_val
 
     allocate(cons(n_a),v_temp(n_a))
 
   ! Compute expected value function EV(a',z) = sum over z' of v(a',z')*z_tran(z,z')
     EV = matmul(v,transpose(z_tran))
 
     ! Step through the state space
     !$omp parallel if (par_fortran==1) default(shared) private(z_c,a_c,a_val,z_val,&
     !$ ap_ind,ap_c,cons,v_temp)
     !$omp do collapse(2)
     do z_c=1,n_z
         do a_c=1,n_a
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
             forall (ap_c=1:n_a, cons(ap_c)>0.0d0)
                v_temp(ap_c) = f_util(cons(ap_c))+beta*EV(ap_c,z_c)
             end forall
             ap_ind = maxloc(v_temp,dim=1)
             v_new(a_c,z_c) = v_temp(ap_ind)
             pol_ap(a_c,z_c) = a_grid(ap_ind)
         enddo
     enddo
     !$omp enddo
     !$omp end parallel
 
     end subroutine bellman_op_vec

and finally bellman_op_vec2:

subroutine bellman_op_vec2(v_new,pol_ap, v)
     ! Purpose: Bellman operator, it does one step of the VFI
     ! This version is partially vectorized (eliminated loops over a and a')
     ! For a loop-based version, see subroutine bellman_op
     implicit none
     ! Inputs:
     real(8), intent(in) :: v(:,:)
     ! Outputs:
     real(8), intent(out) :: v_new(:,:), pol_ap(:,:)
     ! Locals:
     real(8), allocatable :: EV(:,:), cons(:,:), v_temp(:,:)
     integer, allocatable :: ap_ind(:)
     integer :: z_c,ap_c,a_c
     real(8) :: z_val
 
     allocate(cons(n_ap,n_a),v_temp(n_ap,n_a))
 
     ! Compute expected value function EV(a',z) = sum over z' of v(a',z')*z_tran(z,z')
     EV = matmul(v,transpose(z_tran))
 
     ! Step through the state space
     do z_c=1,n_z
         ! Current state
         z_val = z_grid(z_c) 
         ! Compute consumption
         cons = R*a_grid2 + z_val - ap_grid2 ! (n_ap,n_a)
         v_temp = large_negative
         forall (a_c=1:n_a,ap_c=1:n_ap, cons(ap_c,a_c)>0.0d0)
             v_temp(ap_c,a_c) = f_util(cons(ap_c,a_c))+beta*EV(ap_c,z_c)
         end forall
         ! v_temp is a 2-dim array(a',a), we need to find the max along a' 
         ap_ind = maxloc(v_temp,dim=1)
         v_new(:,z_c) = maxval(v_temp,dim=1)
         pol_ap(:,z_c) = a_grid(ap_ind)
     enddo
     
     end subroutine bellman_op_vec2

What I have found, to my partial surprise, is that the vectorized versions are either the same or worse than the loops-based version and that ifx is much slower than ifort (these might be two totally separate issues, but I didn’t want to write two posts). Here are the timings (I report the median time across 10 runs for each case)

|Column 1 | Column 2 | Column 3 | Column 4|

|--- | --- | --- | ---|

|Version | Time in secs (IFORT) | Time in secs (IFX) | |

|All loops (bellman_op) | 0.78 | 1.63 | |
|Vectorized loop over a' (bellman_op_vec) | 0.75 | 1.8 | |
|Vectorized loops over a' and a (bellman_op_vec2) | 1.3 | 2.18 | |

Other observations:

    When working with arrays in the vectorized code, I found that the constructs where and merge are slower than forall or do concurrent. In ifort, forall and do concurrent give similar performance, whereas in ifx do concurrent is super slow. All the timings reported above are based on forall for the vectorized versions.
    As compilation flags I use the following
    /O3 /fast /Qmkl /Qopenmp /Qparallel
    but ifx returns the warning that /Qparallel is not supported. This might be a problem because can the compiler optimize the forall and do concurrent without the parallel option?

If someone is interested in replicating my code, I attach here below the full version. The above code is also available on github here

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
    
    contains
    
    function linspace(my_start, my_stop, n)
        ! Purpose: replicate Matlab function <linspace>
        implicit none
        ! Inputs:
        integer :: n
        real(8) :: my_start, my_stop
        ! Function result:
        real(8) :: linspace(n)
        ! Locals:
        integer :: i
        real(8) :: step, grid(n)

        if (n == 1) then
            grid(n) = (my_start+my_stop)/2.d0
        elseif (n>=2) then
            grid(1) = my_start
            if (n>2) then
                step = (my_stop-my_start)/(real(n-1,8))
                do i = 2, n-1
                    grid(i) = grid(i-1) + step
                end do
            end if
            grid(n) = my_stop
        endif
        linspace = grid
    
    end function linspace

end module mod_numerical
!===============================================================================!

module mod_globals
implicit none

real(8), parameter :: large_negative = -1.0d10
character(len=*), parameter :: savedir = "output\"
real(8) :: R, beta, gamma, a_min, a_max, rho, sig_z
integer :: n_a, n_z, n_ap, par_fortran

real(8), allocatable :: a_grid(:), ap_grid(:), z_grid(:), z_tran(:,:)
real(8), allocatable :: a_grid2(:,:), ap_grid2(:,:) !needed for vectorized code

contains
    
    elemental function f_util(c) result(util)
        ! Purpose: utility function
        ! Assumption: c is positive, this is enforced
        ! elsewhere in the code
        real(8), intent(in) :: c
        real(8) :: util
        
        if (abs(gamma-1.0d0)<1.0d-6) then
            util = log(c)
        else
            util = c**(1.0d0-gamma)/(1.0d0-gamma)
        endif
    
    end function f_util
    
end module mod_globals
!===============================================================================!

module mod_vfi 
    use mod_globals
    use omp_lib
    implicit none
    
    private
    public :: bellman_op, bellman_op_vec, bellman_op_vec2
    
    contains
    
    subroutine bellman_op(v_new,pol_ap, v)
    ! Purpose: Bellman operator, it does one step of the VFI
    ! This version is loop-based and can be parallelized with OpenMP
    ! For vectorized versions, see subroutine bellman_op_vec, bellman_op_vec2
    implicit none
    ! Inputs:
    real(8), intent(in) :: v(:,:)
    ! Outputs:
    real(8), intent(out) :: v_new(:,:), pol_ap(:,:)
    ! Locals:
    real(8), allocatable :: EV(:,:)
    integer :: a_c, z_c, ap_c, ap_ind
    real(8) :: a_val, z_val, aprime_val, cons, v_max, v_temp

    ! Compute expected value function EV(a',z) = sum over z' of v(a',z')*z_tran(z,z')
    EV = matmul(v,transpose(z_tran))

    ! Step through the state space
    !$omp parallel if (par_fortran==1) default(shared) private(z_c,a_c,a_val,z_val,&
    !$ v_max,ap_ind,ap_c,aprime_val,cons,v_temp)
    !$omp do collapse(2)
    do z_c=1,n_z
        do a_c=1,n_a
            ! Current states
            a_val = a_grid(a_c)
            z_val = z_grid(z_c)
            ! Initialize v_new
            v_max = large_negative
            ap_ind = 0
            ! Choose a' optimally by stepping through all possible values
            do ap_c=1,n_a
                aprime_val = ap_grid(ap_c)
                cons = R*a_val + z_val - aprime_val
                if (cons>0.0d0) then
                    v_temp = f_util(cons) + beta*EV(ap_c,z_c)
                    !v_temp = f_util(cons) + beta*sum(v(ap_c,:)*z_tran(z_c,:))
                    if (v_temp>v_max) then
                        v_max = v_temp
                        ap_ind = ap_c
                    end if
                endif
            enddo !end a'
            v_new(a_c,z_c)  = v_max
            pol_ap(a_c,z_c) = ap_grid(ap_ind)
        enddo
    enddo
    !$omp enddo
    !$omp end parallel
    
    end subroutine bellman_op

    subroutine bellman_op_vec(v_new,pol_ap, v)
    ! Purpose: Bellman operator, it does one step of the VFI
    ! This version is partially vectorized (eliminated loop over a')
    ! For a loop-based version, see subroutine bellman_op
    ! For a fully vectorized version, see subroutine bellman_op_vec2
    implicit none
    ! Inputs:
    real(8), intent(in) :: v(:,:)
    ! Outputs:
    real(8), intent(out) :: v_new(:,:), pol_ap(:,:)
    ! Locals:
    real(8), allocatable :: EV(:,:), cons(:), v_temp(:)
    integer :: a_c,z_c,ap_ind,ap_c
    real(8) :: a_val,z_val

    allocate(cons(n_a),v_temp(n_a))

    ! Compute expected value function EV(a',z) = sum over z' of v(a',z')*z_tran(z,z')
    EV = matmul(v,transpose(z_tran))

    ! Step through the state space
    !$omp parallel if (par_fortran==1) default(shared) private(z_c,a_c,a_val,z_val,&
    !$ ap_ind,ap_c,cons,v_temp)
    !$omp do collapse(2)
    do z_c=1,n_z
        do a_c=1,n_a
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
            forall (ap_c=1:n_a, cons(ap_c)>0.0d0)
               v_temp(ap_c) = f_util(cons(ap_c))+beta*EV(ap_c,z_c)
            end forall
            ap_ind = maxloc(v_temp,dim=1)
            v_new(a_c,z_c) = v_temp(ap_ind)
            pol_ap(a_c,z_c) = a_grid(ap_ind)
        enddo
    enddo
    !$omp enddo
    !$omp end parallel

    end subroutine bellman_op_vec

subroutine bellman_op_vec2(v_new,pol_ap, v)
    ! Purpose: Bellman operator, it does one step of the VFI
    ! This version is partially vectorized (eliminated loops over a and a')
    ! For a loop-based version, see subroutine bellman_op
    implicit none
    ! Inputs:
    real(8), intent(in) :: v(:,:)
    ! Outputs:
    real(8), intent(out) :: v_new(:,:), pol_ap(:,:)
    ! Locals:
    real(8), allocatable :: EV(:,:), cons(:,:), v_temp(:,:)
    integer, allocatable :: ap_ind(:)
    integer :: z_c,ap_c,a_c
    real(8) :: z_val

    allocate(cons(n_ap,n_a),v_temp(n_ap,n_a))

    ! Compute expected value function EV(a',z) = sum over z' of v(a',z')*z_tran(z,z')
    EV = matmul(v,transpose(z_tran))

    ! Step through the state space
    do z_c=1,n_z
        ! Current state
        z_val = z_grid(z_c) 
        ! Compute consumption
        cons = R*a_grid2 + z_val - ap_grid2 ! (n_ap,n_a)
        v_temp = large_negative
        forall (a_c=1:n_a,ap_c=1:n_ap, cons(ap_c,a_c)>0.0d0)
            v_temp(ap_c,a_c) = f_util(cons(ap_c,a_c))+beta*EV(ap_c,z_c)
        end forall
        ! v_temp is a 2-dim array(a',a), we need to find the max along a' 
        ap_ind = maxloc(v_temp,dim=1)
        v_new(:,z_c) = maxval(v_temp,dim=1)
        pol_ap(:,z_c) = a_grid(ap_ind)
    enddo
    
    end subroutine bellman_op_vec2

end module mod_vfi
!===============================================================================!

program main
    use mod_globals
    use mod_vfi
    use mod_utilities, only: disp,writescalar,writedim
    use mod_numerical, only: linspace
    use omp_lib 

    implicit none
    
    ! Declarations
    integer :: istat,it,verbose,maxit,z_c,a_c,method 
    real(8) :: t1,t2,err,tol
    !integer, allocatable :: pol_ap_ind(:,:)
    real(8), allocatable :: V(:,:),V_new(:,:),pol_ap(:,:),pol_c(:,:)
    
    write(*,*) "Starting program..."
    
    ! ---------------------- Set numerical parameters ----------------------!
    method      = 1       ! 0: loop-based, 1: vectorized, 2: vectorized2
    verbose     = 1       ! If 1, print iteration info
    par_fortran = 0       ! If 1, parallelize with OpenMP
    tol         = 1.0d-6  ! Tolerance for VFI
    maxit       = 1000    ! Maximum number of iterations

    ! ---------------------- Set economic parameter values ----------------------!
    R     = 1.03d0 ! Gross interest rate, R = 1+r
    beta  = 0.96d0 ! Discount factor
    gamma = 1.0d0  ! Coeff. rel. risk aversion (if 1, log utility)
    
    ! ---------------------- Set grids and probabilities ------------------------!
    a_min = 0.0d0 ! Minimum value of asset grid
    a_max = 4.0d0 ! Maximum value of asset grid
    n_a   = 1000  ! Number of grid points for a (state variable)
    n_ap  = n_a   ! Number of grid points for a' (choice variable)
    a_grid  = linspace(a_min,a_max,n_a)
    ap_grid = a_grid
    !call disp("a_grid = ",a_grid(1:100))
    ! Expand a_grid to 2D
    allocate(a_grid2(n_a,n_a),ap_grid2(n_a,n_a),stat=istat)
    if (istat/=0) error stop "Allocation of a_grid2 and ap_grid2 failed!"
    
    do a_c=1,n_a
        a_grid2(a_c,:)  = a_grid !varies along rows
        ap_grid2(:,a_c) = ap_grid !varies along columns
    enddo

    !call disp("a_grid2  = ",a_grid2)
    !call disp("ap_grid2 = ",ap_grid2)
    !pause 
    
    n_z   = 2 ! Number of grid points for z (shock)
    allocate(z_grid(n_z),z_tran(n_z,n_z),stat=istat)
    if (istat/=0) error stop "Allocation of z_grid and z_tran failed!"
    
    z_grid = [0.5d0,1.0d0]
    ! Transition matrix for z (rows sum to 1)
    z_tran(1,:) = [0.6d0,0.4d0]
    z_tran(2,:) = [0.05d0,0.95d0]
    
    !call disp("z_grid = ",z_grid)
    !call disp("z_tran = ",z_tran)

    ! ---------------------- Initialize value function ----------------------!
    allocate(V(n_a,n_z),V_new(n_a,n_z),pol_ap(n_a,n_z),pol_c(n_a,n_z),stat=istat)
    if (istat/=0) error stop "Allocation of V and pol_ap failed!"

    V = 0.0d0 ! Initial guess for value function
    it = 1
    err = tol+1.0d0
    
    write(*,*) "Starting VFI..."
    write(*,*) "Method = ",method
    t1 = omp_get_wtime()

    do while (err>tol .and. it<maxit)

        if (method==0) then
            call bellman_op(V_new,pol_ap, V)
        elseif (method==1) then
            call bellman_op_vec(V_new,pol_ap, V)
        elseif (method==2) then
            call bellman_op_vec2(V_new,pol_ap, V)
        else
            error stop "Invalid method!"
        endif
        err = maxval(abs(V_new-V))

        if (verbose==1) then
            write(*,*) "Iteration ",it," error = ",err
        endif

        ! Update
        V = V_new
        it = it + 1
    enddo
    
    t2 = omp_get_wtime()
    
    write(*,*) "============================================================="
    write(*,*) 'VFI runs for',real(t2-t1),'seconds.'
    write(*,*) "============================================================="
    !pause 
    
    ! Compute policy function for consumption
    do z_c=1,n_z
        pol_c(:,z_c) = R*a_grid + z_grid(z_c) - pol_ap(:,z_c)
    enddo
    
    write(*,*) "Writing results to file..."
    call sub_write()
    
    
    contains
    
    subroutine sub_write()
    ! NOTE: The subroutines writescalar and writedim are defined in mod_utilities
    ! If you don't have mod_utilities, just comment out the calls to these subroutines
        call writescalar(n_a,trim(savedir)//'n_a.txt')
        call writescalar(n_z,trim(savedir)//'n_z.txt')
        
        ! Write grids and probs to files
        call writedim(z_grid,trim(savedir)//'z_grid.txt')
        call writedim(a_grid,trim(savedir)//'a_grid.txt')
        call writedim(V,trim(savedir)//'V.txt')
        call writedim(pol_ap,trim(savedir)//'pol_ap.txt')
        call writedim(pol_c,trim(savedir)//'pol_c.txt')
        !call writedim(mu,trim(savedir)//'mu.txt')
    
    end subroutine sub_write
    
end program main

Note that in the code there are some openMP directives but they are not used at the moment (since the variable par_fortran is set to zero).

Any comment or suggestion on how to improve the code is welcome! (Including on how to get rid of the smile face in some array definitions :slight_smile:
Overall, I find it quite odd that ifx gives a code that is significantly slower compared to ifort and that now Intel says that ifort is deprecated.

    created
    1h
    last reply
    11m
    3
    replies
    11
    views
    2
    users
    2
    likes

urbanjost
1h

    ```text
    put text here like the table to get fixed-font text
    ``
    ```fortran
    put fortran code here
    ```

start the line with three grave and an optional type like “text” or “bash” or “fortran” and then follow that line with your text or code and then close the text block with a line that is just three grave characters
Markdown Guide:

https://markdown-it.github.io/

You can edit what you already entered and add the block delimiters to make the post text appear properly
1
aledinola
24m

Thanks for your suggestions about typesetting. I hope now the code is easier to read. I did not manage to improve much the layout of the table, though.
ivanpribec
Regular
11m

Two comments:

    mod_utilities are missing, so the internal subroutine sub_write can be removed
    The OpenMP directives are not accepted by gfortran, due to the way the lines are broken. gfortran expects full clauses on a single line, e.g.

    !$omp parallel if (par_fortran==1) default(shared) &
    !$omp private(z_c,a_c,a_val,z_val, ap_ind,ap_c,cons,v_temp)
    !$omp do collapse(2)

1

You will see a count of new replies because you posted a reply to this topic.

 Replies  Views  Activity
A patch was sent to Voyager 2 yesterday 7
 
5.7k  22h
Idioms for exception handling 1
 
170  23h
The counter-intuitive rise of Python in scientific computing 122
 
19.9k  5d
Research articles using Fortran 1
 
13.8k  6d
Introspection in Fortran for generic file I/O libraries 1
 
955  6d
There are 439 unread topics remaining, or view latest topics
Powered by Discourse
