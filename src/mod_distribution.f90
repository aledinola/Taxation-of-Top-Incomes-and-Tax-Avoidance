module mod_distribution
! #VC# V44
!===========================================================================
! INVARIANT DISTRIBUTION MODULE 
! DESCRIPTION:
!   This module stores procedures to compute invariant distribution
!   The main procedure is get_distribution, to be called by the main prog.
! PROGRAMS:
!   get_distribution:  distribution code, with OpenMP
!   get_KN, get_Aggregates
!===========================================================================

use mod_globals
use mod_numerical, only: myerror,outerprod2,locate,mtimes_matmat
use omp_lib
use mod_utilities

implicit none

    
    contains 
!=======================================================================!
    
subroutine find_loc(jstar,omega, a_grid,a_opt)
    implicit none 
    !-----------------------------------------------------------------------!
    ! Purpose: Find location and interp weights of a point on a grid
    ! Declare input and output:
    real(8), intent(in) :: a_grid(:)
    real(8), intent(in) :: a_opt
    integer, intent(out) :: jstar
    real(8), intent(out) :: omega
    ! Declare locals:
    integer :: n
    !-----------------------------------------------------------------------!
    ! Body of find_loc:
    n = size(a_grid)
    jstar = max(min(locate(a_grid,a_opt),n-1),1)
    ! Weight on a_grid(j)
    omega = (a_grid(jstar+1)-a_opt)/(a_grid(jstar+1)-a_grid(jstar))
    omega = max(min(omega,1.0d0),0.0d0)

end subroutine find_loc  
!=======================================================================!

subroutine distrib_one_step(muy_update,mur_update,aopt_ind,omega,aoptr_ind,omegar, &
	temp_ret,distrib,pol)  
implicit none
! Inputs and outputs
real(8),intent(out) :: muy_update(:,:,:,:),mur_update(:)
real(8), intent(in) :: omega(:,:,:,:),omegar(:),temp_ret(:,:)
integer, intent(in) :: aopt_ind(:,:,:,:),aoptr_ind(:)
type(distribfun_m), intent(in) :: distrib 
type(polfun), intent(in) :: pol    
! Local variables
integer :: it,ie,ia,io,iap,iz_min
real(8) :: io_prob,omega_val,check
real(8),allocatable :: muy_z(:,:), muy_z_temp(:,:,:)
   

	muy_update = 0.0d0
	mur_update = 0.0d0
	! Transitions from young given omega and aopt_ind
	
	! We use Gordon's trick to update the distribution of young agents. In particular, 
    ! we first update the endogenous state a, keeping the shocks (eps,theta) at their current 
    ! value. Second, we update the distribution using the transition (eps,theta,eps',theta')
    ! See also Matlab file mu_one_step.m
	
	!$omp parallel if (par_fortran_mu==1) default(shared) private(it,ie,ia,io,iz_min, &
	!$ io_prob,iap,omega_val)  reduction(+:muy_update,mur_update)
    !$omp do collapse(2)
    do it = 1,ntheta ! current shock theta
        do ie = 1,neps ! current shock epsilon
            do ia = 1,na ! current asset holdings
            	do iz_min = 1,nz
				do io = 1,no ! current occupation
					iap       = aopt_ind(ia,ie,it,io)
					omega_val = omega(ia,ie,it,io)
		            io_prob   = pol%occpol(ia,ie,it,iz_min,io)
					!With prob=1-pr_ret, he does not retire
					! with prob=omega_val, next period's asset is on the grid point iap
					
					muy_update(iap,ie,it,io) = muy_update(iap,ie,it,io) + &
						(1.0d0-pr_ret)*omega_val*io_prob*distrib%muy(ia,ie,it,iz_min)
					! with prob=1-omega_val, next period's asset is on the grid point iap+1
					muy_update(iap+1,ie,it,io) = muy_update(iap+1,ie,it,io) + &
						(1.0d0-pr_ret)*(1.0d0-omega_val)*io_prob*distrib%muy(ia,ie,it,iz_min)
						
					!With prob=pr_ret, he retires
					mur_update(iap) =    mur_update(iap)+ &
						pr_ret*omega_val*io_prob*distrib%muy(ia,ie,it,iz_min)
					mur_update(iap+1) =    mur_update(iap+1)+ &
						pr_ret*(1.0d0-omega_val)*io_prob*distrib%muy(ia,ie,it,iz_min)	
						
				enddo !end io		
				enddo !end iz_min

            enddo !END loop over ia (current assets)   
        enddo !ie
    enddo !it
	!$omp enddo
  	!$omp end parallel
  	
  	! Second step of Gordon's trick (updates the distribution mu using the transition 
  	! matrix of eps and theta)
    ! Equivalent formulation to avoid temporary array
    do iz_min = 1,no
        muy_z_temp = muy_update(:,:,:,iz_min)
        !write(*,*) "before reshape"
  		muy_z = reshape(muy_z_temp,[na,neps*ntheta]) !(a',eps*theta)
        !write(*,*) "before matmul"
  		!muy_z = matmul(muy_z,P_eps_theta) ! (a',eps'*theta')
  		!write(*,*) "before mtimes_matmat"
        muy_z = mtimes_matmat(muy_z,P_eps_theta)
        !write(*,*) "before 2nd reshape"
  		muy_update(:,:,:,iz_min)   = reshape(muy_z,[na,neps,ntheta])
  	enddo

 	! Transitions from retired	
	 do ia = 1,na
		 !What assets does guy with choose for tomorrow?
		 iap       = aoptr_ind(ia)
		 omega_val = omegar(ia)

		 !With prob=1-pr_die, agent stays old
		 mur_update(iap)   = mur_update(iap)+(1.0d0-pr_die)*omega_val*distrib%mur(ia)
		 mur_update(iap+1) = mur_update(iap+1)+(1.0d0-pr_die)*(1.0d0-omega_val)*distrib%mur(ia)
		 !With prob=pr_die, agent dies and is replaced by young
		 
		 ! with prob=omega_val, next period's asset is on the grid point iap
		 muy_update(iap,:,:,nz) = muy_update(iap,:,:,nz)+temp_ret*pr_die*omega_val*distrib%mur(ia)
		 ! with prob=1-omega_val, next period's asset is on the grid point iap+1
		 muy_update(iap+1,:,:,nz) = muy_update(iap+1,:,:,nz)+temp_ret*pr_die*(1.0d0-omega_val)*distrib%mur(ia)
		 
	 enddo
  	
    !write(*,*) "Now check that mu_new sum to 1"
    !Check that mu_new sum to 1
    check = sum(muy_update) + sum(mur_update)
    if (abs(check-1.0d0)>0.0000001d0 ) then
    	write(*,*) "sum(muy_update) + sum(mur_update) = ",check
    	write(*,*) "sum(distrib%muy) + sum(distrib%mur) = ",sum(distrib%muy) + sum(distrib%mur)
    	write(*,*) "sum(muy_update),sum(distrib%muy) = ",sum(muy_update),sum(distrib%muy)
    	write(*,*) "sum(mur_update),sum(distrib%mur) = ",sum(mur_update),sum(distrib%mur)
        call myerror("get_distribution: sum(muy_update) + sum(mur_update) does NOT sum to one!")
    endif
    if (abs(sum(muy_update)-sum(distrib%muy))>0.0000001d0 ) then
        write(*,*) "mass of young agents, old vs new: ",sum(distrib%muy),sum(muy_update)
        call myerror("get_distribution: mass of young agents changes!")
    endif
    if (abs(sum(mur_update)-sum(distrib%mur))>0.0000001d0 ) then
        write(*,*) "mass of retirees, old vs new: ",sum(distrib%mur),sum(mur_update)
        call myerror("get_distribution: mass of retired agents changes!")
    endif

end subroutine distrib_one_step
!=======================================================================!
subroutine precompute_distrib(aopt_ind,omega,aoptr_ind,omegar, &
	temp_ret,pol)  
	! DESCRIPTION: Precompute weights and indexes for distribution iteration
implicit none
! Inputs and outputs
real(8), intent(out) :: omega(:,:,:,:),omegar(:),temp_ret(:,:)
integer, intent(out) :: aopt_ind(:,:,:,:),aoptr_ind(:)
type(polfun),intent(in) :: pol
! Local variables
integer :: it, ie, ia, io
real(8) :: aopt_val

!  Pre-compute aopt_ind, and omega
	
	!$omp parallel if (par_fortran==1) default(shared) private(it,ie,ia,io, &
	!$ aopt_val)
    !$omp do collapse(4)
    do it = 1,ntheta ! current shock theta
        do ie = 1,neps ! current shock epsilon
            do ia = 1,na ! current asset holdings
            	! First, update distribution of the young mu1_young
                do io = 1,no ! current occupation
					!What assets does guy with (a,eps,theta,z) choose for tomorrow?
					aopt_val = pol%apoly(ia,ie,it,io)
					call find_loc(aopt_ind(ia,ie,it,io),omega(ia,ie,it,io), a_grid,aopt_val)
                enddo ! io  
            enddo ! ia 
        enddo !ie
    enddo !it
	!$omp enddo
  	!$omp end parallel
  	
  	! Pre-compute aoptr_ind and omegar for retirees
  	 do ia = 1,na
		 io = nz
		 !What assets does guy with (a,eps,theta,nz) choose for tomorrow?
		 aopt_val = pol%apolr(ia)
		 call find_loc(aoptr_ind(ia),omegar(ia), a_grid,aopt_val)
	 enddo
  	
  	! Pre-compute temp_ret is (neps,ntheta) matrix
    temp_ret = outerprod2(prob_eps,prob_theta)

end subroutine precompute_distrib
!=======================================================================!

subroutine get_distribution(distrib,pol,exit_flag)
    
implicit none
!-----------------------------------------------------------!
!PURPOSE
! get_distribution computes the steady-state distribution
! It calls the follwoing subroutines:
! - precompute_distrib
! - distrib_one_step
!Declare inputs
type(distribfun_m), intent(inout) :: distrib 
type(polfun), intent(in) :: pol
!Declare outputs
logical, intent(out) :: exit_flag          ! If .flase., iteration failed

!Declare locals
real(8), allocatable :: muy_update(:,:,:,:),mur_update(:),omega(:,:,:,:),omegar(:)
integer, allocatable :: aopt_ind(:,:,:,:),aoptr_ind(:)
real(8) :: diff, diffy,diffr,t1, t2,check
real(8) :: temp(neps,ntheta),temp_ret(neps,ntheta),io_prob,aopt_val,omega_val
integer :: iter, io,it,ie,ia, istat,iap
!-----------------------------------------------------------!

if (display_mu == 1) then 
	write(*,*) "============================================================="
	write(*,*) "STATIONARY DISTRIBUTION ITERATION..."
	write(*,*) "============================================================="
endif
exit_flag = .true. !This will turn to FALSE if iteration scheme fails
diff      = 10.0d0
iter      = 0

! Check if the sum of mu0_young and mu0_old is 1
!check = sum(distrib%muy) + sum(distrib%mur)
!if (abs(check-1.0d0)>0.00001d0 ) then
!    call myerror("get_distribution: sum(muy) + sum(mur) does NOT sum to one!")
!endif

! Allocate muy_update and mur_update
allocate(muy_update(na,neps,ntheta,nz),mur_update(na))
allocate(omega(na,neps,ntheta,no),aopt_ind(na,neps,ntheta,no),omegar(na),aoptr_ind(na))

t1 = omp_get_wtime()
!call CPU_TIME(t1)

call precompute_distrib(aopt_ind,omega,aoptr_ind,omegar, &
	temp_ret,pol)  

!if (verbose>=1 .and. do_GE==0) write(*,*) "precompute_distrib done!"	

do while (diff > tol_dist .and. iter <= maxit_dist)
    
    iter = iter + 1
    
    !write(*,*) "Before distrib_one_step"
    call distrib_one_step(muy_update,mur_update,aopt_ind,omega,aoptr_ind,omegar, &
	temp_ret,distrib,pol)  
    !write(*,*) "After distrib_one_step"
	
    !Error for updating
    diffy = maxval( ( abs( distrib%muy-muy_update )  ) )
    diffr = maxval( ( abs( distrib%mur-mur_update )  ) )
    diff = max(diffy,diffr)
    !Updating
    distrib%muy = muy_update
    distrib%mur = mur_update
    
    if (display_mu==1) then
        !write(*,*) "-----------------------------------"
        write(*,*) "Iteration: ", iter
        write(*,"(' diff: ',f12.9)") diff
        write(*,"(' diffy: ',f12.9)") diffy
        write(*,"(' diffr: ',f12.9)") diffr
        write(*,*) "-----------------------------------"
    endif
        
enddo !while loop

t2 = omp_get_wtime()
!call CPU_TIME(t2)

if (verbose>=1) then
	!write(*,*) "============================================================="
	write(*,*) "STATIONARY DISTRIBUTION FOUND"
	write(*,*) "============================================================="
	write(*,*) 'DISTRIBUTION runs for',real(t2-t1),'seconds.'
	write(*,*) "============================================================="
endif
if (diff > tol_dist) then
    exit_flag = .false.
    write(*,*) 'WARNING in get_distribution: MU did not converge!'
endif

!Write Distribution mu
! call writedim(distrib%muy,trim(savedir)//trim('muy.txt'))
! call writedim(distrib%mur,trim(savedir)//trim('mur.txt'))

end subroutine get_distribution
!=======================================================================!

subroutine sub_distrib_map(distrib,distrib_m,occpol)
implicit none
! Inputs:
type(distribfun_m), intent(in) :: distrib_m
real(8), intent(in) :: occpol(:,:,:,:,:)
! Outputs:
type(distribfun), intent(inout) :: distrib
! Local variables
integer :: io,iz_min

! Initialize distrib%muy
distrib%muy = 0d0
if (.not. allocated(distrib%muy)) then
	write(*,*) "sub_distrib_map: distrib%muy not allocated!"
	stop
elseif (.not. allocated(distrib_m%muy)) then
	write(*,*) "sub_distrib_map: distrib_m%muy not allocated!"
	stop
endif 

do io = 1,no
	do iz_min=1,nz
		distrib%muy(:,:,:,io) = distrib%muy(:,:,:,io)+occpol(:,:,:,iz_min,io)*distrib_m%muy(:,:,:,iz_min)
	enddo
enddo 

distrib%mur = distrib_m%mur

end subroutine sub_distrib_map
!=======================================================================!

subroutine get_KN(KN_supply,distrib_m,pol)

implicit none
! Inputs/Outputs
real(8),intent(out)            :: KN_supply
type(distribfun_m), intent(in) :: distrib_m
type(polfun), intent(in)       :: pol
! Local variables
type(distribfun) :: distrib
real(8) :: AA, K_entre_tot,labsup,K_C_new,eps_val,n_entre,N_C_new
integer :: ia,ie,it,io

! Convert distrib_m to distrib
allocate(distrib%muy(na,neps,ntheta,no),distrib%mur(na))
call sub_distrib_map(distrib,distrib_m,pol%occpol)

AA = 0.0d0 ! Aggregate savings
!CC = 0.0d0 ! Aggregate consumption

do ia = 1,na
	AA = AA + a_grid(ia)*(sum(distrib%muy(ia,:,:,:)) + distrib%mur(ia))
enddo

! Total capital demand by entrepreneurial sector k(a,theta,o)
K_entre_tot = 0.0d0
do it = 1,ntheta
do ie = 1,neps
do ia = 1,na
	do io = 2,4 !Sum over entrepreneurs only
		K_entre_tot = K_entre_tot + pol%kpol(ia,it,io)*distrib%muy(ia,ie,it,io)
	enddo ! io
enddo ! ia
enddo ! ie
enddo ! it

K_C_new = AA-K_entre_tot

!--------------------------------------------------------------------------------!
! Compute N_C implied by the market clearing condition
! Market clearing for labor:
! N_C + tot labor demand by entre = Aggregate labor supply
! ==> N_C_new = Aggregate labor supply - (tot labor demand by entre)
!--------------------------------------------------------------------------------!
 
! Aggregate (effective) labor supply by workers
labsup = 0.0d0
do it = 1,ntheta
    do ie = 1,neps
    	eps_val = eps_grid(ie)
        do ia = 1,na
            ! sum only if the guy is a worker
            io = 1
            labsup = labsup + eps_val*pol%lpoly_w(ia,ie,it)*distrib%muy(ia,ie,it,io)
        enddo ! ia
    enddo ! ie
enddo ! it

! Total hired labor of entrepreneurs
n_entre = 0.0d0
do it = 1,ntheta
    do ie = 1,neps
        do ia = 1,na
            do io = 2,4
                n_entre = n_entre + pol%npol(ia,it,io)*distrib%muy(ia,ie,it,io)
            enddo !io
        enddo !ia
    enddo !eps
enddo !theta
!
!!N_C_new = labsup - n_entre ! what if N_C_new<0 ??
N_C_new = max(labsup - n_entre, 1.0d-3) 
!write(*,*) "N_C_new: ",N_C_new
!! capital-labor ratio implied by market clearing
KN_supply = K_C_new/N_C_new

end subroutine get_KN
!=======================================================================!

subroutine get_aggregates(agg,prices,distrib_m,pol)

use mod_functions

!--------------------------- LEGEND --------------------------------------%
! DESCRIPTION:
! Given policy function and distribution, compute aggregate moments
! needed for the general equilibrium loop. 
! ALL OTHER MOMENTS ARE COMPUTED IN MOD_TARGETS.
!--------------------------------------------------------------------------%

implicit none
!Declare inputs and outputs:
type(modelAgg),intent(out)    :: agg
type(modelPrices), intent(in) :: prices   ! prices
type(distribfun_m),intent(in) :: distrib_m
type(polfun),intent(in)       :: pol

! Local variables
type(distribfun) :: distrib
integer :: ia,ie,it,io
real(8) :: k_val,nbar_val,y_val,P_occ,wage,rbar,yw,profit_val,n_val,pen,phi_val,yw_deduct,yw_r
real(8) :: eps_val,theta_val,a_val,ss_cap,y_H,labcost_ES,labcost_EC,VA_ES,VA_EC
real(8) :: muy_occ(no),y_vec(no)

! Convert distrib_m to distrib
allocate(distrib%muy(na,neps,ntheta,no),distrib%mur(na))
call sub_distrib_map(distrib,distrib_m,pol%occpol)

wage = prices%w
rbar = prices%r
pen  = prices%pen
y_H  = prices%y_H
ss_cap = prices%ss_cap

! - Aggregate savings
agg%A = 0.0d0 
do ia = 1,na
	agg%A = agg%A + a_grid(ia)*(sum(distrib%muy(ia,:,:,:)) + distrib%mur(ia))
enddo

! - Capital in the entrepreneurial sector
agg%K_entre = 0.0d0
agg%K_share = 0.0d0
do it = 1,ntheta
do ie = 1,neps
do ia = 1,na
	do io = 2,no !Sum over entrepreneurs only
        agg%K_share(io) = agg%K_share(io) + pol%kpol(ia,it,io)*distrib%muy(ia,ie,it,io)
		agg%K_entre = agg%K_entre + pol%kpol(ia,it,io)*distrib%muy(ia,ie,it,io)
	enddo ! io
enddo ! ia
enddo ! ie
enddo ! it

! - Capital in the corporate sector
agg%K_C = agg%A-agg%K_entre
agg%K_share(1) = agg%K_C
! Scale K_share to get shares
agg%K_share = agg%K_share/sum(agg%K_share)

!--------------------------------------------------------------------------------!
! Compute N_C implied by the market clearing condition
! Market clearing for labor:
! N_C + tot labor demand by entre = Aggregate labor supply
! ==> N_C_new = Aggregate labor supply - (tot labor demand by entre)
!--------------------------------------------------------------------------------!
 
! - Aggregate (effective) labor supply by workers
agg%labsup = 0.0d0
do it = 1,ntheta
    do ie = 1,neps
    	eps_val = eps_grid(ie)
        do ia = 1,na
            ! sum only if the guy is a worker
            io = 1
            agg%labsup = agg%labsup + eps_val*pol%lpoly_w(ia,ie,it)*distrib%muy(ia,ie,it,io)
        enddo ! ia
    enddo ! ie
enddo ! it

! - Total hired labor of entrepreneurs and by each legal form
agg%N_entre      = 0.0d0
agg%N_share      = 0.0d0

do it = 1,ntheta
    do ie = 1,neps
        do ia = 1,na
            do io = 2,4
                agg%N_share(io) = agg%N_share(io) + pol%npol(ia,it,io)*distrib%muy(ia,ie,it,io)
                agg%N_entre = agg%N_entre + pol%npol(ia,it,io)*distrib%muy(ia,ie,it,io)
            enddo !io
        enddo !ia
    enddo !eps
enddo !theta

! - Labor in the corporate sector
! Need to make sure that L_C>0!
agg%N_C = max(agg%labsup - agg%N_entre, 1.0d-3) 
agg%N_share(1) = agg%N_C
! Scale agg%N_share to get shares
agg%N_share = agg%N_share/sum(agg%N_share)

!--------------------------------------------------------------------------------!
! Compute K_C implied by the market clearing condition
! Market clearing for capital:
! K_C + (K_EP+K_ES+K_EC) = Aggregate savings
! ==> K_C_new  = Aggregate savings - (K_EP+K_ES+K_EC)
!--------------------------------------------------------------------------------!
!! capital-labor ratio implied by market clearing
agg%KN = agg%K_C/agg%N_C

if (agg%KN<0.0d0) then
   write(*,*) "K-N ratio is negative"
   !pause
!stop
endif

! Excess demand
agg%ED = prices%KN-agg%KN

! - Output in the entrepreneurial sector
agg%Y_entre = 0.0d0
agg%Y_lfo = 0.0d0
do it = 1,ntheta
    theta_val = theta_grid(it)
	do ie = 1,neps
		do ia = 1,na
            do io = 2,4
                k_val       = pol%kpol(ia,it,io)          ! capital
                nbar_val    = pol%npol(ia,it,io) + le     ! total labor used 
                y_val       = theta_val*prodfun(k_val,nbar_val)
                agg%Y_lfo(io) = agg%Y_lfo(io) + y_val*distrib%muy(ia,ie,it,io)
		    	agg%Y_entre = agg%Y_entre+y_val*distrib%muy(ia,ie,it,io)
            enddo ! io
		enddo !ia
	enddo !ie
enddo !it

! - Output in the corporate sector
agg%Y_C = cobb_douglas(agg%K_C,agg%N_C)

! - Aggregate output in both sectors
agg%Y = agg%Y_entre + agg%Y_C



! - Capital-to-output ratio
agg%KY = agg%A/agg%Y 

! Consumption of young agents
agg%C_y = 0.0d0
do io = 1,4
	agg%C_y = agg%C_y + sum(pol%cpoly(:,:,:,io)*distrib%muy(:,:,:,io))
enddo

! - Consumption of retirees
agg%C_r = sum(pol%cpolr*distrib%mur)

! - Aggregate Consumption
agg%C = agg%C_y + agg%C_r 

! - Consumption of entrepreneurs
agg%C_entre = 0.0d0
do io = 2,4
	agg%C_entre = agg%C_entre + sum(pol%cpoly(:,:,:,io)*distrib%muy(:,:,:,io))
enddo

! - Consumption of workers
io = 1
agg%C_worker  = sum(pol%cpoly(:,:,:,io)*distrib%muy(:,:,:,io))


!-------------------------------------------------------------
! Calculate tax revenues, total operating costs, and total tax avoidance costs
!-------------------------------------------------------------
! TODO: tax revenue paid by the top 1% earners


! Aggregate taxes collected by the government
agg%taxes_ss   = 0.0d0 !social security
agg%taxes_inc  = 0.0d0 !income taxes
agg%taxes_corp = 0.0d0 !corporate
agg%taxes_div  = 0.0d0 !dividend
! Tax contribution by each io:
agg%taxes_occ  = 0.0d0 ! dim: 4
agg%taxes_r    = 0.0d0

! Total operating cost
agg%op_cost_tot   = 0.0d0
agg%avoidance_tot = 0.0d0

! Auxiliar terms to update hsv_0 (lambda_i)
agg%calY1 = 0.0d0
agg%calY2 = 0.0d0

! - Total labor cost = wages paid to workers + labor income of entrepreneurs
labcost_ES = 0.0d0
labcost_EC = 0.0d0
! - Value added (profit + labor cost of hired labor)
VA_ES = 0.0d0
VA_EC = 0.0d0

! Loop over state space: (a,eps,theta) x (young or old)
do ia = 1,na
	do it = 1,ntheta
    	do ie = 1,neps
            eps_val    = eps_grid(ie)
            theta_val  = theta_grid(it)
            a_val      = a_grid(ia)
            
            ! - Workers
            	io = 1
                ! Workers pay income and ss taxes
                !   - Social security taxes
                agg%taxes_ss =  agg%taxes_ss+fun_tax_pay(wage*eps_val*pol%lpoly_w(ia,ie,it),ss_cap) &
                	*distrib%muy(ia,ie,it,io)
                !   - Income taxes
                yw = wage*eps_val*pol%lpoly_w(ia,ie,it) - fun_tax_pay(wage*eps_val*pol%lpoly_w(ia,ie,it),ss_cap) + rbar*a_val
                agg%taxes_inc = agg%taxes_inc+ fun_tax_inc(yw,y_H)*distrib%muy(ia,ie,it,io)
                ! Total taxes paid by workers
                agg%taxes_occ(io) = agg%taxes_occ(io) &
                	+ fun_tax_pay(wage*eps_val*pol%lpoly_w(ia,ie,it),ss_cap)*distrib%muy(ia,ie,it,io) &
                	+ fun_tax_inc(yw,y_H)*distrib%muy(ia,ie,it,io)
            	! Auxiliar terms to update hsv_0 (lambda_i)
            	agg%calY1 = agg%calY1 + (min(yw,y_H) + tau_h*max(yw-y_H,0.0d0))*distrib%muy(ia,ie,it,io) 
            	agg%calY2 = agg%calY2 + min(yw,y_H)**(1.0d0-hsv_1)*distrib%muy(ia,ie,it,io) 

            ! - Sole Proprietors
            	io = 2
                ! EP: pay income and ss taxes 
                k_val   = pol%kpol(ia,it,io)
                n_val   = pol%npol(ia,it,io) 
                profit_val = fun_profit(theta_val,k_val,le,n_val,rbar,wage)
                yw = profit_val - fun_tax_pay(profit_val,ss_cap) + rbar*a_val
                !   - Social security taxes
                agg%taxes_ss = agg%taxes_ss + fun_tax_pay(profit_val,ss_cap)*distrib%muy(ia,ie,it,io)
                !   - Income taxes
                agg%taxes_inc = agg%taxes_inc + fun_tax_inc(yw,y_H)*distrib%muy(ia,ie,it,io)
                ! Update taxes paid by EP
                agg%taxes_occ(io) = agg%taxes_occ(io) + fun_tax_pay(profit_val,ss_cap)*distrib%muy(ia,ie,it,io) &
                	+ fun_tax_inc(yw,y_H)*distrib%muy(ia,ie,it,io)
                ! Auxiliar terms to update hsv_0 (lambda_i)
            	agg%calY1 = agg%calY1 + (min(yw,y_H) + tau_h*max(yw-y_H,0.0d0))*distrib%muy(ia,ie,it,io) 
            	agg%calY2 = agg%calY2 + min(yw,y_H)**(1.0d0-hsv_1)*distrib%muy(ia,ie,it,io) 
            	
            ! - S-Corps
            	io = 3
                ! ES: pay income and ss taxes 
                k_val   = pol%kpol(ia,it,io)
                n_val   = pol%npol(ia,it,io) 
                profit_val = fun_profit(theta_val,k_val,le,n_val,rbar,wage)
               	phi_val  = pol%phipol(ia,it,io)
                yw = profit_val -fun_tax_pay(phi_val*profit_val,ss_cap) + rbar*a_val
                !   - Social security taxes
                agg%taxes_ss = agg%taxes_ss + fun_tax_pay(phi_val*profit_val,ss_cap)*distrib%muy(ia,ie,it,io)
                !   - Income taxes
                agg%taxes_inc = agg%taxes_inc + fun_tax_inc(yw-cost_evasion_es(1.0d0-phi_val)-op_cost_es,y_H) &
                						*distrib%muy(ia,ie,it,io)
                 ! Update taxes paid by ES
                agg%taxes_occ(io) = agg%taxes_occ(io) &
                			+ fun_tax_pay(phi_val*profit_val,ss_cap)*distrib%muy(ia,ie,it,io) &
                			+ fun_tax_inc(yw-cost_evasion_es(1d0-phi_val)-op_cost_es,y_H) &
                				*distrib%muy(ia,ie,it,io)
                ! Update operating cost
                agg%op_cost_tot = agg%op_cost_tot + op_cost_es*distrib%muy(ia,ie,it,io)
                ! Update avoidance cost
                agg%avoidance_tot = agg%avoidance_tot + cost_evasion_es(1.0d0-phi_val)*distrib%muy(ia,ie,it,io)
                ! Auxiliar terms to update hsv_0 (lambda_i)
                yw_deduct = max(0.0d0,yw-cost_evasion_es(1d0-phi_val)-op_cost_es)
            	agg%calY1 = agg%calY1 + (min(yw_deduct,y_H) + tau_h*max(yw_deduct-y_H,0.0d0))*distrib%muy(ia,ie,it,io) 
            	agg%calY2 = agg%calY2 + min(yw_deduct,y_H)**(1.0d0-hsv_1)*distrib%muy(ia,ie,it,io) 
                ! Labor income and value added
                labcost_ES = labcost_ES + (phi_val*profit_val + n_val*wage) *distrib%muy(ia,ie,it,io)
                VA_ES = VA_ES + theta_val*prodfun(k_val,n_val)*distrib%muy(ia,ie,it,io)
            ! - C-Corps
            	io = 4
                ! EC: pay ss, income, coporate and dividend taxes
                k_val   = pol%kpol(ia,it,io)
                n_val   = pol%npol(ia,it,io) 
                profit_val = fun_profit(theta_val,k_val,le,n_val,rbar,wage)
               	phi_val  = pol%phipol(ia,it,io)
                yw =phi_val*profit_val-fun_tax_pay(phi_val*profit_val,ss_cap) + rbar*a_val                
                !   - Social security taxes
                agg%taxes_ss = agg%taxes_ss + fun_tax_pay(phi_val*profit_val,ss_cap)*distrib%muy(ia,ie,it,io)
                !   - Income taxes
                agg%taxes_inc = agg%taxes_inc + fun_tax_inc(yw-cost_evasion_ec(phi_val)-op_cost_ec,y_H) &
                				*distrib%muy(ia,ie,it,io)
                !   - Coporate and dividend tax
                agg%taxes_corp = agg%taxes_corp + tau_c*(1.0d0-phi_val)*profit_val*distrib%muy(ia,ie,it,io)
                agg%taxes_div  = agg%taxes_div + tau_d*(1.0d0-tau_c)*(1.0d0-phi_val)*profit_val*distrib%muy(ia,ie,it,io)
                ! Update taxes paid by EC
                agg%taxes_occ(io) = agg%taxes_occ(io) + fun_tax_pay(phi_val*profit_val,ss_cap)*distrib%muy(ia,ie,it,io) &
                                              + fun_tax_inc(yw-cost_evasion_ec(phi_val)-op_cost_ec,y_H) &
                									*distrib%muy(ia,ie,it,io) &
                                              + tau_c*(1.0d0-phi_val)*profit_val*distrib%muy(ia,ie,it,io) &
                                              + tau_d*(1.0d0-tau_c)*(1.0d0-phi_val)*profit_val*distrib%muy(ia,ie,it,io)
                ! Update operating cost
                agg%op_cost_tot = agg%op_cost_tot + op_cost_ec*distrib%muy(ia,ie,it,io)
                ! Update avoidance cost
                agg%avoidance_tot = agg%avoidance_tot + cost_evasion_ec(phi_val)*distrib%muy(ia,ie,it,io)
                 ! Auxiliar terms to update hsv_0 (lambda_i)
                yw_deduct = max(0.0d0,yw-cost_evasion_ec(phi_val)-op_cost_ec)
            	agg%calY1 = agg%calY1 + (min(yw_deduct,y_H) + tau_h*max(yw_deduct-y_H,0.0d0))*distrib%muy(ia,ie,it,io) 
            	agg%calY2 = agg%calY2 + min(yw_deduct,y_H)**(1.0d0-hsv_1)*distrib%muy(ia,ie,it,io) 
                ! Labor income and value added
                labcost_EC = labcost_EC + (phi_val*profit_val + n_val*wage) *distrib%muy(ia,ie,it,io)
                VA_EC = VA_EC + theta_val*prodfun(k_val,n_val)*distrib%muy(ia,ie,it,io)

        enddo ! ie
    enddo ! it
    ! Tax paid by the old
	! Retirees pay income taxes
	agg%taxes_inc = agg%taxes_inc + fun_tax_inc(pen+rbar*a_val,y_H)*distrib%mur(ia)
	! Update taxes paid by retirees
	agg%taxes_r   = agg%taxes_r + fun_tax_inc(pen+rbar*a_val,y_H)*distrib%mur(ia)
	! Auxiliar terms to update hsv_0 (lambda_i)
	yw_r = pen+rbar*a_val
    agg%calY1 = agg%calY1 + (min(yw_r,y_H) + tau_h*max(yw_r-y_H,0.0d0))*distrib%mur(ia)
    agg%calY2 = agg%calY2 + min(yw_r,y_H)**(1.0d0-hsv_1)*distrib%mur(ia)
            	
enddo ! ia

! - Labor share = total labor cost / total output

agg%labshare_ES = labcost_ES/VA_ES
agg%labshare_EC = labcost_EC/VA_EC
agg%labshare_corp = (labcost_ES+labcost_EC)/(VA_ES+VA_EC)
! - Labor share in the corporate sector (including entrepreneurial and non-entrepreneurial)
agg%labshare_allcorp = (wage*agg%N_C+labcost_EC + labcost_ES)/(agg%Y_C + VA_EC + VA_ES)


! - Total tax revenue
agg%taxes_tot = agg%taxes_ss+agg%taxes_inc+agg%taxes_corp+agg%taxes_div
 
! - Government spending
!agg%G = agg%Y * G_frac
agg%G = agg%taxes_inc + agg%taxes_corp + agg%taxes_div




! - Aggregate pension benefits paid by the government
agg%pen_tot   = pen*sum(distrib%mur)       

! - Government budget constraint residual
agg%res_gov = (agg%taxes_inc + agg%taxes_corp + agg%taxes_div)/agg%Y - G_frac ! In eq. this should be zero

agg%bal_gov =  (agg%taxes_inc + agg%taxes_corp + agg%taxes_div) - G_bench

! - Residual of social security budget constrain
agg%bal_ss =agg%taxes_ss - agg%pen_tot

! - Average theta, average k, and average k-l ratio by LFO and for all entre
! avek, avetheta, etc are vectors with 4 elements where
! 1 = all entre
! 2 = sole prop
! 3 = S-corp
! 4 = C-corp
agg%avetheta   = 0.0d0
agg%avek       = 0.0d0 
agg%avey       = 0.0d0 
agg%aveklratio = 0.0d0
agg%avekyratio = 0.0d0
agg%avel       = 0.0d0

! theta_val*prodfun(k_val,nbar_val)

muy_occ  = 0.0d0
do it = 1, ntheta
do ie = 1, neps
do ia = 1, na
	
	agg%avetheta(1) = agg%avetheta(1) +  theta_grid(it)*sum(distrib%muy(ia,ie,it,2:4))
	agg%avek(1) = agg%avek(1) +  sum(pol%kpol(ia,it,2:4)*distrib%muy(ia,ie,it,2:4))
	agg%avel(1) = agg%avel(1) +  sum((pol%npol(ia,it,2:4)+le)*distrib%muy(ia,ie,it,2:4))

	muy_occ(1) = muy_occ(1)  + sum(distrib%muy(ia,ie,it,2:4))
	do io = 2,4
		k_val     = pol%kpol(ia,it,io)
		nbar_val  = pol%npol(ia,it,io)+le
		y_vec(io) = theta_grid(it)*prodfun(k_val,nbar_val)
	
		agg%avetheta(io) = agg%avetheta(io) +  theta_grid(it)*distrib%muy(ia,ie,it,io)
		agg%avek(io) = agg%avek(io) +  pol%kpol(ia,it,io)*distrib%muy(ia,ie,it,io)
		agg%avel(io) = agg%avel(io) +  (pol%npol(ia,it,io)+le)*distrib%muy(ia,ie,it,io)
		agg%avey(io) = agg%avey(io) +  y_vec(io)*distrib%muy(ia,ie,it,io)
		muy_occ(io) = muy_occ(io)  + distrib%muy(ia,ie,it,io)
		
	enddo
	agg%avey(1) = agg%avey(1) +  sum(y_vec(2:4)*distrib%muy(ia,ie,it,2:4))
	
enddo
enddo
enddo
agg%avetheta = agg%avetheta/muy_occ
agg%avek = agg%avek/muy_occ
agg%avel = agg%avel/muy_occ
agg%avey = agg%avey/muy_occ
agg%aveklratio = agg%avek/agg%avel
agg%avekyratio = agg%avek/agg%avey

! Check goods market clearing condition
agg%ED_goods = agg%C + delta*agg%A + agg%G + agg%op_cost_tot + agg%avoidance_tot -agg%Y 
!write(*,*) "ED_goods = ",agg%ED_goods

end subroutine get_aggregates
!=======================================================================!
    
end module mod_distribution