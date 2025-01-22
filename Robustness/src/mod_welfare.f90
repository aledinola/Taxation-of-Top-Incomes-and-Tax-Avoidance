module mod_welfare
    use omp_lib
    use mod_globals
    use mod_utilities
    use mod_functions, only: util_cons,util
    use mod_distribution, only: precompute_distrib
    use mod_numerical, only: myerror, outerprod2,quantili,linspace

    implicit none

! Subroutine to compute V_cons (value from consumption in the initial steady-state
! In this sub, load policy functions from folder ss_init
! Then call a one-step update such that, given an initial V_cons you get an updated V_cons

contains
!=======================================================================!

function compute_cev_agg(val_init,val_final,valc_init,mu_init) result(cev_agg)
! PURPOSE:
	! Compute the aggregate CEV \omega(s,\tau).
	implicit none  
	! Inputs:
	type(valfun), intent(in)  :: val_init,val_final,valc_init
	type(distribfun_m), intent(in) :: mu_init
	! Outputs:
	real(8) :: cev_agg
	! Local variables
	! Function result:
	! cev%Vy has dim (a,eps,theta)
	! cev%Vr has dim (a)
    type(valfun) :: cev
    
	! First compute conditional CEV 
	cev = compute_cev(val_init,val_final,valc_init) 
	! Then, integrate using distribution in the benchmark economy
	cev_agg = 0.0d0

	cev_agg = sum(cev%Vy * mu_init%muy) +  sum(cev%Vr * mu_init%mur) 
	
end function compute_cev_agg
!=======================================================================!

subroutine decomp_cev(cev_agg,cev_aggcomp_agg,cev_distcomp_agg, agg_C_init,agg_labsup_init,pol_init,val_init,valc_init,mu_init,&
	agg_C_final,agg_labsup_final,pol_final,val_final) 
! PURPOSE:
	! Compute the decomposition of CEV into aggregate and distributional components
	implicit none  
	! Inputs:
	real(8), intent(in) :: agg_C_init,agg_labsup_init,agg_C_final,agg_labsup_final
	type(valfun), intent(in) :: val_init,val_final,valc_init
	type(polfun), intent(in) :: pol_init,pol_final
	type(distribfun_m), intent(in) :: mu_init
	! Outputs:
	! Note: Type valfun has fields Vy (4D array) and Vr (1D array)
	real(8),intent(out) :: cev_agg           ! CEV, total
	real(8),intent(out) :: cev_aggcomp_agg   ! CEV, aggregate component
	real(8),intent(out) :: cev_distcomp_agg  ! CEV, distributional component
	
	! Local variables
	type(valfun) :: valhat
	type(valfun) :: cev           ! CEV, total
	type(valfun) :: cev_aggcomp   ! CEV, aggregate component
	!type(valfun) :: cev_distcomp  ! CEV, distributional component

	! First compute conditional CEV, total
	cev = compute_cev(val_init,val_final,valc_init) 
	
	! CEV aggregate component
	call compute_value_hat(valhat, agg_C_init,agg_labsup_init,pol_init,agg_C_final,agg_labsup_final,pol_final,val_final)
    cev_aggcomp = compute_cev(val_init,valhat,valc_init) 

	! CEV, distributional component
	!cev_distcomp%Vy = (cev%Vy - cev_aggcomp%Vy)/(1d0+cev_aggcomp%Vy)
	!cev_distcomp%Vr = (cev%Vr - cev_aggcomp%Vr)/(1d0+cev_aggcomp%Vr)

	! Then, integrate using distribution in the benchmark economy
	cev_agg          = sum(cev%Vy * mu_init%muy) +  sum(cev%Vr * mu_init%mur)
	cev_aggcomp_agg  = sum(cev_aggcomp%Vy * mu_init%muy) +  sum(cev_aggcomp%Vr * mu_init%mur)
	cev_distcomp_agg = (cev_agg - cev_aggcomp_agg)/(1d0+cev_aggcomp_agg)

end subroutine decomp_cev
!=======================================================================!

subroutine decomp_cev_tran(cev_agg,cev_aggcomp_agg,cev_distcomp_agg, agg_C_init,agg_labsup_init, &
	pol_init,val_init,valc_init,mu_init,val_t0,agg_C_final,agg_labsup_final,pol_final,val_final, &
	pol_path,agg_path) 
! PURPOSE:
	! Compute the decomposition of CEV into aggregate and distributional components, considering transition
	implicit none  
	! Inputs:
	real(8), intent(in) :: agg_C_init,agg_labsup_init,agg_C_final,agg_labsup_final
	type(valfun), intent(in) :: val_init,val_final,valc_init,val_t0
	type(polfun), intent(in) :: pol_init,pol_final
	type(polfun), intent(in) :: pol_path(:)
	type(modelAgg), intent(in) :: agg_path(:)
	type(distribfun_m), intent(in) :: mu_init
	! Outputs:
	! Note: Type valfun has fields Vy (4D array) and Vr (1D array)
	real(8),intent(out) :: cev_agg           ! CEV, total
	real(8),intent(out) :: cev_aggcomp_agg   ! CEV, aggregate component
	real(8),intent(out) :: cev_distcomp_agg  ! CEV, distributional component
	
	! Local variables
	type(valfun) :: valhat_t0
	type(valfun) :: cev           ! CEV, total
	type(valfun) :: cev_aggcomp   ! CEV, aggregate component
	!type(valfun) :: cev_distcomp  ! CEV, distributional component

	! First compute conditional CEV, total
	cev = compute_cev(val_init,val_t0,valc_init) 
	
	! CEV aggregate component
	call compute_value_hat_tran(valhat_t0, agg_C_init,agg_labsup_init,agg_C_final,agg_labsup_final, &
			pol_init,pol_final, val_final, pol_path,agg_path)
    cev_aggcomp = compute_cev(val_init,valhat_t0,valc_init) 

	! CEV, distributional component
	!cev_distcomp%Vy = (cev%Vy - cev_aggcomp%Vy)/(1d0+cev_aggcomp%Vy)
	!cev_distcomp%Vr = (cev%Vr - cev_aggcomp%Vr)/(1d0+cev_aggcomp%Vr)

	! Then, integrate using distribution in the benchmark economy
	cev_agg          = sum(cev%Vy * mu_init%muy) +  sum(cev%Vr * mu_init%mur)
	cev_aggcomp_agg  = sum(cev_aggcomp%Vy * mu_init%muy) +  sum(cev_aggcomp%Vr * mu_init%mur)
	cev_distcomp_agg = (cev_agg - cev_aggcomp_agg)/(1d0+cev_aggcomp_agg)

end subroutine decomp_cev_tran
!=======================================================================!

subroutine compute_cev_groups(cev_q,cev_qo,cev_z, val_init,val_final,valc_init,mu_init,occpol,nq) 
! PURPOSE:
	! Compute the CEV by groups (wealth,occupational choice, etc.)
	! We focus only on young agents
	implicit none  
	! Input arguments:
	type(valfun), intent(in)  :: val_init,val_final,valc_init
	type(distribfun_m), intent(in) :: mu_init !muy(a,eps,theta,z_min),mur(a)
	real(8),intent(in) :: occpol(:,:,:,:,:)
	integer, intent(in) :: nq
	! Output arguments:
	real(8), intent(out) :: cev_q(:)    ! CEV by wealth percentiles
	real(8), intent(out) :: cev_qo(:,:) ! CEV by wealth percentiles and occupation (worker = 1, entre = 2)
	real(8), intent(out) :: cev_z(:)    ! CEV by occupation and lfo (worker = 1, sole-prop = 2, S-corp = 3, C-corp = 4)
	! Local variables
	integer :: ia,iq,iz
	real(8) :: a_val
	real(8) :: q_vec(nq),w_q(nq),mass_q(nq),mass_qo(nq,2),mass_z(4)
	type(valfun) :: cev
	real(8), allocatable :: muy_a(:),entre_prob(:,:,:)
	
	write(*,*) "compute_cev_groups starts..."
	
	! First, compute conditional CEV for each point on the state space.
	cev = compute_cev(val_init,val_final,valc_init)
	
	! - Second, compute the CEV by wealth deciles and by occupation (work vs entre)
	! Compute deciles for wealth
	q_vec = linspace(1d0/real(nq,8),1d0,nq)
	
	allocate(muy_a(na))
	do ia = 1,na
		muy_a(ia) = sum(mu_init%muy(ia,:,:,:)) !Marginal distrib of assets for young
	enddo
	
	if (size(a_grid)/=size(muy_a) ) then
		call myerror("a_grid does not match muy_a!")
	endif
	w_q = quantili(a_grid,muy_a,q_vec)
	
	cev_q  = 0d0 !dim: (nq)
	cev_qo = 0.0d0 !dim: (nq,2)
	mass_q = 0d0 !dim: (nq)
	mass_qo = 0.0d0 !dim: (nq,2)
	! First interval
	do ia = 1,na
		a_val = a_grid(ia)
		if (a_val<=w_q(1) ) then
			! all young
			cev_q(1)  = cev_q(1)+ sum(cev%Vy(ia,:,:,:)*mu_init%muy(ia,:,:,:))
			mass_q(1) = mass_q(1)+sum(mu_init%muy(ia,:,:,:))
			! worker
			cev_qo(1,1)  = cev_qo(1,1) + sum(cev%Vy(ia,:,:,:)*occpol(ia,:,:,:,1)*mu_init%muy(ia,:,:,:))
			mass_qo(1,1) = mass_qo(1,1) +sum(occpol(ia,:,:,:,1)*mu_init%muy(ia,:,:,:))
			! entre
			entre_prob = occpol(ia,:,:,:,2) + occpol(ia,:,:,:,3) + occpol(ia,:,:,:,4) 
			cev_qo(1,2)  = cev_qo(1,2) + sum(cev%Vy(ia,:,:,:)*entre_prob*mu_init%muy(ia,:,:,:))
			mass_qo(1,2) = mass_qo(1,2) +sum(entre_prob*mu_init%muy(ia,:,:,:))	 
		else
			exit !since a_grid is strictly increasing
		endif 
	enddo
	! All other intervals
	do ia = 1,na
		a_val = a_grid(ia)
		do iq = 2,nq
			if (a_val>w_q(iq-1) .and. a_val<=w_q(iq)) then
				! all young
				cev_q(iq)  = cev_q(iq)+ sum(cev%Vy(ia,:,:,:)*mu_init%muy(ia,:,:,:))
				mass_q(iq) = mass_q(iq)+sum(mu_init%muy(ia,:,:,:))
				! worker
				cev_qo(iq,1)  = cev_qo(iq,1) + sum(cev%Vy(ia,:,:,:)*occpol(ia,:,:,:,1)*mu_init%muy(ia,:,:,:))
				mass_qo(iq,1) = mass_qo(iq,1) +sum(occpol(ia,:,:,:,1)*mu_init%muy(ia,:,:,:))
				! entre
				entre_prob = occpol(ia,:,:,:,2) + occpol(ia,:,:,:,3) + occpol(ia,:,:,:,4) 
				cev_qo(iq,2)  = cev_qo(iq,2) + sum(cev%Vy(ia,:,:,:)*entre_prob*mu_init%muy(ia,:,:,:))
				mass_qo(iq,2) = mass_qo(iq,2) +sum(entre_prob*mu_init%muy(ia,:,:,:))	
			endif 
		enddo
	enddo
	cev_q = cev_q/mass_q ! dim: (nq)
	cev_qo = cev_qo/mass_qo ! dim: (nq,2)
	
	! CEV by occupation and lfo
	cev_z=0.0d0 ! dim: (no)
	mass_z = 0.0d0 ! dim: (no)
	do iz = 1,no
		cev_z(iz)  = cev_z(iz) + sum(cev%Vy*occpol(:,:,:,:,iz)*mu_init%muy)
		mass_z(iz) = mass_z(iz) +sum(occpol(:,:,:,:,iz)*mu_init%muy)
	enddo
	cev_z = cev_z / mass_z
	! call disp("cev_q",cev_q)
	! call disp("cev_qo",cev_qo)
	! call disp("cev_z",cev_z)
	! pause
	
end subroutine compute_cev_groups
!=======================================================================!

function compute_cev(val_init,val_final,valc_init) result(cev)
	! PURPOSE:
	! Compute the conditional CEV \omega(s,\tau).
	implicit none    
	! Inputs:
	type(valfun), intent(in)  :: val_init,val_final,valc_init
	! Function result:
	! cev%Vy has dim (a,eps,theta)
	! cev%Vr has dim (a)
    type(valfun) :: cev
    ! Local variables:
    real(8) :: aux
    
    aux    = 1.0d0/(1.0d0-sigma1)
    
    ! Compute CEV separately for young and retired agents
    cev%Vy = ((val_final%Vy-val_init%Vy)/valc_init%Vy + 1.0d0)**aux - 1.0d0
    cev%Vr = ((val_final%Vr-val_init%Vr)/valc_init%Vr + 1.0d0)**aux - 1.0d0
    
end function compute_cev
!=======================================================================!
subroutine compute_value_cons(valc,dir)
    implicit none
    ! Input variables
    character(len=*),intent(in) :: dir
    ! Output variables
    type(valfun),intent(out) :: valc
    ! Local variables
    type(valfun)         :: valc0
    ! valc%Vy(a,eps,theta)
    ! valc%Vr(a)
    character(len=300)   :: pwd
    type(polfun)         :: pol
    integer              :: iter, istat
    real(8)              :: err, err_y, err_r
    
    real(8), allocatable:: omega(:,:,:,:),omegar(:),temp_ret(:,:)
    integer, allocatable :: aopt_ind(:,:,:,:),aoptr_ind(:)
    
    
! Body of compute_value_cons

if (verbose>=1) then
    call getcwd(pwd)
    write(*,*) "The current working directory is: ",trim(pwd)
    write(*,*) "We load policies from folder ", dir
    write(*,*) " "
endif


! Load policy functions from initial steady-state
write(*,*) "I am Inside compute_value_cons"
write(*,*) "Loading policy functions from initial steady-state"

!pause

allocate(pol%cpoly(na,neps,ntheta,no),pol%cpolr(na),pol%apoly(na,neps,ntheta,no),pol%apolr(na),&
    pol%occpol(na,neps,ntheta,nz,no),stat=istat)
if (istat/=0) then
	call myerror("compute_value_cons: Allocation failed!")
endif

call readdim(pol%apoly, trim(dir)//trim('apoly.txt'))  !(a,eps,theta,o), o=W,EP,ES,EC
call readdim(pol%cpoly, trim(dir)//trim('cpoly.txt'))  !(a,eps,theta,o), o=W,EP,ES,EC
call readdim(pol%apolr, trim(dir)//trim('apolr.txt'))  !(a), retired
call readdim(pol%cpolr, trim(dir)//trim('cpolr.txt'))  !(a), retired
call readdim(pol%occpol,trim(dir)//trim('occpol.txt')) !(a,eps,theta,o)

write(*,*) "Policy functions read!"

! Precompute some stuff
allocate(omega(na,neps,ntheta,no),aopt_ind(na,neps,ntheta,no),omegar(na),aoptr_ind(na),temp_ret(neps,ntheta))

write(*,*) "Calling precompute_distrib.."

call  precompute_distrib(aopt_ind,omega,aoptr_ind,omegar, &
	temp_ret,pol)  


allocate(valc0%Vy(na,neps,ntheta,nz),valc0%Vr(na))

valc0%Vy = 0d0
valc0%Vr = 0d0
iter     = 1
err      = 10d0

write(*,*) "Start fixed point iteration to find valc.."

do while (iter<=maxitValc .and. err>tolValc)
    
    ! get updated valc (output)
	write(*,*) "before value_cons_one_step"
    call  value_cons_one_step(valc,valc0,aopt_ind,omega,aoptr_ind,omegar, &
	temp_ret,pol)  
	write(*,*) "after value_cons_one_step"
    
    ! Compute distance
	err_y = maxval(abs(valc%Vy - valc0%Vy))/max(maxval(abs(valc0%Vy)),1d0)
	err_r = maxval(abs(valc%Vr - valc0%Vr))/max(maxval(abs(valc0%Vr)),1d0)
    err = max(err_y,err_r)
    
    if (verbose>=1) then
		write(*,*) "======================================"
		write(*,*) "iter = ", iter
		write(*,*) "err_y = ", err_y
		write(*,*) "err_r = ", err_r
		write(*,*) "err   = ", err
		write(*,*) "======================================"
	endif
    
    ! Increase counter and update
    iter  = iter+1
    valc0 = valc
    
enddo

! save resulting valc to external file in folder 
call writedim(valc%Vy, trim(dir)//trim('valc_Vy.txt'))
call writedim(valc%Vr, trim(dir)//trim('valc_Vr.txt'))

end subroutine compute_value_cons

!=======================================================================!
subroutine value_cons_one_step(valc_update,valc,aopt_ind,omega,aoptr_ind,omegar, &
	temp_ret,pol)  
implicit none
! Inputs and outputs
type(valfun),intent(out) :: valc_update
type(valfun),intent(in) :: valc
type(polfun), intent(in) :: pol   
real(8), intent(in):: omega(:,:,:,:),omegar(:),temp_ret(:,:)
integer, intent(in) :: aopt_ind(:,:,:,:),aoptr_ind(:)
! Local variables
integer :: it,ie,ia,io,iap,iz_min
real(8) :: util_vec(4),io_prob,util_val,omega_val,EV_val
real(8),allocatable :: prob_mat(:,:,:),prob_mat_r(:),temp(:,:)
                       !valc_Vy(:,:,:),valc_Vr(:)
	
	valc_update = valc
	valc_update%Vy = 0.0d0
	valc_update%Vr = 0.0d0
    
    allocate(prob_mat(na,neps,ntheta),prob_mat_r(na),temp(neps,ntheta))
    
	! Update valc_update%Vy (young)
	
	!$omp parallel if (par_fortran_valc==1) default(shared) private(it,ie,ia,io, &
	!$ temp,io_prob,iap,iz_min,omega_val,util_val,prob_mat,prob_mat_r,EV_val)
    !$omp do collapse(2)
    do it = 1,ntheta ! current shock theta
        do ie = 1,neps ! current shock epsilon
            do ia = 1,na ! current asset holdings
            	do iz_min = 1,nz ! previous state (W,EP,ES,EC,R)
            		!write(*,*) "it,ie,ia,iz_min",it,ie,ia,iz_min
					!temp is (neps,ntheta) matrix
            		temp = outerprod2(P_eps(ie,:),P_theta(it,:))
                	do io = 1,no ! current period's occupation
                    !What is the probability of choosing io today?
                    io_prob = pol%occpol(ia,ie,it,iz_min,io)
                    ! Flow utility from consumption today
                    util_val= util_cons(pol%cpoly(ia,ie,it,io))
                    ! Next period's asset
					iap       = aopt_ind(ia,ie,it,io)
					omega_val = omega(ia,ie,it,io)
  
                    ! Prob_mat (young->young): dim na',neps',ntheta',no
                    prob_mat = 0.0d0
                    prob_mat(iap,:,:) = omega_val * temp
                    prob_mat(iap+1,:,:) = (1.0d0-omega_val) * temp
					if (abs(sum(prob_mat)-1.0d0)>1d-10) then
                    	call myerror("value_cons_one_step: prob_mat not equal to 1")
					endif
                    
					! Prob_mat (young->retired): dim na'
                    prob_mat_r = 0.0d0
                    prob_mat_r(iap)   = omega_val
                    prob_mat_r(iap+1) = (1.0d0-omega_val)

                    if (abs(sum(prob_mat_r)-1.0d0)>1d-10) then
                    	call myerror("value_cons_one_step: prob_mat_r not equal to 1")
					endif
					
                    ! EV_val: continuation value given ia,ie,it and occ choice io.
					EV_val = (1.0d0-pr_ret) * sum(prob_mat*valc%Vy(:,:,:,io)) + pr_ret* sum(prob_mat_r * valc%Vr)

                    ! Update valc
                    valc_update%Vy(ia,ie,it,iz_min) = valc_update%Vy(ia,ie,it,iz_min) + io_prob*(util_val + beta* EV_val)
                	enddo ! io  
                enddo !iz_min

            enddo !ia (current assets)
        enddo !ie
    enddo !it
	!$omp enddo
  	!$omp end parallel

	! Update valc_update%Vr (old)
	!$omp parallel if (par_fortran_valc==1) default(shared) private(ia, &
	!$ iap,omega_val,util_val,prob_mat,prob_mat_r,EV_val)
    !$omp do
	 do ia = 1,na
         util_val= util_cons(pol%cpolr(ia))
		 !What assets does guy with (a,eps,theta,5) choose for tomorrow?
		 iap      = aoptr_ind(ia)
		 omega_val= omegar(ia)
        
        ! Prob_mat (retired): dim na'
         prob_mat_r = 0.0d0
         prob_mat_r(iap)   = omega_val
         prob_mat_r(iap+1) = (1.0d0-omega_val)
         
        ! Prob_mat (young, after replacement): dim na',neps',ntheta' 
        prob_mat = 0.0d0
        prob_mat(iap,:,:) = omega_val * temp_ret
        prob_mat(iap+1,:,:) = (1.0d0-omega_val) * temp_ret
         
         ! EV_val: continuation value given ia 
         EV_val = (1.0d0-pr_die) * sum(prob_mat_r * valc%Vr) + pr_die* sum(prob_mat*valc%Vy(:,:,:,nz))
         
          ! Update valc
          valc_update%Vr(ia) =  util_val + beta* EV_val
	 enddo ! ia
	!$omp enddo
  	!$omp end parallel

   
end subroutine value_cons_one_step

!=======================================================================!
subroutine compute_value_hat(valhat,agg_C_NR,agg_labsup_NR,pol_NR,agg_C_R,agg_labsup_R,pol_R,val_R,dir_in)
    implicit none
	! This subroutine computes the hypothetical value function valhat holding the share of consumption and labor supply constant at the Reform economy levels.
    ! Note that valhat is the hypothetical value function for the Reform steady state.

	! Input variables
    real(8),intent(in) :: agg_C_NR,agg_labsup_NR,agg_C_R,agg_labsup_R ! aggregate consumption and labor supply in NR (no-reform) and R (reform)
    type(polfun),intent(in) :: pol_NR,pol_R ! policy functions in NR (no-reform) and R (reform)
	type(valfun), intent(in) :: val_R ! value function in R (reform)
	character(len=*),intent(in),optional :: dir_in
	! Output variables
    type(valfun),intent(out) :: valhat
    ! Local variables
	type(polfun) 		 :: polhat 
    type(valfun)         :: valhat0
    character(len=300)   :: pwd
    type(polfun)         :: pol
    integer              :: iter, istat
    real(8)              :: agg_C, agg_labsup,err, err_y, err_ret
    
    real(8), allocatable:: omega(:,:,:,:),omega_ret(:),temp_ret(:,:)
    integer, allocatable :: aopt_ind(:,:,:,:),aopt_ret_ind(:)
    
    
! Body of compute_value_hat

! Load policy functions from initial steady-state
write(*,*) "I am Inside compute_value_hat"
write(*,*) "Computing hypothetical policy functions..."

! Initialize polhat to fix sizes
polhat%cpoly = pol_NR%cpoly
polhat%cpolr = pol_NR%cpolr
polhat%lpoly_w = pol_NR%lpoly_w
! Compute hypothetical values
polhat%cpoly = pol_NR%cpoly / agg_C_NR * agg_C_R
polhat%cpolr = pol_NR%cpolr / agg_C_NR * agg_C_R
polhat%lpoly_w = pol_NR%lpoly_w / agg_labsup_NR * agg_labsup_R

! Initialize valhat0 and valhat to fix sizes
valhat0  = val_R
valhat   = val_R

! Precompute some stuff
write(*,*) "Calling precompute_distrib.."
allocate(omega(na,neps,ntheta,no),aopt_ind(na,neps,ntheta,no),omega_ret(na),aopt_ret_ind(na),temp_ret(neps,ntheta),stat=istat)
if (istat/=0) then
	call myerror("compute_value_hat: Allocation of failed!")
endif
call  precompute_distrib(aopt_ind,omega,aopt_ret_ind,omega_ret,temp_ret, pol_R)  


iter     = 1
err      = 10d0

write(*,*) "Start fixed point iteration to find valhat.."

do while (iter<=maxitValc .and. err>tolValc)
    !write(*,*) "iter = ", iter
    ! get updated valc (output)
    call  value_hat_one_step(valhat, valhat0,aopt_ind,omega,aopt_ret_ind,omega_ret, &
	temp_ret,polhat,pol_R)  
    
    ! Compute distance
    err_y = maxval(abs(valhat%Vy - valhat0%Vy))
    err_ret = maxval(abs(valhat%Vr - valhat0%Vr))
    
    err = max(err_y,err_ret)
    
    if (verbose>=1) then
		write(*,*) "======================================"
		write(*,*) "iter = ", iter
		write(*,*) "err_y = ", err_y
		write(*,*) "err_ret = ", err_ret
		write(*,*) "err   = ", err
		write(*,*) "======================================"
	endif
    
    ! Increase counter and update
    iter  = iter+1
    valhat0 = valhat
    
enddo

! save resulting valc to external file in folder 
if (present(dir_in)) then 
	call writedim(valhat%Vy, trim(dir_in)//trim('valhat_Vy.txt'))
	call writedim(valhat%Vr, trim(dir_in)//trim('valhat_Vr.txt'))
endif 

end subroutine compute_value_hat
!=======================================================================!
subroutine value_hat_one_step(valhat_update,valhat,aopt_ind,omega,aopt_ret_ind,omega_ret, &
	temp_ret,polhat,pol_R)  
implicit none
! Inputs and outputs
! - Suffix _ret refers to retired agents
! - Suffix _R (_NR) refers to reform (no reform)
type(valfun),intent(out) :: valhat_update
type(valfun),intent(in) :: valhat
type(polfun), intent(in) :: polhat,pol_R 
real(8), intent(in):: omega(:,:,:,:),omega_ret(:),temp_ret(:,:)
integer, intent(in) :: aopt_ind(:,:,:,:),aopt_ret_ind(:)
! Local variables
integer :: it,ie,ia,io,iap,iz_min
real(8) :: util_vec(4),temp(neps,ntheta),io_prob,util_val,omega_val,EV_val,lval
real(8),allocatable :: prob_mat(:,:,:),prob_mat_ret(:)
	
	valhat_update = valhat
	valhat_update%Vy = 0.0d0
	valhat_update%Vr = 0.0d0
    
    allocate(prob_mat(na,neps,ntheta),prob_mat_ret(na))
    
	! Update valhat_update%Vy (young)
	
	!$omp parallel if (par_fortran==1) default(shared) private(it,ie,ia,io, &
	!$ temp,io_prob,iap,iz_min,omega_val,util_val,prob_mat,prob_mat_ret,EV_val,lval)
    !$omp do collapse(4)
    do it = 1,ntheta ! current shock theta
        do ie = 1,neps ! current shock epsilon
            do ia = 1,na ! current asset holdings
            	do iz_min = 1,nz ! previous state (W,EP,ES,EC,R)
            		!temp is (neps,ntheta) matrix
            		temp = outerprod2(P_eps(ie,:),P_theta(it,:))
                	do io = 1,no ! current period's occupation
                    !What is the probability of choosing io today?
                    io_prob = pol_R%occpol(ia,ie,it,iz_min,io)
                    ! Flow utility from hypothetical consumption and labor supply today
					if (io == 1) then
						lval = polhat%lpoly_w(ia,ie,it)
					else
						lval = 0.0d0
					endif
                    util_val= util(polhat%cpoly(ia,ie,it,io),lval)
                    ! Next period's asset
					iap       = aopt_ind(ia,ie,it,io)
					omega_val = omega(ia,ie,it,io)
  
                    ! Prob_mat (young->young): dim na',neps',ntheta',no
                    prob_mat = 0.0d0
                    prob_mat(iap,:,:) = omega_val * temp
                    prob_mat(iap+1,:,:) = (1.0d0-omega_val) * temp
					if (abs(sum(prob_mat)-1.0d0)>1d-10) then
                    	call myerror("value_hat_one_step: prob_mat not equal to 1")
					endif
                    ! Prob_mat (young->retired): dim na'
                    prob_mat_ret = 0.0d0
                    prob_mat_ret(iap)   = omega_val
                    prob_mat_ret(iap+1) = (1.0d0-omega_val)
                    if (abs(sum(prob_mat_ret)-1.0d0)>1d-10) then
                    	call myerror("value_hat_one_step: prob_mat_ret not equal to 1")
					endif
					
                    ! EV_val: continuation value given ia,ie,it and occ choice io.
                    EV_val = (1.0d0-pr_ret) * sum(prob_mat * valhat%Vy(:,:,:,io)) + pr_ret* sum(prob_mat_ret * valhat%Vr)
                    
                    ! Update valhat
                    valhat_update%Vy(ia,ie,it,iz_min) = valhat_update%Vy(ia,ie,it,iz_min) + io_prob*(util_val + beta* EV_val)
                	enddo ! io  
                enddo !iz_min

            enddo !ia (current assets)
        enddo !ie
    enddo !it
	!$omp enddo
  	!$omp end parallel

	! Update valc_update%Vr (old)
	!$omp parallel if (par_fortran==1) default(shared) private(ia, &
	!$ iap,omega_val,util_val,prob_mat,prob_mat_ret,EV_val)
    !$omp do
	 do ia = 1,na
         util_val= util(polhat%cpolr(ia),0.0d0)
		 !What assets does guy with (a,eps,theta,nz) choose for tomorrow?
		 iap      = aopt_ret_ind(ia)
		 omega_val= omega_ret(ia)
        
        ! Prob_mat (retired): dim na'
         prob_mat_ret = 0.0d0
         prob_mat_ret(iap)   = omega_val
         prob_mat_ret(iap+1) = (1.0d0-omega_val)
         
        ! Prob_mat (young, after replacement): dim na',neps',ntheta' 
        prob_mat = 0.0d0
        prob_mat(iap,:,:) = omega_val * temp_ret
        prob_mat(iap+1,:,:) = (1.0d0-omega_val) * temp_ret
         
         ! EV_val: continuation value given ia 
         EV_val = (1.0d0-pr_die) * sum(prob_mat_ret * valhat%Vr) + pr_die* sum(prob_mat * valhat%Vy(:,:,:,nz))
         
          ! Update valc
          valhat_update%Vr(ia) =  util_val + beta* EV_val
	 enddo ! ia
	!$omp enddo
  	!$omp end parallel

   
end subroutine value_hat_one_step
!=======================================================================!
subroutine compute_value_hat_tran(valhat_t0, agg_C_init,agg_labsup_init,agg_C_final,agg_labsup_final, &
			pol_init,pol_final, val_final, pol_path,agg_path,dir_in)
    implicit none
	! This subroutine computes the hypothetical value function valhat at t=0 on the transition path
	! holding the share of consumption and labor supply constant at the Reform economy levels.

	! Inputs: 
	real(8),intent(in) :: agg_C_init,agg_labsup_init,agg_C_final,agg_labsup_final
	type(polfun),intent(in) :: pol_init,pol_final
	type(valfun),intent(in) :: val_final
	type(polfun),intent(in) :: pol_path(:)
	type(modelAgg),intent(in) :: agg_path(:)
	character(len=*),intent(in),optional :: dir_in
	! Output variables
    type(valfun),intent(out) :: valhat_t0
    ! Local variables
	type(polfun) 		 :: polhat 
	type(valfun)         :: valhat_t1,valhat_t
    integer              :: t,TT, istat
    real(8), allocatable :: omega(:,:,:,:),omega_ret(:),temp_ret(:,:)
    integer, allocatable :: aopt_ind(:,:,:,:),aopt_ret_ind(:)
    
    
! Body of compute_value_hat_tran

TT = size(agg_path)

! Load policy functions from initial steady-state
write(*,*) "I am Inside compute_value_hat_tran"

write(*,*) "Step 1: Compute valhat at the new steady state..."
call compute_value_hat(valhat_t1, agg_C_init,agg_labsup_init,pol_init, &
	agg_C_final,agg_labsup_final,pol_final,val_final)

write(*,*) "Step 2: Backward induction to compute valhat for t = T-1,T-2,...,1"

allocate(omega(na,neps,ntheta,no),aopt_ind(na,neps,ntheta,no),omega_ret(na),aopt_ret_ind(na), &
	temp_ret(neps,ntheta),stat=istat)
if (istat/=0) then
	call myerror("compute_value_hat_tran: Allocation of failed!")
endif

do t = TT-1,1,-1
	!write(*,*) "t = ", t	
	! Precompute some stuff
	! Compute hypothetical values
	polhat%cpoly = pol_init%cpoly / agg_C_init * agg_path(t)%C
	polhat%cpolr = pol_init%cpolr / agg_C_init * agg_path(t)%C
	polhat%lpoly_w = pol_init%lpoly_w / agg_labsup_init * agg_path(t)%labsup

	!write(*,*) "Calling precompute_distrib.."
	call  precompute_distrib(aopt_ind,omega,aopt_ret_ind,omega_ret,temp_ret, pol_path(t))  

	! Compute valhat for period t
	call  value_hat_one_step(valhat_t,valhat_t1,aopt_ind,omega,aopt_ret_ind,omega_ret, &
		temp_ret,polhat,pol_path(t))

	! Update valhat_t1  
	valhat_t1 = valhat_t

enddo ! t
valhat_t0 = valhat_t1


! save resulting valc to external file in folder 
if (present(dir_in)) then 
	call writedim(valhat_t0%Vy, trim(dir_in)//trim('valhat_t0_Vy.txt'))
	call writedim(valhat_t0%Vr, trim(dir_in)//trim('valhat_t0_Vr.txt'))
endif 

end subroutine compute_value_hat_tran
!=======================================================================!
! End of module mod_welfare
end module mod_welfare