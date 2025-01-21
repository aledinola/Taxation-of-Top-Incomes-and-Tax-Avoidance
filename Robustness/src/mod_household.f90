module mod_household
use mod_numerical, only: myerror,sub_logsum,linint,max_nonconvex
use mod_functions
use mod_globals
!use toolbox, only: plot,execplot
use omp_lib
use mod_utilities
implicit none 

private
public :: sub_vfi, precompute_entre, sub_vfi_onestep

contains

subroutine precompute_entre(cash_e,cash_r,kpol,npol,phipol,rbar,wage,pen,tr,y_ss_cap,y_H)

	implicit none
	!----------------------------------------------------!
	! Inputs and outputs
	real(8), intent(in)  :: rbar,wage,pen,tr,y_ss_cap,y_H  ! Input "prices"
    real(8), intent(out) :: cash_e(:,:,:)  ! Cash-on-hand (a,theta,o) for o=EP,ES,EC
    real(8), intent(out) :: cash_r(:)      ! Cash-on-hand (a) for retired
    real(8), intent(out) :: kpol(:,:,:)    ! k(a,theta,o)  for o=EP,ES,EC
    real(8), intent(out) :: npol(:,:,:)    ! n(a,theta,o)  for o=EP,ES,EC
    real(8), intent(out) :: phipol(:,:,:)  ! phi(a,theta,o)  for o=ES,EC
	! Local variables
	integer :: ia, ith, io
    real(8) :: term1, front1,aux1,aux2,aux3,aux4,lambda_vec(4),x1,f1
    real(8) :: theta_val,a_val,k_unc1,k_unc2,k_con,tot_labor,profit_val,ytaxe,taxe
    real(8) :: ytaxr,taxr
	!----------------------------------------------------!

! Compute some auxiliary terms
term1   = 1.0d0/((1.0d0-gamma)*vi-1.0d0)
front1  = (wage/(vi*(1.0d0-gamma)))**term1	
aux2    = (gamma/(rbar+delta))**(1d0+gamma*vi/(1d0-vi))
aux3    = ((1.0d0-gamma)/wage)**(vi*(1d0-gamma)/(1d0-vi))
aux4    = 1d0/(1d0-(1d0-gamma)*vi)

lambda_vec = [0.0d0,lambda,lambda_es,lambda_ec]
phipol = 0.0d0
!write(*,*) "Start pre-computation static payoffs.."
! call disp("tau_p",tau_p)
! call disp("rbar",rbar)
! call disp("wage",wage)
! call disp("tr",tr)
! call disp("term1",term1)
! call disp("front1",front1)
! call disp("aux2",aux2)
! call disp("aux3",aux3)
! call disp("aux4",aux4)
! pause
do ith = 1,ntheta ! Current theta shock
    do ia = 1,na ! Current assets
    do io = 2,4
        theta_val = theta_grid(ith)
        a_val = a_grid(ia)
        aux1  = (theta_val*vi)**(1d0/(1d0-vi))
        ! STEP 1 - Assume n>=0
        ! Unconstrained capital
        k_unc1 = aux1*aux2*aux3
        k_con = min(k_unc1,lambda_vec(io)*a_val)
        ! tot_labor is nbar
        tot_labor = (theta_val*(1d0-gamma)*vi/wage)**aux4*k_con**(gamma*vi*aux4)
        ! STEP 2 - Check if tot_labor>=le i.e. n>=0
        if (tot_labor>=le) then
            kpol(ia,ith,io) = k_con
            npol(ia,ith,io) = max(0.0d0, tot_labor-le) !redundant
        elseif (tot_labor<le) then
            tot_labor = le
            k_unc2 = (theta_val*vi)**(1d0/(1d0-gamma*vi))*(gamma/(rbar+delta))**(1d0/(1d0-gamma*vi))*&
                tot_labor**(vi*(1d0-gamma)/(1d0-gamma*vi))
            kpol(ia,ith,io) = min(k_unc2,lambda_vec(io)*a_val)
            npol(ia,ith,io) = 0.0d0 ! hired labor
        else
            call myerror("precompute_entre: tot_labor is undefined")
        endif
        ! STEP 3 - Compute optimized profit and hence taxable income
        
    enddo ! end io 
    enddo ! end a
enddo ! end theta

! if (do_plots) then
! 	! Plots to check
! 	call plot(a_grid, lambda_vec(2)*a_grid, legend='lambda*a')
! 	call plot(a_grid,kpol(:,1,2),legend='theta1')
! 	call plot(a_grid,kpol(:,3,2),legend='theta3')
! 	call plot(a_grid,kpol(:,4,2),legend='theta4')
! 	call execplot(title='kpol, EP')

! 	! Plots to check
! 	call plot(a_grid,npol(:,1,2),legend='theta1')
! 	call plot(a_grid,npol(:,3,2),legend='theta3')
! 	call plot(a_grid,npol(:,4,2),legend='theta4')
! 	call execplot(title='npol, EP')
! endif

! - SOLE PROP
do ith = 1,ntheta ! Current theta shock
	theta_val = theta_grid(ith)
    do ia = 1,na ! Current assets
    	io = 2 ! EP is 2
    	a_val = a_grid(ia)
    	! Taxable income of sole-prop
    	profit_val = fun_profit(theta_val,kpol(ia,ith,io),le,npol(ia,ith,io),rbar,wage)
        ytaxe  = profit_val-fun_tax_pay(profit_val,y_ss_cap) + rbar*a_val 
        ! Income taxes paid by entre
        taxe   = fun_tax_inc(ytaxe,y_H)
        cash_e(ia,ith,io) = ytaxe-taxe+a_val+tr
!         if (ia == 450 .and. ith == 8) then
!         	write(*,*) "io",io
!         	write(*,*) "cash", cash_e(ia,ith,io)
!         	write(*,*) "a_val",a_val
!         	write(*,*) "profit_val", profit_val
!         endif
	enddo !end ia
enddo !end i theta

! S-CORPS
do ith = 1,ntheta ! Current theta shock
	theta_val = theta_grid(ith)
    do ia = 1,na ! Current assets
    	io = 3 ! ES is 3
    	a_val = a_grid(ia)
    	profit_val = fun_profit(theta_val,kpol(ia,ith,io),le,npol(ia,ith,io),rbar,wage)
        ! Determine phipol: share of income declared as wage
        if (cf_avoidance == 0) then
        	! baseline model environment
            call max_nonconvex(fun_cash_es, 0.0d0, 1.0d0, x1, f1,mytol=tol_golden,mynx=nx_max_nonconvex_phi)
        	phipol(ia,ith,io) = min(1.0d0,max(0.0d0,x1))
        	cash_e(ia,ith,io) = f1+tr
        elseif (cf_avoidance == 1 .or. cf_avoidance == 2 .or. cf_avoidance == 3) then 
        	! counterfactual environment (no tax avoidance)
        	phipol(ia,ith,io) = 1.0d0 ! all income must be declared as wages
        	f1 = fun_cash_es(1.0d0)
        	cash_e(ia,ith,io) = f1+tr
        else
        	call myerror("precompute_entre: cf_avoidance ill-defined!")
        endif
!         if (ia == 450 .and. ith == 8) then
!         	write(*,*) "io",io
!         	write(*,*) "cash", cash_e(ia,ith,io)
!         	write(*,*) "a_val",a_val
!         	write(*,*) "profit_val", profit_val
!         endif
	enddo !end ia
enddo !end i theta

! C-CORPS
do ith = 1,ntheta ! Current theta shock
	theta_val = theta_grid(ith)
    do ia = 1,na ! Current assets
    	io = 4 ! EC is 4
    	a_val = a_grid(ia)
    	profit_val = fun_profit(theta_val,kpol(ia,ith,io),le,npol(ia,ith,io),rbar,wage)
        ! Determine phipol: share of income declared as wage
        if (cf_avoidance == 0) then
        ! Baseline model environment
			call max_nonconvex(fun_cash_ec, 0.0d0, 1.0d0, x1, f1,mytol=tol_golden,mynx=nx_max_nonconvex_phi)
        	phipol(ia,ith,io) = min(1.0d0,max(0.0d0,x1))
        	cash_e(ia,ith,io) = f1+tr
        elseif (cf_avoidance == 1) then 
        ! counterfactual: no intensive margin of tax avoidance
        	phipol(ia,ith,io) = 0.0d0 ! force to declare all income as profits, so phi (wage share)=0
        	f1 = fun_cash_ec(0.0d0)
        	cash_e(ia,ith,io) = f1+tr
        elseif (cf_avoidance == 2 .or. cf_avoidance == 3) then 
        ! counterfactual: no tax avoidance on either margin
        	phipol(ia,ith,io) = 1.0d0 ! force to declare all income as wages
        	f1 = fun_cash_ec(1.0d0)
        	cash_e(ia,ith,io) = f1+tr
        else
        	call myerror("precompute_entre: cf_avoidance ill-defined!")
        endif
!         if (ia == 450 .and. ith == 8) then
!         	write(*,*) "io",io
!         	write(*,*) "cash", cash_e(ia,ith,io)
!         	write(*,*) "a_val",a_val
!         	write(*,*) "profit_val", profit_val
!         endif
	enddo !end ia
enddo !end i theta

!***Precompute old retirees problem
do ia = 1,na ! Current assets
    a_val = a_grid(ia)
    ytaxr = pen+rbar*a_val
    taxr  = fun_tax_inc(ytaxr,y_H)
    cash_r(ia) = ytaxr-taxr+a_val+tr
enddo

if (display_vfi == 1) then 
	write(*,*) "Pre-computation static payoffs done!"
	write(*,*) " "
endif

contains

	function fun_cash_es(phi) result(F)
    implicit none
    real(8), intent(in) :: phi ! must lie in [0,1]
    real(8) :: F
    ! Local variables
    real(8) :: wage_hat,profit_hat,y_es
    
    wage_hat   = phi*profit_val
    profit_hat = (1.0d0-phi)*profit_val
    y_es = profit_hat + wage_hat-fun_tax_pay(wage_hat,y_ss_cap)+rbar*a_val-cost_evasion_es(1.0d0-phi)-op_cost_es
    
    F = y_es - fun_tax_inc(y_es,y_H) + a_val
    
	end function fun_cash_es
	!---------------------------------------------!
	
	function fun_cash_ec(phi) result(F)
    implicit none
    real(8), intent(in) :: phi ! must lie in [0,1]
    real(8) :: F
    ! Local variables
    real(8) :: wage_hat,profit_hat,y_ec
    
    wage_hat   = phi*profit_val
    profit_hat = (1.0d0-phi)*profit_val
    y_ec = wage_hat-fun_tax_pay(wage_hat,y_ss_cap)+rbar*a_val-cost_evasion_ec(phi)-op_cost_ec
    
    F = (1.0d0-tau_c)*(1.0d0-tau_d)*profit_hat+y_ec-fun_tax_inc(y_ec,y_H) + a_val
    

	end function fun_cash_ec

end subroutine precompute_entre
!==============================================================================!

subroutine sub_expval(EVy,EVnewbw,Vy)
    ! Compute expected value of a young person and 
    ! value of NEWBORNS
    implicit none
    real(8), intent(in) :: Vy(:,:,:,:)
    real(8), intent(out) :: EVy(:,:,:,:), EVnewbw(:)
    ! Locals
    integer :: i,j,ip,jp
    
    ! Expected value of a young person            
    EVy=0.0d0 !dim: (na,neps,ntheta)
    !today
    do i = 1,neps                   
        do j = 1,ntheta
            !next period
            do ip = 1,neps         
                do jp = 1,ntheta
                    EVy(:,i,j,1:no) = EVy(:,i,j,1:no) + Vy(:,ip,jp,1:no)*P_eps(i,ip)*P_theta(j,jp)
                enddo
            enddo
        enddo
    enddo
        
    ! Expected value of V of NEWBORNS (they draw abilities from invariant distribution of eps,theta)
    EVnewbw=0.0d0 !dim: (na)
    do ip=1,neps ! tomorrow
        do jp=1,ntheta
            EVnewbw = EVnewbw + Vy(:,ip,jp,nz)*prob_eps(ip)*prob_theta(jp)
        enddo
    enddo


end subroutine sub_expval

!===============================================================================!

subroutine sub_vfi_onestep(newVy,newVr,cpoly,cpolr,lpoly_w,apoly,apolr,occpol,Vy,Vr,cash_e,cash_r,prices)

	implicit none
	!----------------------------------------------------!
	! Outputs:
	real(8), intent(out)   :: newVy(:,:,:,:), newVr(:)
	real(8), intent(out)   :: cpoly(:,:,:,:),cpolr(:),lpoly_w(:,:,:)
	real(8), intent(out)   :: apoly(:,:,:,:),apolr(:)
	real(8), intent(out)   :: occpol(:,:,:,:,:) !(a,eps,theta,z_min,o)
	! Inputs:
	real(8), intent(in)    :: Vy(:,:,:,:), Vr(:)
	real(8), intent(in)    :: cash_e(:,:,:), cash_r(:)
	type(modelPrices),intent(in) :: prices
	!type(polfun), intent(out) :: pol
	! Local variables
	real(8), allocatable :: EVy(:,:,:,:), EVnewbw(:)
	integer :: it,ie,ia,io,imax(1),iz_min
    real(8) :: err,wage,rbar,tr,l_H,eps_val,a_val,fmax,xmax,yw,yw_opt,cmax,cash,ss_cap
    real(8) :: lab_inc,y_H
    real(8) :: V_aux(4),prob_aux(4)

	!----------------------------------------------------!
	
	wage   = prices%w
	rbar   = prices%r
	tr     = prices%tr
	ss_cap = prices%ss_cap
	y_H    = prices%y_H
	
	! Initialize policies
	cpoly   = -1.0d0
	cpolr   = -1.0d0
	lpoly_w = -1.0d0
	apoly   = -1.0d0
	apolr   = -1.0d0
	occpol  = -1.0d0

	! Given Vy, compute expected values
	! Compute Expected value of a young person, EVy(a',eps,theta)
    ! and Expected value of EVnewbw(a') of NEWBORNS  
    ! (they draw abilities from invariant distribution of eps,theta)
    
    allocate(EVy(na,neps,ntheta,no),EVnewbw(na)) 
    ! Vy is input, Evy and EVnewbw are outputs
    call sub_expval(EVy,EVnewbw,Vy)
	
	!write(*,*) "sub_vfi_onestep: entering loop" 
	
	!$omp parallel if (par_fortran==1) default(shared) private(it,ie,ia,eps_val,a_val,l_H, &
    !$ V_aux,prob_aux,io,xmax,fmax,lab_inc,yw_opt,cmax,cash,iz_min)
    !$omp do collapse(3)
    do it = 1,ntheta
	do ie = 1,neps  ! current state - epsilon
		do ia = 1,na
			eps_val = eps_grid(ie)
			a_val = a_grid(ia)	
			! Get l_H such that taxable income is equal to y_H
			l_H  = fun_l_H(a_val,eps_val,wage,rbar,y_H,ss_cap)
			
			! Reset V_aux and prob_aux
			V_aux    = large_negative
			prob_aux = 0.0d0
			
			! - WORKERS DYNAMIC PROBLEM
    		io = 1 ! worker is 1
			! First, solve for the optimal labor and consumption when labor /= l_H
            ! We max over labor l. Outputs: xmax,fmax     
            call max_fun_rhs_w(xmax,fmax,wage,rbar,tr,y_H,ss_cap,eps_val,a_val,l_H,EVy(:,ie,it,io),Vr)
            V_aux(io)          = fmax
            cpoly(ia,ie,it,io) = fun_c_lab(xmax,wage,rbar,ss_cap,eps_val,a_val,l_H)
            lpoly_w(ia,ie,it)  = xmax
            lab_inc            = wage*eps_val*lpoly_w(ia,ie,it)
            yw_opt             = lab_inc-fun_tax_pay(lab_inc,ss_cap)+rbar*a_val
			apoly(ia,ie,it,io) = yw_opt - fun_tax_inc(yw_opt,y_H) + a_val - cpoly(ia,ie,it,io) +tr
			
			! Then, solve for the optimal labor and consumption when labor == l_H (kink)
			! This requires max over consumption
			if (l_H>=l_min .and. l_H<=l_max) then
				!yw = (1.0d0-tau_p)*wage*eps_val*l_H+rbar*a_val is equal to y_H
				cmax  = y_H - fun_tax_inc(y_H,y_H) + a_val - a_min +tr
				call max_fun_rhs_wH(xmax,fmax,cmax,l_H,y_H,wage,rbar,tr,y_H,eps_val,a_val,EVy(:,ie,it,io),Vr)
		
				! Check if the kink (l_H) is better. If so, update decisions to the kink.
				if (fmax>V_aux(io)) then
					V_aux(io)          = fmax
					cpoly(ia,ie,it,io) = xmax
					lpoly_w(ia,ie,it)  = l_H 
					lab_inc            = wage*eps_val*lpoly_w(ia,ie,it)
					yw_opt             = lab_inc-fun_tax_pay(lab_inc,ss_cap)+rbar*a_val
					apoly(ia,ie,it,io) = yw_opt - fun_tax_inc(yw_opt,y_H) &
					 + a_val - cpoly(ia,ie,it,io)+tr
				endif   
            endif
            
            ! - ENTRE DYNAMIC PROBLEM
			do io = 2,4 ! 2=EP, 3=ES, 4=EC
				cash = cash_e(ia,it,io)
				apoly(ia,ie,it,io)     = a_min
				if (cash>=a_min) then
					! Solve for the optimal savings a'
					call max_fun_rhs_e(xmax,fmax,cash,EVy(:,ie,it,io),Vr)
					V_aux(io)          = fmax
					apoly(ia,ie,it,io) = xmax
					cpoly(ia,ie,it,io) = cash - apoly(ia,ie,it,io)
				endif
			enddo !io
			
			! - OCCUPATIONAL CHOICE PROBABILITY
			if (cf_avoidance == 0) then
				! Baseline model environment
				! Loop over z_minus
				do iz_min = 1,nz
					call sub_logsum(newVy(ia,ie,it,iz_min),prob_aux, V_aux-switch_cost(iz_min,:),sig_e) 
    				occpol(ia,ie,it,iz_min,1:no) = prob_aux ! dim: (no)
    			enddo 
    		elseif (cf_avoidance == 1 .or. cf_avoidance == 2 .or. cf_avoidance == 3) then
    			! Counterfactual: no tax avoidance
    			! Set value of S-corps to a huge negative number since S-corps are strictly dominated by Sole-prop.
				if (cf_occpol == 0) then
    				V_aux(3) = large_negative
				endif
				if (cf_avoidance == 3) then
					V_aux(3) = large_negative
					V_aux(4) = large_negative
				endif
    			do iz_min = 1,nz
					call sub_logsum(newVy(ia,ie,it,iz_min),prob_aux, V_aux-switch_cost(iz_min,:),sig_e) 
    				occpol(ia,ie,it,iz_min,1:no) = prob_aux ! dim: (no)
					if (cf_occpol==1) then
						! Fixed occpol counterfactual: occpol set to benchmark (S-corp become EP)
						occpol(ia,ie,it,iz_min,1) = occpol_bench(ia,ie,it,iz_min,1)
						occpol(ia,ie,it,iz_min,2) = occpol_bench(ia,ie,it,iz_min,2)
						occpol(ia,ie,it,iz_min,3) = occpol_bench(ia,ie,it,iz_min,3)
						occpol(ia,ie,it,iz_min,4) = occpol_bench(ia,ie,it,iz_min,4)
						!write(*,*) "occpol", occpol(ia,ie,it,iz_min,:)
						!write(*,*) "updating newVy"
						!pause
						! Update newVy
						newVy(ia,ie,it,iz_min) = sum((V_aux-switch_cost(iz_min,:))*occpol(ia,ie,it,iz_min,1:no))
						!write(*,*) "newVy updated",newVy(ia,ie,it,iz_min)
						!pause
					endif


    			enddo 
				

    		else
    			call myerror("sub_vfi_onestep: cf_avoidance is ill-defined!")
    		endif
		enddo !ia
	enddo ! ie
    enddo ! it 
    !$omp enddo
    !$omp end parallel
    
    ! Retiree problem
    do ia = 1,na
		cash = cash_r(ia)
		! Solve for the optimal savings a'
		call max_fun_rhs_r(xmax,fmax,cash,EVnewbw,Vr)
		newVr(ia) = fmax
		apolr(ia) = xmax
		cpolr(ia) = cash - apolr(ia) 
	enddo !ia
    
! Check intermediate results

end subroutine sub_vfi_onestep
!==============================================================================!

 function fun_c_lab(lab,wage,rbar,ss_cap,eps_val,a_val,l_H) result(c_lab)
	 implicit none
	 ! Compute optimal consumption "c_lab" conditional on labor supply "lab"
	 real(8), intent(in) :: lab,wage,rbar,ss_cap,eps_val,a_val,l_H
	 real(8) :: c_lab
	 ! Local variables
	 real(8) :: l_ss_cap,taxable_inc,num,den
	 
	 ! Cutoff for social security
	 l_ss_cap = ss_cap/(wage*eps_val)
	 
	 if (lab<=l_H .and. lab<l_ss_cap) then
	 	taxable_inc = (1.0d0-tau_p)*wage*eps_val*lab+rbar*a_val
	 	num  = (1.0d0-tau_p)*wage*eps_val*hsv_0*(1.0d0-hsv_1)*taxable_inc**(-hsv_1)
	 elseif (lab>l_H .and. lab<l_ss_cap) then
		 num = (1.0d0-tau_h)*(1-tau_p)*wage*eps_val
	 elseif (lab<=l_H .and. lab>=l_ss_cap) then
	 	taxable_inc = (1.0d0-tau_p)*ss_cap+wage*eps_val*lab-ss_cap+rbar*a_val
	 	num  = wage*eps_val*hsv_0*(1.0d0-hsv_1)*taxable_inc**(-hsv_1)
	 elseif (lab>l_H .and. lab>=l_ss_cap) then
		 num = (1.0d0-tau_h)*wage*eps_val
	 else
		 call myerror("fun_c_lab: lab == undefined. This should not happen!")
	 endif
	 
	 den   = chi * lab**(sigma2)
	 c_lab = (num/den)**(1.0d0/sigma1)

 end function fun_c_lab
 !==============================================================================!
 function fun_rhs_w(lab,wage,rbar,tr,y_H,ss_cap,eps_val,a_val,l_H,EVy_vec,Vr) result(rhs_w)
    ! RHS of worker's value function given labor supply lab
    ! consumption is found using first order conditions (fun_c_lab)
	 implicit none
	 ! Inputs/outputs
	 real(8), intent(in) :: lab
	 real(8),intent(in)  :: wage,rbar,tr,y_H,ss_cap,eps_val,a_val,l_H
	 real(8),intent(in)  :: EVy_vec(:),Vr(:)
	 real(8) :: rhs_w
	 ! Local variables
	 real(8) :: EVy_int, Vr_int
	 real(8) :: c_lab, taxable_inc, aprime, aprime_lim, lab_inc
	 ! Inherited from host
	 ! ia, ieps, ith, EVy, Vo
 
     rhs_w = large_negative
 
	 ! Optimal consumption conditional on labor
	 c_lab       = fun_c_lab(lab,wage,rbar,ss_cap,eps_val,a_val,l_H)
	 lab_inc     = wage*eps_val*lab
	 taxable_inc = lab_inc-fun_tax_pay(lab_inc,ss_cap) + rbar*a_val
	 aprime      = taxable_inc-fun_tax_inc(taxable_inc,y_H) + a_val - c_lab + tr
	 aprime      = min(aprime,a_max)
	 !aprime_lim = max(aprime,a_min)
 
 	if (aprime>=a_min .and. c_lab>=0.0d0 .and.  lab>=0.0d0) then
	 	EVy_int = linint(a_grid,EVy_vec,aprime)
	 	Vr_int  = linint(a_grid,Vr,aprime)
	 	rhs_w   = util(c_lab,lab)+&
		 	beta*(1.0d0-pr_ret)*EVy_int+ beta*pr_ret*Vr_int
     endif
     
	 ! Add penalty for negative a'
	 !rhs_w = rhs_w+large_negative*abs(aprime-aprime_lim)

 end function fun_rhs_w
 !==============================================================================!
 subroutine max_fun_rhs_w(xmax,fmax,wage,rbar,tr,y_H,ss_cap,eps_val,a_val,l_H,EVy_vec,Vr)
 	implicit none
 	! Outputs/inputs
 	real(8),intent(out) :: xmax,fmax
 	real(8),intent(in)  :: wage,rbar,tr,y_H,ss_cap,eps_val,a_val,l_H
 	real(8),intent(in)  :: EVy_vec(:),Vr(:)
 	
 	! bounds for l must be chosen such that aprime is between a_grid(1) and a_grid(n_a)
    call max_nonconvex(wrapper_fun_rhs_w, l_min, l_max, xmax, fmax,mytol=tol_golden,mynx=nx_max_nonconvex)
	!call golden_method(wrapper_fun_rhs_w, l_min, l_max, xmax, fmax,mytol=tol_golden)
 
 contains  
    !------------------------------------------------
	 function wrapper_fun_rhs_w(lab)  result(rhs_w)
	 	implicit none
	 	real(8), intent(in) :: lab
	 	real(8) :: rhs_w
	 	rhs_w = fun_rhs_w(lab,wage,rbar,tr,y_H,ss_cap,eps_val,a_val,l_H,EVy_vec,Vr)
	 end function wrapper_fun_rhs_w
	!------------------------------------------------
 end subroutine max_fun_rhs_w
 !==============================================================================!
 
 function fun_rhs_w_cl(cons,lab,yw,wage,rbar,tr,y_H,eps_val,a_val,EVy_vec,Vr) result(rhs)
	 ! RHS of worker's value function given consumption cons and labor supply lab
	 implicit none
	 real(8), intent(in) :: cons,lab
	 real(8), intent(in) :: yw,wage,rbar,tr,y_H
	 real(8), intent(in) :: eps_val,a_val
	 real(8), intent(in) :: EVy_vec(:),Vr(:)
	 real(8) :: rhs
	 ! Local variables
	 real(8) :: EVy_int, Vr_int
	 real(8) :: c_lab, aprime, aprime_lim
	
	 if (cons<0d0) then
		 write(*,*) "cons = ", cons
		 call myerror("fun_rhs_w_cl: cons<=0")
	 endif
     
     rhs = large_negative
     
     aprime = yw - fun_tax_inc(yw,y_H) + a_val - cons + tr
	 aprime = min(aprime,a_max)
	
     if (cons>=0.0d0 .and. aprime>=a_min .and. lab>=0.0d0) then
		

		EVy_int = linint(a_grid,EVy_vec,aprime)
		Vr_int  = linint(a_grid,Vr,aprime)
	 	rhs   = util(cons,lab)+&
		 	beta*(1.0d0-pr_ret)*EVy_int+ beta*pr_ret*Vr_int
     endif
      
	 ! Add penalty for negative a'
	 !rhs = rhs+large_negative*abs(aprime-aprime_lim)

 end function fun_rhs_w_cl
!==============================================================================!
 subroutine max_fun_rhs_wH(xmax,fmax,cmax,l_H,yw,wage,rbar,tr,y_H,eps_val,a_val,EVy_vec,Vr)
 	implicit none
 	! Outputs/inputs
 	real(8),intent(out) :: xmax,fmax
 	real(8), intent(in) :: cmax,l_H,yw,wage,rbar,tr,y_H,eps_val,a_val
	real(8), intent(in) :: EVy_vec(:),Vr(:)
 	
 	call max_nonconvex(wrapper_fun_rhs_wH,0.0d0,cmax,xmax,fmax,mytol=tol_golden,mynx=nx_max_nonconvex)
	!call golden_method(wrapper_fun_rhs_wH,0.0d0,cmax,xmax,fmax,mytol=tol_golden)
 
 contains  
    !------------------------------------------------
	 function wrapper_fun_rhs_wH(cons)  result(rhs_wH)
	 	implicit none
	 	real(8), intent(in) :: cons
	 	real(8) :: rhs_wH
	 	rhs_wH = fun_rhs_w_cl(cons,l_H,yw,wage,rbar,tr,y_H,eps_val,a_val,EVy_vec,Vr) 
	 end function wrapper_fun_rhs_wH
	!------------------------------------------------
 end subroutine max_fun_rhs_wH
 !==============================================================================!
 function fun_rhs_e(aprime,cash,EVy_vec,Vr) result(rhs)
    ! RHS of entrepreneur's value function given a'
    ! Labor supply is fixed at parameter le
    ! Consumption follows from the budget constraint
	 implicit none
	 ! Inputs/outputs
	 real(8), intent(in) :: aprime
	 real(8),intent(in)  :: cash
	 real(8),intent(in)  :: EVy_vec(:),Vr(:)
	 real(8) :: rhs
	 ! Local variables
	 real(8) :: cons, cons_lim
	 real(8) :: EVy_int, Vr_int
 
 	 rhs = large_negative
 
	 ! Optimal labor conditional on c: requires root-finding
	 cons     = cash - aprime
	 !cons_lim = max(cons,1d-10)
 	 if (cons>=0.0d0 .and. aprime>=a_min ) then
	 	EVy_int = linint(a_grid,EVy_vec,aprime)
	 	Vr_int  = linint(a_grid,Vr,aprime)
	 	rhs   = util(cons,le)+&
		 	beta*(1.0d0-pr_ret)*EVy_int+ beta*pr_ret*Vr_int
     endif
     
	 ! Add penalty for negative c
	 !rhs = rhs+large_negative*abs(cons-cons_lim)

 end function fun_rhs_e
 !==============================================================================!
 subroutine max_fun_rhs_e(xmax,fmax,cash,EVy_vec,Vr)
 	implicit none
 	! Outputs/inputs
 	real(8),intent(out) :: xmax,fmax
	real(8),intent(in)  :: cash
	real(8),intent(in)  :: EVy_vec(:),Vr(:)
	
 	call max_nonconvex(wrapper_fun_rhs_e,a_min,min(a_max,cash),xmax,fmax,mytol=tol_golden,mynx=nx_max_nonconvex)
	!call golden_method(wrapper_fun_rhs_e,a_min,min(a_max,cash),xmax,fmax,mytol=tol_golden)
 
 contains  
    !------------------------------------------------
	 function wrapper_fun_rhs_e(aprime)  result(rhs)
	 	implicit none
	 	real(8), intent(in) :: aprime
	 	real(8) :: rhs
	 	rhs = fun_rhs_e(aprime,cash,EVy_vec,Vr)
	 end function wrapper_fun_rhs_e
	!------------------------------------------------
 end subroutine max_fun_rhs_e
 
 !===============================================================================!

subroutine howard_operator(Vy_inout,Vr_inout,pol,wage,rbar,tr,y_H,ss_cap,cash_e,cash_r)
    ! Perform one step of Howard acceleration
    ! Update once the value functions Vy(a,eps,theta) and Vo(a)
    implicit none
    ! Inputs/outputs
    real(8), intent(inout)  :: Vy_inout(:,:,:,:), Vr_inout(:)
    type(polfun),intent(in) :: pol
    real(8),intent(in)      :: wage,rbar,tr,y_H,ss_cap
    real(8),intent(in)      :: cash_e(:,:,:),cash_r(:)
    ! Local variables
    real(8),allocatable :: EVy(:,:,:,:), EVnewbw(:),Vy_update(:,:,:,:),Vr_update(:)
    integer :: ia,ie,it,io,iz_min
    real(8) :: eps_val,a_val,cons,lab,yw,aprime,cash,prob_aux(no),V_aux(no)
    
	! Given Vy_inout, compute expected values
    allocate(EVy(na,neps,ntheta,no),EVnewbw(na)) 
    call sub_expval(EVy,EVnewbw,Vy_inout)

	allocate(Vy_update(na,neps,ntheta,nz),Vr_update(na)) 
	
    ! Young Agents:
    
    !$omp parallel if (par_fortran==1) default(shared) private(it,ie,ia,eps_val,a_val, &
    !$ io,lab,cons,yw,aprime,V_aux,cash,prob_aux,iz_min)
    !$omp do collapse(3)
    do it=1,ntheta ! Current theta shock
        do ie=1,neps ! Current eps shock
            do ia = 1,na ! Current assets
            	eps_val = eps_grid(ie)
            	a_val   = a_grid(ia)
            	
            	! Reset V_aux and prob_aux
            	prob_aux = 0.0d0
            	V_aux    = large_negative
            	
            	! - WORKERS
            	io = 1
            	lab  = pol%lpoly_w(ia,ie,it)
            	cons = pol%cpoly(ia,ie,it,io)
            	yw   = wage*eps_val*lab - fun_tax_pay(wage*eps_val*lab,ss_cap)+rbar*a_val
                V_aux(io) = fun_rhs_w_cl(cons,lab,yw,wage,rbar,tr,y_H,eps_val,a_val,EVy(:,ie,it,io),Vr_inout) 

                ! - ENTREPRENEURS
                do io = 2,4
                	aprime = pol%apoly(ia,ie,it,io)
                	cash   = cash_e(ia,it,io)
                	if (cash>=a_min) then
                		V_aux(io) = fun_rhs_e(aprime,cash,EVy(:,ie,it,io),Vr_inout)
                	endif	
                enddo !io
                
                ! - UPDATE Vy
                if (cf_avoidance == 0) then
                	! Baseline model environment
                	! Loop over z_minus
					do iz_min = 1,nz
						call sub_logsum(Vy_update(ia,ie,it,iz_min),prob_aux, V_aux-switch_cost(iz_min,:),sig_e) 
    					!occpol(ia,ie,it,iz_min,1:no) = prob_aux ! dim: (no)
    				enddo                
                elseif (cf_avoidance == 1 .or. cf_avoidance == 2 .or. cf_avoidance ==3) then
    				! Counterfactual: no tax avoidance
    				! Set value of S-corps to a huge negative number since S-corps are strictly dominated by Sole-prop.
    				if (cf_occpol == 0) then
						V_aux(3) = large_negative
					endif 
					if (cf_avoidance == 3) then
						V_aux(3) = large_negative
						V_aux(4) = large_negative
					endif
                	do iz_min = 1,nz

						call sub_logsum(Vy_update(ia,ie,it,iz_min),prob_aux, V_aux-switch_cost(iz_min,:),sig_e) 
    					!occpol(ia,ie,it,iz_min,1:no) = prob_aux ! dim: (no)
						
    				enddo 
    			else
    				call myerror("howard_operator: cf_avoidance is ill-defined!")
    			endif
                
            enddo ! ia
        enddo ! ie
    enddo ! it
  	!$omp enddo
  	!$omp end parallel
    
    ! Retirees:
    do ia = 1,na ! Current assets
        aprime = pol%apolr(ia)
        cash   = cash_r(ia)
        Vr_update(ia)  =  fun_rhs_r(aprime,cash,EVnewbw,Vr_inout) 
    enddo ! ia
    
    Vy_inout = Vy_update
    Vr_inout = Vr_update
    
end subroutine howard_operator
 
 !==============================================================================!
 function fun_rhs_r(aprime,cash,EVnewbw,Vr) result(rhs)
    ! RHS of retiree's value function given consumption
	 implicit none
	 ! Inputs/outputs
	 real(8), intent(in) :: aprime
	 real(8),intent(in)  :: cash
	 real(8),intent(in)  :: EVnewbw(:),Vr(:)
	 real(8) :: rhs
	 ! Local variables
	 real(8) :: cons, cons_lim
	 real(8) :: EVy_int, Vr_int
     
     rhs = -huge(0.0d0)
     
	 cons     = cash - aprime
	 cons_lim = max(cons,1d-10)
 
	 Vr_int  = linint(a_grid,Vr,aprime)
 	 EVy_int = linint(a_grid,EVnewbw,aprime)

     if (cons>=1d-10) then
	 	rhs   = util(cons_lim,0.0d0)+&
		 	beta*(1.0d0-pr_die)*Vr_int+ beta*pr_die*EVy_int
     endif
     
	 ! Add penalty for negative consumption
	 !rhs = rhs+large_negative*abs(cons-cons_lim)

 end function fun_rhs_r
  !==============================================================================!
 subroutine max_fun_rhs_r(xmax,fmax,cash,EVnewbw,Vr)
 	implicit none
 	! Outputs/inputs
 	real(8),intent(out) :: xmax,fmax
	real(8),intent(in)  :: cash
	real(8),intent(in)  :: EVnewbw(:),Vr(:)
 	
 	call max_nonconvex(wrapper_fun_rhs_r,a_grid(1),min(a_grid(na),cash),xmax,fmax,mytol=tol_golden,mynx=nx_max_nonconvex)
	!call golden_method(wrapper_fun_rhs_r,a_grid(1),min(a_grid(na),cash),xmax,fmax,mytol=tol_golden)
 
 contains  
    !------------------------------------------------
	 function wrapper_fun_rhs_r(aprime)  result(rhs)
	 	implicit none
	 	real(8), intent(in) :: aprime
	 	real(8) :: rhs
	 	rhs = fun_rhs_r(aprime,cash,EVnewbw,Vr)
	 end function wrapper_fun_rhs_r
	!------------------------------------------------
 end subroutine max_fun_rhs_r
!==============================================================================!

subroutine sub_vfi(val,pol,exit_flag,prices)
	implicit none
	!----------------------------------------------------!
	! Inputs and outputs
	type(valfun), intent(inout)   :: val
	type(polfun), intent(out)     :: pol
	logical, intent(out)          :: exit_flag
	type(modelPrices), intent(in) :: prices
	! Local variables
	integer :: istat,iter,i_howard
    real(8) :: err,err_y,err_r
    real(8) :: t1, t2
    real(8), allocatable :: cash_e(:,:,:), cash_r(:), newVy(:,:,:,:), newVr(:)
	!----------------------------------------------------!

exit_flag = .true.

t1 = omp_get_wtime()
!call CPU_TIME(t1)

! - STEP 1: Solve entrepreneurs profit maximization 
allocate(cash_e(na,ntheta,4),cash_r(na),pol%kpol(na,ntheta,4),pol%npol(na,ntheta,4),&
	pol%phipol(na,ntheta,4),stat=istat)
if (istat/=0) then
	call myerror("sub_vfi: Allocation failed!")
endif

call precompute_entre(cash_e,cash_r,pol%kpol,pol%npol,pol%phipol,prices%r,prices%w,prices%pen,&
	prices%tr,prices%ss_cap,prices%y_H)
!call writedim(cash_e,trim(savedir)//trim('ss_bench')//trim(slash)//trim('cash_e.txt'))

! STEP 2 - Value function iteration
if (display_vfi == 1) then 
	write(*,*) "Iterations.."
endif
iter = 1
err  = 10.0d0

allocate(newVy(na,neps,ntheta,nz),newVr(na),pol%cpoly(na,neps,ntheta,no),pol%cpolr(na),pol%lpoly_w(na,neps,ntheta),&
	pol%apoly(na,neps,ntheta,no),pol%apolr(na),pol%occpol(na,neps,ntheta,nz,no),stat=istat)
if (istat/=0) then
	call myerror("sub_vfi: Allocation failed!")
endif

do while (iter<=maxitV .and. err>tolV)

	call sub_vfi_onestep(newVy,newVr,pol%cpoly,pol%cpolr,pol%lpoly_w,pol%apoly,pol%apolr,&
		pol%occpol,val%Vy,val%Vr,cash_e,cash_r,prices)
	
	! Howard accelerator
	if (do_howard == 1 .and. cf_occpol==0) then
		do i_howard = 1,n_howard		
			call howard_operator(newVy, newVr,pol,prices%w,prices%r,prices%tr,prices%y_H,prices%ss_cap,cash_e,cash_r)
		enddo
	endif
	
	err_y = maxval(abs(newVy-val%Vy))/max(maxval(abs(val%Vy)),1d0)
	err_r = maxval(abs(newVr-val%Vr))/max(maxval(abs(val%Vr)),1d0)
	err   = max(err_y,err_r)
	
	if (display_vfi==1) then
		write(*,*) "======================================"
		write(*,*) "iter = ", iter
		write(*,*) "err_y = ", err_y
		write(*,*) "err_r = ", err_r
		write(*,*) "err   = ", err
		write(*,*) "======================================"
	endif
	
	! Update
	val%Vy = newVy
	val%Vr = newVr
	iter   = iter +1
	
enddo

t2 = omp_get_wtime()
!call CPU_TIME(t2)
if (verbose>=1) then
	write(*,*) "============================================================="
	write(*,*) "VALUE FUNCTION ITER FINISHED"
	write(*,*) 'CASH+VFI runs for',real(t2-t1),'seconds.'
	write(*,*) "============================================================="
endif
if (err>tolV) then
    exit_flag = .false.
    write(*,*) 'WARNING: VFI did NOT converge!'
    write(*,*) 'final VFI distance     = ', err
    !pause
endif

! Check policy functions
! if (any(pol%occpol<-0.00001d0) .or. any(pol%occpol>1.00001d0)  ) then
! 	write(*,*) "max occpol", maxval(pol%occpol)
! 	write(*,*) "min occpol", minval(pol%occpol)
! 	call myerror("sub_vfi: occpol<0 or occpol>1")
! endif
! if (maxval(abs(sum(pol%occpol,dim=5)-1.0d0))>1d-8  ) then
! 	call myerror("sub_vfi: occpol does not sum to one")
! endif


end subroutine sub_vfi
!==============================================================================!

end module mod_household