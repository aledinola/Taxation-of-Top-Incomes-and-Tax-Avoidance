module mod_steady_state

	use mod_globals
	use mod_functions
	use mod_targets, only: targets_compute,txt_export,mom_export
	use mod_household, only: sub_vfi
	use mod_initialize,   only: initialize_model
	use mod_distribution, only: get_distribution,get_KN,get_aggregates
	use mod_numerical, only: zbrent,myerror,mywarning
	use omp_lib
	implicit none
	
	
	contains
	
	subroutine compute_steady_state(tau_p0,hsv_0_0,prices,model_targets,val,pol,distrib,agg,flag,dir_save_mom)
	implicit none
	! Declare inputs/outputs
	real(8),intent(inout)		     :: tau_p0,hsv_0_0
	type(modelPrices), intent(inout) :: prices
	type(modelResults), intent(out)  :: model_targets
	type(valfun), intent(out)        :: val
	type(polfun),   intent(out)      :: pol
	type(distribfun_m), intent(out)  :: distrib
	type(modelAgg), intent(out)      :: agg
	integer,intent(out)              :: flag ! convergence flag for VFI
	! Optional inputs:
	character(len=*), intent(in), optional :: dir_save_mom
	! Declare locals
	real(8) :: r0, ED,r_min,r_max,hsv_0_new,r_new,hsv_0_diff,r_diff,tau_p_diff,tau_p_new
	real(8) :: t1,t2
	integer :: iter,unitno,ierr
	character(len=256)  :: file_save_targets,filename

	t1 = omp_get_wtime()

	flag = 0
	r0   = prices%r
	
	
	if (do_GE==0) then
		tau_p = tau_p0
		hsv_0 = hsv_0_0	
	    ED = excess_demand(r0)
	elseif (do_GE==1) then
		! We clear the following (2) contraints:
		!	i. adjust r to clear: capital excess demand = 0
		!	ii. adjust tau_p to clear: total pension = total payroll taxes
		
		! Reset file to record results
    	filename = trim(savedir)//trim("GE_loop.txt")
    	open(newunit=unitno, file=trim(filename), status='replace',  iostat=ierr)
    	write(unitno,*) "iter		r		tau_p		r_diff		tau_p_diff"
		close(unitno)
				
		iter       = 0
		tau_p_diff = 1.0d0
		r_diff     = 1.0d0
		
		do while (iter <= maxit_ge .and. (tau_p_diff>tol_ge_tau_p .or. r_diff>tol_ge_r))
			
			! excess_demand updates distrib and pol
			tau_p = tau_p0
			hsv_0 = hsv_0_0

			if (mydebug==1) then      
				write(*,*) "before excess demand, r0 = ", r0	
				pause
			endif  
		    ED    = excess_demand(r0)    
			! Get new bal_gov and new KN
			call get_aggregates(agg,prices,distrib,pol)
			! call fun_prices_KN to update interest rate given KN (prices is inout)		
			call fun_prices_KN(prices,agg%KN)
			
			! Update r and tau_p
			r_new     = prices%r
			!write(*,*) "after excess demand, r_new = ", r_new
			!pause 
			tau_p_new = agg%pen_tot / (agg%taxes_ss / tau_p)
			
			! Compute difference
		    tau_p_diff = abs(tau_p_new-tau_p0)
		    r_diff     = abs(r_new - r0)
		    		
			! Update tr and r using a convex combination of old guess and implied values
			tau_p0 = old_weight_tau_p*tau_p0 + (1.0d0-old_weight_tau_p) * tau_p_new
			r0     = old_weight_r *r0        + (1.0d0-old_weight_r)  * r_new
			! Communicate the updated tr and r to structure with prices
			prices%r  = r0

			if (mydebug==1) then 
				write(*,*) "Updated r0 = ", r0
				pause
			endif 

			if (verbose>=1) then
				write(*,*) "==========================="
				write(*,"(A,I)")     "iter       = ", iter
				write(*,"(A,F12.8)") "r          = ", prices%r
				write(*,"(A,F12.8)") "tau_p      = ", tau_p0
				write(*,"(A,F12.8)") "r_diff     = ", r_diff
				write(*,"(A,F12.8)") "tau_p_diff = ", tau_p_diff
			endif
			open(newunit=unitno,file=trim(filename),form='formatted',status='old',position='append',iostat=ierr)
			write(unitno,*) iter,prices%r,tau_p0,r_diff,tau_p_diff
			close(unitno)
			
			iter = iter + 1

		enddo !iter
		
	elseif (do_GE==2) then
        ! We clear the following (3) contraints:
		!	i. adjust r to clear: capital excess demand = 0
		!	ii. adjust tau_p to clear: total pension = total payroll taxes
		!	iii. adjust hsv_0 (lambda_i) to clear: taxes_inc+taxes_corp+taxes_div = G_bench
		
		! Reset file to record results
    	filename = trim(savedir)//trim("GE_loop.txt")
    	open(newunit=unitno, file=trim(filename), status='replace',  iostat=ierr)
    	write(unitno,*) "iter		r		tau_p		hsv_0		r_diff		tau_p_diff		hsv_0_diff"
		close(unitno)
		
		iter = 0
		hsv_0_diff = 1.0d0
		r_diff  = 1.0d0
		tau_p_diff = 1.0d0
		
		do while (iter <= maxit_ge2 .and. (tau_p_diff>tol_ge_tau_p .or. hsv_0_diff>tol_ge_hsv_0 .or. r_diff>tol_ge_r))
			
			! excess_demand updates distrib and pol
			tau_p = tau_p0
			hsv_0 = hsv_0_0
		    ED = excess_demand(r0)    
			! Get new bal_gov and new KN
			call get_aggregates(agg,prices,distrib,pol)
			prices%KN = agg%KN
			! call fun_prices_KN to update interest rate given KN (prices is inout)		
			call fun_prices_KN(prices,prices%KN)
			
			! Update r and tr
			r_new  = prices%r
			hsv_0_new = fun_update_hsv_0(agg%calY1,agg%calY2,G_bench,agg%taxes_corp,agg%taxes_div)
			tau_p_new = agg%pen_tot / (agg%taxes_ss / tau_p)

			! Compute difference
		    hsv_0_diff = abs(hsv_0_new-hsv_0_0)
		    r_diff  = abs(r_new - r0)
		    tau_p_diff = abs(tau_p_new-tau_p0)
				
			! Update tr and r using a convex combination of old guess and implied values
			hsv_0_0 = old_weight_hsv_0*hsv_0_0 + (1.0d0-old_weight_hsv_0) * hsv_0_new
			r0  = old_weight_r *r0  + (1.0d0-old_weight_r)  * r_new
			tau_p0 = old_weight_tau_p*tau_p0 + (1.0d0-old_weight_tau_p) * tau_p_new
			
			! Communicate the updated r to structure with prices
			prices%r  = r0
			
			if (verbose>=1) then
				write(*,*) "==========================="
				write(*,"(A,I)")     "iter       = ", iter
				write(*,"(A,F12.8)") "r          = ", prices%r
				write(*,"(A,F12.8)") "tau_p      = ", tau_p0				
				write(*,"(A,F12.8)") "hsv_0         = ", hsv_0_0
				write(*,"(A,F12.8)") "r_diff     = ", r_diff
				write(*,"(A,F12.8)") "tau_p_diff = ", tau_p_diff
				write(*,"(A,F12.8)") "hsv_0_diff    = ", hsv_0_diff
			endif
			open(newunit=unitno,file=trim(filename),form='formatted',status='old',position='append',iostat=ierr)
			write(unitno,*) iter,prices%r,tau_p0,hsv_0_0,r_diff,hsv_0_diff,tau_p_diff
			close(unitno)
			
			iter = iter + 1

		enddo !iter

	elseif (do_GE==3) then


		! We clear the following (2) contraints:
		!	i. adjust tau_p to clear: total pension = total payroll taxes
		!	ii. adjust hsv_0 (lambda_i) to clear: taxes_inc+taxes_corp+taxes_div = G_bench
		
		! Reset file to record results
    	filename = trim(savedir)//trim("GE_loop.txt")
    	open(newunit=unitno, file=trim(filename), status='replace',  iostat=ierr)
    	write(unitno,*) "iter		tau_p		hsv_0		tau_p_diff		hsv_0_diff"
		close(unitno)
		
		iter = 0
		hsv_0_diff = 1.0d0
		tau_p_diff = 1.0d0
		
		do while (iter <= maxit_ge2 .and. (tau_p_diff>tol_ge_tau_p .or. hsv_0_diff>tol_ge_hsv_0))
			
			! excess_demand updates distrib and pol
			tau_p = tau_p0
			hsv_0 = hsv_0_0
		    ED = excess_demand(r0)    
			! Get new bal_gov and new KN
			call get_aggregates(agg,prices,distrib,pol)
			prices%KN = agg%KN
			! We do not call fun_prices_KN since we don't want to update r and w		
			!call fun_prices_KN(prices)
			
			! Update hsv_0 and tau_p
			hsv_0_new = fun_update_hsv_0(agg%calY1,agg%calY2,G_bench,agg%taxes_corp,agg%taxes_div)
			tau_p_new = agg%pen_tot / (agg%taxes_ss / tau_p)

			! Compute difference
		    hsv_0_diff = abs(hsv_0_new-hsv_0_0)
		    tau_p_diff = abs(tau_p_new-tau_p0)
				
			! Update tr and r using a convex combination of old guess and implied values
			hsv_0_0 = old_weight_hsv_0*hsv_0_0 + (1.0d0-old_weight_hsv_0) * hsv_0_new
			tau_p0 = old_weight_tau_p*tau_p0 + (1.0d0-old_weight_tau_p) * tau_p_new
			
			if (verbose>=1) then
				write(*,*) "==========================="
				write(*,"(A,I)")     "iter       = ", iter
				write(*,"(A,F12.8)") "r          = ", prices%r
				write(*,"(A,F12.8)") "tau_p      = ", tau_p0				
				write(*,"(A,F12.8)") "hsv_0         = ", hsv_0_0
				write(*,"(A,F12.8)") "tau_p_diff = ", tau_p_diff
				write(*,"(A,F12.8)") "hsv_0_diff    = ", hsv_0_diff
			endif
			open(newunit=unitno,file=trim(filename),form='formatted',status='old',position='append',iostat=ierr)
			write(unitno,*) iter,tau_p0,hsv_0_0,hsv_0_diff,tau_p_diff
			close(unitno)
			
			iter = iter + 1

		enddo !iter


	else 

		call myerror("compute_steady_state: do_GE must be 0,1,2,3")
		
	endif
	
	call get_aggregates(agg,prices,distrib,pol)
	if (verbose>=1) write(*,*) "aggregates computed."
	
	if (verbose>=1) write(*,*) "Computing targets..."
	call targets_compute(model_targets,prices, distrib,pol)
	
	file_save_targets =  trim(savedir)//trim("targets_model_manual.txt")
	call txt_export(file_save_targets,prices, agg,model_targets)
	! Write certain moments to txt files to be plotted
	if (present(dir_save_mom)) then
		call mom_export(dir_save_mom,model_targets)
	endif


	t2 = omp_get_wtime()
	
	if (verbose>=1) then
		write(*,*) "============================================================="
		write(*,*) " STEADY STATE COMPUTED ", "do_GE = ",do_GE
		write(*,*) "============================================================="
		write(*,'(a,f15.10)') " Interest rate       = ", r0
		write(*,'(a,f15.10)') " Excess demand       = ", ED
		write(*,*) 'Total time',real(t2-t1),'seconds.'
	endif
	

	contains
		!------------------------------------
		function excess_demand(r) result(ED)
			implicit none
			! Inputs/outputs
			real(8), intent(in) :: r
			real(8) :: ED
			! Local variables
			logical :: flag_vfi,flag_distrib
			real(8) :: KN_supply
		
			! Update r and prices
			prices%r = r
			! Given r, update other fields in prices, including KN, w, y_ave,pen,ss_cap,y_H
			call fun_prices(prices)
			
			! Compute income cutoff 
			if (do_benchmark) then
            	prices%y_H = fun_y_H(hsv_0,hsv_1,tau_H)
            else
				! Here we fix y_H if tau_H>tau_H_bench=39.6%
                prices%y_H = min(y_H_bench, fun_y_H(hsv_0,hsv_1,tau_H) )
            endif
		
			! Reset value function, distribution
			call initialize_model(val,distrib)
		
			! Value function 
			call sub_vfi(val,pol,flag_vfi,prices)
 			if (flag_vfi==.false.) then
 				call mywarning("excess_demand: VFI did not converge!")
				flag = -1
 			endif
		
			if (verbose>=1 .and. do_GE==0) write(*,*) "VFI completed!"

			! Distribution
			if (verbose>=1 .and. do_GE==0) write(*,*) "Starting distribution..."
			call get_distribution(distrib,pol,flag_distrib)
 			if (flag_distrib==.false.) then
 				call myerror("excess_demand: Distrib. did not converge!")
 			endif
			if (verbose>=1 .and. do_GE==0) write(*,*) "Distribution done!"
			
			! Aggregate
			call get_KN(KN_supply,distrib,pol)
		
			ED = prices%KN-KN_supply
		
			if (verbose>=1) then
				write(*,'(a,f15.10)') " Current r       = ", r
				write(*,'(a,f15.10)') " Current ED      = ", ED
				write(*,*) "============================================================="
			endif
		
			end function excess_demand
			!----------------------------
	end subroutine compute_steady_state
	
!==========================================================
	
end module mod_steady_state