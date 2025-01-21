module mod_transition
use mod_globals
use mod_numerical,only: myerror
use mod_utilities
use mod_functions
use mod_household, only: precompute_entre, sub_vfi_onestep
use mod_distribution, only: precompute_distrib, distrib_one_step, get_aggregates 
use mod_targets, only: targets_compute
!use toolbox, only: plot,execplot
use omp_lib
implicit none
contains

subroutine compute_transition(tau_h_path,KN_path,hsv_0_path,tau_p_path,pol_path_out,agg_path_out,val_t1,val_final, &
	 distrib_init,agg_final,dir_name_in)
	implicit none
	! Inputs/Outputs
	real(8),intent(in) :: tau_h_path(:)
	real(8),intent(inout) :: KN_path(:),hsv_0_path(:),tau_p_path(:)
	type(modelAgg), intent(in) :: agg_final
	type(distribfun_m),intent(in) :: distrib_init
	type(valfun),intent(in) :: val_final
	type(valfun),intent(out) :: val_t1
	type(polfun),intent(out) :: pol_path_out(:)
	type(modelAgg),intent(out) :: agg_path_out(:)
	! Optional inputs:
	character(len=*), intent(in), optional :: dir_name_in
	! Locals
	real(8) :: t1,t2
	real(8) :: err_tran_KN,err_tran_hsv_0,err_tran_tau_p
	real(8),allocatable :: KN_path_old(:),hsv_0_path_old(:),tau_p_path_old(:),T_vec(:), avek_path(:,:)
	real(8),allocatable :: share_entre_path(:,:),inc_share_EP_path(:), inc_share_ES_path(:), inc_share_EC_path(:), share_wage_ES_path(:),netinc_share_EP_path(:), netinc_share_ES_path(:), netinc_share_EC_path(:)
	integer :: T,pathcounter,t_c,unitno,ierr
	character(len=256) :: filename,dir_name
	type(polfun), allocatable :: pol_path(:)
	type(modelAgg), allocatable :: agg_path(:)
	type(modelPrices), allocatable :: prices_path(:)

	! Initialize optional argument
	if (present(dir_name_in)) then
		dir_name = dir_name_in
	else
		dir_name = trim(savedir)//trim("transition")//trim(slash)
		
	endif

	! Initialize convergence parameters	
	pathcounter         = 0
	err_tran_KN         = 10.0d0
	err_tran_hsv_0      = 10.0d0
	
	! Initialize paths
	T = size(tau_h_path)
	allocate(KN_path_old(T),hsv_0_path_old(T),tau_p_path_old(T),pol_path(T),T_vec(T),agg_path(T), &
		prices_path(T),avek_path(T,4),share_entre_path(T,4),inc_share_EP_path(T), inc_share_ES_path(T), inc_share_EC_path(T), share_wage_ES_path(T),netinc_share_EP_path(T), netinc_share_ES_path(T), netinc_share_EC_path(T))	
	! Create T_vec for the plots
	do t_c = 1,T
		T_vec(t_c) = real(t_c,8)
	enddo
	! T_vec = [ ( real(t_c,8), t_c=1,T )  ]
	
	agg_path(T)       = agg_final
	prices_path(T)%KN = KN_path(T)
	prices_path(T)%tr = hsv_0_path(T)
	call fun_prices_KN(prices_path(T),prices_path(T)%KN)
	
	! Initialize KN_path_old, tr_path_old and tau_p_path_old for the fixed point iteration
	KN_path_old = KN_path
	hsv_0_path_old = hsv_0_path
	tau_p_path_old = tau_p_path
	
	filename = trim(savedir)//trim("transition_loop.txt")
    open(newunit=unitno, file=trim(filename), status='replace',  iostat=ierr)
    write(unitno,*) "pathcounter  err_tran_KN		     err_tran_tau_p			err_tran_hsv_0"
	close(unitno)
	
	do while (pathcounter<maxiterations .and. (err_tran_KN >tol_KN .or. err_tran_hsv_0 >tol_hsv_0 .or. err_tran_tau_p >tol_tau_p ))
		
		t1 = omp_get_wtime()
		!call CPU_TIME(t1)
		pathcounter = pathcounter+1
	
		! Backward iteration on Value Function
		!First, go from T-1 to 1 calculating the Value function and Optimal
		!policy function at each step. Since we won't need to keep the value
		!functions for anything later we just store the next period one in
		!Vnext, and the current period one to be calculated in V
		!NOTE: for welfare we ned also the value function
		call sub_vfi_transition(pol_path,val_t1)
		!write(*,*) "sub_distrib_transition starting..."
		! Compute next-period distribution
		call sub_distrib_transition(KN_path,hsv_0_path,tau_p_path)

		! Compute distance
		! See how far apart the price paths are and then update
    	err_tran_KN=maxval(abs(KN_path-KN_path_old))
		err_tran_hsv_0=maxval(abs(hsv_0_path-hsv_0_path_old))
		err_tran_tau_p=maxval(abs(tau_p_path-tau_p_path_old))
   		! update
   		KN_path_old(1:T-1) = oldpathweight_KN*KN_path_old(1:T-1) + (1.0d0-oldpathweight_KN)*KN_path(1:T-1)
   		hsv_0_path_old(1:T-1) = oldpathweight_hsv_0*hsv_0_path_old(1:T-1) + (1.0d0-oldpathweight_hsv_0)*hsv_0_path(1:T-1)
   		tau_p_path_old(1:T-1) = oldpathweight_tau_p*tau_p_path_old(1:T-1) + (1.0d0-oldpathweight_tau_p)*tau_p_path(1:T-1)

		t2 = omp_get_wtime()
		!call CPU_TIME(t2)
		! Plots to check
		
		write(*,*) "======================================"
		write(*,"(A,I2)")    "pathcounter    = ", pathcounter
		write(*,"(A,F12.8)") "err_tran_KN 	 = ", err_tran_KN
		write(*,"(A,F12.8)") "err_tran_hsv_0 = ", err_tran_hsv_0
		write(*,"(A,F12.8)") "err_tran_tau_p = ", err_tran_tau_p
		write(*,"(A,F12.8)") "time (sec.)	 = ", t2-t1
		write(*,*) "======================================"
		
		open(newunit=unitno,file=trim(filename),form='formatted',status='old',position='append',iostat=ierr)
		write(unitno,*) pathcounter,err_tran_KN,err_tran_tau_p,err_tran_hsv_0
		close(unitno)
		
		
	enddo
	
	! if (do_plots) then
	! 	call plot(T_vec, KN_path_old)
	! 	call execplot(title='KN path')
	
	! 	call plot(T_vec, tau_p_path_old)
	! 	call execplot(title='Payroll tax path')
	
	! 	call plot(T_vec, hsv_0_path_old)
	! 	call execplot(title='hsv_0 (lambda_i) path')
	! endif
	write(*,*) "Saving some policy functions to output pol_path_out.."

	do t_c = 1,T
		pol_path_out(t_c)%apoly = pol_path(t_c)%apoly
		pol_path_out(t_c)%apolr = pol_path(t_c)%apolr
		pol_path_out(t_c)%occpol = pol_path(t_c)%occpol
		agg_path_out(t_c) = agg_path(t_c)
	enddo

	write(*,*) "Writing objects to separate text files.."
			
	call sub_write_transition(trim(dir_name))
	
	
contains
!============================================================================
		subroutine sub_vfi_transition(pol_path,val_t1)
		implicit none
		! Inputs/outputs
		type(polfun),intent(out)  :: pol_path(:)
		type(valfun), intent(out) :: val_t1
		! Local variables
		integer :: t_c,istat
		type(modelPrices) :: prices_t
		type(polfun) :: pol_t
		type(valfun) :: val_t,val_new
		real(8),allocatable :: cash_e(:,:,:),cash_r(:)
		
		allocate(cash_e(na,ntheta,no),cash_r(na),pol_t%kpol(na,ntheta,no),pol_t%npol(na,ntheta,no),&
				pol_t%phipol(na,ntheta,no),pol_t%cpoly(na,neps,ntheta,no),pol_t%cpolr(na), &
				pol_t%lpoly_w(na,neps,ntheta),pol_t%apoly(na,neps,ntheta,no),pol_t%apolr(na), &
				pol_t%occpol(na,neps,ntheta,nz,no),val_t%Vy(na,neps,ntheta,nz),val_t%Vr(na), &
				val_new%Vy(na,neps,ntheta,nz),val_new%Vr(na), stat=istat)
		if (istat/=0) then
			call myerror("sub_vfi_transition: Allocation failed!")
		endif
		
		val_t   = val_final
		val_new = val_final
		do t_c = T-1,1,-1
    		!write(*,*) "t_c",t_c
    		! Assign values to globals according to transition path in t
    		tau_H = tau_h_path(t_c)
    		tau_p = tau_p_path_old(t_c)
    		hsv_0 = hsv_0_path_old(t_c)
    		
    		! Prices implied by KN
    		prices_t%KN = KN_path_old(t_c)
    		call fun_prices_KN(prices_t,prices_t%KN)
    		
    		! Compute income cutoff 
			if (do_benchmark) then
            	prices_t%y_H = fun_y_H(hsv_0,hsv_1,tau_H)
            else
                prices_t%y_H = min(y_H_bench, fun_y_H(hsv_0,hsv_1,tau_H) )
            endif
    		
! 			write(*,*) "KN",prices_t%KN
! 			write(*,*) "r",prices_t%r
! 			write(*,*) "w",prices_t%w
! 			write(*,*) "pen",prices_t%pen

			call precompute_entre(cash_e,cash_r,pol_t%kpol,pol_t%npol,pol_t%phipol, &
				prices_t%r,prices_t%w,prices_t%pen,prices_t%tr,prices_t%ss_cap,prices_t%y_H)
				
			! One-step vfi update
			call sub_vfi_onestep(val_new%Vy,val_new%Vr,pol_t%cpoly,pol_t%cpolr,pol_t%lpoly_w, &
				pol_t%apoly,pol_t%apolr,pol_t%occpol,val_t%Vy,val_t%Vr,cash_e,cash_r,prices_t)
			! Update
			val_t            = val_new
			pol_path(t_c)    = pol_t
			prices_path(t_c) = prices_t
   			
		enddo ! t_c
		! Assign the value function in the first period of the transition as output
		val_t1 = val_t
		
		end subroutine sub_vfi_transition
		!===============================================================!
		
		subroutine sub_distrib_transition(KN_path,hsv_0_path,tau_p_path)
		implicit none
		! Inputs/outputs
		real(8),intent(out) :: KN_path(:),hsv_0_path(:),tau_p_path(:)
		! Local variables
		real(8),allocatable :: muy_update(:,:,:,:),mur_update(:)
		real(8),allocatable :: omega(:,:,:,:),omegar(:)
		real(8) :: temp_ret(neps,ntheta),share_entre,share_EP,share_ES,share_EC
		integer,allocatable :: aopt_ind(:,:,:,:),aoptr_ind(:)
		integer :: t_c
		type(polfun) :: pol_t
		type(distribfun_m) :: distrib_t
		type(modelPrices) :: prices_t
		type(modelAgg) :: agg_t
		type(modelResults) :: model_targets_t
		! Allocate muy_update and mur_update
		allocate(muy_update(na,neps,ntheta,nz),mur_update(na))
		allocate(omega(na,neps,ntheta,no),aopt_ind(na,neps,ntheta,no),omegar(na),aoptr_ind(na))
		
		distrib_t = distrib_init
		do t_c = 1,T-1
       		!write(*,*) "t_c",t_c
       		pol_t = pol_path(t_c)  		
			call precompute_distrib(aopt_ind,omega,aoptr_ind,omegar, &
				temp_ret,pol_t)  
			call distrib_one_step(muy_update,mur_update,aopt_ind,omega,aoptr_ind,omegar, &
				temp_ret,distrib_t,pol_t)  
			
			! Update
			distrib_t%muy = muy_update
			distrib_t%mur = mur_update
			! Calculate entrepreneur and LFO shares of young population
			share_EP   = sum(distrib_t%muy(:,:,:,2))
   	 		share_ES   = sum(distrib_t%muy(:,:,:,3))
    		share_EC   = sum(distrib_t%muy(:,:,:,4))
			share_entre = (share_EP + share_ES + share_EC)/sum(distrib_t%muy)
			share_EP = share_EP / sum(distrib_t%muy(:,:,:,2:4))
			share_ES = share_ES / sum(distrib_t%muy(:,:,:,2:4))
			share_EC = share_EC / sum(distrib_t%muy(:,:,:,2:4))
			! Calculate prices and aggregates
			prices_t = prices_path(t_c)
			tau_p    = tau_p_path_old(t_c)
			call get_aggregates(agg_t,prices_t,distrib_t,pol_t)
			call targets_compute(model_targets_t,prices_t,distrib_t,pol_t)
			KN_path(t_c)    = agg_t%KN
			hsv_0_path(t_c) = fun_update_hsv_0(agg_t%calY1,agg_t%calY2,G_bench,agg_t%taxes_corp,agg_t%taxes_div)
			agg_path(t_c)   = agg_t
			tau_p_path(t_c) = agg_t%pen_tot / (agg_t%taxes_ss / tau_p)
			avek_path(t_c,:) = agg_t%avek
			share_entre_path(t_c,1) = share_entre
			share_entre_path(t_c,2) = share_EP
			share_entre_path(t_c,3) = share_ES
			share_entre_path(t_c,4) = share_EC
			! Path of selected model moments: inc_share_EP_path, inc_share_ES_path, inc_share_EC_path, share_wage_ES_path 
			inc_share_EP_path(t_c) = model_targets_t%inc_share_EP
			inc_share_ES_path(t_c) = model_targets_t%inc_share_ES
			inc_share_EC_path(t_c) = model_targets_t%inc_share_EC
			netinc_share_EP_path(t_c) = model_targets_t%netinc_share_EP
			netinc_share_ES_path(t_c) = model_targets_t%netinc_share_ES
			netinc_share_EC_path(t_c) = model_targets_t%netinc_share_EC
			share_wage_ES_path(t_c)       = model_targets_t%share_wage_ES
		enddo ! t_c		
		
		end subroutine sub_distrib_transition
		!===============================================================!
		
	subroutine sub_write_transition(dir)
		implicit none
		character(len=*), intent(in) :: dir
		
		! Paths of equilibrium objects
		call writedim(KN_path_old,  trim(dir)//trim('KN_path.txt')) 
		call writedim(hsv_0_path_old,  trim(dir)//trim('hsv_0_path.txt'))
		call writedim(tau_p_path_old,  trim(dir)//trim('tau_p_path.txt'))
		! Paths of aggregates and tax revenues
		call writedim(agg_path%Y,  trim(dir)//trim('Y_path.txt'))
		call writedim(agg_path%Y_entre,  trim(dir)//trim('Y_entre_path.txt'))

		! A_path is aggregate capital = K_entre+K_C
		call writedim(agg_path%A,  trim(dir)//trim('A_path.txt'))
		call writedim(agg_path%C,  trim(dir)//trim('C_path.txt'))
		call writedim(agg_path%taxes_ss,  trim(dir)//trim('taxes_ss_path.txt'))
		call writedim(agg_path%taxes_inc,  trim(dir)//trim('taxes_inc_path.txt'))
		call writedim(agg_path%taxes_div,  trim(dir)//trim('taxes_div_path.txt'))
		call writedim(agg_path%taxes_corp,  trim(dir)//trim('taxes_corp_path.txt'))
		call writedim(agg_path%taxes_tot,  trim(dir)//trim('taxes_tot_path.txt'))
		call writedim(agg_path%labshare_ES,  trim(dir)//trim('labshare_ES_path.txt'))
		call writedim(agg_path%labshare_EC,  trim(dir)//trim('labshare_EC_path.txt'))
		call writedim(agg_path%labshare_corp,  trim(dir)//trim('labshare_corp_path.txt'))
		call writedim(agg_path%labshare_allcorp,  trim(dir)//trim('labshare_allcorp_path.txt'))

		call writedim(avek_path, trim(dir)//trim('avek_path.txt'))
		
		! Paths of prices
		call writedim(prices_path%r, trim(dir)//trim('r_path.txt'))
		call writedim(prices_path%w, trim(dir)//trim('wage_path.txt'))
		
		! Paths of entre. and LFO shares
		call writedim(share_entre_path, trim(dir)//trim('share_entre_path.txt'))
		
		! Paths of selected model moments 
		call writedim(inc_share_EP_path,  trim(dir)//trim('inc_share_EP_path.txt'))
		call writedim(inc_share_ES_path,  trim(dir)//trim('inc_share_ES_path.txt'))
		call writedim(inc_share_EC_path,  trim(dir)//trim('inc_share_EC_path.txt'))
		call writedim(share_wage_ES_path,  trim(dir)//trim('share_wage_ES_path.txt'))
		call writedim(netinc_share_EP_path,  trim(dir)//trim('netinc_share_EP_path.txt'))
		call writedim(netinc_share_ES_path,  trim(dir)//trim('netinc_share_ES_path.txt'))
		call writedim(netinc_share_EC_path,  trim(dir)//trim('netinc_share_EC_path.txt'))
		! Paths of some policy functions
		!call writedim(pol_path%apoly, trim(dir)//trim('apoly_path.txt'))
		!call writedim(pol_path%apolr, trim(dir)//trim('apolr_path.txt'))
		!call writedim(pol_path%occpol, trim(dir)//trim('occpol_path.txt'))


		end subroutine sub_write_transition
		!===============================================================!
		
end subroutine compute_transition

end module mod_transition