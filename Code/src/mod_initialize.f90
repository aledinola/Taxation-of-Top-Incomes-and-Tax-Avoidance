module mod_initialize
	
    use mod_globals
    use mod_utilities
    use mod_numerical, only: myerror,discretize_pareto,linspace,grid,my_ss,tauchen,discretize_AR,kron
    use mod_functions
    implicit none
    
contains
    !======================================================
	! Initialize prices
    subroutine initialize_prices(prices0)
    implicit none
    ! Declare inputs/outputs
    type(modelPrices),intent(out) :: prices0
    ! Local variables
    real(8) :: r_min,r_max
    
    ! Initial values for r and tau_p are read from file estim_params
    !prices0%r = 0.021228699d0 !0.01901953d0! 0.0190195268d0  
    !prices0%tau_p = 0.12056381d0 
    prices0%r = r_glob 
    
    ! call fun_prices to update the other fields in "prices", including KN,w,y_ave,pen
	call fun_prices(prices0)

    ! The bracket is relevant only for bisection, not for fzero
    r_min = 0.01d0
    r_max = 0.05d0 !0.06d0
	prices0%r_brac = [r_min, r_max]
	
	! Lump-sum transfer to balance gov budget constraint
	prices0%tr = 0.0d0
    
    end subroutine initialize_prices
    !======================================================
    
    ! For initializing invariant model objects
    subroutine initialize_params()
    implicit none
    
    !Declare local variables:
    integer :: i,istat
    real(8) :: che1, res_eps(neps),res_theta(ntheta), ergo_theta(ntheta)
    real(8) :: ave_eps,sig_esp_stdev
    real(8) :: P_eps_help(neps-1,neps-1), eps_grid_help(neps-1)
    real(8), allocatable :: res_eps_theta(:)

    if (verbose>0) then
        write(*,*) "============================================================="
        write(*,*) "NEW MODEL RUN"
        write(*,*) "============================================================="
        write(*,'(a,f10.4)') ' a_max    =', a_max
        write(*,'(a,f10.2)') ' a_space  =', a_space
        write(*,'(a,I5)')    ' na       =', na
        write(*,'(a,I5)')    ' neps     =', neps
        write(*,'(a,I5)')    ' ntheta   =', ntheta
        !write(*,'(a,I5)')    ' nx       =', nx_max_nonconvex
        write(*,'(a,f10.4)') ' Taste shock st.dev.= ', sig_e
    endif
    

    !=============================================================
    ! Initialize income tax threshold y_H (such that marginal tax function is continuous)
    !============================================================= 
    tau_h = top_rate ! benchmark SS, from US tax code
    
    !=============================================================
    ! Initialize G_bench: government spending (in level) in the benchmark model
    !============================================================= 
    if (do_benchmark) then
    	G_bench   = 0.0d0
    	y_H_bench = 0.0d0
    elseif (.not. do_benchmark) then
    	call readscalar(G_bench,trim(dir_bench)//trim("agg_G.txt"))
    	call readscalar(y_H_bench,trim(dir_bench)//trim("y_H.txt"))
    endif
    
    
    !=============================================================
    ! EPSILON SHOCK - WORKERS
    !============================================================= 
    ! Discretize the shock process for epsilon (worker ability)
    
    if (super_shock == .true.) then
		if (verbose>=1) write(*,*) "Epsilon: Discretize AR(1) process (with superstar shock)"
		! Call discretize_AR to get the grid and trans mat for the non-superstar shock epsilons.
		! eps_grid_help has dimension neps-1
		! P_eps_help has dimension neps-1 by neps-1
    	call discretize_AR(rho_eps, 0d0, sig_eps, eps_grid_help, P_eps_help)	
        ! NOTE: Tauchen needs standard deviation as input, while discretize_AR wants the variance!!
        !sig_esp_stdev = SQRT(sig_eps)
		!call tauchen(rho_eps,sig_esp_stdev,0d0,3d0,neps-1, eps_grid_help, P_eps_help)
    
		! Add superstar epsilon to eps grid
        eps_grid(1:neps-1) = exp(eps_grid_help)
        eps_grid(neps) = eps_super
		
		! Construct the transition prob P_eps, including transitions to superstar shock
		P_eps = 0.0d0
		P_eps(:,neps) = prob_super
		do i = 1,(neps-1)
			P_eps(i,1:neps-1) = P_eps_help(i,:)
			! Rescale so that each row sums to 1
			P_eps(i,:) = P_eps(i,:)/sum(P_eps(i,:))		
        enddo
        
        ! superstar returns to the median
        !P_eps(neps,(neps/2+1)) = min(1.0d0,prob_super_back)
        !P_eps(neps,neps) = 1.0d0-min(1.0d0,prob_super_back) 
        
        ! superstar returns to all points with equal probability
        P_eps(neps,1:neps-1) = prob_super_back/real(neps-1,8)      
        P_eps(neps,neps) = 1.0d0 - sum(P_eps(neps,1:neps-1))
         
    else
    	! No superstar shock
    	if (verbose>=1) write(*,*) "Epsilon: Discretize AR(1) process (no superstar shock)"
    	call discretize_AR(rho_eps, 0d0, sig_eps, eps_grid, P_eps)
    	!call tauchen(rho_eps,sig_eps,0d0,3d0,neps, eps_grid, P_eps)
    	eps_grid = exp(eps_grid)
    endif
    
    !=============================================================
    ! THETA SHOCK - ENTREPRENEURS
    !=============================================================
    ! Discretize the shock process for theta (entre ability)
    ! We have three options for doing so.
    if (theta_calib==1) then
        
        if (verbose>=1) write(*,*) "Theta: Discretize AR(1) process"
        ! IMPORTANT NOTE: 
        ! Subroutine "discretize_AR" from toolbox uses Rouwenhorst method. IMPORTANT: sig_theta is the variance of the innovation.
        ! Subroutine "tauchen" from mod_numerical uses Tauchen's method. IMPORTANT: sig_theta is the stdev of the innovation.
        ! Note: Important difference between "discretize_AR" and "tauchen": 
        ! in "discretize_AR", if you increase the no. of grid points, the grid becomes wider. 
        ! In "tauchen" instead, the lower and upper bound of the grid are the same irrespective 
        ! of the number og grid points.
        
        ! CHOOSE discretize_AR or tauchen
        !call discretize_AR(rho_theta,uncmean_theta,sig_theta,theta_grid,P_theta)
        call tauchen(rho_theta,sig_theta,uncmean_theta,3d0,ntheta, theta_grid, P_theta)
    
        theta_grid = exp(theta_grid)
      
    else
    
        call myerror("initialize: theta_calib out of range")
    
    endif
    
	!Check if rows of Markov chains for (eps,theta) sum to 1:
    res_eps   = sum(P_eps,2)-1d0
    res_theta = sum(P_theta,2)-1d0
    if ( any(abs(res_eps)>1d-6) )  then
        call myerror('initialize: rows of P_eps do not sum to 1.')
    endif
    
    if ( any(abs(res_theta)>1d-6) )  then
        call myerror('initialize: rows of P_theta do not sum to 1.')
    endif
    
    !!! Combine eps and theta into a single transition matrix.
    ! Order: eps comes first
    P_eps_theta = kron(P_theta,P_eps)
    
    res_eps_theta = sum(P_eps_theta,2)-1d0
    if ( any(abs(res_eps_theta)>1d-6) )  then
        call myerror('initialize: rows of P_eps_theta do not sum to 1.')
    endif
    
    if (verbose>=1) then
        
        write(*,*) "Epsilon grid:" 
        do i=1,neps
            write(*,'(f10.4)') eps_grid(i)
        enddo
        write(*,*) "P_eps (Markov chain of epsilon):"
        call printMatrix(P_eps)
        
        write(*,*) " "
        
        write(*,*) "Theta grid:" 
        do i=1,ntheta
            write(*,'(f10.6)') theta_grid(i)
        enddo
        write(*,*) "P_theta (Markov chain of theta):"
        call printMatrix(P_theta)
        
        write(*,*) "Size of P_eps_theta", shape(P_eps_theta)
    
    endif 
    
    !Markov chain for age={Y,R}, where Y=young and R=retiree
    P_age(1,:) = [1.0d0-pr_ret, pr_ret]
    P_age(2,:) = [pr_die, 1.0d0-pr_die]
    
    !Compute invariant distrib of P_theta
    prob_theta = my_ss(P_theta,ntheta)
    che1 = sum(prob_theta)
    prob_theta = prob_theta/che1
    
    !Compute invariant distrib of P_eps
    prob_eps = my_ss(P_eps,neps)
    che1 = sum(prob_eps)
    prob_eps = prob_eps/che1
    
    ! Compute unconditional average of epsilon
    ave_eps = sum(prob_eps*eps_grid)
    
    ! Normalization so that E(eps) = 1
    eps_grid = eps_grid/ave_eps
    
    if (verbose>=1) write(*,*) "E(eps) = ", sum(prob_eps*eps_grid)
    
    ! initialize grid for assets
    call grid(a_grid,a_min,a_max,a_space)
    !do i=2,5
    !write(*,*) i,a_grid(i)-a_grid(i-1),a_grid(i)>a_grid(i-1)
    !enddo
    !pause
    
    ! Initialize switching cost
    
    switch_cost = switch_cost_CP ! Uniform switching cost
    
    do i = 1,no
     switch_cost(i,i) = 0.0d0
    enddo
    switch_cost(nz,:) = 0.0d0

    ! Set switch cost from EC (io=4) to Passthroughs (io=2,3)
    !switch_cost(4,1:3)   = switch_cost_CP 
    ! Set switch cost from Passthroughs (io=2,3) to EC (io=4)
    !switch_cost(1:3,4)   = switch_cost_PC 
    
    if (verbose>=1) then
	
		if (par_fortran==0) then
			write(*,*) "OpenMP: OFF"
		elseif (par_fortran==1) then
			write(*,*) "OpenMP: ON"
		endif
	
		if (do_howard==0) then
			write(*,*) "Howard: OFF"
		elseif (do_howard==1) then
			write(*,*) "Howard: ON"
		endif
	
		if (v_old_guess==0) then
			write(*,*) "Reset initial value function"
		elseif (v_old_guess==1) then
			write(*,*) "Load initial value function from file"
		endif
        write(*,*) "checking do_GE = ",do_GE
		if (do_GE==1) then
			write(*,*) "GENERAL EQUILIBRIUM"
		elseif (do_GE==0) then
			write(*,*) "PARTIAL EQUILIBRIUM"
		endif
	endif !verbose    
    if (skew_max_non_convex<1d0) then
        call myerror("mod_initialize: skew_max_non_convex cannot be smaller than 1")
    endif
    
    ! Export some grids
    !call write1dim(a_grid,savedir//'a_grid.txt')
    !write(*,*) "exit initialize_params"
    end subroutine initialize_params
!======================================================
	! Initialize value functions, and distributions
    subroutine initialize_model(val0,distrib0)
    implicit none
    ! Declare inputs/outputs
	type(valfun),     intent(out) :: val0
	type(distribfun_m), intent(out) :: distrib0
	    
    !Declare local variables:
    integer :: i,istat
    real(8) :: sharey

    ! initialize value function 
    allocate(val0%Vy(na,neps,ntheta,nz),val0%Vr(na),stat=istat)
    if (istat/=0) then
		call myerror("initialize: Allocation of val0 failed!")
	endif
    if (v_old_guess==0) then
        val0%Vy = 0.0d0 !dim: (na,neps,ntheta,nz)
        val0%Vr = 0.0d0 !dim: (na)
    elseif (v_old_guess==1) then
        write(*,*) "Use old value function as initial guess"
        call readdim(val0%Vy,trim(savedir)//trim('V0.txt'))
        call readdim(val0%Vr,trim(savedir)//trim('VR_0.txt'))        
    endif
       
    ! initialize distribution
    allocate(distrib0%muy(na,neps,ntheta,nz),distrib0%mur(na),stat=istat)
    if (istat/=0) then
		call myerror("initialize: Allocation of distrib0 failed!")
	endif
	sharey       = pr_die / (pr_ret + pr_die)
    distrib0%muy = sharey/real(na*neps*ntheta*nz,8) !dim: (na,neps,ntheta,no)
    distrib0%mur = (1.0d0-sharey)/real(na,8)     !dim: (na)


    end subroutine initialize_model

!======================================================
	! Initialize data targets for steady state estimation
    subroutine initialize_data()
    implicit none
    ! This subroutine updates data_targets (global var)
	
	data_targets%share_entre_act = 0.15161
	data_targets%share_EP_entre = 0.6736
	data_targets%share_ES_entre = 0.2363
	data_targets%share_EC_entre = 0.09
	data_targets%share_wage_ES = 0.3627
	data_targets%share_wage_EC = 0.1988
	data_targets%ratio_aveinc_entre_worker = 2.601
	data_targets%ratio_medinc_entre_worker = 1.5574
	data_targets%empsize_rel_ES = 1.5024
	data_targets%empsize_rel_EC = 6.7736
	data_targets%share_inc_entre =0.3173
	data_targets%share_wealth_entre =0.5355
	data_targets%wealth_share_EP = 0.4006
	data_targets%wealth_share_ES = 0.4008
	data_targets%wealth_share_EC = 0.1986
	data_targets%inc_share_EP = 0.4618
	data_targets%inc_share_ES = 0.3581
	data_targets%inc_share_EC = 0.1801
	data_targets%gini_wealth_all =0.8422
	data_targets%gini_wealth_entre =0.7808
	data_targets%gini_inc_all =0.5445
	data_targets%gini_inc_entre =0.6219
	data_targets%wealth_share(1)  =0.3345
	data_targets%wealth_share(2)  =0.736
	data_targets%inc_share(1)  =0.1905
	data_targets%inc_share(2)  =0.4489
	data_targets%share_entre_top(1)  =0.8265
	data_targets%share_entre_top(2)  =0.6183
	data_targets%share_entre_top(3)  =0.4688
	data_targets%share_entre_top_inc(1)  =0.6681
	data_targets%share_entre_top_inc(2)  =0.5003
	data_targets%share_entre_top_inc(3)  =0.3771
	data_targets%share_entre_quint(1)  =0.0503
	data_targets%share_entre_quint(2)  =0.04548
	data_targets%share_entre_quint(3)  =0.1164
	data_targets%share_entre_quint(4)  =0.1863
	data_targets%share_entre_quint(5)  =0.3598
	data_targets%share_entre_quint_inc(1)  =0.08612
	data_targets%share_entre_quint_inc(2)  =0.08374
	data_targets%share_entre_quint_inc(3)  =0.14914
	data_targets%share_entre_quint_inc(4)  =0.1517
	data_targets%share_entre_quint_inc(5)  =0.28953
	data_targets%share_lfo_top(2,3) = 0.4387
	data_targets%share_lfo_top(3,3) = 0.3812
	data_targets%share_lfo_top(4,3)  = 0.18
	data_targets%share_lfo_top_inc(2,3)  =0.4574
	data_targets%share_lfo_top_inc(3,3)  =0.3645
	data_targets%share_lfo_top_inc(4,3)  =0.1781
	data_targets%share_lfo_quint(2,5)  =0.54029
	data_targets%share_lfo_quint(3,5)  =0.32606
	data_targets%share_lfo_quint(4,5)  =0.13365
	data_targets%share_lfo_quint_inc(2,5)  =0.49747
	data_targets%share_lfo_quint_inc(3,5)  =0.33923
	data_targets%share_lfo_quint_inc(4,5)  =0.1633
	data_targets%trans_entre_work  =0.2456

	
	end subroutine initialize_data
end module mod_initialize