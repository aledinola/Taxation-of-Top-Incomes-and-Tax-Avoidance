module mod_globals

	use mod_numerical, only: myerror
	implicit none
    !Here define type of parameters and arrays (as allocatable)
    !Arrays are then assigned a numerical value (i.e. they are initialized)
    !in the subroutine initialize()
	
    !Folder from where we read inputs
    character(len=10) :: inputdir
    !Folder where we save outputs
    character(len=10) :: savedir
    character(len=10) :: slash
    character(len=256) :: dir_bench
    ! '..\output\' means: go up one folder and then enter folder output
    integer, parameter :: par_fortran    = 1   ! Flag 0-1
    integer, parameter :: par_fortran_mu = 0   ! Flag 0-1
    integer, parameter :: par_fortran_valc = 1 ! Flag 0-1
    integer :: mydebug 
    
    !=============================================================
    ! DECLARING ECONOMIC PARAMETERS
    !=============================================================    
    ! Number of estimated parameters
    integer, parameter :: N_params = 23
    
    ! - SUPER STAR SHOCK
    logical, parameter :: super_shock     = .true. ! T if super-star shock for epsilon is added
    real(8) :: eps_super      ! this is in levels, not in log
    real(8) :: prob_super     !0.02d0 /8.0d0*9.0d0
    real(8) :: prob_super_back  !Bettina
    ! - PREFERENCES
    integer, parameter :: endo_ls = 1       ! 0=exogenous, 1= endogenous labor supply
    real(8)            :: beta              ! discount factor 
    real(8), parameter :: sigma1  = 1.5d0   ! CRRA utility function u(c,l), best at 2!!
    real(8), parameter :: sigma2  = 1.7d0   ! Based on Bettina (2021)
    real(8) :: chi      ! utility weight of leisure in u(c,l)
    real(8), parameter :: sig_e   = 0.01d0  ! Standard deviation of Gumbel shock (if =0d0, no smoothing, take max)
    
    ! - PRODUCTION TECHNOLOGY and FINANCIAL FRICTIONS
    real(8), parameter :: alpha     = 0.33d0  ! Cobb-Douglas coeff. on capital
    real(8), parameter :: delta     = 0.06d0  ! Depreciation rate for capital 
    real(8) :: vi        ! Span of control prodfun 0.88
    real(8) :: gamma     ! alpha/vi ! Capital share in entre sector
    real(8) :: lambda    ! Collateral constraint parameter for EP (Sole-prop.)
    real(8) :: lambda_es ! Collateral constraint parameter for ES (S-corp.)
    real(8) :: lambda_ec !2.1d0  ! Collateral constraint parameter for EC (C-corp.)
    
    ! - TAX AVOIDANCE TECHNOLOGY
    integer, parameter :: no_avoidance = 0       ! Flag 0/1 IMPLEMENTED??
    real(8) :: c0_es  ! cost of avoidance function, S-corp
    real(8), parameter :: c1_es        = 2.0d0   ! cost of avoidance function, S-corp
    real(8) :: c0_ec  ! cost of avoidance function, C-corp
    real(8), parameter :: c1_ec        = 2.0d0   ! cost of avoidance function, C-corp
    real(8) :: op_cost_es  ! cost of operating an S-corp (\kappa)
    real(8) :: op_cost_ec    !10000d0  ! cost of operating a C-corp (\kappa)
    
    ! - SHOCKS DISTRIBUTION - Epsilon (worker ability)
    real(8), parameter :: rho_eps   = 0.94d0  ! Persistence epsilon
    real(8), parameter :: sig_eps   = 0.02d0  ! Dispersion epsilon (*variance*)
    
    ! - SHOCKS DISTRIBUTION - Theta (entrepren ability)
    integer, parameter :: theta_calib   = 1     ! 1=Tauchen/Rouwen (AR1), 2=Brueggeman (2020),3=Pareto (Buera-Moll)
    !The params below are relevant only if theta_calib = 1
	real(8) :: rho_theta     ! Persistence theta
	real(8) :: sig_theta    ! Dispersion theta (standard deviation of the innovation)
	real(8) :: uncmean_theta  !-0.085d0  ! Mean theta

    ! - STOCHASTIC LIFE CYCLE
    real(8), parameter :: pr_ret = 0.022d0     !Prob of retiring (source: Bettina)
    real(8), parameter :: pr_die = 0.089d0     !Prob of dying (source: Bettina)
    
    ! - CORPORATE TAX T_C(y)
    real(8), parameter :: tau_c = 0.35d0    !0.21d0 !   !Corporate tax rate
    
    ! - DIVIDEND TAX T_d(y)
    real(8), parameter :: tau_d = 0.1806d0 ! average marginal dividend tax
    !real(8), parameter :: tau_d =  0.4d0*0.238d0 !0.2d0  !Maximum tax rate for dividends
    
    ! - PERSONAL INCOME TAX T_i(y)
    integer, parameter :: taxfunc = 4     ! 1= HSV, 2= flat rate, 3= step-wise. 4 = HSV + top income tax rate
    !Relevant only if taxfunc = 1:
    real(8) :: hsv_0 != 0.72d0    ! HSV lambda, level. 
    !real(8), parameter :: hsv_1 = 0.06d0  ! HSV tau, progressivity
    real(8) :: hsv_1 != 0.13d0  ! HSV tau, progressivity. 0.17 (Bakis, Kayma, Poschke), 
    
    ! Parameters for top income tax rate, only used when taxfunc=4
    real(8) :: tau_h   ! top income tax rate								
    real(8) :: y_H_bench 

    !proportional income tax rate (relevant only if taxfunc=2)
    real(8), parameter :: tau_lin = 0.15d0       
    real(8), parameter :: deduc   = 0.23d0 ! as a frac of GDP pc (in 2013: 12200/53117)
    !Number of tax brackets, relevant only if taxfunc = 3:
    integer, parameter :: n_brackets = 7
    !Marginal tax rates, relevant only if taxfunc = 3:
    real(8), parameter :: top_rate = 0.396d0 ! U.S. tax code 2013
    real(8), dimension(n_brackets) :: tau_m    = [0.1d0,0.15d0,0.25d0,0.28d0,0.33d0,0.35d0,top_rate]
    !Tax brackets (as frac of ave inc), relevant only if taxfunc = 3:
    real(8), dimension(n_brackets) :: tax_brac = [0.0d0,0.336d0,1.365d0,2.756d0,4.199d0,7.499d0,8.472d0] !norm by GDP pc.
    ! Number of average tax rates to compute
    integer, parameter :: n_atr = 6

    ! - SOCIAL SECURITY TAX AND PENSION BENEFITS
    real(8), parameter :: G_frac    = 0.17d0  !government spending (source: Kinderman and Krueger)
    real(8) :: G_bench   ! To be set to G_frac times equilibrium YY (GDP)
    real(8), parameter :: repl = 0.4d0  !replacement rate (source: Bettina)    
    real(8), parameter :: B_tot    = 0.400d0  !pension benefits (source: Bettina)
    real(8) :: r_glob    ! global interest rate
    real(8) :: tau_p    != 0.124d0  !Payroll tax rate, endogenous
    ! ss_cap_ratio: social security income cap divided by annual wage income
    ! - ss cap in 2013 is 113,700. Ref:https://www.ssa.gov/oact/COLA/cbb.html#Series
    ! - annual wage in 2013 is 49,808 from Quarterly Census of Employment and Wages. 
    !		Ref: 2013 excel file from https://www.bls.gov/cew/downloadable-data-files.htm. 
    real(8), parameter :: ss_cap_ratio = 2.283d0 !113,700/49,808 social security cap/ave annual wage
    										
    ! Switch cost (occupational/LFO) parameters
    real(8) :: switch_cost_CP, switch_cost_PC
    
    !=============================================================
    ! DECLARING COMPUTATIONAL PARAMETERS
    !=============================================================
    !integer, parameter :: dp = kind(0.0d0)
    real(8), parameter :: large_negative = -1d100 
    real(8), parameter :: small          = 1.0d-8
    
    logical :: do_plots,do_benchmark
    integer :: myOS,cf_avoidance,cf_occpol
    ! - GENERAL EQUILIBRIUM
    integer ::            do_GE               ! Flag 0/1/2
    integer ::			  do_run              ! Flag, set in main
    integer ::            do_fix_ss_cap       ! Flag 0/1, set in main
    integer ::            do_fix_pen          ! Flag 0/1, set in main

    real(8),parameter :: old_weight_hsv_0 = 0.5d0 ! Dampening parameter for hsv_0 (lambda_i) update
    real(8),parameter :: old_weight_r     = 0.9d0 ! Dampening parameter for r update
    real(8),parameter :: old_weight_tau_p = 0.5d0 ! Dampening parameter for tau_p update

	real(8),parameter  :: tol_ge_hsv_0 = 1.0d-6	! tolerance for hsv_0 (lambda_i) update (v41: 1d-6) 
	real(8),parameter  :: tol_ge_r     = 1.0d-6 ! tolerance for r update (v41: 1d-6)
	real(8),parameter  :: tol_ge_tau_p = 1.0d-6	! tolerance for tol_ge_tau_p update (v41: 1d-6)
	
    integer, parameter :: maxit_ge      = 30  ! for do_GE = 1
    integer, parameter :: maxit_ge2     = 100 ! for do_GE = 2

    ! - TRANSITION
    real(8),parameter  :: oldpathweight_KN    = 0.9d0
	real(8),parameter  :: tol_KN              = 1.0d-3 !2.0d-3
	real(8),parameter  :: oldpathweight_hsv_0 = 0.5d0
	real(8),parameter  :: tol_hsv_0           = 1.0d-3
	real(8),parameter  :: oldpathweight_tau_p = 0.5d0
	real(8),parameter  :: tol_tau_p           = 1.0d-3
	integer,parameter  :: maxiterations       = 100

    ! - VALUE FUNCTION ITERATION
    integer, parameter :: v_old_guess = 0     ! fFlag 0/1, if 1 it recycles old Value as initial guess
    real(8), parameter :: tolV        = 1d-5  ! Tolerance for value function (v41: 1d-5)
    integer, parameter :: maxitV      = 200    ! Maximum number of iterations for the value function (with howard algo, you never need more than 15-20 iterations)
    real(8), parameter :: tol_golden  = 1e-7  ! Tolerance crit for golden method (default is 1d-6)
    integer, parameter :: monotone    = 0     ! Flag 0/1/2 to use monotonicity trick in VFI
    integer, parameter :: concave     = 0     ! Flag 0/1/2 to use concavity trick in VFI
    integer, parameter :: do_howard   = 1     ! Flag 0/1 to use Howard algo in VFI
    integer, parameter :: n_howard    = 30    ! Num of iterations for Howard algo
    integer, parameter :: nx_max_nonconvex = 150 ! nx in the subroutine max_nonconvex (v41: 80)
    										   ! 1 = pure golden search; 
    integer, parameter :: nx_max_nonconvex_phi = 150  ! same as previous one but for computation of phi (v41: 80)                                        
    real(8), parameter :: skew_max_non_convex = 3d0 ! Skewness internal grid in max_non_convex
    real(8), parameter :: tolValc = 1d-5       ! Tolerance for Valc (see compute_value_cons in mod_welfare) (v41: 1d-5)
    integer, parameter :: maxitValc = 1000     ! Max no. of iterations for Valc
    
    ! - DISTRIBUTION
    real(8), parameter :: tol_dist    = 1d-11  ! Tolerance for the distribution (v41: 1d-9)
    integer, parameter :: maxit_dist  = 10000 ! Maximum number of iterations for the distribution mu
    
    ! - DISPLAY FLAGS
    integer            :: display_mu          ! Flag 0/1 to display MU iterns 
    integer            :: display_vfi         ! Flag 0/1 to display VF iterns
    integer            :: display_howard      ! Flag 0/1 to display Howard iterns
    
    ! - GRID SPECIFICATIONS
    integer, parameter :: nz   = 5  ! number of occ status {1=W, 2=EP, 3=ES, 4=EC, 5=R}
    integer, parameter :: no   = 4  !number of occ status when young {1=W, 2=EP, 3=ES, 4=EC}
    integer, parameter :: neps = 11 ! no. gridpoints for espilon shock
    !If we use Bettina's method (theta_calib=2), we need to 
    !set ntheta equal to 4.
    integer, parameter :: ntheta = 15    ! no. gridpoints theta shock

    ! Grid for assets:
    integer, parameter :: na      = 500     ! assets grid (v41: 500)
    real(8), parameter :: a_min   = 1d-6    ! Lower bound asset grid
    real(8), parameter :: a_max   = 1000.0d0 ! Upper bound asset grid
    real(8), parameter :: a_space = 3.5d0    ! grid spacing (if 1, evenly spaced)
    
    ! Grid for labor supply (hours worked):
    real(8), parameter :: l_min = 0.00d0   ! minimum hours
    real(8), parameter :: l_max = 1.6d0 ! maximum hours
    real(8), parameter :: le    = 0.0d0 !0.33d0    ! Fixed hours for entre
    
    integer :: nall = na*neps*ntheta*nz  !All dimensions for state vars 
    
    !=============================================================
    !   DECLARING GLOBAL SCALARS,SMALL ARRAYS,TYPES
    !   to be shared across subprograms
    !=============================================================
    
    integer :: verbose
    
    ! - Structure (derived data type) to collect value functions:
    type valfun
    	real(8), allocatable :: Vy(:,:,:,:) !dim: (na,eps,ntheta,nz)
    	real(8), allocatable :: Vr(:)     !dim: (na)
    end type valfun
    
    ! - Structure (derived data type) to collect policy functions:
	type polfun
		real(8), allocatable :: apoly(:,:,:,:)   ! a'(a,eps,theta,o), young, given o=W,EP,ES,EC
		real(8), allocatable :: apolr(:)         ! a'(a), retired
		real(8), allocatable :: cpoly(:,:,:,:)   ! c(a,eps,theta,o), young, given o=W,EP,ES,EC
		real(8), allocatable :: cpolr(:)         ! c(a), retired
		real(8), allocatable :: lpoly_w(:,:,:)   ! Hours worked l(a,eps,theta), given worker
		real(8), allocatable :: occpol(:,:,:,:,:)  ! prob_occ(a,eps,theta,z_,o), young, given o=W,EP,ES,EC
		! Only for entrepreneurs, static problem
		real(8), allocatable :: phipol(:,:,:)    ! phi(a,theta,o), given o=ES,EC
		real(8), allocatable :: kpol(:,:,:)      ! k(a,theta,o), given o=EP,ES,EC
		real(8), allocatable :: npol(:,:,:)      ! n(a,theta,o), given o=EP,ES,EC
	end type polfun
	
    ! For cf_occpol=1, occpol_bench is the occpol from the benchmark model
    real(8), allocatable :: occpol_bench(:,:,:,:,:)  ! prob_occ(a,eps,theta,z_,o), young, given o=W,EP,ES,EC


	! - Structure (derived data type) to collect distributions:
	type distribfun_m
    	real(8), allocatable :: muy(:,:,:,:) ! m(na,neps,ntheta,nz_min), young
    	real(8), allocatable :: mur(:)       ! m(na), retired
    end type distribfun_m
    type distribfun
    	real(8), allocatable :: muy(:,:,:,:) ! mu(na,neps,ntheta,no), young
    	real(8), allocatable :: mur(:)       ! mu(na), retired
    end type distribfun
    
    ! - Structure (derived data type) to collect model targets:
    type modelResults
    	! - Occ/LFO distribution
		real(8) :: share_entre_act,share_EP_entre,share_ES_entre,share_EC_entre,share_active
    	! - Share of income declared as wage (ES and EC)
		real(8) :: share_wage_ES, share_wage_EC, share_wage_qn_ES(3),share_wage_qn_EC(3),share_wage_qk_ES(3),share_wage_qk_EC(3),share_wage_EC_cf
		! - Share facing financial constraints and leverages (EP,ES,EC)
		real(8) :: frac_bor_ep, frac_bor_es, frac_bor_ec
		real(8) :: leverage_ep,leverage_es,leverage_ec
		! - Average eps,income,labor supply,income ratio
		real(8) :: ratio_aveinc_entre_worker, ratio_medinc_entre_worker,ave_l_work
		! - Share of employer
		real(8) :: share_empl_EP, share_empl_ES, share_empl_EC, share_empl_entre
		! - Average firm size relative to EP
		real(8) :: empsize_rel_ES, empsize_rel_EC
        real(8) :: payrollsize_lfo(4)
		! - Share of income and wealth owned by entre
		real(8) :: share_inc_entre, share_wealth_entre,share_toptax
		real(8) :: wealth_share_EP,wealth_share_ES,wealth_share_EC
		real(8) :: inc_share_EP,inc_share_ES,inc_share_EC
        real(8) :: netinc_share_EP,netinc_share_ES,netinc_share_EC
        real(8) :: median_inc_lfo(4),median_netinc_lfo(4),median_emp_lfo(4)
		! - Taxes paid by top earners
		real(8) :: taxes_top(4)
        real(8) :: atr_vec(n_atr) ! (1)Top 0.01% (2) P99.9+, (3) P99-P99.9, (4) P90-P99, (5) P50-P90, (6) up to P50
		! - Inequality measures
		real(8) :: gini_wealth_all, gini_wealth_entre,wealth_share(4)
		real(8) :: gini_inc_all, gini_inc_entre,inc_share(4)
		! - Occ/LFO choice by income and wealth
		real(8) :: share_entre_top(3), share_lfo_top(2:4,3),share_entre_quint(5),share_lfo_quint(2:4,5)    
		real(8) :: share_entre_top_inc(3), share_lfo_top_inc(2:4,3),share_entre_quint_inc(5),share_lfo_quint_inc(2:4,5)    
		! - Employment share by firm size
		real(8) :: share_emp(4)
		! - Transition rates
		real(8) :: trans_entre_work,trans_work_entre,trans_CP,trans_PC,trans_SC,trans_mat(4,4)
    
    end type modelResults
    
    type(modelResults) :: data_targets
    
    type modelPrices
    	real(8) :: r
    	real(8) :: r_brac(2)
    	real(8) :: w
    	real(8) :: KN
    	real(8) :: y_ave
    	real(8) :: pen
    	real(8) :: ss_cap
    	real(8) :: tr
    	real(8) :: y_H
    end type modelPrices
    
	 type modelAgg
	  	real(8) :: A,K_entre,K_C,labsup,N_entre,N_C,KN
		real(8) :: G,Y,Y_C,Y_entre,Y_lfo(4),KY, C, C_y,C_r,C_entre,C_worker,calY1,calY2
		real(8) :: taxes_ss,taxes_inc,taxes_div,taxes_corp,taxes_tot,taxes_occ(4),taxes_r
		real(8) :: pen_tot,res_gov,ED,bal_gov,bal_ss,op_cost_tot,avoidance_tot,ED_goods
		real(8) :: avetheta(4),avek(4),avey(4),avel(4),aveklratio(4),avekyratio(4)
        real(8) :: N_share(4), K_share(4)
        ! - labor share
        real(8) :: labshare,labshare_ES, labshare_EC, labshare_corp,labshare_allcorp
	 end type modelAgg
         
     type modelShares
         real(8) :: share_work_act
         real(8) :: share_entre_act
         real(8) :: share_EP_entre
         real(8) :: share_ES_entre
         real(8) :: share_EC_entre
         
     end type modelShares
     
     type(modelShares) :: shares
     
    !=============================================================
    ! DECLARING GRIDS, MARKOV TRANSITION MATRICES AND STATIONARY DISTRIBUTION FOR ABILITIES
    !=============================================================
    real(8), allocatable :: a_grid(:), eps_grid(:), theta_grid(:)
    real(8), allocatable :: P_eps(:,:), P_theta(:,:), P_age(:,:) 
    real(8), allocatable :: P_eps_theta(:,:)
    real(8), allocatable :: prob_theta(:), prob_eps(:)  !dim: (ntheta), (neps)
    real(8), allocatable :: switch_cost(:,:) !dim: (nz,no)
 
    !=================================================
    ! CHECKS 
    !=================================================
    real(8) :: wealth_share_ep, wealth_share_es, wealth_share_ec, inc_share_ep, inc_share_es, inc_share_ec = 0.0d0
    real(8) :: maxinc_work, maxinc_ep, maxinc_es, maxinc_ec

	contains
	!====================================================! 
    
    subroutine pack_params(par_vec)
    !! pack_params converts individual parameters (global variables)
    ! into a parameter vector "par_vec", to be used as input
    ! into fun_obj
    implicit none
    ! Declare variables
    real(8), intent(out) :: par_vec(:)
    integer :: i
    
    ! Order of parameters has to be the same as in unpack_params
    
    ! Execution
    i=1
    par_vec(i) = beta  ! beta
    i=i+1
    par_vec(i) = r_glob  ! interest rate 
    i = i+1
    par_vec(i) = tau_p  ! payroll tax rate
    i=i+1
    par_vec(i) = hsv_0 ! lambda_hsv
    i=i+1
    par_vec(i) = hsv_1 ! tau_hsv
    i=i+1
    par_vec(i) = chi ! chi
    i=i+1
    par_vec(i) = vi ! span of control
    i=i+1
    par_vec(i) = gamma ! capital share in the entrepreneurial sector
    i=i+1
    par_vec(i) = rho_theta 
    i=i+1
    par_vec(i) = sig_theta
    i=i+1
    par_vec(i) = uncmean_theta
    i=i+1
    par_vec(i) = lambda
    i=i+1
    par_vec(i) = lambda_es
    i=i+1
    par_vec(i) = lambda_ec
    i=i+1
    par_vec(i) = c0_es
    i=i+1
    par_vec(i) = c0_ec
    i=i+1
    par_vec(i) = op_cost_es
    i=i+1
    par_vec(i) = op_cost_ec
    i=i+1
    par_vec(i) = eps_super
    i=i+1
    par_vec(i) = prob_super
    i=i+1
    par_vec(i) = prob_super_back
     i=i+1
    par_vec(i) = switch_cost_CP
     i=i+1
    par_vec(i) = switch_cost_PC
    
    if (size(par_vec)/=i) then
    	write(*,*) "Size of par_vec = ",size(par_vec)
    	write(*,*) "i = ",i
    	call myerror("Error in pack_params")
    endif 
    
    end subroutine pack_params
!====================================================!
    
    subroutine unpack_params(par_vec)
    !! unpack_params converts parameter vector par_vec
    ! into global variables
    implicit none
    ! Declare variables
    real(8), intent(in) :: par_vec(:)
    integer :: i
        
    ! Order of parameters has to be the same as in pack_params
    
    ! Execution
    !i=1
    !p_offer_a(1)=par_vec(i)
    i=1
    beta  = par_vec(i)         ! beta
    i=i+1
    r_glob  = par_vec(i)       ! global interest rate
    i=i+1
    tau_p  = par_vec(i)         ! payroll tax rate
    i=i+1
    hsv_0 = par_vec(i)         ! lambda_hsv
    i=i+1
    hsv_1 = par_vec(i)         ! tau_hsv
    i=i+1
    chi = par_vec(i)            ! marginal utility
    i=i+1
    vi = par_vec(i)            ! span of control
    i=i+1
    gamma = par_vec(i)         ! capital share in the entrepreneurial sector
    i=i+1
    rho_theta = par_vec(i)     ! Persistence theta shock
    i=i+1 
    sig_theta = par_vec(i)     ! Dispersion theta shock 
    i=i+1
    uncmean_theta = par_vec(i) ! Unconditional mean of theta shock
    i=i+1
    lambda = par_vec(i)        ! Collateral constraint of EP
    i=i+1
    lambda_es = par_vec(i)     ! Collateral constraint of ES
    i=i+1
    lambda_ec = par_vec(i)     ! Collateral constraint of EC
    i=i+1
    c0_es = par_vec(i)         ! Tax avoidance cost ES
    i=i+1
    c0_ec = par_vec(i)         ! Tax avoidance cost EC
    i=i+1
    op_cost_es = par_vec(i)    ! Operating cost ES
    i=i+1
    op_cost_ec = par_vec(i)    ! Operating cost EC
    i=i+1
    eps_super = par_vec(i)     ! Super-star shock for eps, level
    i=i+1
    prob_super = par_vec(i)    ! Super-star shock for eps, inflow rate
    i=i+1
    prob_super_back = par_vec(i)  ! Super-star shock for eps, outflow rate
    i=i+1
    switch_cost_CP = par_vec(i)  ! Switch cost C->Passthrough
    i=i+1
    switch_cost_PC = par_vec(i)  ! Switch cost Passthrough->C
    
	if (size(par_vec)/=i) then
    	call myerror("Error in unpack_params")
    endif  
 
      
    end subroutine unpack_params

!====================================================!


end module mod_globals
