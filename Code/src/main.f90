!##############################################################################
! PROGRAM MAIN
!##############################################################################

program main

use mod_globals
use mod_initialize
use mod_steady_state, only:  compute_steady_state
use mod_transition
use mod_welfare, only: compute_value_cons, compute_cev,compute_cev_agg,compute_cev_groups,decomp_cev,decomp_cev_tran
use mod_numerical,    only: myerror,linspace
use mod_utilities
use omp_lib
implicit none

!============================================================================
! DECLARATION 
!============================================================================
integer, parameter   :: iMaxThreads = 14
integer              :: i,istat,flag_ss,uno,T
real(8)              :: tau_h_init,tau_h_final
character(len=256)   :: pwd, dir_name
logical              :: flag_root
character(len=256)   :: dir_comptran
logical              :: dirExists,valc_exists
integer              :: do_write 
integer   		     :: n_tau_h      
real(8), parameter   :: max_tau_h = 0.6d0 !This is only for comparative statics do_RUN=2 or 3
character(len=:), allocatable :: file_params 

real(8)            :: tau_p_init,tau_p_final,tau_p0,hsv_0_0,hsv_0_init,hsv_0_final
type(modelPrices)  :: prices0,prices_init,prices_final
type(valfun)       :: val0,val_init,val_final,val_t1,valc_init,cev
type(polfun)       :: pol,pol_init,pol_final
type(polfun),allocatable :: pol_path(:)
type(distribfun_m) :: distrib0,distrib_init,distrib_final
type(modelAgg)     :: agg,agg_init,agg_final
type(modelAgg),allocatable :: agg_path(:)
type(modelResults) :: model_targets
integer            :: ind_max,do_GE_save,do_compute_valc,do_CEV_decomp
integer, parameter :: nq = 4 ! number of wealth groups in welfare analysis
real(8) :: cev_agg,cev_aggcomp_agg,cev_distcomp_agg,cev_q(nq),cev_qo(nq,2),cev_z(4)
real(8),allocatable :: tau_h_path(:),KN_path(:),hsv_0_path(:),tau_p_path(:)
real(8),allocatable,dimension(:) :: tau_h_vec,tau_p0_vec,taxes_tot_vec, taxes_ss_vec, &
			taxes_inc_vec,taxes_corp_vec,taxes_div_vec, res_gov_vec,  share_toptax_vec, Y_vec,Y_entre_vec,  &
			r_vec,w_vec,share_entre_act_vec,share_EP_entre_vec,share_ES_entre_vec, &
			share_EC_entre_vec, gini_wealth_all_vec, gini_wealth_entre_vec, &
			gini_inc_all_vec, gini_inc_entre_vec, K_entre_vec,K_C_vec,N_entre_vec,N_C_vec, &
            C_vec,share_wage_ES_vec,share_wage_EC_vec,hsv_0_vec,cev_vec,cev_aggcomp_vec,cev_distcomp_vec,&
			y_H_vec,ss_cap_vec,taxes_r_vec,labshare_ES_vec,labshare_EC_vec,labshare_corp_vec,labshare_allcorp_vec,inc_share_EP_vec,inc_share_ES_vec,inc_share_EC_vec, &
			netinc_share_EP_vec,netinc_share_ES_vec,netinc_share_EC_vec
real(8),allocatable,dimension(:,:) :: wealth_share_vec,avek_vec,avetheta_vec,inc_share_vec,taxes_top_vec, &
			share_lfo_top1,share_lfo_top5,share_lfo_top10,taxes_occ_vec,Y_lfo_vec
real(8),allocatable,dimension(:,:,:) :: trans_mat_vec
integer, ALLOCATABLE :: flag_ss_vec(:)

!============================================================================
! SET SOME FLAGS HERE 
!============================================================================
cf_avoidance   = 0		! 0 = baseline environment with avoidance
						! 1 = counterfactual: no income shifting
						! 2 = counterfactual: no tax avoidance incentives (everyone pays income + social security tax) 
						! 3 = counterfactual: no S- or C-corporations.
						! 4 = no C-corporations. 
cf_occpol      = 0 	    ! 0 = occpol choice is flexible (for calibration etc, choose this)
						! 1 = fix occpol to the benchmark 
do_run         = 0		! 0 = evaluate ss once, 
						! 1 = Model validation: Compare tau_h = 0.28 vs tau_h = 0.5
						! 2 = comparative statics, 
						! 3 = fiscal neutral welfare with transition (loop over tau_h) 
                        ! 4 = transition (once, at welfare max tau_h), 
						! 5 = tax avoidance welfare (with transition, for no-avoidance experiments)
n_tau_h        = 52     ! If do_run = 2 or 3, set n_tau_h to 52. If do_run = 1, set n_tau_h to 3.
do_fix_ss_cap  = 0		! 0 = ss_cap is endogenous, 1 = ss_cap is fixed at the benchmark level
do_fix_pen     = 0		! 0 = pension replacement rate is fixed, but pension level is endogeneous, 1 = pension level is fixed (so that replacement rate is endogeneous). Only relevant if do_benchmark = .false.
do_compute_valc = 0		! 1 = compute valc of benchmark (only relevant if do_run = 0 and do_benchmark = .true.) 
do_GE          = 1     ! Option only relevant for steady state or comparative statics. (do_run = 0 or 2)
						! 0 =PE
						! 1 =GE r and tau_p   (for baseline model calibration)
					    ! 2 =GE with r, tau_p and hsv_0 (for counterfactuals with fiscal neutrality)
						! 3 =PE with fiscal neutrality (r is fixed)
do_CEV_decomp  = 0		! 0 = do not decompose CEV, 1 = decompose CEV (slow)
do_benchmark   = .true. ! true = benchmark steady state, false = else.
						 ! do_benchmark should be .false., unless running steady state (do_run = 0)
						 ! in the benchmark model.
verbose        = 1       ! Flag 0/1/2 to control output on the screen in the whole code.
display_mu     = 0        ! Flag 0/1 to control output on the screen in the distrib.
display_vfi    = 0        ! Flag 0/1 to control output on the screen in the vfi.
display_howard = 0        ! Flag 0/1 to control output on the screen in vfi/howard
mydebug        = 0
do_write       = 1       ! 1= write results of comparative statics to txt

!============================================================================
! EXECUTION
!============================================================================
write(*,*) "Happy version v47!"

call OMP_SET_NUM_THREADS(iMaxThreads)

!$omp parallel if (par_fortran==1)
if (verbose>=1) write(*,*) "Parallel hello to you!"
!$omp end parallel

! Identify OS
! "#ifdef" is preprocessor directive. Must use "-fpp" in compilation. 
! No indention for preprocessor directive.
#ifdef __APPLE__
	myOS = 1
#endif
#ifdef _WIN64
	myOS = 2
#endif
    
if (myOS == 1) then
	! MAC
	slash    = '/'
	savedir  = 'output/'
	inputdir = 'input/'
elseif (myOS == 2) then
	! WINDOWS
	slash   = '\'
	savedir = 'output\'
	inputdir = 'input\'
else
	! UNIX and everything else
	slash   = '/'
    savedir='output/'
	inputdir = 'input/'
endif

file_params = "estim_params_benchmark.txt" ! File with parameters to be read
dir_bench   = trim(savedir)//trim("ss_bench")//trim(slash) ! Directory for benchmark results

if (verbose>=1) then
    call getcwd(pwd)
    write(*,*) "The current working directory is: ",trim(pwd)
	write(*,*) "Reading model inputs from dir  :",trim(inputdir)
	write(*,*) "Reading parameters from file   : ", file_params
	write(*,*) "Saving benchmark results to dir: ", trim(dir_bench)
    write(*,*) " "
endif

!============================================================================
! INITIALIZATION
!============================================================================

!Generates grids and other invariant model objects
write(*,*) "calling initialize_params..."
allocate(a_grid(na),eps_grid(neps),theta_grid(ntheta), &
	P_eps(neps,neps),P_theta(ntheta,ntheta),P_age(2,2),prob_theta(ntheta),prob_eps(neps),&
	switch_cost(nz,no),stat=istat)
if (istat/=0) then
	call myerror("main: Allocation failed!")
endif

if (verbose>=1) write(*,*) "Reading parameters from file.."
call read_params(trim(file_params))
if (verbose>=1) write(*,*) "Parameters read!"

call initialize_params()

!=============================================================
! SOLVE Steady-State MODEL or Transition
!=============================================================    

if (do_run == 0) then
	write(*,*) "============================================================="
	write(*,*) "SOLVING FOR STEADY STATE EQUILIBRIUM"
	write(*,*) "Do_GE = ",do_GE
	call initialize_prices(prices0)
	tau_p0 = tau_p
	hsv_0_0 = hsv_0
	write(*,*) "prices0%r = ",prices0%r  !Set equal to r_glob
	write(*,*) "tau_p0    = ",tau_p0
	write(*,*) "hsv_0     = ",hsv_0
	
	if (verbose>=1) write(*,*) "Prices initialization done!"

	if (cf_occpol == 1) then 
		allocate(occpol_bench(na,neps,ntheta,nz,no), stat=istat)
		if (istat/=0) then
			call myerror("occpol_bench: Allocation failed!")
		endif
		call readdim(occpol_bench,  trim(dir_bench)//trim('occpol.txt')) !(a,eps,theta,z_min,o)
	endif
	
	! Export results to external files 
	if (do_benchmark) then
		if (verbose>=1) then
			!pause
			write(*,*) " "
			write(*,*) "STARTING STEADY-STATE COMPUTATION FOR BENCHMARK"
			write(*,*) " "
		endif
		!pause 
		call compute_steady_state(tau_p0,hsv_0_0,prices0,model_targets,val0,pol,distrib0,agg,flag_ss,dir_bench)
		write(*,*) "Writing objects to separate text files.."
		call sub_write(trim(dir_bench),val0,pol,distrib0,prices0,agg,tau_p0,hsv_0_0)
		if (do_compute_valc==1) then 
			write(*,*) " "
			write(*,*) "COMPUTING valc, VALUE FROM CONSUMPTION IN THE BENCHMARK STEADY STATE"
			! We save a copy of valc in dir_bench
			call compute_value_cons(valc_init,trim(dir_bench) ) 
		endif

	else
		call compute_steady_state(tau_p0,hsv_0_0,prices0,model_targets,val0,pol,distrib0,agg,flag_ss)
		if ((do_GE == 2 .or. do_GE==3) .and. (cf_avoidance == 1 .or. cf_avoidance ==2 .or. cf_avoidance == 3 .or. cf_avoidance ==4)) then ! Fiscal neutrality tax avoidance
			! Save value functions to calculate welfare
			! Determine directory for output files
			if (do_GE == 2 .and. cf_avoidance==1 .and. cf_occpol == 0) then
				! Experiment 1: income shifting, flexible occpol
				dir_name = trim(savedir)//trim("exp1_fn")//trim(slash)
			elseif (do_GE == 2 .and. cf_avoidance==2 .and. cf_occpol == 0) then
				! Experiment 2: no tax avoidance at all, flexible occpol
				dir_name = trim(savedir)//trim("exp2_fn")//trim(slash)
			elseif (do_GE == 2 .and. cf_avoidance==1 .and. cf_occpol == 1) then
				! Experiment 3: income shifting, fixed occpol
				dir_name = trim(savedir)//trim("exp3_fn")//trim(slash)
			elseif (do_GE == 2 .and. cf_avoidance==2 .and. cf_occpol == 1) then
				! Experiment 4: no tax avoidance at all, fixed occpol
				dir_name = trim(savedir)//trim("exp4_fn")//trim(slash)
			elseif (do_GE == 2 .and. cf_avoidance==3 .and. cf_occpol == 0) then
				! Experiment 5: no S- or C-corporations, flexible occpol
				dir_name = trim(savedir)//trim("exp5_fn")//trim(slash)
			elseif (do_GE == 2 .and. cf_avoidance==4 .and. cf_occpol == 0) then
				! Experiment 6: no C-corporations, S-corp income shifting allowed
				dir_name = trim(savedir)//trim("exp6_fn")//trim(slash)
			endif
			call sub_write(trim(dir_name),val0,pol,distrib0,prices0,agg,tau_p0,hsv_0_0)
		endif ! Fiscal neutrality tax avoidance
	endif ! do_benchmark
elseif (do_run == 1) then
	! Model validation: Compare tau_h = 0.28 vs tau_h = 0.5
	write(*,*) "============================================================="
	write(*,*) "Comparing tau_h = 0.28 vs tau_h = 0.5 (steady state with fiscal neutrality)"
	if (cf_avoidance   /= 0) then
		call myerror("do_run = 1: cf_avoidance must be set to 0!")
	endif
	do_GE = 2
	n_tau_h = 3
	! Initial tax rate
	tau_h_init = tau_h
	
	write(*,*) "============================================================="
	write(*,*) "Comparative Statics"
	write(*,*) "do_GE = ",do_GE
    
    dir_name = trim(savedir)//trim("validation")//trim(slash)
   
	allocate(tau_h_vec(n_tau_h),tau_p0_vec(n_tau_h),hsv_0_vec(n_tau_h),taxes_tot_vec(n_tau_h), taxes_ss_vec(n_tau_h), &
			taxes_inc_vec(n_tau_h),taxes_corp_vec(n_tau_h),taxes_div_vec(n_tau_h), res_gov_vec(n_tau_h), share_toptax_vec(n_tau_h),  Y_vec(n_tau_h),  &
			r_vec(n_tau_h),w_vec(n_tau_h),share_entre_act_vec(n_tau_h),share_EP_entre_vec(n_tau_h),share_ES_entre_vec(n_tau_h), &
			share_EC_entre_vec(n_tau_h), gini_wealth_all_vec(n_tau_h), gini_wealth_entre_vec(n_tau_h), &
			gini_inc_all_vec(n_tau_h), gini_inc_entre_vec(n_tau_h), K_entre_vec(n_tau_h),K_C_vec(n_tau_h),N_entre_vec(n_tau_h),N_C_vec(n_tau_h), &
            C_vec(n_tau_h),share_wage_ES_vec(n_tau_h),share_wage_EC_vec(n_tau_h),cev_vec(n_tau_h),cev_aggcomp_vec(n_tau_h),cev_distcomp_vec(n_tau_h),&
			ss_cap_vec(n_tau_h),flag_ss_vec(n_tau_h),taxes_r_vec(n_tau_h),labshare_ES_vec(n_tau_h),labshare_EC_vec(n_tau_h),labshare_corp_vec(n_tau_h),labshare_allcorp_vec(n_tau_h), &
			inc_share_EP_vec(n_tau_h),inc_share_ES_vec(n_tau_h),inc_share_EC_vec(n_tau_h),&
			netinc_share_EP_vec(n_tau_h),netinc_share_ES_vec(n_tau_h),netinc_share_EC_vec(n_tau_h))
	allocate(wealth_share_vec(n_tau_h,4),avek_vec(n_tau_h,4),avetheta_vec(n_tau_h,4), &
		inc_share_vec(n_tau_h,4),taxes_top_vec(n_tau_h,4),share_lfo_top1(n_tau_h,3), &
		share_lfo_top5(n_tau_h,3),share_lfo_top10(n_tau_h,3),taxes_occ_vec(n_tau_h,no),&
		Y_lfo_vec(n_tau_h,no),trans_mat_vec(n_tau_h,4,4))
	! Compute initial steady state (in GE)
	write(*,*) " "
	write(*,*) "COMPUTING INITIAL STEADY STATE..."
	tau_h = tau_h_init
	do_GE_save = do_GE ! Save initial do_GE setting
    do_GE          = 1 ! No fiscal requirement in benchmark model
    do_benchmark   = .true.
	call initialize_prices(prices_init)
	write(*,*) "tau_h",tau_h
	tau_p_init = tau_p
	hsv_0_init = hsv_0
	call compute_steady_state(tau_p_init,hsv_0_init,prices_init,model_targets,val_init,pol_init,distrib_init,agg_init,flag_ss)
	! Export results to external files
	!call sub_write(trim(savedir)//trim("ss_init")//trim(slash),val_init,pol_init,distrib_init,prices_init,agg_init,tau_p_init)
	!call sub_write(trim(dir_bench),val_init,pol_init,distrib_init,prices_init,agg_init,tau_p_init)

	! Update do_GE
	do_GE = do_GE_save
	do_benchmark = .false.
	call initialize_params() ! Reset G_bench
	! Construct tax vector
    tau_h_vec = [0.28d0,tau_h,0.5d0]

	call initialize_prices(prices0)
    
	do i = 1,n_tau_h
		
		! Update tau_h
		tau_h = tau_h_vec(i)
		! Reset initial prices
		call initialize_prices(prices0)

		write(*,*) "============================================================="
		write(*,*) "Computing steady state, tau_h = ", tau_h
		write(*,*) i," out of ",n_tau_h
        write(*,*) "do_GE = ", do_GE
		write(*,*) "============================================================="
        
		tau_p0 = tau_p
		hsv_0_0 = hsv_0
		
		!write(*,*) "Initial condition for r = ", prices0%r
		!pause 
		call compute_steady_state(tau_p0,hsv_0_0,prices0,model_targets,val0,pol,distrib0,agg,flag_ss)
		
		tau_p0_vec(i)    = tau_p0
		hsv_0_vec(i)     = hsv_0_0
		taxes_tot_vec(i) =  agg%taxes_tot
		taxes_ss_vec(i)  = agg%taxes_ss
		taxes_inc_vec(i) = agg%taxes_inc
		taxes_corp_vec(i)= agg%taxes_corp
		taxes_div_vec(i) = agg%taxes_div
        taxes_occ_vec(i,:) = agg%taxes_occ ! (no)
        taxes_r_vec(i) = agg%taxes_r ! scalar
		labshare_ES_vec(i) = agg%labshare_ES
		labshare_EC_vec(i) = agg%labshare_EC
		labshare_corp_vec(i) = agg%labshare_corp
		labshare_allcorp_vec(i) = agg%labshare_allcorp
		res_gov_vec(i)= agg%res_gov
		share_toptax_vec(i) = model_targets%share_toptax 
		Y_vec(i)= agg%Y
		r_vec(i)= prices0%r
		w_vec(i)= prices0%w
		ss_cap_vec(i) = prices0%ss_cap
		share_entre_act_vec(i)= model_targets%share_entre_act
		share_EP_entre_vec(i)= model_targets%share_EP_entre
		share_ES_entre_vec(i)= model_targets%share_ES_entre
		share_EC_entre_vec(i)= model_targets%share_EC_entre
		inc_share_EP_vec(i)= model_targets%inc_share_EP
		inc_share_ES_vec(i)= model_targets%inc_share_ES
		inc_share_EC_vec(i)= model_targets%inc_share_EC
		netinc_share_EP_vec(i)= model_targets%netinc_share_EP
		netinc_share_ES_vec(i)= model_targets%netinc_share_ES
		netinc_share_EC_vec(i)= model_targets%netinc_share_EC
		gini_wealth_all_vec(i)= model_targets%gini_wealth_all
		gini_wealth_entre_vec(i)= model_targets%gini_wealth_entre
		gini_inc_all_vec(i)= model_targets%gini_inc_all
		gini_inc_entre_vec(i)=  model_targets%gini_inc_entre
		K_entre_vec(i)= agg%K_entre
		K_C_vec(i)= agg%K_C
		N_entre_vec(i)= agg%N_entre
		N_C_vec(i)= agg%N_C     
		C_vec(i)= agg%C
		share_wage_ES_vec(i)= model_targets%share_wage_ES
		share_wage_EC_vec(i)= model_targets%share_wage_EC
		flag_ss_vec(i) = flag_ss
		
		wealth_share_vec(i,:)=model_targets%wealth_share
		avek_vec(i,:)=agg%avek
		avetheta_vec(i,:)=agg%avetheta
		inc_share_vec(i,:)=model_targets%inc_share
		taxes_top_vec(i,:)=model_targets%taxes_top
		share_lfo_top1(i,:) = model_targets%share_lfo_top(2:4,1)
		share_lfo_top5(i,:) = model_targets%share_lfo_top(2:4,2)
		share_lfo_top10(i,:) = model_targets%share_lfo_top(2:4,3)
		trans_mat_vec(i,:,:) = model_targets%trans_mat
		Y_lfo_vec(i,:) = agg%Y_lfo

	enddo !end tau_h
	
	!Export
	write(*,*) "Validation: Exporting comparative static results to external files.."
	call writedim(tau_h_vec,trim(dir_name)//trim('tau_h_vec.txt'))
	call writedim(tau_p0_vec,trim(dir_name)//trim('tau_p0_vec.txt'))
	call writedim(hsv_0_vec,trim(dir_name)//trim('hsv_0_vec.txt'))
	call writedim(taxes_tot_vec,trim(dir_name)//trim('taxes_tot_vec.txt'))
	call writedim(taxes_ss_vec,trim(dir_name)//trim('taxes_ss_vec.txt'))
	call writedim(taxes_inc_vec,trim(dir_name)//trim('taxes_inc_vec.txt'))
	call writedim(taxes_corp_vec,trim(dir_name)//trim('taxes_corp_vec.txt'))
	call writedim(taxes_div_vec,trim(dir_name)//trim('taxes_div_vec.txt'))
    call writedim(taxes_occ_vec,trim(dir_name)//trim('taxes_occ_vec.txt'))
    call writedim(taxes_r_vec,trim(dir_name)//trim('taxes_r_vec.txt'))
	call writedim(labshare_ES_vec,trim(dir_name)//trim('labshare_ES_vec.txt'))
	call writedim(labshare_EC_vec,trim(dir_name)//trim('labshare_EC_vec.txt'))
	call writedim(labshare_corp_vec,trim(dir_name)//trim('labshare_corp_vec.txt'))
	call writedim(labshare_allcorp_vec,trim(dir_name)//trim('labshare_allcorp_vec.txt'))
	call writedim(inc_share_EP_vec,trim(dir_name)//trim('inc_share_EP_vec.txt'))
	call writedim(inc_share_ES_vec,trim(dir_name)//trim('inc_share_ES_vec.txt'))
	call writedim(inc_share_EC_vec,trim(dir_name)//trim('inc_share_EC_vec.txt'))
	call writedim(netinc_share_EP_vec,trim(dir_name)//trim('netinc_share_EP_vec.txt'))
	call writedim(netinc_share_ES_vec,trim(dir_name)//trim('netinc_share_ES_vec.txt'))
	call writedim(netinc_share_EC_vec,trim(dir_name)//trim('netinc_share_EC_vec.txt'))
	call writedim(res_gov_vec,trim(dir_name)//trim('res_gov_vec.txt'))
	call writedim(share_toptax_vec,trim(dir_name)//trim('share_toptax_vec.txt'))
	call writedim(Y_vec,trim(dir_name)//trim('Y_vec.txt'))
	call writedim(r_vec,trim(dir_name)//trim('r_vec.txt'))
	call writedim(w_vec,trim(dir_name)//trim('w_vec.txt'))
	call writedim(ss_cap_vec,trim(dir_name)//trim('ss_cap_vec.txt'))
	call writedim(share_entre_act_vec,trim(dir_name)//trim('share_entre_act_vec.txt'))
	call writedim(share_EP_entre_vec,trim(dir_name)//trim('share_EP_entre_vec.txt'))
	call writedim(share_ES_entre_vec,trim(dir_name)//trim('share_ES_entre_vec.txt'))
	call writedim(share_EC_entre_vec,trim(dir_name)//trim('share_EC_entre_vec.txt'))
	call writedim(gini_wealth_all_vec,trim(dir_name)//trim('gini_wealth_all_vec.txt'))
	call writedim(gini_wealth_entre_vec,trim(dir_name)//trim('gini_wealth_entre_vec.txt'))
	call writedim(gini_inc_all_vec,trim(dir_name)//trim('gini_inc_all_vec.txt'))
	call writedim(gini_inc_entre_vec,trim(dir_name)//trim('gini_inc_entre_vec.txt'))
	call writedim(K_entre_vec,trim(dir_name)//trim('K_entre_vec.txt'))
	call writedim(K_C_vec,trim(dir_name)//trim('K_C_vec.txt'))
	call writedim(N_entre_vec,trim(dir_name)//trim('N_entre_vec.txt'))
	call writedim(N_C_vec,trim(dir_name)//trim('N_C_vec.txt'))
	call writedim(C_vec,trim(dir_name)//trim('C_vec.txt'))
	call writedim(share_wage_ES_vec,trim(dir_name)//trim('share_wage_ES_vec.txt'))
	call writedim(share_wage_EC_vec,trim(dir_name)//trim('share_wage_EC_vec.txt'))
	call writedim(flag_ss_vec,trim(dir_name)//trim('flag_ss_vec.txt'))
	
	call writedim(wealth_share_vec,trim(dir_name)//trim('wealth_share_vec.txt'))
	call writedim(avek_vec,trim(dir_name)//trim('avek_vec.txt'))
	call writedim(avetheta_vec,trim(dir_name)//trim('avetheta_vec.txt'))
	call writedim(inc_share_vec,trim(dir_name)//trim('inc_share_vec.txt'))
	call writedim(taxes_top_vec,trim(dir_name)//trim('taxes_top_vec.txt'))
	call writedim(share_lfo_top1,trim(dir_name)//trim('share_lfo_top1.txt'))
	call writedim(share_lfo_top5,trim(dir_name)//trim('share_lfo_top5.txt'))
	call writedim(share_lfo_top10,trim(dir_name)//trim('share_lfo_top10.txt'))
	call writedim(Y_lfo_vec,trim(dir_name)//trim('Y_lfo_vec.txt'))
	call writedim(trans_mat_vec,trim(dir_name)//trim('trans_mat_vec.txt'))



	write(*,*) "============================================================="
	write(*,*) "Validation: Transition from tau_h = 0.5 to tau_h=0.28"
	
	! Final tax rate
	tau_h_final = tau_h_vec(1) ! 
	! Initial tax rate
	tau_h_init = tau_h_vec(n_tau_h) ! 
	! Number of time periods for the transition
	T = 100
	!call read_params()

	! Compute initial steady state (in GE)
	write(*,*) " "
	write(*,*) "COMPUTING INITIAL STEADY STATE..."
	write(*,*) "tau_h = ", tau_h_init
	tau_h = tau_h_init
	do_GE          = 2
	do_benchmark   = .false.
	call initialize_params() ! Reset G_bench
	call initialize_prices(prices_init) 
	tau_p_init = tau_p
	hsv_0_init = hsv_0
	call compute_steady_state(tau_p_init,hsv_0_init,prices_init,model_targets,val_init,pol_init,distrib_init,agg_init,flag_ss)

	write(*,*) " "
	 write(*,*) "LOADING valc, VALUE FROM CONSUMPTION IN THE INITIAL STEADY STATE"
	 allocate(valc_init%Vy(na,neps,ntheta,nz),valc_init%Vr(na))
	call readdim(valc_init%Vy,trim(dir_bench)//trim('valc_Vy.txt'))
	call readdim(valc_init%Vr,trim(dir_bench)//trim('valc_Vr.txt')) 

	! Compute final steady state in fiscal neutrality GE
	do_GE = 2 ! fiscal neutrality GE
	do_benchmark   = .false.
	call initialize_params() ! Reset G_bench
	write(*,*) " "
	write(*,*) "COMPUTING FINAL STEADY STATE"
	write(*,*) "tau_h = ", tau_h_final
	! Set the global variable tau_h to the current tau_h_final
	tau_h = tau_h_final
	! Initialize prices_final and tau_p_final. They will be updated to equilibrium values in compute_steady_state.
	call initialize_prices(prices_final)
	tau_p_final = tau_p
	hsv_0_final = hsv_0
	call compute_steady_state(tau_p_final,hsv_0_final,prices_final,model_targets,val_final,pol_final,distrib_final,agg_final,flag_ss)

	write(*,*) " "
	write(*,*) "SOLVING FOR TRANSITION"

	! Permanent increase in tau_h
	allocate(tau_h_path(T),pol_path(T),agg_path(T))
	tau_h_path = tau_h_final

	! Initial guess for the "price" path
	allocate(KN_path(T),hsv_0_path(T),tau_p_path(T))
	KN_path(1:T/2)      = linspace(prices_init%KN, prices_final%KN, T/2)
	KN_path(T/2+1:T)    = prices_final%KN
	hsv_0_path(1:T/2)      = linspace(hsv_0_init, hsv_0_final, T/2)
	hsv_0_path(T/2+1:T)    = hsv_0_final
	tau_p_path(1:T/2)   = linspace(tau_p_init, tau_p_final, T/2)
	tau_p_path(T/2+1:T) = tau_p_final


	write(*,*) "KN_init,KN_final",prices_init%KN, prices_final%KN
	write(*,*) "r_init,r_final",prices_init%r, prices_final%r
	write(*,*) "hsv_0_init,hsv_0_final",hsv_0_init,hsv_0_final
	write(*,*) "tau_p_init,tau_p_final",tau_p_init,tau_p_final
	!pause

	! Compute transition
	call compute_transition(tau_h_path,KN_path,hsv_0_path,tau_p_path,pol_path,agg_path, &
		val_t1,val_final, distrib_init,agg_final,dir_name)


elseif (do_run == 2) then
	! Initial tax rate
	tau_h_init = tau_h
	
	write(*,*) "============================================================="
	write(*,*) "Comparative Statics"
	write(*,*) "do_GE = ",do_GE
    
    ! do_GE should be set properly at the top of this program.
    if (cf_avoidance==0 .and. do_GE == 1) then
    	! Doing comparative statics in benchmark model (with avoidance)
    	dir_name = trim(savedir)//trim("compstat")//trim(slash)
        write(*,*) "Saving results in following directory:"
        write(*,*) dir_name
        !pause
    elseif (cf_avoidance==0 .and. do_GE == 0) then
	    dir_name = trim(savedir)//trim("compstat_PE")//trim(slash)
	elseif (cf_avoidance==0 .and. do_GE == 2) then
	    dir_name = trim(savedir)//trim("compstat_fn")//trim(slash)
	elseif (cf_avoidance==0 .and. do_GE == 10) then
	    dir_name = trim(savedir)//trim("compstat_GEr")//trim(slash)
    elseif (cf_avoidance==2 .and. do_GE == 1) then
	    dir_name = trim(savedir)//trim("compstat_CF2")//trim(slash)
	elseif (cf_avoidance==2 .and. do_GE == 0) then
	    dir_name = trim(savedir)//trim("compstat_CF2_PE")//trim(slash)
	elseif (cf_avoidance==2 .and. do_GE == 10) then
	    dir_name = trim(savedir)//trim("compstat_CF2_GEr")//trim(slash)
	elseif (cf_avoidance==2 .and. do_GE == 2) then
	    dir_name = trim(savedir)//trim("compstat_CF2_fn")//trim(slash)
	endif
	
	allocate(tau_h_vec(n_tau_h),tau_p0_vec(n_tau_h),hsv_0_vec(n_tau_h),taxes_tot_vec(n_tau_h), taxes_ss_vec(n_tau_h), &
			taxes_inc_vec(n_tau_h),taxes_corp_vec(n_tau_h),taxes_div_vec(n_tau_h), res_gov_vec(n_tau_h), share_toptax_vec(n_tau_h),  Y_vec(n_tau_h),  &
			 Y_entre_vec(n_tau_h),r_vec(n_tau_h),w_vec(n_tau_h),share_entre_act_vec(n_tau_h),share_EP_entre_vec(n_tau_h),share_ES_entre_vec(n_tau_h), &
			share_EC_entre_vec(n_tau_h), gini_wealth_all_vec(n_tau_h), gini_wealth_entre_vec(n_tau_h), &
			gini_inc_all_vec(n_tau_h), gini_inc_entre_vec(n_tau_h), K_entre_vec(n_tau_h),K_C_vec(n_tau_h),N_entre_vec(n_tau_h),N_C_vec(n_tau_h), &
            C_vec(n_tau_h),share_wage_ES_vec(n_tau_h),share_wage_EC_vec(n_tau_h),cev_vec(n_tau_h),cev_aggcomp_vec(n_tau_h),cev_distcomp_vec(n_tau_h),&
			ss_cap_vec(n_tau_h),flag_ss_vec(n_tau_h),taxes_r_vec(n_tau_h))
	allocate(wealth_share_vec(n_tau_h,4),avek_vec(n_tau_h,4),avetheta_vec(n_tau_h,4), &
		inc_share_vec(n_tau_h,4),taxes_top_vec(n_tau_h,4),share_lfo_top1(n_tau_h,3), &
		share_lfo_top5(n_tau_h,3),share_lfo_top10(n_tau_h,3),taxes_occ_vec(n_tau_h,no))

	! Compute initial steady state (in GE)
	write(*,*) " "
	write(*,*) "COMPUTING INITIAL STEADY STATE..."
	tau_h = tau_h_init
	do_GE_save = do_GE ! Save initial do_GE setting
    do_GE          = 1 ! No fiscal requirement in benchmark model
    do_benchmark   = .true.
	call initialize_prices(prices_init)
	write(*,*) "tau_h",tau_h
	tau_p_init = tau_p
	hsv_0_init = hsv_0
	call compute_steady_state(tau_p_init,hsv_0_init,prices_init,model_targets,val_init,pol_init,distrib_init,agg_init,flag_ss)
	! Export results to external files
	call sub_write(trim(savedir)//trim("ss_init")//trim(slash),val_init,pol_init,distrib_init,prices_init,agg_init,tau_p_init,hsv_0_init)
	call sub_write(trim(dir_bench),val_init,pol_init,distrib_init,prices_init,agg_init,tau_p_init,hsv_0_init)

	! Update do_GE
	do_GE = do_GE_save
	do_benchmark = .false.
	call initialize_params() ! Reset G_bench, y_H_bench
	! Construct tax vector
	tau_h_vec = linspace(top_rate,max_tau_h,n_tau_h)

	write(*,*) "tau_h_vec = "
	write(*,*) (tau_h_vec(i), i=1,n_tau_h)
	!pause

	!tau_h_vec = linspace(0.32d0,0.52d0,n_tau_h)
	
    
    !Load value functions in the benchmark
    if (do_GE == 2) then ! fiscal neutral case
    		
		allocate(valc_init%Vy(na,neps,ntheta,nz),valc_init%Vr(na))
		call readdim(valc_init%Vy,trim(dir_bench)//trim('valc_Vy.txt'))
	 	call readdim(valc_init%Vr,trim(dir_bench)//trim('valc_Vr.txt'))

		! Compute value from consumption in the initial steady state
		! When we do the loop over tax rates, we compute valc only once
		!write(*,*) " "
		!write(*,*) "COMPUTING valc, VALUE FROM CONSUMPTION IN THE INITIAL STEADY STATE"
		!call compute_value_cons(valc_init,trim(savedir)//trim("ss_init")//trim(slash)) 
		!allocate(valc_init%Vy(na,neps,ntheta,nz),valc_init%Vr(na))
	 	!call readdim(valc_init%Vy,trim(savedir)//trim("ss_init")//trim(slash)//trim('valc_Vy.txt'))
	 	!call readdim(valc_init%Vr,trim(savedir)//trim("ss_init")//trim(slash)//trim('valc_Vr.txt')) 
    endif 
    
	do i = 1,n_tau_h
		! Update tau_h
		tau_h = tau_h_vec(i)
		! Reset initial prices
		call initialize_prices(prices0)

		write(*,*) "============================================================="
		write(*,*) "Computing steady state, tau_h = ", tau_h
		write(*,*) i," out of ",n_tau_h
        write(*,*) "do_GE = ", do_GE
		write(*,*) "============================================================="
        
		tau_p0 = tau_p
		hsv_0_0 = hsv_0
		
		!write(*,*) "Initial condition for r = ", prices0%r
		!pause 
		call compute_steady_state(tau_p0,hsv_0_0,prices0,model_targets,val0,pol,distrib0,agg,flag_ss)
		!write(*,*) "Equil. value for r = ", prices0%r

		tau_p0_vec(i)    = tau_p0
		hsv_0_vec(i)     = hsv_0_0
		taxes_tot_vec(i) =  agg%taxes_tot
		taxes_ss_vec(i)  = agg%taxes_ss
		taxes_inc_vec(i) = agg%taxes_inc
		taxes_corp_vec(i)= agg%taxes_corp
		taxes_div_vec(i) = agg%taxes_div
        taxes_occ_vec(i,:) = agg%taxes_occ ! (no)
        taxes_r_vec(i) = agg%taxes_r ! scalar
		res_gov_vec(i)= agg%res_gov
		share_toptax_vec(i) = model_targets%share_toptax 
		Y_vec(i)= agg%Y
		Y_entre_vec(i)= agg%Y_entre
		r_vec(i)= prices0%r
		w_vec(i)= prices0%w
		ss_cap_vec(i) = prices0%ss_cap
		share_entre_act_vec(i)= model_targets%share_entre_act
		share_EP_entre_vec(i)= model_targets%share_EP_entre
		share_ES_entre_vec(i)= model_targets%share_ES_entre
		share_EC_entre_vec(i)= model_targets%share_EC_entre
		gini_wealth_all_vec(i)= model_targets%gini_wealth_all
		gini_wealth_entre_vec(i)= model_targets%gini_wealth_entre
		gini_inc_all_vec(i)= model_targets%gini_inc_all
		gini_inc_entre_vec(i)=  model_targets%gini_inc_entre
		K_entre_vec(i)= agg%K_entre
		K_C_vec(i)= agg%K_C
		N_entre_vec(i)= agg%N_entre
		N_C_vec(i)= agg%N_C     
		C_vec(i)= agg%C
		share_wage_ES_vec(i)= model_targets%share_wage_ES
		share_wage_EC_vec(i)= model_targets%share_wage_EC
		flag_ss_vec(i) = flag_ss
		
		wealth_share_vec(i,:)=model_targets%wealth_share
		avek_vec(i,:)=agg%avek
		avetheta_vec(i,:)=agg%avetheta
		inc_share_vec(i,:)=model_targets%inc_share
		taxes_top_vec(i,:)=model_targets%taxes_top
		share_lfo_top1(i,:) = model_targets%share_lfo_top(2:4,1)
		share_lfo_top5(i,:) = model_targets%share_lfo_top(2:4,2)
		share_lfo_top10(i,:) = model_targets%share_lfo_top(2:4,3)
		
		! CEV (only in the fiscal neutrality case)
		if (do_GE == 2) then
			if (do_CEV_decomp == 1) then
				call decomp_cev(cev_vec(i),cev_aggcomp_vec(i),cev_distcomp_vec(i), agg_init%C,agg_init%labsup,pol_init,val_init,valc_init,distrib_init,&
					agg%C,agg%labsup,pol,val0)
			else
				cev_vec(i) = compute_cev_agg(val_init,val0,valc_init,distrib_init) 
			endif
		endif
	enddo !end tau_h
	
	!Export
	if (do_write==1) then
	call writedim(tau_h_vec,trim(dir_name)//trim('tau_h_vec.txt'))
	call writedim(tau_p0_vec,trim(dir_name)//trim('tau_p0_vec.txt'))
	call writedim(hsv_0_vec,trim(dir_name)//trim('hsv_0_vec.txt'))
	call writedim(taxes_tot_vec,trim(dir_name)//trim('taxes_tot_vec.txt'))
	call writedim(taxes_ss_vec,trim(dir_name)//trim('taxes_ss_vec.txt'))
	call writedim(taxes_inc_vec,trim(dir_name)//trim('taxes_inc_vec.txt'))
	call writedim(taxes_corp_vec,trim(dir_name)//trim('taxes_corp_vec.txt'))
	call writedim(taxes_div_vec,trim(dir_name)//trim('taxes_div_vec.txt'))
    call writedim(taxes_occ_vec,trim(dir_name)//trim('taxes_occ_vec.txt'))
    call writedim(taxes_r_vec,trim(dir_name)//trim('taxes_r_vec.txt'))
	call writedim(res_gov_vec,trim(dir_name)//trim('res_gov_vec.txt'))
	call writedim(share_toptax_vec,trim(dir_name)//trim('share_toptax_vec.txt'))
	call writedim(Y_vec,trim(dir_name)//trim('Y_vec.txt'))
	call writedim(Y_entre_vec,trim(dir_name)//trim('Y_entre_vec.txt'))
	call writedim(r_vec,trim(dir_name)//trim('r_vec.txt'))
	call writedim(w_vec,trim(dir_name)//trim('w_vec.txt'))
	call writedim(ss_cap_vec,trim(dir_name)//trim('ss_cap_vec.txt'))
	call writedim(share_entre_act_vec,trim(dir_name)//trim('share_entre_act_vec.txt'))
	call writedim(share_EP_entre_vec,trim(dir_name)//trim('share_EP_entre_vec.txt'))
	call writedim(share_ES_entre_vec,trim(dir_name)//trim('share_ES_entre_vec.txt'))
	call writedim(share_EC_entre_vec,trim(dir_name)//trim('share_EC_entre_vec.txt'))
	call writedim(gini_wealth_all_vec,trim(dir_name)//trim('gini_wealth_all_vec.txt'))
	call writedim(gini_wealth_entre_vec,trim(dir_name)//trim('gini_wealth_entre_vec.txt'))
	call writedim(gini_inc_all_vec,trim(dir_name)//trim('gini_inc_all_vec.txt'))
	call writedim(gini_inc_entre_vec,trim(dir_name)//trim('gini_inc_entre_vec.txt'))
	call writedim(K_entre_vec,trim(dir_name)//trim('K_entre_vec.txt'))
	call writedim(K_C_vec,trim(dir_name)//trim('K_C_vec.txt'))
	call writedim(N_entre_vec,trim(dir_name)//trim('N_entre_vec.txt'))
	call writedim(N_C_vec,trim(dir_name)//trim('N_C_vec.txt'))
	call writedim(C_vec,trim(dir_name)//trim('C_vec.txt'))
	call writedim(share_wage_ES_vec,trim(dir_name)//trim('share_wage_ES_vec.txt'))
	call writedim(share_wage_EC_vec,trim(dir_name)//trim('share_wage_EC_vec.txt'))
	call writedim(flag_ss_vec,trim(dir_name)//trim('flag_ss_vec.txt'))
	
	call writedim(wealth_share_vec,trim(dir_name)//trim('wealth_share_vec.txt'))
	call writedim(avek_vec,trim(dir_name)//trim('avek_vec.txt'))
	call writedim(avetheta_vec,trim(dir_name)//trim('avetheta_vec.txt'))
	call writedim(inc_share_vec,trim(dir_name)//trim('inc_share_vec.txt'))
	call writedim(taxes_top_vec,trim(dir_name)//trim('taxes_top_vec.txt'))
	call writedim(share_lfo_top1,trim(dir_name)//trim('share_lfo_top1.txt'))
	call writedim(share_lfo_top5,trim(dir_name)//trim('share_lfo_top5.txt'))
	call writedim(share_lfo_top10,trim(dir_name)//trim('share_lfo_top10.txt'))
	endif ! IF do_write

	if (do_GE == 2) then ! Fiscal neutrality case, export welfare measure
		if (do_write==1) then
			call writedim(cev_vec,trim(dir_name)//trim('cev_vec.txt'))
			if (do_CEV_decomp == 1) then
				call writedim(cev_aggcomp_vec,trim(dir_name)//trim('cev_aggcomp_vec.txt'))
				call writedim(cev_distcomp_vec,trim(dir_name)//trim('cev_distcomp_vec.txt'))
			endif
		endif 
		! Calculate cev by groups for the optimal tau_h
		write(*,*) "Computing welfare by groups..."
		! Find the tau_h that maximizes cev on the grid
		ind_max     = maxloc(cev_vec,dim=1)
		tau_h = tau_h_vec(ind_max)
		! Recompute steady state at the optimal tau_h
		tau_p0 = tau_p
		hsv_0_0 = hsv_0
		call compute_steady_state(tau_p0,hsv_0_0,prices0,model_targets,val0,pol,distrib0,agg,flag_ss)
		! Write results (value functions, etc) to external files in dir_name
		if (do_write) call sub_write(trim(dir_name),val0,pol,distrib0,prices0,agg,tau_p0,hsv_0_0)
		! Compute CEV by groups 
		call sub_compute_cev(trim(dir_name))

	endif

    
elseif (do_run==3) then
    ! 3 = fiscal neutral welfare with transition (loop over tau_h) 
    write(*,*) "============================================================="
	write(*,*) "Optimal Tax Exercise with Transition: Loop over tau_h"
	write(*,*) "do_GE = ",do_GE

    ! do_GE should be set properly at the top of this program.
	dir_comptran = find_dir_comptran(cf_avoidance,do_GE,slash,savedir)

    ! Initial and final tax rate
	tau_h_init = tau_h ! The initial tau_h is set in initialize_params
	! Number of time periods for the transition
	T = 100

	! Compute initial steady state (in GE)
	write(*,*) "Computing initial steady state..."
	tau_h = tau_h_init
    do_GE          = 1
    do_benchmark   = .true.
	call initialize_prices(prices_init)
	write(*,*) "tau_h",tau_h
	tau_p_init = tau_p
	hsv_0_init = hsv_0
	
	call compute_steady_state(tau_p_init,hsv_0_init,prices_init,model_targets,val_init,pol_init,distrib_init,agg_init,flag_ss)
	
	! Export results to external files
	call sub_write(trim(savedir)//trim("ss_init")//trim(slash),val_init,pol_init,distrib_init,prices_init,agg_init,tau_p_init,hsv_0_init)
	!call sub_write(trim(dir_bench),val_init,pol_init,distrib_init,prices_init,agg_init,tau_p_init,hsv_0_init)
	
	write(*,*) "sub_write in do_run=3 done!"
	write(*,*) "dir_bench = ", dir_bench

	! Read value from consumption for the initial steady-state
	! The value from consumption does not depend on the tax rate
	allocate(valc_init%Vy(na,neps,ntheta,nz),valc_init%Vr(na))
	call readdim(valc_init%Vy,trim(dir_bench)//trim('valc_Vy.txt'))
	call readdim(valc_init%Vr,trim(dir_bench)//trim('valc_Vr.txt')) 

	! Loop over final tax rates
	do_GE = 2 ! fiscal neutrality GE
	do_benchmark   = .false.
	call initialize_params() ! Reset G_bench and y_H_bench
	!write(*,*) "initializing params done!"

	allocate(tau_h_path(T),KN_path(T),hsv_0_path(T),tau_p_path(T),pol_path(T),agg_path(T), &
		cev_vec(n_tau_h),cev_aggcomp_vec(n_tau_h),cev_distcomp_vec(n_tau_h),y_H_vec(n_tau_h))
	tau_h_vec = linspace(top_rate,max_tau_h,n_tau_h)
	
	! overwrite/create files to store outputs
	open(uno, file=trim(dir_comptran)//trim("cev_vec.txt"), status="replace", action="write")
	close(uno)
	open(uno, file=trim(dir_comptran)//trim("cev_aggcomp_vec.txt"),  status="replace", action="write")
	close(uno)
	open(uno, file=trim(dir_comptran)//trim("cev_distcomp_vec.txt"),  status="replace", action="write")
	close(uno)
	open(uno, file=trim(dir_comptran)//trim("tau_h_vec.txt"),  status="replace", action="write")
	close(uno)
	open(uno, file=trim(dir_comptran)//trim("y_H_vec.txt"),  status="replace", action="write")
	close(uno)

	do i = 1,n_tau_h
		tau_h_final = tau_h_vec(i)
		! Compute final steady state
		write(*,*) "============================================================="
		write(*,*) "Computing final steady state, tau_h = ", tau_h_final
		write(*,*) i," out of ",n_tau_h
		! Set the global variable tau_h to the current tau_h_final
		tau_h = tau_h_final
	
		! Initialize prices_final and tau_p_final. They will be updated to equilibrium values in compute_steady_state.
		call initialize_prices(prices_final)
		tau_p_final = tau_p
		hsv_0_final = hsv_0
		call compute_steady_state(tau_p_final,hsv_0_final,prices_final,model_targets,val_final,pol_final,distrib_final,agg_final,flag_ss)
		y_H_vec(i) = prices_final%y_H
		! Export results to external files
		!call sub_write(trim(savedir)//trim("ss_final")//trim(slash),val_final,pol_final,distrib_final,prices_final,agg_final,tau_p_final)
 
		write(*,*) "SOLVING FOR TRANSITION"
		! Load initial and final steady state
		!  initial steady state
		!call load_ss(trim(savedir)//trim("ss_init")//trim(slash),val_init,pol_init,distrib_init,prices_init,agg_init,tau_p_init)
		!  final steady state
		!call load_ss(trim(savedir)//trim("ss_final")//trim(slash),val_final,pol_final,distrib_final,prices_final,agg_final,tau_p_final)
	
		! Permanent increase in tau_h
		tau_h_path = tau_h_final
	
		! Initial guess for the "price" path
		KN_path(1:T/2)      = linspace(prices_init%KN, prices_final%KN, T/2)
		KN_path(T/2+1:T)    = prices_final%KN
		hsv_0_path(1:T/2)   = linspace(hsv_0_init, hsv_0_final, T/2)
		hsv_0_path(T/2+1:T) = hsv_0_final
		tau_p_path(1:T/2)   = linspace(tau_p_init, tau_p_final, T/2)
		tau_p_path(T/2+1:T) = tau_p_final
	
		write(*,*) "KN_init,KN_final",prices_init%KN, prices_final%KN
		write(*,*) "r_init,r_final",prices_init%r, prices_final%r
		write(*,*) "hsv_0_init, hsv_0_final",hsv_0_init, hsv_0_final
		write(*,*) "tau_p_init,tau_p_final",tau_p_init,tau_p_final
		!pause
		! Compute transition
		call compute_transition(tau_h_path,KN_path,hsv_0_path,tau_p_path,pol_path,agg_path, &
			val_t1,val_final, distrib_init,agg_final)
	
		! Compute CEV aggregate
		write(*,*) "============================================================="
		if (do_CEV_decomp==0) then
			write(*,*) "Calling compute_cev_agg for tau_h_vec(",i,")=", tau_h_final	
			cev_vec(i) = compute_cev_agg(val_init,val_t1,valc_init,distrib_init) 
			open(uno, file=trim(dir_comptran)//trim("cev_vec.txt"), status="old", position="append", action="write")
			write(uno,*) cev_vec(i)
			close(uno)
			open(uno, file=trim(dir_comptran)//trim("tau_h_vec.txt"), status="old", position="append", action="write")
			write(uno,*) tau_h_vec(i)
			close(uno)
			open(uno, file=trim(dir_comptran)//trim("y_H_vec.txt"), status="old", position="append", action="write")
			write(uno,*) y_H_vec(i)
			close(uno)
		elseif (do_CEV_decomp==1) then
			write(*,*) "Calling decomp_cev_tran for tau_h_vec(",i,")=", tau_h_final	
			call decomp_cev_tran(cev_vec(i),cev_aggcomp_vec(i),cev_distcomp_vec(i), &
				agg_init%C,agg_init%labsup, pol_init,val_init,valc_init,distrib_init, &
				val_t1,agg_final%C,agg_final%labsup,pol_final,val_final,pol_path,agg_path) 
			! write(*,*) "cev = ", cev_vec(i)
			! write(*,*) "cev_aggcomp = ", cev_aggcomp_vec(i)
			! write(*,*) "cev_distcomp = ", cev_distcomp_vec(i)
			! Write results to external files
			open(uno, file=trim(dir_comptran)//trim("cev_vec.txt"), status="old", position="append", action="write")
			write(uno,*) cev_vec(i)
			close(uno)
			open(uno, file=trim(dir_comptran)//trim("cev_aggcomp_vec.txt"), status="old", position="append", action="write")
			write(uno,*) cev_aggcomp_vec(i)
			close(uno)
			open(uno, file=trim(dir_comptran)//trim("cev_distcomp_vec.txt"), status="old", position="append", action="write")
			write(uno,*) cev_distcomp_vec(i)
			close(uno)
			open(uno, file=trim(dir_comptran)//trim("tau_h_vec.txt"), status="old", position="append", action="write")
			write(uno,*) tau_h_vec(i)
			close(uno)
			open(uno, file=trim(dir_comptran)//trim("y_H_vec.txt"), status="old", position="append", action="write")
			write(uno,*) y_H_vec(i)
			close(uno)
		endif ! do_CEV_decomp
! 		pause
	enddo ! i
elseif (do_run == 4) then
    ! Do transition (once, at welfare max tau_h)
     ! 4 = transition to the optimal tax tau_h (fiscal neutral) 
    write(*,*) "============================================================="
	write(*,*) "Transition to Optimal Tax with Fiscal Neutrality"

	! do_GE should be set properly at the top of this program.
	! Find the directory for the optimal tax
	dir_comptran = find_dir_comptran(cf_avoidance,do_GE,slash,savedir)

    ! Find the optimal tau_h 
	allocate(tau_h_vec(n_tau_h),cev_vec(n_tau_h),y_H_vec(n_tau_h))
	! Read vector cev_tracev_vecn_agg in folder "dir_comptran"
	call readdim(cev_vec,trim(dir_comptran)//trim("cev_vec.txt"))
	! Read tau_h_vec, the vector of top tax rates
	call readdim(tau_h_vec,trim(dir_comptran)//trim("tau_h_vec.txt"))
	! Find the tau_h that maximizes cev on the grid
	ind_max     = maxloc(cev_vec,dim=1)
	! Final tax rate
	tau_h_final = tau_h_vec(ind_max) ! optimal top marginal tax rate

	! Initial tax rate
	tau_h_init = top_rate ! Benchmark model tax rate
	! Number of time periods for the transition
	T = 100
	!call read_params()

	! Compute initial steady state (in GE)
	write(*,*) " "
	write(*,*) "COMPUTING INITIAL STEADY STATE..."
	tau_h = tau_h_init
    do_GE          = 1
    do_benchmark   = .true.
	call initialize_prices(prices_init)
	write(*,*) "tau_h",tau_h
	tau_p_init = tau_p
	hsv_0_init = hsv_0
	call compute_steady_state(tau_p_init,hsv_0_init,prices_init,model_targets,val_init,pol_init,distrib_init,agg_init,flag_ss)
	! Export results to external files
	call sub_write(trim(savedir)//trim("ss_init")//trim(slash),val_init,pol_init,distrib_init,prices_init,agg_init,tau_p_init,hsv_0_init)
	!call sub_write(trim(dir_bench),val_init,pol_init,distrib_init,prices_init,agg_init,tau_p_init,hsv_0_init)
			
	! Compute value from consumption in the initial steady state
	! When we do the loop over tax rates, we compute valc only once
	
	! Compute value from consumption in the initial steady state
	! When we do the loop over tax rates, we compute valc only once
	write(*,*) " "
	!write(*,*) "COMPUTING valc, VALUE FROM CONSUMPTION IN THE INITIAL STEADY STATE"
	!call compute_value_cons(valc_init,trim(savedir)//trim("ss_init")//trim(slash))
	! If we skip compute_value_cons, we need to load valc0 from file
	 write(*,*) "LOADING valc, VALUE FROM CONSUMPTION IN THE INITIAL STEADY STATE"
	 allocate(valc_init%Vy(na,neps,ntheta,nz),valc_init%Vr(na))
	call readdim(valc_init%Vy,trim(dir_bench)//trim('valc_Vy.txt'))
	call readdim(valc_init%Vr,trim(dir_bench)//trim('valc_Vr.txt')) 

	! Compute final steady state in fiscal neutrality GE
	do_GE = 2 ! fiscal neutrality GE
	do_benchmark   = .false.
	call initialize_params() ! Reset G_bench
	write(*,*) " "
	write(*,*) "COMPUTING FINAL STEADY STATE"
	write(*,*) "tau_h = ", tau_h_final
	! Set the global variable tau_h to the current tau_h_final
	tau_h = tau_h_final
	! Initialize prices_final and tau_p_final. They will be updated to equilibrium values in compute_steady_state.
	call initialize_prices(prices_final)
	tau_p_final = tau_p
	hsv_0_final = hsv_0
	call compute_steady_state(tau_p_final,hsv_0_final,prices_final,model_targets,val_final,pol_final,distrib_final,agg_final,flag_ss)
	
	
	! Export results to external files
	call sub_write(trim(savedir)//trim("ss_final")//trim(slash),val_final,pol_final,distrib_final,prices_final,agg_final,tau_p_final,hsv_0_final)
 	write(*,*) " "
	write(*,*) "SOLVING FOR TRANSITION"

	! Permanent increase in tau_h
	allocate(tau_h_path(T),pol_path(T),agg_path(T))
	tau_h_path = tau_h_final
	
	! Initial guess for the "price" path
	allocate(KN_path(T),hsv_0_path(T),tau_p_path(T))
	KN_path(1:T/2)      = linspace(prices_init%KN, prices_final%KN, T/2)
	KN_path(T/2+1:T)    = prices_final%KN
	hsv_0_path(1:T/2)      = linspace(hsv_0_init, hsv_0_final, T/2)
	hsv_0_path(T/2+1:T)    = hsv_0_final
	tau_p_path(1:T/2)   = linspace(tau_p_init, tau_p_final, T/2)
	tau_p_path(T/2+1:T) = tau_p_final
	
	
	write(*,*) "KN_init,KN_final",prices_init%KN, prices_final%KN
	write(*,*) "r_init,r_final",prices_init%r, prices_final%r
	write(*,*) "hsv_0_init,hsv_0_final",hsv_0_init,hsv_0_final
	write(*,*) "tau_p_init,tau_p_final",tau_p_init,tau_p_final
	!pause
	
	! Compute transition
	call compute_transition(tau_h_path,KN_path,hsv_0_path,tau_p_path,pol_path,agg_path, &
		val_t1,val_final, distrib_init,agg_final,dir_comptran)
    
    write(*,*) " "
	write(*,*) "COMPUTING CEV'S"
	! Aggregate CEV
	! write(*,*) "Calling decomp_cev_tran"
 	! call decomp_cev_tran(cev_agg,cev_aggcomp_agg,cev_distcomp_agg, &
	! 		agg_init%C,agg_init%labsup, pol_init,val_init,valc_init,distrib_init, &
	! 		agg_final%C,agg_final%labsup,pol_final,val_final,pol_path,agg_path) 

    ! Write cev_agg to a file
    ! call writescalar(cev_agg,trim(dir_comptran)//trim("cev_agg.txt"))
	! call writescalar(cev_aggcomp_agg,trim(dir_comptran)//trim("cev_aggcomp_agg.txt"))
	! call writescalar(cev_distcomp_agg,trim(dir_comptran)//trim("cev_distcomp_agg.txt"))
    ! CEV by groups 
    call compute_cev_groups(cev_q,cev_qo,cev_z, val_init,val_t1,valc_init,distrib_init,pol_init%occpol,nq) 
    call writedim(cev_q,trim(dir_comptran)//trim("cev_q.txt"))
    call writedim(cev_qo,trim(dir_comptran)//trim("cev_qo.txt"))
    call writedim(cev_z,trim(dir_comptran)//trim("cev_z.txt"))
elseif (do_run==5) then
	write(*,*) "============================================================="
	write(*,*) "WELFARE OF TAX AVOIDANCE with transitional dynamics"

	if (cf_avoidance /=0) then
		write(*,*) "cf_avoidance /=0, cf_avoidance = ", cf_avoidance
		call myerror('do_run=5: cf_avoidance must be 0 ')
	endif
	! Number of time periods for the transition
	T = 100
	allocate(tau_h_path(T),pol_path(T),agg_path(T))
	allocate(KN_path(T),hsv_0_path(T),tau_p_path(T))
	! Compute initial steady state (in do_GE = 1)
	write(*,*) " "
	write(*,*) "COMPUTING INITIAL STEADY STATE..."
	tau_h = top_rate
    do_GE          = 1
    do_benchmark   = .true.
	call initialize_prices(prices_init)
	write(*,*) "tau_h",tau_h
	tau_p_init = tau_p
	hsv_0_init = hsv_0
	call compute_steady_state(tau_p_init,hsv_0_init,prices_init,model_targets,val_init,pol_init,distrib_init,agg_init,flag_ss)
	! Export results to external files
	! call sub_write(trim(savedir)//trim("ss_init")//trim(slash),val_init,pol_init,distrib_init,prices_init,agg_init,tau_p_init,hsv_0_init)
	!call sub_write(trim(dir_bench),val_init,pol_init,distrib_init,prices_init,agg_init,tau_p_init,hsv_0_init)
	
	! If we skip compute_value_cons, we need to load valc0 from file
	 write(*,*) "LOADING valc, VALUE FROM CONSUMPTION IN THE INITIAL STEADY STATE"
	 allocate(valc_init%Vy(na,neps,ntheta,nz),valc_init%Vr(na))
	call readdim(valc_init%Vy,trim(dir_bench)//trim('valc_Vy.txt'))
	call readdim(valc_init%Vr,trim(dir_bench)//trim('valc_Vr.txt')) 

! Experiment 6
	write(*,*) "Experiment 6"
	dir_name = trim(savedir)//trim("exp6_fn")//trim(slash)
	! Compute final steady state in fiscal neutrality GE
	do_GE = 2 ! fiscal neutrality GE
	do_benchmark   = .false.
	call initialize_params() ! Reset G_bench
	write(*,*) " "
	write(*,*) "COMPUTING FINAL STEADY STATE"
	cf_avoidance = 4 
	cf_occpol = 0
	! Initialize prices_final and tau_p_final. They will be updated to equilibrium values in compute_steady_state.
	call initialize_prices(prices_final)
	tau_p_final = tau_p
	hsv_0_final = hsv_0
	call compute_steady_state(tau_p_final,hsv_0_final,prices_final,model_targets,val_final,pol_final,distrib_final,agg_final,flag_ss)
 	write(*,*) " "
	write(*,*) "SOLVING FOR TRANSITION"

	! 
	tau_h_path = tau_h
	
	! Initial guess for the "price" path
	KN_path(1:T/2)      = linspace(prices_init%KN, prices_final%KN, T/2)
	KN_path(T/2+1:T)    = prices_final%KN
	hsv_0_path(1:T/2)      = linspace(hsv_0_init, hsv_0_final, T/2)
	hsv_0_path(T/2+1:T)    = hsv_0_final
	tau_p_path(1:T/2)   = linspace(tau_p_init, tau_p_final, T/2)
	tau_p_path(T/2+1:T) = tau_p_final
	
	write(*,*) "KN_init,KN_final",prices_init%KN, prices_final%KN
	write(*,*) "r_init,r_final",prices_init%r, prices_final%r
	write(*,*) "hsv_0_init,hsv_0_final",hsv_0_init,hsv_0_final
	write(*,*) "tau_p_init,tau_p_final",tau_p_init,tau_p_final
	
	! Compute transition
	call compute_transition(tau_h_path,KN_path,hsv_0_path,tau_p_path,pol_path,agg_path, &
		val_t1,val_final, distrib_init,agg_final,dir_name)
    
    write(*,*) " "
	write(*,*) "COMPUTING CEV'S"
	cev_agg = compute_cev_agg(val_init,val_t1,valc_init,distrib_init)
    call compute_cev_groups(cev_q,cev_qo,cev_z, val_init,val_t1,valc_init,distrib_init,pol_init%occpol,nq) 
    ! Write cev results to file
    call writescalar(cev_agg,trim(dir_name)//trim("cev_agg.txt"))
	call writedim(cev_q,trim(dir_name)//trim("cev_q.txt"))
    call writedim(cev_qo,trim(dir_name)//trim("cev_qo.txt"))
    call writedim(cev_z,trim(dir_name)//trim("cev_z.txt"))

else
    
    call myerror('do_run not defined!')
    
endif ! do_run


write(*,*) "Program completed successfully!"

contains

subroutine sub_write(dir,val,pol,distrib,prices,agg,tau_p0,hsv_0_0)
    use mod_utilities 
	implicit none
	! Declare inputs
	character(len=*),intent(in) :: dir
	type(valfun),intent(in) :: val
	type(polfun),intent(in) :: pol
	type(distribfun_m),intent(in) :: distrib
	type(modelPrices),intent(in) :: prices
	type(modelAgg),intent(in) :: agg
	real(8),intent(in) :: tau_p0,hsv_0_0
    !Declare local variables
    integer :: constants(5)
    
    ! Write the outputs as txt files (to be imported into matlab) 
    write(*,*) 'Saving txt files in the following directory:'
    write(*,*) dir
    write(*,*) '...Done!'
    
    ! Export dimensions of arrays 
    constants(1) = na      ! assets
    constants(2) = neps    ! worker ability
    constants(3) = ntheta  ! entre ability
    constants(4) = nz      ! occupational states x young/old (W,EP,ES,EC,R)
	constants(5) = no      ! occupational states (W,EP,ES,EC)
	
    call writedim(constants,trim(dir)//'constants.txt')

    !Write interest rate, transfer, and tau_p
    call writescalar(prices%r,trim(dir)//trim('r.txt'))
    call writescalar(prices%KN,trim(dir)//trim('KN.txt'))
    call writescalar(prices%w,trim(dir)//trim('w.txt'))
    call writescalar(prices%ss_cap,trim(dir)//trim('ss_cap.txt'))
    call writescalar(prices%y_H,trim(dir)//trim('y_H.txt'))
    call writescalar(prices%pen,trim(dir)//trim('pen.txt'))
    call writescalar(tau_p0,trim(dir)//trim('tau_p.txt'))
    call writescalar(hsv_0_0,trim(dir)//trim('hsv_0.txt'))
    
    !Write government spending G (in levels)
    call writescalar(agg%G,trim(dir)//trim('agg_G.txt'))
    call writescalar(agg%A,trim(dir)//trim('agg_A.txt'))
    call writescalar(agg%C,trim(dir)//trim('agg_C.txt'))
    call writescalar(agg%Y,trim(dir)//trim('agg_Y.txt'))
	call writescalar(agg%labsup,trim(dir)//trim('labsup.txt'))
    call writescalar(agg%taxes_tot,trim(dir)//trim('taxes_tot.txt'))
    call writescalar(agg%taxes_inc,trim(dir)//trim('taxes_inc.txt'))
    call writescalar(agg%avek(1),trim(dir)//trim('avek_all.txt'))

    
    !Write grids and transition matrices
    call writedim(a_grid,trim(dir)//trim('a_grid.txt'))
    call writedim(eps_grid,trim(dir)//trim('eps_grid.txt'))
    call writedim(theta_grid,trim(dir)//trim('theta_grid.txt'))
    call writedim(P_eps,trim(dir)//trim('P_eps.txt'))
    call writedim(P_theta,trim(dir)//trim('P_theta.txt'))
    call writedim(P_age,trim(dir)//trim('P_age.txt'))
    call writedim(prob_theta,trim(dir)//trim('prob_theta.txt'))
    call writedim(prob_eps,trim(dir)//trim('prob_eps.txt'))
    
    !Write policy functions
    call writedim(pol%occpol,  trim(dir)//trim('occpol.txt')) !(a,eps,theta,o)
    call writedim(pol%kpol,    trim(dir)//trim('kpol.txt')) !(a,theta,o)
    call writedim(pol%npol,    trim(dir)//trim('npol.txt')) !(a,theta,o)
    call writedim(pol%phipol,  trim(dir)//trim('phipol.txt')) !(a,theta,o)
    
    call writedim(pol%apoly, trim(dir)//trim('apoly.txt')) !(a,eps,theta,o), o=W,EP,ES,EC
    call writedim(pol%cpoly, trim(dir)//trim('cpoly.txt')) !(a,eps,theta,o), o=W,EP,ES,EC
    call writedim(pol%apolr, trim(dir)//trim('apolr.txt')) !(a), retired
    call writedim(pol%cpolr, trim(dir)//trim('cpolr.txt')) !(a), retired
    
    call writedim(pol%lpoly_w, trim(dir)//trim('lpoly_w.txt')) !(a,eps,theta), workers
    
    !Write Distribution 
    call writedim(distrib%muy,trim(dir)//trim('muy.txt')) ! (a,eps,theta,z_min)
    call writedim(distrib%mur,trim(dir)//trim('mur.txt')) ! (a)
    
    !Write value function
    call writedim(val%Vy,trim(dir)//trim('Vy.txt')) ! (a,eps,theta,nz)
    call writedim(val%Vr,trim(dir)//trim('Vr.txt')) ! (a)

end subroutine sub_write
 !===============================================================  
    subroutine read_params(file_params)
    ! read_params reads parameters from a file into vector par_vec and
    ! then assigns global parameter
    
		implicit none
		character(len=*), intent(in) :: file_params
		real(8) :: par_vec(N_params)
	
		call readdim(par_vec,trim(inputdir)//trim(file_params))
		call unpack_params(par_vec)
    end subroutine read_params
 !===============================================================  
    subroutine load_ss(dir,val,pol,distrib,prices,agg,tau_p0,hsv_0_0)
    use mod_utilities 
	implicit none
	! Declare inputs/outputs
	character(len=*),intent(in) :: dir
	type(valfun),intent(out) :: val
	type(polfun),intent(out) :: pol
	type(distribfun_m),intent(out) :: distrib
	type(modelPrices),intent(out) :: prices
	type(modelAgg),intent(out) :: agg
	real(8),intent(out) :: tau_p0,hsv_0_0
	! Locals
	integer :: istat
	real(8) :: check
	
	
	! Read interest rate and transfer
	call readscalar(prices%r,trim(dir)//trim('r.txt'))
	call readscalar(tau_p0,trim(dir)//trim('tau_p.txt'))
	call readscalar(hsv_0_0,trim(dir)//trim('hsv_0.txt'))
	call fun_prices(prices)
	
	
	! Read aggregates 
    call readscalar(agg%A,trim(dir)//trim('agg_A.txt'))
    call readscalar(agg%C,trim(dir)//trim('agg_C.txt'))
    call readscalar(agg%Y,trim(dir)//trim('agg_Y.txt'))
    call readscalar(agg%taxes_tot,trim(dir)//trim('taxes_tot.txt'))
    call readscalar(agg%taxes_inc,trim(dir)//trim('taxes_inc.txt'))
    call readscalar(agg%avek(1),trim(dir)//trim('avek_all.txt'))

	
	! Read pol
	allocate(pol%kpol(na,ntheta,no),pol%npol(na,ntheta,no),&
			pol%phipol(na,ntheta,no),pol%cpoly(na,neps,ntheta,no),pol%cpolr(na), &
			pol%lpoly_w(na,neps,ntheta),pol%apoly(na,neps,ntheta,no),pol%apolr(na), &
			pol%occpol(na,neps,ntheta,nz,no),distrib%muy(na,neps,ntheta,nz),distrib%mur(na), &
			val%Vy(na,neps,ntheta,nz),val%Vr(na), stat=istat)
	if (istat/=0) then
		call myerror("load_ss: Allocation failed!")
	endif
	
    !Read policy functions
    call readdim(pol%occpol,  trim(dir)//trim('occpol.txt')) !(a,eps,theta,o)
    call readdim(pol%kpol,    trim(dir)//trim('kpol.txt')) !(a,theta,o)
    call readdim(pol%npol,    trim(dir)//trim('npol.txt')) !(a,theta,o)
    call readdim(pol%phipol,  trim(dir)//trim('phipol.txt')) !(a,theta,o)
    
    call readdim(pol%apoly, trim(dir)//trim('apoly.txt')) !(a,eps,theta,o), o=W,EP,ES,EC
    call readdim(pol%cpoly, trim(dir)//trim('cpoly.txt')) !(a,eps,theta,o), o=W,EP,ES,EC
    call readdim(pol%apolr, trim(dir)//trim('apolr.txt')) !(a), retired
    call readdim(pol%cpolr, trim(dir)//trim('cpolr.txt')) !(a), retired
    
    call readdim(pol%lpoly_w, trim(dir)//trim('lpoly_w.txt')) !(a,eps,theta), workers
    
    !Read Distribution
    call readdim(distrib%muy,trim(dir)//trim('muy.txt')) ! (a,eps,theta,z_min)
    call readdim(distrib%mur,trim(dir)//trim('mur.txt')) ! (a)
	
	!Read value function
    call readdim(val%Vy,trim(dir)//trim('Vy.txt')) ! (a,eps,theta)
    call readdim(val%Vr,trim(dir)//trim('Vr.txt')) ! (a)
    
    ! Check that the loaded distribution sums to one
    check = sum(distrib%muy) + sum(distrib%mur)
    if (abs(check-1.0d0)>0.0000001d0 ) then
    	write(*,*) "sum(distrib%muy) + sum(distrib%mur) = ",check
        call myerror("main:load_ss: sum(distrib%muy) + sum(distrib%mur) does NOT sum to one!")
    endif
	
    end subroutine load_ss
!=============================================================
    subroutine sub_compute_cev(dir)
    use mod_utilities 
	implicit none
	! Declare inputs/outputs
	character(len=*),intent(in) :: dir
	! Declare local variable
	type(valfun)       :: val_exp

	! Load from dir_cf: val
	call load_ss(dir,val_exp,pol_final,distrib_final,prices_final,agg_final,tau_p0,hsv_0_0)
	! Compute aggregate CEV: 
	cev_agg = compute_cev_agg(val_init,val_exp,valc_init,distrib_init)
    ! Write cev_agg to a file
    call writescalar(cev_agg,trim(dir)//trim("cev_agg.txt"))
    ! Compute CEV by groups 
    call compute_cev_groups(cev_q,cev_qo,cev_z, val_init,val_exp,valc_init,distrib_init,pol_init%occpol,nq) 
    call writedim(cev_q,trim(dir)//trim("cev_q.txt"))
    call writedim(cev_qo,trim(dir)//trim("cev_qo.txt"))
    call writedim(cev_z,trim(dir)//trim("cev_z.txt"))

	end subroutine sub_compute_cev
!=============================================================

function find_dir_comptran(cf_avoidance,do_GE,slash,savedir)  result(dir_comptran)
	implicit none 
	! Declare inputs
	integer, intent(in) :: cf_avoidance,do_GE
	character(len=*), intent(in) :: slash, savedir
	! Declare function result
	character(len=256) :: dir_comptran

	if (cf_avoidance==0 .and. do_GE == 1) then
    	! Doing comparative statics in benchmark model (with avoidance)
    	dir_comptran = trim(savedir)//trim("comptran")//trim(slash)
	elseif (cf_avoidance==0 .and. do_GE == 2) then
	    dir_comptran = trim(savedir)//trim("comptran_fn")//trim(slash)
	elseif (cf_avoidance==0 .and. do_GE == 10) then
	    dir_comptran = trim(savedir)//trim("comptran_GEr")//trim(slash)
    elseif (cf_avoidance==2 .and. do_GE == 1) then
	    dir_comptran = trim(savedir)//trim("comptran_CF2")//trim(slash)
	elseif (cf_avoidance==2 .and. do_GE == 10) then
	    dir_comptran = trim(savedir)//trim("comptran_CF2_GEr")//trim(slash)
	elseif (cf_avoidance==2 .and. do_GE == 2) then
	    dir_comptran = trim(savedir)//trim("comptran_CF2_fn")//trim(slash)
	endif

end function find_dir_comptran

!=============================================================

end program main

