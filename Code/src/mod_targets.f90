module mod_targets
    use mod_globals
    use mod_functions
    use mod_numerical, only: myerror, my_ss,lrzcurve,outerprod,quantili, my_closest,locate,linint,unique,cumsum,quick_sort
!    use toolbox, only: plot, execplot, sort
    use mod_utilities
    use mod_distribution, only: sub_distrib_map
    implicit none
    
	! Variables for checking wealth and income distributions 
    real(8) :: mass_q(5),inc_mass_q(5),mass_a1,mass_aN,mass_a1_a3
    
    
contains
    
    !=================================================================!
    subroutine targets_compute(model_targets,prices,distrib_m,pol)
    implicit none
    ! Declare output and inputs 
    type(modelResults),intent(out) :: model_targets
    type(modelPrices),intent(in)   :: prices
    type(distribfun_m), intent(in)   :: distrib_m
    type(polfun),intent(in)        :: pol    

    ! Local utility variables
    integer :: it,ie,ia,io,j,istat,iq,iop,itp,iep
    real(8) :: t1,t2
    real(8) :: theta_val,eps_val,a_val,l_val,k_val,nbar_val,phi_val,profit_val,n_val,P_occ,yw,inc_temp
	type(distribfun) :: distrib
	! Useful scalars and arrays
	real(8), allocatable :: pregov_inc(:,:,:,:),pregov_incr(:),taxableinc(:,:,:,:),inctax(:,:,:,:), &
		inctaxr(:),taxableincr(:),net_inc(:,:,:,:),net_incr(:),payrollsize(:,:,:,:)
	real(8) :: opcost_vec(2:4),wealth_entre,wealth_all,inc_entre,inc_all
	real(8) :: inc_p(1),qw(1),r,w,pen,ss_cap,y_H
    real(8), allocatable :: mu_vec(:),inc_vec(:),mu(:,:,:,:),mur(:)
	! Marginal distributions and occ shares
	real(8) :: mu_a(na), mu_eps(neps), mu_theta(ntheta), mu_age(2)
	real(8) :: mu_age1(2),res_eps(neps),res_theta(ntheta),res_age(2)
	real(8) :: share_work,share_EP,share_ES,share_EC,share_ret,share_entre
	real(8) :: share_all,share_work_act,netinc_entre
	real(8) :: mean_eps, cond_mean_eps_work, aveinc_work,aveinc_entre,medinc_work,medinc_entre
	
	! Moments:	
	! - Occ/LFO distribution
	real(8) :: share_entre_act,share_EP_entre,share_ES_entre,share_EC_entre,share_active

	! - Share of income declared as wage (ES and EC)
	real(8) :: share_wage_ES, share_wage_EC,share_wage_qn_ES(3), share_wage_qn_EC(3),share_wage_qk_ES(3), share_wage_qk_EC(3),share_wage_EC_cf
    ! - Share facing financial constraints and leverages (EP,ES,EC)
    real(8) :: frac_bor_ep, frac_bor_es, frac_bor_ec
    real(8) :: leverage_ep,leverage_es,leverage_ec
    ! - Average income ratio, ave hours worked by workers
    real(8) :: ratio_aveinc_entre_worker, ave_l_work
    real(8) :: ratio_medinc_entre_worker
    ! - Share of employer
    real(8) :: share_empl_EP, share_empl_ES, share_empl_EC, share_empl_entre
    ! - Firm size (# employees) by LFO, and size relative to EP
    real(8) :: empsize_EP,empsize_ES,empsize_EC,empsize_rel_ES,empsize_rel_EC,payrollsize_lfo(4)
    ! - Share of income and wealth owned by entre
    real(8) :: share_inc_entre, share_wealth_entre
    ! - Share of net income by lfo
    real(8) :: netinc_share_EP, netinc_share_ES, netinc_share_EC
    ! - Taxes paid by top earners
    real(8) :: taxes_top(4),atr_vec(n_atr)
    ! - Inequality measures
    real(8) :: gini_wealth_all, gini_wealth_entre,wealth_share(4)
    real(8) :: gini_inc_all, gini_inc_entre,inc_share(4),share_toptax
    ! - Occ/LFO choice by income and wealth
    real(8) :: share_entre_top(3), share_lfo_top(2:4,3),share_entre_quint(5),share_lfo_quint(2:4,5)    
    real(8) :: share_entre_top_inc(3), share_lfo_top_inc(2:4,3),share_entre_quint_inc(5),share_lfo_quint_inc(2:4,5)    
	! - Employment share by firm size
    real(8) :: share_emp(4)
    real(8) :: median_inc_lfo(4),median_netinc_lfo(4),median_emp_lfo(4)
	! - Transition rates
    real(8) :: trans_entre_work,trans_work_entre,trans_CP,trans_PC,trans_SC,trans_mat(4,4)
    ! - Ave entrep productivity
    real(8) :: avetheta(4),avek(4),aveklratio(4),avekyratio(4)


    call CPU_TIME(t1)


	! Convert distrib_m to distrib
	allocate(distrib%muy(na,neps,ntheta,no),distrib%mur(na))
	call sub_distrib_map(distrib,distrib_m,pol%occpol)

	w = prices%w
	r = prices%r
	pen  = prices%pen
	ss_cap = prices%ss_cap
	y_H  = prices%y_H
    !-------------------------------------------------------------
    ! Marginal distributions
    !-------------------------------------------------------------
    
    !Given mu, compute marginal distributions
    mu_a = 0.0d0
    do ia = 1,na
    	mu_a(ia) = sum(distrib%muy(ia,:,:,:)) + distrib%mur(ia)
    enddo

    mu_eps = 0.0d0
    do ie = 1,neps
        mu_eps(ie) = sum(distrib%muy(:,ie,:,:))
    enddo
    mu_eps = mu_eps/sum(mu_eps)

    mu_theta = 0.0d0
    do it = 1,ntheta
        mu_theta(it) = sum(distrib%muy(:,:,it,:))
    enddo
    mu_theta = mu_theta/sum(mu_theta)

    ! Share of young and retirees
    mu_age(1) = sum(distrib%muy)
    mu_age(2) = sum(distrib%mur)
    mu_age1 = my_ss(P_age,2)
        
    !The residuals below must be close to 0
	res_eps = abs(mu_eps-prob_eps)
    if ( any(res_eps>1d-6) ) then
 		write(*,*) ([mu_eps(ie),prob_eps(ie)],ie=1,neps)
        call myerror('distrib. residual eps must be close to 0')
    endif
    res_theta = abs(mu_theta-prob_theta)
    if ( any(res_theta>1d-6) ) then
    	write(*,*) ([mu_theta(it),prob_theta(it)],it=1,ntheta)
        call myerror('distrib. residual theta must be close to 0')
    endif
    res_age = abs(mu_age-mu_age1)
    if ( any(res_age>1d-6) ) then
    	write(*,*) ([mu_age(ia),mu_age1(ia)],ia=1,2)
        call myerror('distrib. residual age must be close to 0')
    endif
    
    ! Define distribution of young agents over a, eps,theta, and occupation 
    mu = distrib%muy
    mur = distrib%mur
    !-------------------------------------------------------------
    ! Calculate pre-gov, taxable, net incomes, and payroll by lfo (pregov_inc, taxableinc, net_inc,payrollsize_EP/ES/EC)
    ! * Definition for pregov_inc: before taxes but after deductions
    ! * Definition for taxableinc: (gross) income subject to income taxation, or dividend taxation
    ! * Definition for inctax: income tax and dividend tax
    ! * Definition for net_inc: after tax income
    ! * Definition for pregov_incr: before taxes but after deductions for retirees
    ! * Definition for payroll: size of payroll in entreneurial business subject to payroll tax.
    !-------------------------------------------------------------
    
    allocate(pregov_inc(na,neps,ntheta,no),pregov_incr(na),taxableinc(na,neps,ntheta,no), &
    	taxableincr(na),inctax(na,neps,ntheta,no),payrollsize(na,neps,ntheta,no),inctaxr(na), &
        net_inc(na,neps,ntheta,no), net_incr(na),        stat=istat)
    if (istat/=0) then
        call myerror("mod_targets: allocation of pregov_inc etc. failed")
    endif
    pregov_inc = 0.0d0
    net_inc    = 0.0d0
    taxableinc = 0.0d0
    inctax     = 0.0d0
    payrollsize = 0.0d0
    do it = 1,ntheta
        theta_val  = theta_grid(it)           
        do ie = 1,neps
            eps_val    = eps_grid(ie)
	        do ia = 1,na
                a_val      = a_grid(ia)
                ! Young
                ! - Worker
                	io = 1
                    l_val    = pol%lpoly_w(ia,ie,it)
                    pregov_inc(ia,ie,it,io) = w*l_val*eps_val + r*a_val
                    taxableinc(ia,ie,it,io) = w*l_val*eps_val - fun_tax_pay(w*l_val*eps_val,ss_cap) + r*a_val
                    inctax(ia,ie,it,io)     = fun_tax_inc(taxableinc(ia,ie,it,io),y_H)
                    !net_inc(ia,ie,it,io)    = taxableinc(ia,ie,it,io) - inctax(ia,ie,it,io)
                ! - Sole Proprietor
                	io = 2
                    k_val = pol%kpol(ia,it,io)
                    n_val = pol%npol(ia,it,io)
                    profit_val = fun_profit(theta_val,k_val,le,n_val,r,w)
                	pregov_inc(ia,ie,it,io) = profit_val + r*a_val
                	taxableinc(ia,ie,it,io) = profit_val -fun_tax_pay(profit_val,ss_cap) + r*a_val
                	inctax(ia,ie,it,io)     = fun_tax_inc(taxableinc(ia,ie,it,io),y_H)
                	net_inc(ia,ie,it,io)    = profit_val
                    ! We use the SUSB payroll definition: for unincorporated businesses, payroll does not include officer compensation
                    payrollsize(ia,ie,it,io)= (w*n_val)
                ! - S-Corps
                	io = 3
                    k_val = pol%kpol(ia,it,io)
                    n_val = pol%npol(ia,it,io)
                    phi_val  = pol%phipol(ia,it,io)
                    profit_val = fun_profit(theta_val,k_val,le,n_val,r,w) &
                    	 -  cost_evasion_es(1.0d0-phi_val) - op_cost_es
                    pregov_inc(ia,ie,it,io) = profit_val + r*a_val
                    taxableinc(ia,ie,it,io) = phi_val*profit_val - fun_tax_pay(phi_val*profit_val,ss_cap) &
                    	 						+ (1.0d0-phi_val)*profit_val + r*a_val &
                    	 						-  cost_evasion_es(1.0d0-phi_val) - op_cost_es
                    inctax(ia,ie,it,io)     = fun_tax_inc(taxableinc(ia,ie,it,io),y_H)
                    net_inc(ia,ie,it,io)    = (1.0d0-phi_val)*profit_val
                    payrollsize(ia,ie,it,io)     = w*n_val+phi_val*profit_val
                ! - C-corps
                	io = 4
                	k_val = pol%kpol(ia,it,io)
                    n_val = pol%npol(ia,it,io)
                    phi_val  = pol%phipol(ia,it,io)
                    profit_val = fun_profit(theta_val,k_val,le,n_val,r,w) &
                    	 -  cost_evasion_ec(phi_val) - op_cost_ec
                    pregov_inc(ia,ie,it,io) = profit_val + r*a_val
                    ! inc_temp: income subject to the personal income tax
                    inc_temp = phi_val*profit_val - fun_tax_pay(phi_val*profit_val,ss_cap) &
                    	 						+ r*a_val -  cost_evasion_ec(phi_val) - op_cost_ec
         			taxableinc(ia,ie,it,io) = inc_temp + (1.0d0-phi_val)*(1.0d0-tau_c)*profit_val 
                    inctax(ia,ie,it,io)     = fun_tax_inc(inc_temp,y_H) &
                                                + tau_d*(1.0d0-phi_val)*(1.0d0-tau_c)*profit_val
                    net_inc(ia,ie,it,io)    = (1.0d0-tau_c)*(1.0d0-tau_d)*(1.0d0-phi_val)*profit_val 
                    payrollsize(ia,ie,it,io)= w*n_val+phi_val*profit_val
	        enddo ! ia
        enddo !ie
    enddo ! it    
	! Retirees
	pregov_incr = pen + r*a_grid
	taxableincr = pregov_incr
	do ia = 1,na
		inctaxr(ia)     = fun_tax_inc(taxableincr(ia),y_H)
        net_incr(ia)    = taxableincr(ia) - inctaxr(ia)
	enddo
	
    !! Some checks ---------------------
    maxinc_work = maxval(pregov_inc(:,:,:,1))
    maxinc_ep = maxval(pregov_inc(:,:,:,2))
    maxinc_es = maxval(pregov_inc(:,:,:,3))
    maxinc_ec = maxval(pregov_inc(:,:,:,4))

    !write(*,*) "targets_compute: OK until here"

    ! Check if each quintile of wealth and income contains 20% of mass
    ! Outputs: mass_q, inc_mass_q
    call checkdist()
    
    !-------------------------------------------------------------
    ! Calculate the average tax rate by income percentiles
    !-------------------------------------------------------------
    if (verbose>=2) write(*,*) "targets_compute: before atr_by_inc"
    call atr_by_inc()
    ! Outputs: atr_vec(5)
    
    !-------------------------------------------------------------
    ! Calculate shares of work/sole-prop/S-corp/C-corp/Retirees
    !-------------------------------------------------------------
    
	share_work = sum(distrib%muy(:,:,:,1))
    share_EP   = sum(distrib%muy(:,:,:,2))
    share_ES   = sum(distrib%muy(:,:,:,3))
    share_EC   = sum(distrib%muy(:,:,:,4))
    share_ret  = sum(distrib%mur)

    share_active = share_work + share_EP + share_ES + share_EC
    share_all    = share_active + share_ret

    if ( abs(share_all-1.0d0)>1d-6 ) call myerror("targets_compute: all does not sum to 1")

    share_work_act  = share_work/share_active
    share_entre     = share_EP + share_ES + share_EC
    share_entre_act = share_entre/share_active
    
    ! Relative share of LFO given entrepreneurship 
    share_EP_entre = share_EP/share_entre
    share_ES_entre = share_ES/share_entre
    share_EC_entre = share_EC/share_entre
! 	write(*,*) "Calculate shares of work/sole-prop/S-corp/C-corp/Retirees...done!"
      
    !-------------------------------------------------------------
    ! Calculate shares of incomes reported as wage by S and C-corps 
    ! 	by profit (bottom 40%, middle 40%, top 20% by LFO)
    ! Note: shares = agg wages / agg incomes
    !-------------------------------------------------------------    
    ! Outputs: share_wage_ES,share_wage_EC
    call income_declare()
     if (verbose>=1) write(*,*) "Calculate income_declare...done!"
	!-------------------------------------------------------------
    ! Average eps, labor supply, and income, average income ratio (entre to worker)
    !-------------------------------------------------------------    
	! Outputs from this part: 
	!	mean_eps,cond_mean_eps_work, aveinc_work, aveinc_entre,ave_l_work, ratio_aveinc_entre_worker

	
    ! Average epsilon in population and among workers
    mean_eps = 0.0d0
    cond_mean_eps_work = 0.0d0
    ! Average income of workers and entrepreneurs
    aveinc_work = 0.0d0
    aveinc_entre = 0.0d0
    ! Average labor supply of workers
    ave_l_work = 0.0d0
    do it=1,ntheta
        do ie=1,neps
            do ia=1,na
                mean_eps     = mean_eps + eps_grid(ie)*sum(mu(ia,ie,it,:))
                ! Worker
                io = 1
                cond_mean_eps_work = cond_mean_eps_work + eps_grid(ie)*mu(ia,ie,it,io)
                aveinc_work= aveinc_work + pregov_inc(ia,ie,it,io)*mu(ia,ie,it,io)
                ave_l_work = ave_l_work + pol%lpoly_w(ia,ie,it)*mu(ia,ie,it,io)
                ! Entrepreneur
                do io = 2,4
                	aveinc_entre = aveinc_entre + pregov_inc(ia,ie,it,io)*mu(ia,ie,it,io)
                enddo ! io
            enddo ! ia
        enddo ! ie
    enddo ! it
    mean_eps            = mean_eps/share_active
    cond_mean_eps_work  = cond_mean_eps_work/share_work
    aveinc_work   = aveinc_work / share_work
	aveinc_entre = aveinc_entre / share_entre
	ave_l_work = ave_l_work/share_work
	! Ave inc of entrepreneurs to ave inc of workers
	ratio_aveinc_entre_worker =  aveinc_entre/aveinc_work
	
	! Median inc of entrepreneurs to median inc of workers
	
	! Median income of entrepreneurs
	
	mu_vec = reshape(mu(:,:,:,2:4),[na*neps*ntheta*3])
	inc_vec = reshape(pregov_inc(:,:,:,2:4),[na*neps*ntheta*3])
    ! All percentiles we want to look at
    qw = [0.5d0] 
	! Income percentiles corresponding to qw
    inc_p = quantili(inc_vec,mu_vec,qw)
    medinc_entre = inc_p(1)
    
    ! Median income of workers
	mu_vec = reshape(mu(:,:,:,1),[na*neps*ntheta])
	inc_vec = reshape(pregov_inc(:,:,:,1),[na*neps*ntheta])
    ! All percentiles we want to look at
    qw = [0.5d0] 
	! Income percentiles corresponding to qw
    inc_p = quantili(inc_vec,mu_vec,qw)
    medinc_work = inc_p(1)
    ! The ratio
    ratio_medinc_entre_worker = medinc_entre/medinc_work
    deallocate(mu_vec,inc_vec)

    ! Median inc, netinc, and size by LFO: median_inc_lfo(no), median_netinc_lfo(no), median_emp_lfo(no)
	
	! Median income by lfo
	do io = 2,4
        mu_vec = reshape(mu(:,:,:,io),[na*neps*ntheta])
        inc_vec = reshape(pregov_inc(:,:,:,io),[na*neps*ntheta])
        ! All percentiles we want to look at
        qw = [0.5d0] 
        ! Income percentiles corresponding to qw
        inc_p = quantili(inc_vec,mu_vec,qw)
        median_inc_lfo(io) = inc_p(1)
        deallocate(mu_vec,inc_vec)
    enddo
	! Median net income by lfo
	do io = 2,4
        mu_vec = reshape(mu(:,:,:,io),[na*neps*ntheta])
        inc_vec = reshape(net_inc(:,:,:,io),[na*neps*ntheta])
        ! All percentiles we want to look at
        qw = [0.5d0] 
        ! Income percentiles corresponding to qw
        inc_p = quantili(inc_vec,mu_vec,qw)
        median_netinc_lfo(io) = inc_p(1)
        deallocate(mu_vec,inc_vec)
    enddo
    ! Median emp by lfo
	do io = 2,4
        mu_vec = reshape(sum(mu(:,:,:,io),dim=2),[na*ntheta])
        inc_vec = reshape(pol%npol(:,:,io),[na*ntheta])
        ! All percentiles we want to look at
        qw = [0.5d0] 
        ! Income percentiles corresponding to qw
        inc_p = quantili(inc_vec,mu_vec,qw)
        median_emp_lfo(io) = inc_p(1)
        deallocate(mu_vec,inc_vec)
    enddo

     if (verbose>=1) write(*,*) "Calculate Average eps, labor supply, and income, average and median income ratios...done!"
    !-------------------------------------------------------------
    ! Fraction of borrowing-constrained entre
    ! Leverage ratio of entre by LFO
    !-------------------------------------------------------------    
    ! Outputs: frac_bor_ep, frac_bor_es, frac_bor_ec
    !			leverage_ep,leverage_es,leverage_ec
    call frac_borcon()
     if (verbose>=1) write(*,*) "Calculating frac_borcon...done!"
    !-------------------------------------------------------------
    ! Share of employers (hiring>0), average firm size by LFO
    !-------------------------------------------------------------
    ! Outputs: share_empl_EP, share_empl_ES, share_empl_EC, share_empl_entre
    !			empsize_EP, empsize_ES,empsize_EC,empsize_rel_ES,empsize_rel_EC,
    !           payrollsize_lfo
    call employment_stats()
	 if (verbose>=1) write(*,*) "Calculating employment_stats...done!"
    !-------------------------------------------------------------
    ! Calculate share of income and wealth owned by entre
    !-------------------------------------------------------------
    
    ! Share of wealth and income held by entrepreneurs in active labor force
    ! Note: income is pre-gov income (use pregov_inc)
    wealth_entre = 0d0
    wealth_all   = 0d0
    share_wealth_entre = 0d0
	
	inc_entre = 0.0d0
	inc_all = 0.0d0
	share_inc_entre = 0.0d0
	
    do io = 1,4 ! active population only
        do it = 1,ntheta
            do ie = 1,neps
            	inc_all = inc_all + sum(pregov_inc(:,ie,it,io)*mu(:,ie,it,io))
                wealth_all = wealth_all + sum(a_grid*mu(:,ie,it,io))
                if (io>=2 .and. io<=4) then !entre
                    wealth_entre = wealth_entre + sum(a_grid*mu(:,ie,it,io))
                    inc_entre = inc_entre + sum(pregov_inc(:,ie,it,io)*mu(:,ie,it,io)) 
                endif
            enddo
        enddo
    enddo
    
    share_wealth_entre = wealth_entre/wealth_all
    share_inc_entre = inc_entre/inc_all
    if (share_wealth_entre<0.0d0 .or. share_wealth_entre>1.0d0) then
        write(*,*) "share_wealth_entre not correct"
    endif
    if (share_inc_entre<0.0d0 .or. share_inc_entre>1.0d0) then
        write(*,*) "share_inc_entre not correct"
    endif
     if (verbose>=1) write(*,*) "Calculating share_wealth_entre and share_inc_entre...done!"  
    
    !Calculate shares of wealth and income held by LFO
    wealth_share_ep = 0.0d0
    wealth_share_es = 0.0d0
    wealth_share_ec = 0.0d0
    inc_share_ep = 0.0d0
    inc_share_es = 0.0d0
    inc_share_ec = 0.0d0
    netinc_share_ep = 0.0d0
    netinc_share_es = 0.0d0 
    netinc_share_ec = 0.0d0

    do io = 2,4 
        do it = 1,ntheta
            do ie = 1,neps
                if (io==2) then !EP
                    wealth_share_ep = wealth_share_ep + sum(a_grid*mu(:,ie,it,io))
                    inc_share_ep    = inc_share_ep + sum(pregov_inc(:,ie,it,io)*mu(:,ie,it,io)) 
                    netinc_share_ep = netinc_share_ep + sum(net_inc(:,ie,it,io)*mu(:,ie,it,io)) 

                endif
                if (io==3) then !ES
                    wealth_share_es = wealth_share_es + sum(a_grid*mu(:,ie,it,io))
                    inc_share_es    = inc_share_es + sum(pregov_inc(:,ie,it,io)*mu(:,ie,it,io)) 
                    netinc_share_es = netinc_share_es + sum(net_inc(:,ie,it,io)*mu(:,ie,it,io))
                endif
                if (io==4) then !EC
                    wealth_share_ec = wealth_share_ec + sum(a_grid*mu(:,ie,it,io))
                    inc_share_ec    = inc_share_ec + sum(pregov_inc(:,ie,it,io)*mu(:,ie,it,io)) 
                    netinc_share_ec = netinc_share_ec + sum(net_inc(:,ie,it,io)*mu(:,ie,it,io))
                endif
            enddo
        enddo
    enddo
    
    wealth_share_ep = wealth_share_ep/wealth_entre
    wealth_share_es = wealth_share_es/wealth_entre
    wealth_share_ec = wealth_share_ec/wealth_entre
    inc_share_ep = inc_share_ep/inc_entre
    inc_share_es = inc_share_es/inc_entre
    inc_share_ec = inc_share_ec/inc_entre
    netinc_entre = netinc_share_ep+netinc_share_es+netinc_share_ec
    netinc_share_ep = netinc_share_ep/netinc_entre
    netinc_share_es = netinc_share_es/netinc_entre
    netinc_share_ec = netinc_share_ec/netinc_entre

    !-------------------------------------------------------------
    ! Calculate taxes paid by top earners
    !-------------------------------------------------------------
    ! Output: taxes_top(1:4)
    call taxes_by_inc()
     if (verbose>=1) write(*,*) "Calculating taxes paid by top earners...done!"    

    !-------------------------------------------------------------
    ! Calculate occupational shares by wealth and income
    !-------------------------------------------------------------
    ! Outputs:   share_entre_top, share_lfo_top,share_entre_quint,share_lfo_quint 
    call share_entre_by_wealth()
    
    ! Outputs: share_entre_top_inc, share_lfo_top_inc,share_entre_quint_inc,share_lfo_quint_inc
	call share_entre_by_inc()
	
 	 if (verbose>=1) write(*,*) "Calculating occupational shares by wealth and income...done!"    
    !-------------------------------------------------------------
    ! Calculate share of young population paying top marginal tax
    !-------------------------------------------------------------
     ! Outputs: share_toptax
    call sub_share_toptax()
     if (verbose>=1) write(*,*) "Calculating share of young population paying top marginal tax...done!"    

    
    !-------------------------------------------------------------
    ! Calculate inequality statistics
    !-------------------------------------------------------------

    
    ! Wealth inequality 
    ! Outputs: gini_wealth_all,gini_wealth_entre,wealth_share
    call wealth_ineq()
    
    ! Income inequality (gini, top1, etc.)
    ! Outputs: gini_inc_all, gini_inc_entre,inc_share
    call income_ineq()
     if (verbose>=1) write(*,*) "Calculating inequality statistics...done!"    
    !-------------------------------------------------------------
    ! Employment share by firm size bins 
    !-------------------------------------------------------------
    ! Outputs: share_emp
    call share_emp_firmsize()
  	 if (verbose>=1) write(*,*) "Calculating share_emp_firmsize...done!"    
    !-------------------------------------------------------------
    ! Entrepreneurship exit rate
    !-------------------------------------------------------------
    ! Outputs:
    ! trans_entre_work,trans_work_entre,trans_CP,trans_PC,trans_SC,trans_mat
    call occ_trans()
     if (verbose>=1) write(*,*) "Calculating occ_trans...done!"    

    !----------------------------------------------------------------
    ! Pack results into a structure (called derived data type in F)
    ! The derived type model_targets is defined in mod_gloabls
    !----------------------------------------------------------------

    ! - Occ/LFO distribution
    model_targets%share_entre_act  = share_entre_act
    model_targets%share_EP_entre   = share_EP_entre
    model_targets%share_ES_entre   = share_ES_entre
    model_targets%share_EC_entre   = share_EC_entre
    model_targets%share_active     = share_active
    ! - Share of income declared as wage (ES and EC)
    model_targets%share_wage_ES = share_wage_ES
    model_targets%share_wage_EC = share_wage_EC
    model_targets%share_wage_EC_cf = share_wage_EC_cf
    model_targets%share_wage_qn_ES = share_wage_qn_ES
    model_targets%share_wage_qn_EC = share_wage_qn_EC
    model_targets%share_wage_qk_ES = share_wage_qk_ES
    model_targets%share_wage_qk_EC = share_wage_qk_EC
    ! - Share facing financial constraints and leverages (EP,ES,EC)
    model_targets%frac_bor_ep = frac_bor_ep
    model_targets%frac_bor_es = frac_bor_es
    model_targets%frac_bor_ec = frac_bor_ec
    model_targets%leverage_ep = leverage_ep
    model_targets%leverage_es = leverage_es
    model_targets%leverage_ec = leverage_ec
	! - Average eps,income,labor supply,income ratio
	model_targets%ratio_aveinc_entre_worker = ratio_aveinc_entre_worker
	model_targets%ratio_medinc_entre_worker = ratio_medinc_entre_worker
	model_targets%ave_l_work = ave_l_work
	! - Share of employer
	model_targets%share_empl_EP = share_empl_EP
	model_targets%share_empl_ES = share_empl_ES
	model_targets%share_empl_EC = share_empl_EC
	model_targets%share_empl_entre = share_empl_entre
	! - Average firm size relative to EP
	model_targets%empsize_rel_ES = empsize_rel_ES
	model_targets%empsize_rel_EC = empsize_rel_EC
    ! - payroll size by lfo 
    model_targets%payrollsize_lfo = payrollsize_lfo
	! - Share of income and wealth owned by entre
	model_targets%share_inc_entre = share_inc_entre
	model_targets%share_wealth_entre = share_wealth_entre
	! - Inequality measures
	model_targets%gini_wealth_all = gini_wealth_all
	model_targets%gini_wealth_entre = gini_wealth_entre
	model_targets%wealth_share = wealth_share
	model_targets%gini_inc_all = gini_inc_all
	model_targets%gini_inc_entre = gini_inc_entre
	model_targets%inc_share = inc_share
	! - Occ/LFO choice by income and wealth
	model_targets%share_entre_top = share_entre_top
	model_targets%share_lfo_top = share_lfo_top
	model_targets%share_entre_quint = share_entre_quint
	model_targets%share_lfo_quint = share_lfo_quint
	model_targets%share_entre_top_inc=share_entre_top_inc
	model_targets%share_lfo_top_inc=share_lfo_top_inc
	model_targets%share_entre_quint_inc=share_entre_quint_inc
	model_targets%share_lfo_quint_inc=share_lfo_quint_inc
	model_targets%share_toptax=share_toptax
	model_targets%taxes_top = taxes_top
    ! - median size by lfo
    model_targets%median_emp_lfo = median_emp_lfo
    model_targets%median_inc_lfo = median_inc_lfo
    model_targets%median_netinc_lfo = median_netinc_lfo
	! - Employment share by firm size
	model_targets%share_emp = share_emp
	model_targets%trans_entre_work = trans_entre_work
	model_targets%trans_work_entre = trans_work_entre
	model_targets%trans_CP = trans_CP
	model_targets%trans_PC = trans_PC
	model_targets%trans_SC = trans_SC

	model_targets%trans_mat = trans_mat
	! - income and wealth share by LFO
	model_targets%wealth_share_EP = wealth_share_EP
	model_targets%wealth_share_ES = wealth_share_ES
	model_targets%wealth_share_EC = wealth_share_EC
	model_targets%inc_share_EP = inc_share_EP
	model_targets%inc_share_ES = inc_share_ES
	model_targets%inc_share_EC = inc_share_EC
    model_targets%netinc_share_EP = netinc_share_EP
    model_targets%netinc_share_ES = netinc_share_ES
    model_targets%netinc_share_EC = netinc_share_EC
	! - average tax rate by income percentiles
	model_targets%atr_vec = atr_vec
	
!     call CPU_TIME(t2)
     if (verbose>=1) write(*,*) 'Time to compute targets: ',real(t2-t1),'seconds.'
    
    contains    
    !===============================================================================!
    subroutine checkdist()
    implicit none
    ! Purpose: check if each wealth or income quintile contains 20% of mass.
    !	This may be a problem if there is some grid points that contains more than 20% of mass.
   
    ! Declare locals
    real(8) :: q_vec(4),w_q(4),inc_q(4)
    !real(8) :: inc_vec(na*ntheta*neps*4),mu_vec(na*ntheta*neps*4)
    real(8) :: muy_a(na)
    real(8), allocatable :: inc_vec(:), mu_vec(:)
    
    ! Check wealth distribution of the young
	do ia = 1,na
		muy_a(ia) = sum(distrib%muy(ia,:,:,:))
	enddo
	
	! vectorize income
	inc_vec = reshape(pregov_inc,[na*ntheta*neps*no])
	! vectorize income distribution
	mu_vec = reshape(mu,[na*ntheta*neps*no])
	
	q_vec = [0.2d0,0.4d0,0.6d0,0.8d0]
	w_q = quantili(a_grid,muy_a,q_vec)
	inc_q = quantili(inc_vec,mu_vec,q_vec)
	deallocate(mu_vec,inc_vec)
	mass_q = 0.0d0
	inc_mass_q = 0.0d0
	do ia = 1,na
		a_val = a_grid(ia)
		if (a_val<=w_q(1)) then ! Q1
			mass_q(1) = mass_q(1) + muy_a(ia)
		elseif (a_val>w_q(1) .and. a_val<=w_q(2)) then ! Q2
			mass_q(2) = mass_q(2) + muy_a(ia)
		elseif (a_val>w_q(2) .and. a_val<=w_q(3)) then ! Q3
			mass_q(3) = mass_q(3) + muy_a(ia)
		elseif (a_val>w_q(3) .and. a_val<=w_q(4)) then ! Q4
			mass_q(4) = mass_q(4) + muy_a(ia)
		else! Q5
			mass_q(5) = mass_q(5) + muy_a(ia)
		endif
		do it = 1,ntheta
		do ie = 1,neps
			do io = 1,no
				if (pregov_inc(ia,ie,it,io)<=inc_q(1)) then !Q1
					inc_mass_q(1) = inc_mass_q(1)  + mu(ia,ie,it,io)
				elseif (pregov_inc(ia,ie,it,io)>inc_q(1) .and. pregov_inc(ia,ie,it,io)<=inc_q(2) ) then !Q2
					inc_mass_q(2) = inc_mass_q(2)  + mu(ia,ie,it,io)
				elseif (pregov_inc(ia,ie,it,io)>inc_q(2) .and. pregov_inc(ia,ie,it,io)<=inc_q(3) ) then !Q3
					inc_mass_q(3) = inc_mass_q(3)  + mu(ia,ie,it,io)
				elseif (pregov_inc(ia,ie,it,io)>inc_q(3) .and. pregov_inc(ia,ie,it,io)<=inc_q(4) ) then !Q4
					inc_mass_q(4) = inc_mass_q(4)  + mu(ia,ie,it,io)
				else !Q5
					inc_mass_q(5) = inc_mass_q(5)  + mu(ia,ie,it,io)
				endif
			enddo ! io
		enddo !ie
		enddo !it
	enddo ! ia
	
	! Normalize
	mass_q = mass_q/sum(mass_q)
	inc_mass_q = inc_mass_q/sum(inc_mass_q)
	!call disp('mass_q',mass_q)
	!call disp('inc_mass_q',inc_mass_q)
	! Debug
	!call disp('mu_a(1:10)',mu_a(1:10))
	!pause
	! Calculate mass with ia = 1 and ia = na
	mass_a1    = sum(mu(1,:,:,:))
	mass_a1_a3 = sum(mu(1:3,:,:,:))
	mass_aN = sum(mu(na,:,:,:))
	
    end subroutine checkdist

    !===============================================================================!
    subroutine income_declare()
    implicit none
    ! Purpose: Calculate shares of incomes reported as wage by S and C-corps 
    ! Note: shares = agg wages / agg incomes
    ! All firms and by firm size bins
    ! Note: This is different from average phi (see subroutine commented out above)
    ! Declare local variables
    real(8) :: positive_tiny
    real(8) :: qn_e(2),wage_qn_ES(3),wage_qn_EC(3),profit_qn_ES(3),profit_qn_EC(3)
    real(8) :: qk_e(2),wage_qk_ES(3),wage_qk_EC(3),profit_qk_ES(3),profit_qk_EC(3)
    real(8) :: wage_ES,wage_EC,wage_EC_cf,profit_ES,profit_EC,profit_EC_cf
    
    positive_tiny = 1.0d-10
    ! Find q40 and q80 of firm size (n) among all firms using quantili 
    inc_vec = reshape(pol%npol(:,:,2:4), [na*ntheta*3])
    mu_vec = reshape(sum(mu(:,:,:,2:4),dim=2), [na*ntheta*3])
    qn_e = quantili(inc_vec,mu_vec,[0.4d0,0.8d0])
    ! Find q40 and q80 of firm size (k) among all firms using quantili 
    inc_vec = reshape(pol%kpol(:,:,2:4), [na*ntheta*3])
    mu_vec = reshape(sum(mu(:,:,:,2:4),dim=2), [na*ntheta*3])
    qk_e = quantili(inc_vec,mu_vec,[0.4d0,0.8d0]) 

    share_wage_ES = 0d0
    share_wage_EC = 0d0
    share_wage_EC_cf = 0.0d0
    wage_ES = 0d0
    wage_EC = 0d0
    wage_EC_cf = 0.0d0
    profit_ES = 0d0
    profit_EC = 0d0
    profit_EC_cf = 0.0d0
    share_wage_qn_ES = 0.0d0
    share_wage_qn_EC = 0.0d0
    wage_qn_ES = 0.0d0
    wage_qn_EC = 0.0d0
    profit_qn_ES = 0.0d0
    profit_qn_EC = 0.0d0
    share_wage_qk_ES = 0.0d0
    share_wage_qk_EC = 0.0d0
    wage_qk_ES = 0.0d0
    wage_qk_EC = 0.0d0
    profit_qk_ES = 0.0d0
    profit_qk_EC = 0.0d0
    do it = 1,ntheta
        theta_val  = theta_grid(it)           
        do ie = 1,neps
            eps_val    = eps_grid(ie)
	        do ia = 1,na
                a_val      = a_grid(ia)
                ! ES
                	io = 3
                    k_val    = pol%kpol(ia,it,io)
                    n_val    = pol%npol(ia,it,io)
                    phi_val  = pol%phipol(ia,it,io)
                    profit_val = fun_profit(theta_val,k_val,le,n_val,r,w) 
                    profit_ES = profit_ES + profit_val*mu(ia,ie,it,io)
                    wage_ES  = wage_ES + phi_val*profit_val*mu(ia,ie,it,io)
                    ! By emp size bins (n)
                    if (n_val<=qn_e(1)) then ! bottom 40% 
                    	profit_qn_ES(1) = profit_qn_ES(1)  + profit_val*mu(ia,ie,it,io)
                    	wage_qn_ES(1) = wage_qn_ES(1) + phi_val*profit_val*mu(ia,ie,it,io)
                    elseif (n_val>qn_e(1) .and. n_val<=qn_e(2)) then ! middle 40%
                    	profit_qn_ES(2) = profit_qn_ES(2)  + profit_val*mu(ia,ie,it,io)
                    	wage_qn_ES(2) = wage_qn_ES(2) + phi_val*profit_val*mu(ia,ie,it,io)
                    else ! top 20%
                    	profit_qn_ES(3) = profit_qn_ES(3)  + profit_val*mu(ia,ie,it,io)
                    	wage_qn_ES(3) = wage_qn_ES(3) + phi_val*profit_val*mu(ia,ie,it,io)
                    endif 
                    ! By capital size bins (k)
                    if (k_val<=qk_e(1)) then ! bottom 40% 
                    	profit_qk_ES(1) = profit_qk_ES(1)  + profit_val*mu(ia,ie,it,io)
                    	wage_qk_ES(1) = wage_qk_ES(1) + phi_val*profit_val*mu(ia,ie,it,io)
                    elseif (k_val>qk_e(1) .and. k_val<=qk_e(2)) then ! middle 40%
                    	profit_qk_ES(2) = profit_qk_ES(2)  + profit_val*mu(ia,ie,it,io)
                    	wage_qk_ES(2) = wage_qk_ES(2) + phi_val*profit_val*mu(ia,ie,it,io)
                    else ! top 20%
                    	profit_qk_ES(3) = profit_qk_ES(3)  + profit_val*mu(ia,ie,it,io)
                    	wage_qk_ES(3) = wage_qk_ES(3) + phi_val*profit_val*mu(ia,ie,it,io)
                    endif
                ! EC
                	io = 4
                    k_val    = pol%kpol(ia,it,io)
                    n_val    = pol%npol(ia,it,io)
                    phi_val  = pol%phipol(ia,it,io)
                    profit_val = fun_profit(theta_val,k_val,le,n_val,r,w) 
                    profit_EC = profit_EC + profit_val*mu(ia,ie,it,io)
                    wage_EC  = wage_EC + phi_val*profit_val*mu(ia,ie,it,io)
                    ! Counterfactual: if EC are actually running S-corps
                    k_val    = pol%kpol(ia,it,3)
                    n_val    = pol%npol(ia,it,3)
                    phi_val  = pol%phipol(ia,it,3)
                    profit_val = fun_profit(theta_val,k_val,le,n_val,r,w) 
                    profit_EC_cf = profit_EC_cf + profit_val*mu(ia,ie,it,io)
                    wage_EC_cf = wage_EC_cf + phi_val*profit_val*mu(ia,ie,it,io)
                    ! By emp size bins (n)
                    if (n_val<=qn_e(1)) then ! bottom 40% 
                    	profit_qn_EC(1) = profit_qn_EC(1)  + profit_val*mu(ia,ie,it,io)
                    	wage_qn_EC(1) = wage_qn_EC(1) + phi_val*profit_val*mu(ia,ie,it,io)
                    elseif (n_val>qn_e(1) .and. n_val<=qn_e(2)) then ! middle 40%
                    	profit_qn_EC(2) = profit_qn_EC(2)  + profit_val*mu(ia,ie,it,io)
                    	wage_qn_EC(2) = wage_qn_EC(2) + phi_val*profit_val*mu(ia,ie,it,io)
                    else ! top 20%
                    	profit_qn_EC(3) = profit_qn_EC(3)  + profit_val*mu(ia,ie,it,io)
                    	wage_qn_EC(3) = wage_qn_EC(3) + phi_val*profit_val*mu(ia,ie,it,io)
                    endif 
                    ! By capital size bins (k)
                    if (k_val<=qk_e(1)) then ! bottom 40% 
                    	profit_qk_EC(1) = profit_qk_EC(1)  + profit_val*mu(ia,ie,it,io)
                    	wage_qk_EC(1) = wage_qk_EC(1) + phi_val*profit_val*mu(ia,ie,it,io)
                    elseif (k_val>qk_e(1) .and. k_val<=qk_e(2)) then ! middle 40%
                    	profit_qk_EC(2) = profit_qk_EC(2)  + profit_val*mu(ia,ie,it,io)
                    	wage_qk_EC(2) = wage_qk_EC(2) + phi_val*profit_val*mu(ia,ie,it,io)
                    else ! top 20%
                    	profit_qk_EC(3) = profit_qk_EC(3)  + profit_val*mu(ia,ie,it,io)
                    	wage_qk_EC(3) = wage_qk_EC(3) + phi_val*profit_val*mu(ia,ie,it,io)
                    endif
	        enddo ! ia
        enddo !ie
    enddo ! it

    share_wage_ES = wage_ES/(positive_tiny+profit_ES)
    share_wage_EC = wage_EC/(positive_tiny+profit_EC)
    share_wage_EC_cf = wage_EC_cf/(positive_tiny+profit_EC_cf)
    do iq = 1,3
	    share_wage_qn_ES(iq) = wage_qn_ES(iq)/(positive_tiny + profit_qn_ES(iq))
	    share_wage_qn_EC(iq) = wage_qn_EC(iq)/(positive_tiny + profit_qn_EC(iq))
        share_wage_qk_ES(iq) = wage_qk_ES(iq)/(positive_tiny + profit_qk_ES(iq))
	    share_wage_qk_EC(iq) = wage_qk_EC(iq)/(positive_tiny + profit_qk_EC(iq))
	enddo
    end subroutine income_declare
!===============================================================================!
  
    subroutine frac_borcon()
    implicit none
    ! Purpose: Compute moments related to financial constraints for entre
    !	- fraction of borrowing constraint entre by lfo
    !	- leverage ratio by lfo (k-a)/k
    
    frac_bor_ep = 0.0d0
    frac_bor_es = 0.0d0
    frac_bor_ec = 0.0d0
    leverage_ep = 0.0d0
    leverage_es = 0.0d0
    leverage_ec = 0.0d0
     do it = 1,ntheta
        theta_val  = theta_grid(it)           
        do ie = 1,neps
            eps_val    = eps_grid(ie)
	        do ia = 1,na
                a_val      = a_grid(ia)
                ! EP
                	io = 2
                    k_val    = pol%kpol(ia,it,io)
                    if (abs(k_val-lambda*a_val)<=1d-12) then
                        frac_bor_ep = frac_bor_ep + mu(ia,ie,it,io)
                    endif
                    leverage_ep = leverage_ep + max(k_val-a_val,0.0d0)/k_val*mu(ia,ie,it,io)
                ! ES
                	io = 3
                    k_val    = pol%kpol(ia,it,io)
                    if (abs(k_val-lambda_es*a_val)<=1d-12) then
                        frac_bor_es = frac_bor_es + mu(ia,ie,it,io)
                    endif
                    leverage_es = leverage_es + max(k_val-a_val,0.0d0)/k_val*mu(ia,ie,it,io)
                ! EC
                	io = 4
                    k_val    = pol%kpol(ia,it,io)
                    if (abs(k_val-lambda_ec*a_val)<=1d-12) then
                        frac_bor_ec = frac_bor_ec + mu(ia,ie,it,io)
                    endif
                    leverage_ec = leverage_ec + max(k_val-a_val,0.0d0)/k_val*mu(ia,ie,it,io)
	        enddo ! ia
        enddo !ie
    enddo ! it
    frac_bor_ep = frac_bor_ep/share_EP
    frac_bor_es = frac_bor_es/share_ES
    frac_bor_ec = frac_bor_ec/share_EC
   	leverage_ep = leverage_ep/share_EP
    leverage_es = leverage_es/share_ES
   	leverage_ec = leverage_ec/share_EC

    end subroutine frac_borcon
    !===============================================================================!
    
    subroutine employment_stats()
    implicit none
    ! Purpose: Compute moments related to firm size, share of employers etc.
    
    
    ! Share of employer by LFO
    share_empl_EP = 0d0 
    share_empl_ES = 0d0 
    share_empl_EC = 0d0 
    share_empl_entre = 0d0
    
    ! Firm size by LFO
    empsize_EP = 0.0d0
    empsize_ES = 0.0d0
    empsize_EC = 0.0d0

    ! Payroll by lfo
    payrollsize_lfo = 0.0d0

    do it = 1,ntheta
        theta_val  = theta_grid(it)           
        do ie = 1,neps
            eps_val    = eps_grid(ie)
	        do ia = 1,na
                a_val      = a_grid(ia)
                ! EP
                	io = 2
                    n_val    = pol%npol(ia,it,io)
                    if (n_val>1.0d-6) then
                        share_empl_EP = share_empl_EP + mu(ia,ie,it,io)
                        share_empl_entre = share_empl_entre + mu(ia,ie,it,io)
                    endif 
                    empsize_EP = empsize_EP + (le+n_val)*mu(ia,ie,it,io)
                    payrollsize_lfo(io) = payrollsize_lfo(io) + payrollsize(ia,ie,it,io)*mu(ia,ie,it,io)
                ! ES
                	io = 3
                    n_val    = pol%npol(ia,it,io)
                    if (n_val>1.0d-6) then
                        share_empl_ES = share_empl_ES + mu(ia,ie,it,io)
                        share_empl_entre = share_empl_entre + mu(ia,ie,it,io)
                    endif
                    empsize_ES = empsize_ES + (le+n_val)*mu(ia,ie,it,io)
                    payrollsize_lfo(io) = payrollsize_lfo(io) + payrollsize(ia,ie,it,io)*mu(ia,ie,it,io)
                ! EC
                	io = 4
                    n_val    = pol%npol(ia,it,io)
                    if (n_val>1.0d-6) then
                        share_empl_EC = share_empl_EC + mu(ia,ie,it,io)
                        share_empl_entre = share_empl_entre + mu(ia,ie,it,io)
                    endif
                    empsize_EC = empsize_EC + (le+n_val)*mu(ia,ie,it,io)
                    payrollsize_lfo(io) = payrollsize_lfo(io) + payrollsize(ia,ie,it,io)*mu(ia,ie,it,io)
                
	        enddo ! ia
        enddo !ie
    enddo ! it
    share_empl_entre = share_empl_entre/share_entre
    share_empl_EP = share_empl_EP/share_EP
    share_empl_ES = share_empl_ES/share_ES
    share_empl_EC = share_empl_EC/share_EC
    
    empsize_EP = empsize_EP/share_EP
    empsize_ES = empsize_ES/share_ES
    empsize_EC = empsize_EC/share_EC
    
    empsize_rel_ES = empsize_ES/empsize_EP
    empsize_rel_EC = empsize_EC/empsize_EP 
    
    payrollsize_lfo(2) = payrollsize_lfo(2)/share_EP
    payrollsize_lfo(3) = payrollsize_lfo(3)/share_ES
    payrollsize_lfo(4) = payrollsize_lfo(4)/share_EC


!     !Calculate relative size based on assets/capital
!     assetsize_EP = 0.0d0
!     assetsize_ES = 0.0d0
!     assetsize_EC = 0.0d0
!     
!     do it = 1,ntheta
!         theta_val  = theta_grid(it)           
!         do ie = 1,neps
!             eps_val    = eps_grid(ie)
! 	        do ia = 1,na
!                 a_val      = a_grid(ia)
!                 EP
!                 	iz = 2
!                     assetsize_EP = assetsize_EP + a_val*mu(ia,ie,it,iz)
!                 ES
!                 	iz = 3
!                     assetsize_ES = assetsize_ES + a_val*mu(ia,ie,it,iz)
!                 EC
!                 	iz = 4
!                     assetsize_EC = assetsize_EC + a_val*mu(ia,ie,it,iz)
!             
! 	        enddo ! ia
!         enddo !ie
!     enddo ! it
    
    end subroutine employment_stats
    !===============================================================================!
    
    subroutine wealth_ineq()
    implicit none
    ! Purpose: Compute wealth inequality statistics
    !	- Gini: unconditional and among entrepreneurs
    !	- share of wealth in the top 1, 10, 20% and bottom 40%
    
    ! Local variables
    real(8), allocatable :: mu_a_pmf(:),mu_a_pmf_entre(:)
    real(8), allocatable :: mu_cum(:),wealth_cum(:)

    allocate(mu_a_pmf(na),mu_a_pmf_entre(na))

    ! Marginal dist of mu wrt wealth
	do ia = 1,na
		mu_a_pmf(ia) = sum(mu(ia,:,:,1:4))
		mu_a_pmf_entre(ia) = sum(mu(ia,:,:,2:4))
	enddo
    ! Normalize 
    mu_a_pmf = mu_a_pmf/sum(mu_a_pmf)
    mu_a_pmf_entre = mu_a_pmf_entre/sum(mu_a_pmf_entre)
    
	! Compute Gini of wealth 
	call lrzcurve(mu_a_pmf,a_grid,gini_wealth_all,mu_cum,wealth_cum)
	call lrzcurve(mu_a_pmf_entre,a_grid,gini_wealth_entre)

    !Shares of wealth owned by a certain pct of the population:
	!	[top 1, top 10, top 20, bottom 40%]
	
	wealth_share = 0.0d0
    wealth_share(1) = 1.0d0-linint(mu_cum,wealth_cum,0.99d0)
    wealth_share(2) = 1.0d0-linint(mu_cum,wealth_cum,0.9d0)
    wealth_share(3) = 1.0d0-linint(mu_cum,wealth_cum,0.8d0)
    wealth_share(4) = linint(mu_cum,wealth_cum,0.40d0)
  


    end subroutine wealth_ineq
    !===============================================================================!
    
    subroutine income_ineq()
    implicit none

    ! Purpose: Compute income inequality statistics
    ! Note: for income, use pregov_inc, and mu, both have dim:(na,neps,ntheta,nz).
    !	- Gini: unconditional and among entrepreneurs
    !	- share of income in the top 1, 10, 20% and bottom 40%
    
    ! Local variables
    real(8),allocatable :: mu_cum(:),inc_cum(:)
    real(8),allocatable :: mu_vec(:),inc_vec(:),mu_vec_entre(:),inc_vec_entre(:)
    !real(8) :: mu_vec(na*neps*ntheta*4),inc_vec(na*neps*ntheta*4), &
    !	mu_vec_entre(na*neps*ntheta*3),inc_vec_entre(na*neps*ntheta*3)
    
    
    ! Vectorize pregov_inc and mu
    mu_vec = reshape(mu(:,:,:,1:4),[na*neps*ntheta*4]) 
    inc_vec = reshape(pregov_inc(:,:,:,1:4),[na*neps*ntheta*4]) 
    mu_vec_entre = reshape(mu(:,:,:,2:4),[na*neps*ntheta*3]) 
    inc_vec_entre = reshape(pregov_inc(:,:,:,2:4),[na*neps*ntheta*3]) 
    
    ! Normalize 
    mu_vec = mu_vec/sum(mu_vec)
    mu_vec_entre = mu_vec_entre/sum(mu_vec_entre)

	! Compute Gini of income 
	call lrzcurve(mu_vec,inc_vec,gini_inc_all,mu_cum,inc_cum)
	call lrzcurve(mu_vec_entre,inc_vec_entre,gini_inc_entre)

    !Shares of income owned by a certain pct of the population:
	!	[top 1, top 10, top 20, bottom 40%]
	
	inc_share = 0.0d0
    inc_share(1) = 1.0d0-linint(mu_cum,inc_cum,0.99d0)
    inc_share(2) = 1.0d0-linint(mu_cum,inc_cum,0.9d0)
    inc_share(3) = 1.0d0-linint(mu_cum,inc_cum,0.8d0)
    inc_share(4) = linint(mu_cum,inc_cum,0.40d0)
	
	
    end subroutine income_ineq    
    

!===============================================================================!

    subroutine share_entre_by_inc()
    implicit none
    ! Purpose: Compute the share of entre and LFO distribution in the top percentiles 
    !	of the income distrib (top 1%, 5%, and 10%), and by income quintiles
    ! Note: consider active population only!
    
    ! Local variables:
    real(8) :: qw(7), inc_p(7),iq
    
    !real(8) :: mu_vec(na*neps*ntheta*4),inc_vec(na*neps*ntheta*4)
    real(8), allocatable :: mu_vec(:), inc_vec(:)
	real(8) :: inc_quintiles(0:5),weight0,weight1,mass_young_q,mass_young_toppct(3), &
		mass_between_q,mass_at_q0, inc_toppctile(3), &
    	mass_at_q1,mass_entre_between_q,mass_entre_at_q0,mass_entre_at_q1,mass_entre_q, &
    	mass_lfo_between_q,mass_lfo_at_q0,mass_lfo_at_q1,mass_lfo_q
    	
	
	! Vectorize mu and pregov_inc
	mu_vec = reshape(mu(:,:,:,1:4),[na*neps*ntheta*4])
	inc_vec = reshape(pregov_inc(:,:,:,1:4),[na*neps*ntheta*4])
	
    ! All percentiles we want to look at
    qw = [0.2d0,0.4d0,0.6d0,0.8d0,0.9d0,0.95d0,0.99d0] 
	! Income percentiles corresponding to qw
    inc_p = quantili(inc_vec,mu_vec,qw)
    if (size(qw)/=size(inc_p)) then
        call myerror("share_entre_by_inc: perc. size not correct")
    endif
    
    ! Find share of entre and LFO dist for
    !	- By quintiles: 0-20%, 20-40%,...,80-100%
    !	- Top percentiles: 1%, top 5%, top 10%

    ! By quintiles

    inc_quintiles(0) = minval(inc_vec) - 1.0d0
    inc_quintiles(1:4) = inc_p(1:4)
    inc_quintiles(5) = maxval(inc_vec) + 1.0d0
   
    ! mass of young agent per quintile
    mass_young_q = sum(mu(:,:,:,1:4))/5.0d0
    ! mass of young agent at the top 1, 5, and 10%
    mass_young_toppct(1) = sum(mu(:,:,:,1:4)) * 0.01d0
    mass_young_toppct(2) = sum(mu(:,:,:,1:4)) * 0.05d0
    mass_young_toppct(3) = sum(mu(:,:,:,1:4)) * 0.10d0

    weight0 = 0.0d0
    weight1 = 0.0d0
    do iq = 1,5
		! mass of young agent with income strictly between two quintile values
		mass_between_q = sum(mu(:,:,:,1:4),mask=(pregov_inc(:,:,:,1:4)>inc_quintiles(iq-1) .and. &
			pregov_inc(:,:,:,1:4)<inc_quintiles(iq)))    
		! mass of young agent with income equal to the lower quintile value
		mass_at_q0 = sum(mu(:,:,:,1:4),mask=(pregov_inc(:,:,:,1:4)==inc_quintiles(iq-1)))    
		mass_at_q1 = sum(mu(:,:,:,1:4),mask=(pregov_inc(:,:,:,1:4)==inc_quintiles(iq)))    

		! what fraction of mass_at_q0 should be assigned to quintile iq?
		weight0 = 1.0d0 - weight1
		! what fraction of mass_at_q1 should be assigned to quintile iq?
		if (mass_at_q1 >0.0d0) then
			weight1 = (mass_young_q-mass_between_q - weight0*mass_at_q0)/mass_at_q1
		else
			weight1 = 1.0d0
		endif
		! weight1 must be between 0 and 1
		if (weight1<0.0d0 .or. weight1 > 1.0d0) then
			write(*,*) "weight1 = ",weight1
			write(*,*) "mass_young_q = ",mass_young_q
			write(*,*) "mass_between_q = ",mass_between_q
			write(*,*) "weight0 = ",weight0
			write(*,*) "mass_at_q0 = ",mass_at_q0
			write(*,*) "mass_at_q1 = ",mass_at_q1
			!pause
		endif
		
		! - Share of entrepreneurs by quintile
		! mass of entrepreneurs with income strictly between two quintile values
		mass_entre_between_q = sum(mu(:,:,:,2:4),mask=(pregov_inc(:,:,:,2:4)>inc_quintiles(iq-1) .and. &
			pregov_inc(:,:,:,2:4)<inc_quintiles(iq)))    
		! mass of entrepreneurs with income equal to the lower quintile value
		mass_entre_at_q0 = weight0*sum(mu(:,:,:,2:4),mask=(pregov_inc(:,:,:,2:4)==inc_quintiles(iq-1)))    
		mass_entre_at_q1 = weight1*sum(mu(:,:,:,2:4),mask=(pregov_inc(:,:,:,2:4)==inc_quintiles(iq)))    
		! total mass of entrepreneurs in quintile iq
		mass_entre_q = mass_entre_between_q + mass_entre_at_q0 + mass_entre_at_q1
		! share of entrepreneurs in quintile iq
		share_entre_quint_inc(iq) = mass_entre_q/mass_young_q
		
		! - Share of lfo by quintile
		do io = 2,4
			! mass of lfo iz with income strictly between two quintile values
			mass_lfo_between_q = sum(mu(:,:,:,io),mask=(pregov_inc(:,:,:,io)>inc_quintiles(iq-1) .and. &
				pregov_inc(:,:,:,io)<inc_quintiles(iq)))    
			! mass of lfo iz with income equal to the lower quintile value
			mass_lfo_at_q0 = weight0*sum(mu(:,:,:,io),mask=(pregov_inc(:,:,:,io)==inc_quintiles(iq-1)))    
			mass_lfo_at_q1 = weight1*sum(mu(:,:,:,io),mask=(pregov_inc(:,:,:,io)==inc_quintiles(iq)))    
			! total mass of lfo iz in quintile iq
			mass_lfo_q = mass_lfo_between_q + mass_lfo_at_q0 + mass_lfo_at_q1
			! share of lfo iz in quintile iq
			share_lfo_quint_inc(io,iq) = mass_lfo_q/mass_entre_q
		enddo ! iz
    enddo ! iq

    ! Top percentiles
    !	- [top 1%, top 5%, top 10%]
    
    inc_toppctile(1) = inc_p(7) ! top 1%
    inc_toppctile(2) = inc_p(6) ! top 5%
    inc_toppctile(3) = inc_p(5) ! top 10%
    
    weight0 = 0.0d0
    do iq = 1,3
		! mass of young agent with income strictly above the percentile value
		mass_between_q = sum(mu(:,:,:,1:4),mask=(pregov_inc(:,:,:,1:4)>inc_toppctile(iq)))    
		! mass of young agent with income equal to the percentile value
		mass_at_q0 = sum(mu(:,:,:,1:4),mask=(pregov_inc(:,:,:,1:4)==inc_toppctile(iq)))    

		! what fraction of mass_at_q0 should be assigned to the top percentile iq?
		if (mass_at_q0 >0.0d0) then
			weight0 = (mass_young_toppct(iq)-mass_between_q)/mass_at_q0
		else
			weight0 = 1.0d0
		endif
		! weight0 must be between 0 and 1
		if (weight0<0.0d0 .or. weight0 > 1.0d0) then
			write(*,*) "weight0 = ",weight0
			write(*,*) "mass_young_toppct(iq) = ",mass_young_toppct(iq)
			write(*,*) "mass_between_q = ",mass_between_q
			write(*,*) "mass_at_q0 = ",mass_at_q0
			!pause
		endif
		
		! - Share of entrepreneurs at top percentile
		! mass of entrepreneurs with income strictly above the percentile value
		mass_entre_between_q = sum(mu(:,:,:,2:4),mask=(pregov_inc(:,:,:,2:4)>inc_toppctile(iq)))    
		! mass of entrepreneurs with income equal to percentile value
		mass_entre_at_q0 = weight0*sum(mu(:,:,:,2:4),mask=(pregov_inc(:,:,:,2:4)==inc_toppctile(iq)))    
		! total mass of entrepreneurs at the top percentile iq
		mass_entre_q = mass_entre_between_q + mass_entre_at_q0
		! share of entrepreneurs at the top percentile iq
		share_entre_top_inc(iq) = mass_entre_q/mass_young_toppct(iq)
		
		! - Share of lfo by quintile
		do io = 2,4
			! mass of lfo iz with income strictly above the percentile value
			mass_lfo_between_q = sum(mu(:,:,:,io),mask=(pregov_inc(:,:,:,io)>inc_toppctile(iq)))    
			! mass of lfo iz with income equal to percentile value
			mass_lfo_at_q0 = weight0*sum(mu(:,:,:,io),mask=(pregov_inc(:,:,:,io)==inc_toppctile(iq)))    
			! total mass of lfo iz at the top percentile iq
			mass_lfo_q = mass_lfo_between_q + mass_lfo_at_q0
			! share of lfo iz at the top percentile iq
			share_lfo_top_inc(io,iq) = mass_lfo_q/mass_entre_q
		enddo ! io
    enddo ! iq 

    
    end subroutine share_entre_by_inc
    
!===============================================================================!

    subroutine taxes_by_inc()
    implicit none
    ! Computes tax revenue by income percentile
    ! Outputs: taxes_top = [taxes_top1, taxes_top3, taxes_top5, taxes_top10]
    real(8), allocatable :: mu_vec(:), inc_vec(:)
    real(8) :: qw(4), inc_p(4),iq, mass_top(4)

    ! Vectorize mu and pregov_inc
	mu_vec = reshape(mu(:,:,:,1:4),[na*neps*ntheta*4])
	inc_vec = reshape(pregov_inc(:,:,:,1:4),[na*neps*ntheta*4])
	
    ! All percentiles we want to look at
    qw = [0.99d0,0.97d0,0.95d0,0.9d0] 
	! Income percentiles corresponding to qw
    inc_p = quantili(inc_vec,mu_vec,qw)
    if (size(qw)/=size(inc_p)) then
        call myerror("taxes_by_inc: perc. size not correct")
    endif
    
    taxes_top = 0.0d0
    mass_top  = 0.0d0
    
	! Loop over state space: (a,eps,theta) x (young or old)
	do ia = 1,na
		do it = 1,ntheta
			do ie = 1,neps
				eps_val    = eps_grid(ie)
				theta_val  = theta_grid(it)
				a_val      = a_grid(ia)
				
				do iq = 1,4 
					! - Workers
					io = 1
					if (pregov_inc(ia,ie,it,io)>=inc_p(iq)) then
					! If worker is a top earner
						! (income tax) taxable income
						yw = w*eps_val*pol%lpoly_w(ia,ie,it) - fun_tax_pay(w*eps_val*pol%lpoly_w(ia,ie,it),ss_cap) + r*a_val
						! Workers pay income and ss taxes
						taxes_top(iq) =  taxes_top(iq)  &
							+ fun_tax_pay(w*eps_val*pol%lpoly_w(ia,ie,it),ss_cap)*mu(ia,ie,it,io) &
							+ fun_tax_inc(yw,y_H)*distrib%muy(ia,ie,it,io)
						mass_top(iq)  = mass_top(iq) + distrib%muy(ia,ie,it,io)
					endif
			
					! - Sole Proprietors
					io = 2
					if (pregov_inc(ia,ie,it,io)>=inc_p(iq)) then
					! If EP is a top earner
						k_val   = pol%kpol(ia,it,io)
						n_val   = pol%npol(ia,it,io) 
						profit_val = fun_profit(theta_val,k_val,le,n_val,r,w)
						yw = profit_val - fun_tax_pay(profit_val,ss_cap) + r*a_val
						taxes_top(iq) =  taxes_top(iq) &
							+ fun_tax_pay(profit_val,ss_cap)*mu(ia,ie,it,io) &
							+ fun_tax_inc(yw,y_H)*mu(ia,ie,it,io)
						mass_top(iq)  = mass_top(iq) + mu(ia,ie,it,io)
					endif				
					! - S-Corps
					io = 3
					if (pregov_inc(ia,ie,it,io)>=inc_p(iq)) then
					! If ES is a top earner
						k_val   = pol%kpol(ia,it,io)
						n_val   = pol%npol(ia,it,io) 
						profit_val = fun_profit(theta_val,k_val,le,n_val,r,w)
						phi_val  = pol%phipol(ia,it,io)
						yw = profit_val - fun_tax_pay(phi_val*profit_val,ss_cap) + r*a_val
						taxes_top(iq) =  taxes_top(iq) &
							+ fun_tax_pay(phi_val*profit_val,ss_cap)*mu(ia,ie,it,io) &
							+ fun_tax_inc(yw-cost_evasion_es(1.0d0-phi_val)-op_cost_es,y_H) &
									*mu(ia,ie,it,io)
						mass_top(iq)  = mass_top(iq) + mu(ia,ie,it,io)
					endif				
					! - C-Corps
					io = 4
					if (pregov_inc(ia,ie,it,io)>=inc_p(iq)) then
					! If EC is a top earner
						k_val   = pol%kpol(ia,it,io)
						n_val   = pol%npol(ia,it,io) 
						profit_val = fun_profit(theta_val,k_val,le,n_val,r,w)
						phi_val  = pol%phipol(ia,it,io)
						yw = phi_val*profit_val - fun_tax_pay(phi_val*profit_val,ss_cap) + r*a_val                
						taxes_top(iq) =  taxes_top(iq) & 
							+  fun_tax_pay(phi_val*profit_val,ss_cap)*mu(ia,ie,it,io) &
							+ fun_tax_inc(yw-cost_evasion_ec(phi_val)-op_cost_ec,y_H) &
										*mu(ia,ie,it,io) &
							+ tau_c*(1.0d0-phi_val)*profit_val*mu(ia,ie,it,io) &
							+ tau_d*(1.0d0-tau_c)*(1.0d0-phi_val)*profit_val*mu(ia,ie,it,io)
						mass_top(iq)  = mass_top(iq) + mu(ia,ie,it,io)
					endif
				enddo ! iq
			enddo ! ie
		enddo ! it
	enddo ! ia
    
    ! Adjust taxes_top if the top x percent includes more or less than exactly x percent due to 
    !	model discreteness.
    taxes_top = taxes_top / mass_top * (1.0d0 - qw)*sum(mu)
    

	end subroutine taxes_by_inc
!===============================================================================!

    subroutine atr_by_inc()
    implicit none
    ! Computes average tax rate (atr) by income percentiles
    ! Outputs: atr_vec  (1) P99.99+ (2) P99.9+, (3) P99-P99.9, (4) P90-P99, (5) P50-P90, (6) up to P50
    ! Note: we exclude negative income taxes from the calculation of the total income taxes. 
    real(8), allocatable :: mu_vec(:), taxableinc_vec(:),inctax_vec(:)
    real(8) :: qw(n_atr), inc_p(0:n_atr), tot_taxableinc(n_atr),tot_inctax(n_atr),temp_vec(n_atr)
	integer :: istat,iq
	
    !write(*,*) "I AM IN atr_by_inc"

	allocate(mu_vec(na*neps*ntheta*no+na),taxableinc_vec(na*neps*ntheta*no+na), &
		inctax_vec(na*neps*ntheta*no+na),stat=istat )
    if (istat/=0) then
        call myerror("atr_by_inc: Allocation failed!")
    endif
	
	! Vectorize mu, taxableinc, and inctax for both the young and reitrees
	! Young
    !write(*,*) "atr_by_inc: before reshape"
	mu_vec(1:na*neps*ntheta*no)         = reshape(mu,[na*neps*ntheta*no])
	taxableinc_vec(1:na*neps*ntheta*no) = reshape(taxableinc(:,:,:,:),[na*neps*ntheta*no])
	inctax_vec(1:na*neps*ntheta*no)     = reshape(inctax(:,:,:,:),[na*neps*ntheta*no])
	!write(*,*) "atr_by_inc: after reshape"
    ! Retirees
	mu_vec(na*neps*ntheta*no+1:na*neps*ntheta*no+na)         = mur
	taxableinc_vec(na*neps*ntheta*no+1:na*neps*ntheta*no+na) = taxableincr
	inctax_vec(na*neps*ntheta*no+1:na*neps*ntheta*no+na)     = inctaxr
		
	! All percentiles we want to look at
    qw = [0.9995d0, 0.999d0, 0.99d0,0.90d0,0.5d0,0.0d0] 
	! Income percentiles corresponding to qw
    !write(*,*) "atr_by_inc: now calling function quantili"
    inc_p(1:n_atr) = quantili(taxableinc_vec,mu_vec,qw)
    inc_p(0)     = huge(0.0d0)
    inc_p(n_atr) = -huge(0.0d0)
    if (size(qw)/=size(inc_p)-1) then
        call myerror("atr_by_inc: perc. size not correct")
    endif
    
    ! Loop over state space: (a,eps,theta) x (young or old)
    tot_taxableinc = 0.0d0
    tot_inctax     = 0.0d0
    temp_vec = 0.0d0
	do ia = 1,na
		! Young
		do it = 1,ntheta
			do ie = 1,neps
				do io = 1,no
					do iq = 1,n_atr
						if (taxableinc(ia,ie,it,io)>=inc_p(iq) .and. taxableinc(ia,ie,it,io)<inc_p(iq-1)) then
							tot_taxableinc(iq) = tot_taxableinc(iq) + taxableinc(ia,ie,it,io)*mu(ia,ie,it,io)
							tot_inctax(iq)     = tot_inctax(iq)     + max(0.0d0,inctax(ia,ie,it,io))*mu(ia,ie,it,io)
                            temp_vec(iq)      = temp_vec(iq)       + mu(ia,ie,it,io)
						endif
					enddo ! iq
				enddo ! io
			enddo ! ie
		enddo ! it
		! Retirees
		do iq = 1,n_atr
			if (taxableincr(ia)>=inc_p(iq) .and. taxableincr(ia)<inc_p(iq-1)) then
				tot_taxableinc(iq) = tot_taxableinc(iq) + taxableincr(ia)*mur(ia)
				tot_inctax(iq)     = tot_inctax(iq)     + max(0.0d0,inctaxr(ia))*mur(ia)
                temp_vec(iq)      = temp_vec(iq)       + mur(ia)
			endif
		enddo ! iq
	enddo ! ia
    
    ! Calculate averate tax rate as the sum of total income taxes divided by the sum of total taxable incomes
    !write(*,*) 'inc_p',inc_p
    !write(*,*) 'temp_vec', temp_vec
    !write(*,*) 'tot_inctax', tot_inctax
    !write(*,*) 'tot_taxableinc', tot_taxableinc
    !pause
    atr_vec = tot_inctax/tot_taxableinc
    
	
	end subroutine atr_by_inc
!===============================================================================!

    subroutine sub_share_toptax()
    implicit none
    ! Purpose: Compute the share of young agents paying the top income tax rate (tau_h)
    ! output: share_toptax
    ! Local variables
    real(8),allocatable :: mu_vec(:),inc_vec(:),xs(:),ws(:),cums(:),xs_u(:),cums_u(:)
    integer, allocatable :: ix(:),ind_u(:)
    integer :: n,n_u,istat
    ! Vectorize mu and pregov_inc
	mu_vec = reshape(mu(:,:,:,1:4),[na*neps*ntheta*4])
	inc_vec = reshape(taxableinc(:,:,:,1:4),[na*neps*ntheta*4])
    ! Find the percentile correspond to y_H
    
    n = size(inc_vec)
    allocate( ix(n), xs(n),ws(n),cums(n),stat=istat )
    if (istat/=0) then
        call myerror("sub_share_toptax: Allocation failed!")
    endif
    
    !xs is x sorted, ix is the sorting index
    xs = inc_vec 
    call quick_sort(xs,ix)
    ws=mu_vec(ix)
    ws=ws/sum(ws)
    cums=cumsum(ws)
    
    ![xs_u,ind_u] = unique(xs,'last')
    call unique(xs, xs_u,ind_u)
    cums_u = cums(ind_u)
    n_u = size(xs_u)
   ! share_toptax: share of young population paying top marginal tax rate (tau_h)
	if (y_H<=xs_u(1)) then
		share_toptax = 1.0d0
	elseif (y_H>=xs_u(n_u)) then
		share_toptax = 0.0d0
	else
		share_toptax = 1.0d0-linint(xs_u,cums_u,y_H)  
	endif
  	end subroutine sub_share_toptax
    !===============================================================================!
    
    subroutine share_entre_by_wealth()
    implicit none
    ! Purpose: Compute the share of entre and LFO distribution in the top percentiles 
    !	of the wealth distrib (top 1%, 5%, and 10%), and by wealth quintiles
    ! Note: consider active population only!
    
    ! Local variables:
    real(8) :: qw(7), wealth_p(7),mu_a_vec(na),ia_p(7),positive_tiny
    integer :: ieps,itheta,iq
    real(8) :: wealth_quintiles(0:5),mass_young_q,mass_young_toppct(3),weight0,weight1, &
    	mass_between_q,mass_at_q0,mass_at_q1,mass_entre_between_q,mass_entre_at_q0,mass_entre_at_q1, &
    	mass_entre_q,mass_lfo_between_q,mass_lfo_at_q0,mass_lfo_at_q1,mass_lfo_q,wealth_toppctile(3)
	! Define a tiny positive number
	positive_tiny = 1d-10
	
	! Marginal dist of mu wrt wealth
	do ia = 1,na
		mu_a_vec(ia) = sum(mu(ia,:,:,1:4))
	enddo
    ! All percentiles we want to look at
    qw = [0.2d0,0.4d0,0.6d0,0.8d0,0.9d0,0.95d0,0.99d0] 
	! Wealth percentiles corresponding to qw
    wealth_p = quantili(a_grid,mu_a_vec,qw)
    if (size(qw)/=size(wealth_p)) then
        call myerror("share_entre_by_wealth: perc. size not correct")
    endif
    
    ! Find the closest a_grid index for each wealth_p
    do iq = 1,size(qw)
    	ia_p(iq) = my_closest(a_grid,na,wealth_p(iq))
    enddo
    
    ! Find share of entre and LFO dist for
    !	- By quintiles: 0-20%, 20-40%,...,80-100%
    !	- Top percentiles: 1%, top 5%, top 10%

    
    ! By quintiles

    wealth_quintiles(0) = minval(a_grid) - 1.0d0
    wealth_quintiles(1:4) = wealth_p(1:4)
    wealth_quintiles(5) = maxval(a_grid) + 1.0d0
   
    ! mass of young agent per quintile
    mass_young_q = sum(mu(:,:,:,1:4))/5.0d0
    ! mass of young agent at the top 1, 5, and 10%
    mass_young_toppct(1) = sum(mu(:,:,:,1:4)) * 0.01d0
    mass_young_toppct(2) = sum(mu(:,:,:,1:4)) * 0.05d0
    mass_young_toppct(3) = sum(mu(:,:,:,1:4)) * 0.10d0

    weight0 = 0.0d0
    weight1 = 0.0d0
    do iq = 1,5
		! mass of young agent with wealth strictly between two quintile values
		mass_between_q = 0.0d0
		mass_at_q0 = 0.0d0
		mass_at_q1 = 0.0d0
		do ieps = 1,neps
		do itheta = 1,ntheta
		do io = 1,4
			mass_between_q = mass_between_q + sum(mu(:,ieps,itheta,io),mask=(a_grid>wealth_quintiles(iq-1) .and. &
				a_grid<wealth_quintiles(iq)))    
			! mass of young agent with wealth equal to the lower and higher quintile value
			mass_at_q0 = mass_at_q0 + sum(mu(:,ieps,itheta,io),mask=(a_grid==wealth_quintiles(iq-1)))    
			mass_at_q1 = mass_at_q1 + sum(mu(:,ieps,itheta,io),mask=(a_grid==wealth_quintiles(iq)))    
		enddo
		enddo
		enddo
		! what fraction of mass_at_q0 should be assigned to quintile iq?
		weight0 = 1.0d0 - weight1
		! what fraction of mass_at_q1 should be assigned to quintile iq?
		if (mass_at_q1 >0.0d0) then
			weight1 = (mass_young_q-mass_between_q - weight0*mass_at_q0)/mass_at_q1
		else
			weight1 = 1.0d0
		endif
		! weight1 must be between 0 and 1
		if (weight1<0.0d0 .or. weight1 > 1.0d0) then
			write(*,*) "weight1 = ",weight1
			write(*,*) "mass_young_q = ",mass_young_q
			write(*,*) "mass_between_q = ",mass_between_q
			write(*,*) "weight0 = ",weight0
			write(*,*) "mass_at_q0 = ",mass_at_q0
			write(*,*) "mass_at_q1 = ",mass_at_q1
			!pause
		endif
		
		! - Share of entrepreneurs by quintile
		! mass of entrepreneurs with wealth strictly between two quintile values
		mass_entre_between_q = 0.0d0
		mass_entre_at_q0 = 0.0d0
		mass_entre_at_q1 = 0.0d0
		do ieps = 1,neps
		do itheta = 1,ntheta
		do io = 2,4
			mass_entre_between_q = mass_entre_between_q + sum(mu(:,ieps,itheta,io),mask=(a_grid>wealth_quintiles(iq-1) .and. &
				a_grid<wealth_quintiles(iq)))    
			! mass of entrepreneurs with wealth equal to the lower quintile value
			mass_entre_at_q0 = mass_entre_at_q0 + weight0*sum(mu(:,ieps,itheta,io),mask=(a_grid==wealth_quintiles(iq-1)))    
			mass_entre_at_q1 = mass_entre_at_q1 + weight1*sum(mu(:,ieps,itheta,io),mask=(a_grid==wealth_quintiles(iq)))    
		enddo
		enddo
		enddo 
		! total mass of entrepreneurs in quintile iq
		mass_entre_q = mass_entre_between_q + mass_entre_at_q0 + mass_entre_at_q1
		! share of entrepreneurs in quintile iq
		share_entre_quint(iq) = mass_entre_q/mass_young_q
		
		! - Share of lfo by quintile
		do io = 2,4
			mass_lfo_between_q = 0.0d0
			mass_lfo_at_q0 = 0.0d0
			mass_lfo_at_q1 = 0.0d0
			do ieps = 1,neps
			do itheta = 1,ntheta
				! mass of lfo iz with wealth strictly between two quintile values
				mass_lfo_between_q = mass_lfo_between_q + sum(mu(:,ieps,itheta,io),mask=(a_grid>wealth_quintiles(iq-1) .and. &
					a_grid<wealth_quintiles(iq)))    
				! mass of lfo iz with wealth equal to the lower quintile value
				mass_lfo_at_q0 = mass_lfo_at_q0 + weight0*sum(mu(:,ieps,itheta,io),mask=(a_grid==wealth_quintiles(iq-1)))    
				mass_lfo_at_q1 = mass_lfo_at_q1 + weight1*sum(mu(:,ieps,itheta,io),mask=(a_grid==wealth_quintiles(iq)))    
			enddo ! itheta
			enddo ! ieps
			! total mass of lfo iz in quintile iq
			mass_lfo_q = mass_lfo_between_q + mass_lfo_at_q0 + mass_lfo_at_q1
			! share of lfo iz in quintile iq
			share_lfo_quint(io,iq) = mass_lfo_q/(positive_tiny+mass_entre_q)
		enddo ! iz
    enddo ! iq
    ! up to here!
    
    ! Top percentiles
    !	- [top 1%, top 5%, top 10%]
    
    wealth_toppctile(1) = wealth_p(7) ! top 1%
    wealth_toppctile(2) = wealth_p(6) ! top 5%
    wealth_toppctile(3) = wealth_p(5) ! top 10%
    
    weight0 = 0.0d0
    do iq = 1,3
    	mass_at_q0 = 0.0d0
    	mass_between_q = 0.0d0
    	do ieps = 1,neps
    	do itheta = 1,ntheta
    	do io = 1,4
			! mass of young agent with wealth strictly above the percentile value
			mass_between_q = mass_between_q + sum(mu(:,ieps,itheta,io),mask=(a_grid>wealth_toppctile(iq)))    
			! mass of young agent with wealth equal to the percentile value
			mass_at_q0 = mass_at_q0 + sum(mu(:,ieps,itheta,io),mask=(a_grid==wealth_toppctile(iq)))    
		enddo
		enddo
		enddo
		
		! what fraction of mass_at_q0 should be assigned to the top percentile iq?
		if (mass_at_q0 >0.0d0) then
			weight0 = (mass_young_toppct(iq)-mass_between_q)/mass_at_q0
		else
			weight0 = 1.0d0
		endif
		! weight0 must be between 0 and 1
		if (weight0<0.0d0 .or. weight0 > 1.0d0) then
			write(*,*) "weight0 = ",weight0
			write(*,*) "mass_young_toppct(iq) = ",mass_young_toppct(iq)
			write(*,*) "mass_between_q = ",mass_between_q
			write(*,*) "mass_at_q0 = ",mass_at_q0
			!pause
		endif
		
		! - Share of entrepreneurs at top percentile
		mass_entre_between_q = 0.0d0
		mass_entre_at_q0 = 0.0d0
		do ieps = 1,neps
		do itheta = 1,ntheta
		do io = 2,4
			! mass of entrepreneurs with wealth strictly above the percentile value
			mass_entre_between_q = mass_entre_between_q + sum(mu(:,ieps,itheta,io),mask=(a_grid>wealth_toppctile(iq)))    
			! mass of entrepreneurs with wealth equal to percentile value
			mass_entre_at_q0 = mass_entre_at_q0 + weight0*sum(mu(:,ieps,itheta,io),mask=(a_grid==wealth_toppctile(iq)))    			
		enddo
		enddo
		enddo
		! total mass of entrepreneurs at the top percentile iq
		mass_entre_q = mass_entre_between_q + mass_entre_at_q0
		! share of entrepreneurs at the top percentile iq
		share_entre_top(iq) = mass_entre_q/mass_young_toppct(iq)
		! up to here!
		! - Share of lfo by quintile
		do io = 2,4
			mass_lfo_at_q0 = 0.0d0
			mass_lfo_between_q = 0.0d0
			do ieps = 1,neps
			do itheta = 1,ntheta
				! mass of lfo iz with income strictly above the percentile value
				mass_lfo_between_q = mass_lfo_between_q + sum(mu(:,ieps,itheta,io),mask=(a_grid>wealth_toppctile(iq)))    
				! mass of lfo iz with income equal to percentile value
				mass_lfo_at_q0 = mass_lfo_at_q0 + weight0*sum(mu(:,ieps,itheta,io),mask=(a_grid==wealth_toppctile(iq)))    
			enddo
			enddo
			! total mass of lfo iz at the top percentile iq
			mass_lfo_q = mass_lfo_between_q + mass_lfo_at_q0
			! share of lfo iz at the top percentile iq
			share_lfo_top(io,iq) = mass_lfo_q/(positive_tiny + mass_entre_q)
		enddo ! io
    enddo ! iq 
    
    end subroutine share_entre_by_wealth
!===============================================================================!
    
    subroutine share_emp_firmsize()
    implicit none
    ! Purpose: Computes employment share by firm size bins
	! * Definition: Firm size = nbar
	! Declare locals
	real(8) :: mu_entre(na,ntheta,2:4), qvec(3), n_q(3)
    !nbar_vec(na*ntheta*3),mu_vec(na*ntheta*3)
	real(8), allocatable :: n_vec(:), mu_vec(:)
    
    ! construct mu_entre
    do io = 2,4
    do it = 1,ntheta
    do ia = 1,na
    	mu_entre(ia,it,io) = sum(mu(ia,:,it,io))
    enddo
    enddo
    enddo
    ! vectorize nbarpol and mu_entre
    n_vec = reshape(pol%npol(:,:,2:4),[na*ntheta*3])
    mu_vec = reshape(mu_entre,[na*ntheta*3])
    ! bin cutoffs
    qvec = [0.6211d0,0.7935d0,0.8978d0]
    n_q = quantili(n_vec,mu_vec,qvec)

    ! Employment shares by firm size bins
    
    share_emp = 0.0d0
    
    do io = 2,4
    do it = 1,ntheta
    do ia = 1,na
    	if (pol%npol(ia,it,io)<=n_q(1)) then ! Bin1
    		share_emp(1) = share_emp(1) + pol%npol(ia,it,io)*mu_entre(ia,it,io) 
    	elseif (pol%npol(ia,it,io)>n_q(1) .and. pol%npol(ia,it,io)<=n_q(2) ) then ! Bin2
    		share_emp(2) = share_emp(2) +pol%npol(ia,it,io)*mu_entre(ia,it,io) 
    	elseif (pol%npol(ia,it,io)>n_q(2) .and. pol%npol(ia,it,io)<=n_q(3) ) then ! Bin3
    		share_emp(3) = share_emp(3) + pol%npol(ia,it,io)*mu_entre(ia,it,io) 
    	else
    		share_emp(4) = share_emp(4) + pol%npol(ia,it,io)*mu_entre(ia,it,io) 
    	endif
    enddo
    enddo
    enddo
    share_emp = share_emp/sum(share_emp)
	end subroutine share_emp_firmsize   
!===============================================================================!
    
    subroutine occ_trans()
    implicit none      
    ! This subroutine computes the transition rate between workers and entrepreneurs
    ! Declare locals
    real(8) :: temp(neps,ntheta),omega,aopt_val,iop_prob,trans_mat_weighted(no,no)
    integer :: aopt_ind
    
    ! trans_mat: element (i,j) = transition rate from state i to j, conditional on not retiring. i,j between 1 and 4
    
    trans_mat = 0.0d0
    trans_mat_weighted = 0.0d0
    do it = 1,ntheta
	do ie = 1,neps
		!temp is (neps,ntheta) matrix
        temp = outerprod(P_eps(ie,:),P_theta(it,:))
        
		do ia = 1,na
			do io = 1,4
				! Asset for tomorrow?
				aopt_val = pol%apoly(ia,ie,it,io)
				aopt_ind = max(min(locate(a_grid,aopt_val),na-1),1)
				! Split weight between aopt_ind and aopt_ind+1
				omega = (a_grid(aopt_ind+1)-aopt_val)/(a_grid(aopt_ind+1)-a_grid(aopt_ind))
				omega = max(min(omega,1d0),0d0)
				! Occ choice tomorrow
				do itp = 1,ntheta
				do iep = 1,neps            	
					do iop = 1,4
					! Probability of occupation izp given aopt_ind,iep,itp
					iop_prob = pol%occpol(aopt_ind,iep,itp,io,iop)
					! Update trans_mat
					trans_mat(io,iop) = trans_mat(io,iop) + omega*temp(iep,itp)*iop_prob*mu(ia,ie,it,io)
					trans_mat_weighted(io,iop) = trans_mat_weighted(io,iop) + omega*temp(iep,itp)*iop_prob*mu(ia,ie,it,io)*pol%npol(ia,it,io)
					! Probability of occupation izp given aopt_ind+1,iep,itp
					iop_prob = pol%occpol(aopt_ind+1,iep,itp,io,iop)
					! Update trans_mat
					trans_mat(io,iop) = trans_mat(io,iop) + (1.0d0-omega)*temp(iep,itp)*iop_prob*mu(ia,ie,it,io)
                    trans_mat_weighted(io,iop) = trans_mat_weighted(io,iop) + (1.0d0-omega)*temp(iep,itp)*iop_prob*mu(ia,ie,it,io)*pol%npol(ia,it,io)
					enddo !izp
				enddo ! itp
				enddo ! iep
			enddo ! io
		enddo !ia
    enddo !ie
    enddo !it
    ! Transition rate from entre to worker
    trans_entre_work = sum(trans_mat(2:4,1))/sum(trans_mat(2:4,:))
    ! Transition rate from worker to entre
    trans_work_entre = sum(trans_mat(1,2:4))/sum(trans_mat(1,:))
	! Transition rate from C to Passthroughs (S or P), conditional on remaining an entrepreneur
    trans_CP = sum(trans_mat_weighted(4,2:3))/sum(trans_mat_weighted(4,2:4))
    ! Transition rate from Passthroughs (S or P) to C, conditional on remaining an entrepreneur
    trans_PC = trans_mat_weighted(2,4)/sum(trans_mat_weighted(2,2:4))
    trans_SC = trans_mat_weighted(3,4)/sum(trans_mat_weighted(3,2:4))



	! 4 by 4 transition rate
	do io = 1,4
		trans_mat(io,:) = trans_mat(io,:)/sum(trans_mat(io,:))    
    enddo
    
    end subroutine occ_trans
    

!===============================================================================!

    !-------------------------------------------------------------------------------!
    
    end subroutine targets_compute
    
    
    
    !===============================================================================!
    
    subroutine txt_export(filename,prices, agg,model_targets)
    ! Purpose: WRITES PARAMETERS AND TARGETS TO A FILE

    !Declare inputs:
    character(LEN=*), intent(in)    :: filename
    type(modelAgg), intent(in) :: agg
    type(modelResults), intent(in) :: model_targets
    type(modelPrices),intent(in) :: prices
    !Note: model parameters and flags are defined as globals
    
    !Declare locals:
    integer :: uno,io


    open(newunit=uno,file=trim(filename), status='replace')

    !Two blank lines (align with excel)
    write(uno,'(a)')' '

    write(uno, '(a)')'FLAGS'
    write(uno, '(a,I)')'do_GE        : ',do_GE 
    write(uno, '(a,I)')'cf_avoidance  : ',cf_avoidance 
    write(uno, '(a,I)')'cf_occpol  : ',cf_occpol
    write(uno, '(a,I)')'howard       : ',do_howard
    write(uno, '(a,f12.6)')'a_space : ',a_space
    write(uno, '(a,f12.6)')'G_frac : ',G_frac
    write(uno, '(a,f12.6)')'G_bench : ',G_bench
    
    
    write(uno, '(a,l2)')'super_shock        : ',super_shock 
    write(uno, '(a,f12.6)')'eps_super : ',eps_super
    write(uno, '(a,f12.6)')'prob_super : ',prob_super
    write(uno, '(a,f12.6)')'prob_super_back : ',prob_super_back
    
    write(uno, '(a)')' '
    write(uno, '(a)')'GRIDS'
    write(uno, '(a,I4)')'na     : ',na
    write(uno, '(a,I4)')'neps   : ',neps 
    write(uno, '(a,I4)')'ntheta : ',ntheta
    write(uno, '(a,f12.6)')'a_max : ',a_max
	write(uno, '(a,f12.6)')'l_max : ',l_max
	write(uno, '(a,f12.6)')'le : ',le
    write(uno, '(a,I4)')'nx_max_nonconvex : ',nx_max_nonconvex
    write(uno, '(a,f12.9)')'tolV : ',tolV
	
    write(uno, '(a)')' '
    write(uno, '(a)')'EXTERNAL PARAMETERS'
    write(uno, '(a,f12.6)')'sigma1      : ',sigma1
    write(uno, '(a,f12.6)')'sigma2     : ',sigma2
    write(uno, '(a,f12.6)')'alpha      : ',alpha
    write(uno, '(a,f12.6)')'lambda     : ',lambda
    write(uno, '(a,f12.6)')'lambda_es  : ',lambda_es
    write(uno, '(a,f12.6)')'lambda_ec  : ',lambda_ec
    write(uno, '(a,f12.6)')'rho_eps    : ',rho_eps
    write(uno, '(a,f12.6)')'sig_eps    : ',sig_eps
    write(uno, '(a,f12.6)')'pr_ret     : ',pr_ret
    write(uno, '(a,f12.6)')'pr_die     : ',pr_die
	write(uno, '(a,f12.6)')'sig_e      : ',sig_e
	
    write(uno, '(a)')' '
    write(uno, '(a)')'TAX PARAMETERS'
    write(uno, '(a,I)')'taxfunc          : ',taxfunc
    write(uno, '(a,f12.6)')'tau_c        : ',tau_c
    write(uno, '(a,f12.6)')'tau_p        : ',tau_p
    write(uno, '(a,f12.6)')'ss_cap_ratio : ',ss_cap_ratio
    !write(uno, '(a,f12.6)')'tau_med  : ',tau_med
    !write(uno, '(a,f12.6)')'tau_sur  : ',tau_sur
    write(uno, '(a,f12.6)')'tau_d        : ',tau_d
    !write(uno, '(a,f12.6)')'top_tau  : ',tau_m(n_brackets)
    ! - Gouveia and Strauss (maybe not needed)
    write(uno, '(a,f12.6)')'lambda (HSV) : ',hsv_0
    write(uno, '(a,f12.6)')'tau (HSV)    : ',hsv_1
    write(uno, '(a,f12.6)')'repl         : ',repl
	write(uno, '(a,f12.6)')'tau_h        : ',tau_h
	write(uno, '(a,f12.6)')'y_H          : ',prices%y_H
	write(uno, '(a,f12.6)')'share_toptax : ', model_targets%share_toptax

    write(uno, '(a)')' '
    write(uno, '(a)')'ESTIMATED PARAMETERS'
    write(uno, '(a,f12.6)')'beta      : ',beta
    write(uno, '(a,f12.6)')'chi       : ',chi
    write(uno, '(a,f12.6)')'delta     : ',delta
    write(uno, '(a,f12.6)')'vi        : ',vi
    write(uno, '(a,f12.6)')'gamma     : ',gamma

    ! theta_shock AR1
    write(uno, '(a,f12.6)')'rho_theta :     ',rho_theta
    write(uno, '(a,f12.6)')'sig_theta :     ',sig_theta
    write(uno, '(a,f12.6)')'uncmean_theta : ',uncmean_theta

    ! theta_shock Bettina
    !write(uno, '(a,f12.6)')'theta_bar  :    ',theta_bar
    !write(uno, '(a,f12.6)')'theta_hat  :    ',theta_hat

    ! theta_shock Buera-Moll
    !write(uno, '(a,f12.6)')'psi        :    ',psi
    !write(uno, '(a,f12.6)')'shape      :    ',shape
    !write(uno, '(a,f12.6)')'theta_min  :    ',theta_min
    !write(uno, '(a,f12.6)')'theta_max  :    ',theta_max

    write(uno, '(a,f12.6)')'c0_es      :    ',c0_es
    write(uno, '(a,f12.6)')'c1_es      :    ',c1_es
    write(uno, '(a,f12.6)')'c0_ec      :    ',c0_ec
    write(uno, '(a,f12.6)')'c1_ec      :    ',c1_ec
    write(uno, '(a,f12.6)')'op_cost_es :    ',op_cost_es
    write(uno, '(a,f12.6)')'op_cost_ec :    ',op_cost_ec
	write(uno, '(a,f12.6)')'switch_cost_CP :    ',switch_cost_CP
	write(uno, '(a,f12.6)')'switch_cost_PC :    ',switch_cost_PC

    write(uno, '(a)')' '
    write(uno, '(a)')'Aggregates'
    write(uno, '(a,f20.9)')'ED         : ',agg%ED
    write(uno, '(a,f20.9)')'ED goods (in output units): ',agg%ED_goods/agg%Y
    write(uno, '(a,f20.9)')'bal_ss     : ',agg%bal_ss
    write(uno, '(a,f12.9)')'res_gov    : ',agg%res_gov
    write(uno, '(a,f12.9)')'tax_inc_gdp: ',agg%taxes_inc/agg%Y
    write(uno, '(a,f12.9)')'r          : ',prices%r
    write(uno, '(a,f12.6)')'w          : ',prices%w
    write(uno, '(a,f12.6)')'pen        : ',prices%pen
    write(uno, '(a,f12.6)')'tr         : ',prices%tr
    write(uno, '(a,f12.6)')'ss_cap     : ',prices%ss_cap
    write(uno, '(a,f12.6)')'bal_gov    : ',agg%bal_gov
    write(uno, '(a,f12.6)')'y_ave      : ',prices%y_ave
    write(uno, '(a,f12.9)')'A              : ',agg%A
    write(uno, '(a,f12.9)')'Y              : ',agg%Y
    write(uno, '(a,f12.9)')'Y_C         : ',agg%Y_C
    write(uno, '(a,f12.9)')'Y_entre              : ',agg%Y_entre
    write(uno, '(a,f12.9)')'Y_entre_share              : ',agg%Y_entre/agg%Y
    write(uno, '(a,f12.9)')'KY              : ',agg%KY   
    write(uno, '(a,f12.6)')'ave_l_work      : ', model_targets%ave_l_work
    write(uno, '(a,f12.9)')'taxes_ss/taxes_tot        : ',agg%taxes_ss/agg%taxes_tot   
    write(uno, '(a,f12.9)')'taxes_inctaxes_tot       : ',agg%taxes_inc/agg%taxes_tot  
    write(uno, '(a,f12.9)')'taxes_corptaxes_tot      : ',agg%taxes_corp/agg%taxes_tot   
    write(uno, '(a,f12.9)')'taxes_divtaxes_tot       : ',agg%taxes_div/agg%taxes_tot    
    write(uno, '(a,f12.9)')'taxes_tot       : ',agg%taxes_tot  
    
    write(uno, '(a,f12.9)')'taxes_inc_corp_div       : ',agg%taxes_inc +  agg%taxes_corp + agg%taxes_div
    write(uno, '(a,f12.9)')'taxes_ss       : ',agg%taxes_ss  

    write(uno, '(a,f12.9)')'taxes_tot_gdp   : ',agg%taxes_tot/agg%Y
    write(uno, '(a,f12.9)')'pen_tot       : ',agg%pen_tot  
    write(uno, '(a,f12.9)')'taxes_occ(1)       : ',agg%taxes_occ(1)  
    write(uno, '(a,f12.9)')'taxes_occ(2)       : ',agg%taxes_occ(2)  
    write(uno, '(a,f12.9)')'taxes_occ(3)       : ',agg%taxes_occ(3)  
    write(uno, '(a,f12.9)')'taxes_occ(4)       : ',agg%taxes_occ(4)  
    write(uno, '(a,f12.9)')'taxes_r      : ',agg%taxes_r
    write(uno, '(a,f12.9)')'K_entre            : ',agg%K_entre
    write(uno, '(a,f12.9)')'K_C                : ',agg%K_C
    write(uno, '(a,f12.9)')'K                : ',agg%K_entre + agg%K_C
    write(uno, '(a,f12.9)')'N_C                : ',agg%N_C
    write(uno, '(a,f12.9)')'N_entre            : ',agg%N_entre
    write(uno, '(a,f12.9)')'labsup             : ',agg%labsup
    write(uno, '(a,f12.9)')'C                 : ',agg%C
    write(uno, '(a,f12.9)')'C_entre           : ',agg%C_entre
    write(uno, '(a,f12.9)')'C_worker            : ',agg%C_worker
    write(uno, '(a,f12.9)')'C_r             : ',agg%C_r
    write(uno, '(a,f12.9)')'avetheta_allentre  : ',agg%avetheta(1)
    write(uno, '(a,f12.9)')'avetheta_soleprop  : ',agg%avetheta(2)
    write(uno, '(a,f12.9)')'avetheta_Scorp     : ',agg%avetheta(3)
    write(uno, '(a,f12.9)')'avetheta_Ccorp     : ',agg%avetheta(4)
    write(uno, '(a,f12.9)')'avek_allentre      : ',agg%avek(1)
    write(uno, '(a,f12.9)')'avek_soleprop      : ',agg%avek(2)
    write(uno, '(a,f12.9)')'avek_Scorp         : ',agg%avek(3)
    write(uno, '(a,f12.9)')'avek_Ccorp          : ',agg%avek(4)
    write(uno, '(a,f12.9)')'aveklratio_allentre	 : ',agg%aveklratio(1)
    write(uno, '(a,f12.9)')'aveklratio_soleprop	 : ',agg%aveklratio(2)
    write(uno, '(a,f12.9)')'aveklratio_Scorp    : ',agg%aveklratio(3)
    write(uno, '(a,f12.9)')'aveklratio_Ccorp    : ',agg%aveklratio(4)
    
    write(uno, '(a,f12.9)')'avekyratio_allentre	 : ',agg%avekyratio(1)
    write(uno, '(a,f12.9)')'avekyratio_soleprop	 : ',agg%avekyratio(2)
    write(uno, '(a,f12.9)')'avekyratio_Scorp    : ',agg%avekyratio(3)
    write(uno, '(a,f12.9)')'avekyratio_Ccorp    : ',agg%avekyratio(4)

    write(uno, '(a,f12.9)')'N_share_C	 : ',agg%N_share(1)
    write(uno, '(a,f12.9)')'N_share_EP	 : ',agg%N_share(2)
    write(uno, '(a,f12.9)')'N_share_ES   : ',agg%N_share(3)
    write(uno, '(a,f12.9)')'N_share_EC   : ',agg%N_share(4)
    write(uno, '(a,f12.9)')'N_share_allC   : ',agg%N_share(1)+agg%N_share(4)

    write(uno, '(a,f12.9)')'K_share_C	 : ',agg%K_share(1)
    write(uno, '(a,f12.9)')'K_share_EP	 : ',agg%K_share(2)
    write(uno, '(a,f12.9)')'K_share_ES   : ',agg%K_share(3)
    write(uno, '(a,f12.9)')'K_share_EC   : ',agg%K_share(4)
    write(uno, '(a,f12.9)')'K_share_allC   : ',agg%K_share(1)+agg%K_share(4)

    write(uno, '(a,f12.9)')'labshare_ES	 : ',agg%labshare_ES
    write(uno, '(a,f12.9)')'labshare_EC	 : ',agg%labshare_EC
    write(uno, '(a,f12.9)')'labshare_corp	 : ',agg%labshare_corp
    write(uno, '(a,f12.9)')'labshare_allcorp	 : ',agg%labshare_allcorp

    write(uno, '(a)')' '
    write(uno, '(a)')'TARGETS and RESULTS'
    ! - Occ/LFO distribution
    write(uno, '(a)')'Occ/LFO distribution'
    write(uno, '(a,f12.6)')'share_entre_act : ', model_targets%share_entre_act
    write(uno, '(a,f12.6)')'share_EP_entre : ', model_targets%share_EP_entre
    write(uno, '(a,f12.6)')'share_ES_entre : ', model_targets%share_ES_entre
    write(uno, '(a,f12.6)')'share_EC_entre : ', model_targets%share_EC_entre
    write(uno, '(a,f12.6)')'share_active : ', model_targets%share_active
    ! - Share of income declared as wage (ES and EC)
    write(uno, '(a)')'Share of income declared as wage (ES and EC)'
    write(uno, '(a,f12.6)')'share_wage_ES : ', model_targets%share_wage_ES
    write(uno, '(a,f12.6)')'share_wage_EC : ', model_targets%share_wage_EC
    write(uno, '(a,f12.6)')'share_wage_EC_cf : ', model_targets%share_wage_EC_cf

    ! - Share of income declared as wage (ES and EC), by emp size bins
    write(uno, '(a,f12.6)')'share_wage_ES(bottom 40%), by emp: ', model_targets%share_wage_qn_ES(1)
    write(uno, '(a,f12.6)')'share_wage_ES(middle 40%), by emp : ', model_targets%share_wage_qn_ES(2)
    write(uno, '(a,f12.6)')'share_wage_ES(top 20%), by emp : ', model_targets%share_wage_qn_ES(3)
    write(uno, '(a,f12.6)')'share_wage_EC(bottom 40%), by emp : ', model_targets%share_wage_qn_EC(1)
    write(uno, '(a,f12.6)')'share_wage_EC(middle 40%), by emp : ', model_targets%share_wage_qn_EC(2)
    write(uno, '(a,f12.6)')'share_wage_EC(top 20%), by emp : ', model_targets%share_wage_qn_EC(3)
    ! - Share of income declared as wage (ES and EC), by capital size bins
    write(uno, '(a,f12.6)')'share_wage_ES(bottom 40%), by capital: ', model_targets%share_wage_qk_ES(1)
    write(uno, '(a,f12.6)')'share_wage_ES(middle 40%), by capital : ', model_targets%share_wage_qk_ES(2)
    write(uno, '(a,f12.6)')'share_wage_ES(top 20%), by capital : ', model_targets%share_wage_qk_ES(3)
    write(uno, '(a,f12.6)')'share_wage_EC(bottom 40%), by capital : ', model_targets%share_wage_qk_EC(1)
    write(uno, '(a,f12.6)')'share_wage_EC(middle 40%), by capital : ', model_targets%share_wage_qk_EC(2)
    write(uno, '(a,f12.6)')'share_wage_EC(top 20%), by capital : ', model_targets%share_wage_qk_EC(3)
    ! - Share facing financial constraints and leverages (EP,ES,EC)
    write(uno, '(a)')'Share facing financial constraints and leverages (EP,ES,EC)'
    write(uno, '(a,f12.6)')'frac_bor_ep : ', model_targets%frac_bor_ep
    write(uno, '(a,f12.6)')'frac_bor_es : ', model_targets%frac_bor_es
    write(uno, '(a,f12.6)')'frac_bor_ec : ', model_targets%frac_bor_ec
    write(uno, '(a,f12.6)')'leverage_ep : ', model_targets%leverage_ep
    write(uno, '(a,f12.6)')'leverage_es : ', model_targets%leverage_es
    write(uno, '(a,f12.6)')'leverage_ec : ', model_targets%leverage_ec
	! - Average income ratio
    write(uno, '(a)')'Average and median income ratio'
	write(uno, '(a,f12.6)')'ratio_aveinc_entre_worker : ', model_targets%ratio_aveinc_entre_worker
	write(uno, '(a,f12.6)')'ratio_medinc_entre_worker : ', model_targets%ratio_medinc_entre_worker

	! - Share of employer
    write(uno, '(a)')'Share of employer'
	write(uno, '(a,f12.6)')'share_empl_EP : ', model_targets%share_empl_EP
	write(uno, '(a,f12.6)')'share_empl_ES : ', model_targets%share_empl_ES
	write(uno, '(a,f12.6)')'share_empl_EC : ', model_targets%share_empl_EC
	write(uno, '(a,f12.6)')'share_empl_entre : ', model_targets%share_empl_entre
	! - Average firm size by LFO
	write(uno, '(a)')'Average firm size (# employees) rel to EP'
	write(uno, '(a,f12.6)')'empsize_rel_ES : ', model_targets%empsize_rel_ES
	write(uno, '(a,f12.6)')'empsize_rel_EC : ', model_targets%empsize_rel_EC
    write(uno, '(a,f12.6)')'payrollsize_lfo(EP) : ', model_targets%payrollsize_lfo(2)
    write(uno, '(a,f12.6)')'payrollsize_lfo(ES) : ', model_targets%payrollsize_lfo(3)
    write(uno, '(a,f12.6)')'payrollsize_lfo(EC) : ', model_targets%payrollsize_lfo(4)

    ! - Median size by lfo
    write(uno, '(a)')'Median income by lfo'
    do io = 2,4
        write(uno, '(a,f12.6)')'median_inc_lfo : ', model_targets%median_inc_lfo(io)
    enddo
    write(uno, '(a)')'Median net income by lfo'
    do io = 2,4
        write(uno, '(a,f12.6)')'median_netinc_lfo : ', model_targets%median_netinc_lfo(io)
    enddo
    write(uno, '(a)')'Median emp by lfo'
    do io = 2,4
        write(uno, '(a,f12.6)')'median_emp_lfo : ', model_targets%median_emp_lfo(io)
    enddo
    ! - Wealth and income shares by LFO
	write(uno, '(a)')'wealth and income shares by LFO'
	write(uno, '(a,f12.6)')'wealth_share_ep : ', model_targets%wealth_share_ep
    write(uno, '(a,f12.6)')'wealth_share_es : ', model_targets%wealth_share_es
    write(uno, '(a,f12.6)')'wealth_share_ec : ', model_targets%wealth_share_ec
    write(uno, '(a,f12.6)')'inc_share_ep : ', model_targets%inc_share_ep
    write(uno, '(a,f12.6)')'inc_share_es : ', model_targets%inc_share_es
    write(uno, '(a,f12.6)')'inc_share_ec : ', model_targets%inc_share_ec
    write(uno, '(a,f12.6)')'netinc_share_ep : ', model_targets%netinc_share_ep
    write(uno, '(a,f12.6)')'netinc_share_es : ', model_targets%netinc_share_es
    write(uno, '(a,f12.6)')'netinc_share_ec : ', model_targets%netinc_share_ec
	! - Share of income and wealth owned by entre
    write(uno, '(a)')'Share of income and wealth owned by entre'
	write(uno, '(a,f12.6)')'share_inc_entre : ', model_targets%share_inc_entre
	write(uno, '(a,f12.6)')'share_wealth_entre : ', model_targets%share_wealth_entre
	! - Inequality measures
    write(uno, '(a)')' '
    write(uno, '(a)')'Inequality measures'
	write(uno, '(a,f12.6)')'gini_wealth_all : ', model_targets%gini_wealth_all
	write(uno, '(a,f12.6)')'gini_wealth_entre : ', model_targets%gini_wealth_entre
	write(uno, '(a,f12.6)')'wealth_share(1) : ', model_targets%wealth_share(1)
	write(uno, '(a,f12.6)')'wealth_share(2) : ', model_targets%wealth_share(2)
	write(uno, '(a,f12.6)')'wealth_share(3) : ', model_targets%wealth_share(3)
	write(uno, '(a,f12.6)')'wealth_share(4) : ', model_targets%wealth_share(4)
	write(uno, '(a,f12.6)')'gini_inc_all : ', model_targets%gini_inc_all
	write(uno, '(a,f12.6)')'gini_inc_entre : ', model_targets%gini_inc_entre
	write(uno, '(a,f12.6)')'inc_share(1) : ', model_targets%inc_share(1)
	write(uno, '(a,f12.6)')'inc_share(2) : ', model_targets%inc_share(2)
	write(uno, '(a,f12.6)')'inc_share(3) : ', model_targets%inc_share(3)
	write(uno, '(a,f12.6)')'inc_share(4) : ', model_targets%inc_share(4)
	
	! - Occ/LFO choice by income and wealth
    write(uno, '(a)')' '
    write(uno, '(a)')'Occ/LFO choice by income and wealth'
	write(uno, '(a,f12.6)')'share_entre_top(1) : ', model_targets%share_entre_top(1)
	write(uno, '(a,f12.6)')'share_entre_top(2) : ', model_targets%share_entre_top(2)
	write(uno, '(a,f12.6)')'share_entre_top(3) : ', model_targets%share_entre_top(3)
	write(uno, '(a,f12.6)')'share_entre_top_inc(1) : ', model_targets%share_entre_top_inc(1)
	write(uno, '(a,f12.6)')'share_entre_top_inc(2) : ', model_targets%share_entre_top_inc(2)
	write(uno, '(a,f12.6)')'share_entre_top_inc(3) : ', model_targets%share_entre_top_inc(3)
	write(uno, '(a,f12.6)')'share_entre_quint(1) : ', model_targets%share_entre_quint(1)
	write(uno, '(a,f12.6)')'share_entre_quint(2) : ', model_targets%share_entre_quint(2)
	write(uno, '(a,f12.6)')'share_entre_quint(3) : ', model_targets%share_entre_quint(3)
	write(uno, '(a,f12.6)')'share_entre_quint(4) : ', model_targets%share_entre_quint(4)
	write(uno, '(a,f12.6)')'share_entre_quint(5) : ', model_targets%share_entre_quint(5)
	write(uno, '(a,f12.6)')'share_entre_quint_inc(1) : ', model_targets%share_entre_quint_inc(1)
	write(uno, '(a,f12.6)')'share_entre_quint_inc(2) : ', model_targets%share_entre_quint_inc(2)
	write(uno, '(a,f12.6)')'share_entre_quint_inc(3) : ', model_targets%share_entre_quint_inc(3)
	write(uno, '(a,f12.6)')'share_entre_quint_inc(4) : ', model_targets%share_entre_quint_inc(4)
	write(uno, '(a,f12.6)')'share_entre_quint_inc(5) : ', model_targets%share_entre_quint_inc(5)
	write(uno, '(a,f12.6)')'share_lfo_top(2,1) : ', model_targets%share_lfo_top(2,1)
	write(uno, '(a,f12.6)')'share_lfo_top(3,1) : ', model_targets%share_lfo_top(3,1)
	write(uno, '(a,f12.6)')'share_lfo_top(4,1) : ', model_targets%share_lfo_top(4,1)
	write(uno, '(a,f12.6)')'share_lfo_top(2,2) : ', model_targets%share_lfo_top(2,2)
	write(uno, '(a,f12.6)')'share_lfo_top(3,2) : ', model_targets%share_lfo_top(3,2)
	write(uno, '(a,f12.6)')'share_lfo_top(4,2) : ', model_targets%share_lfo_top(4,2)
	write(uno, '(a,f12.6)')'share_lfo_top(2,3) : ', model_targets%share_lfo_top(2,3)
	write(uno, '(a,f12.6)')'share_lfo_top(3,3) : ', model_targets%share_lfo_top(3,3)
	write(uno, '(a,f12.6)')'share_lfo_top(4,3) : ', model_targets%share_lfo_top(4,3)
	write(uno, '(a,f12.6)')'share_lfo_top_inc(2,1) : ', model_targets%share_lfo_top_inc(2,1)
	write(uno, '(a,f12.6)')'share_lfo_top_inc(3,1) : ', model_targets%share_lfo_top_inc(3,1)
	write(uno, '(a,f12.6)')'share_lfo_top_inc(4,1) : ', model_targets%share_lfo_top_inc(4,1)
	write(uno, '(a,f12.6)')'share_lfo_top_inc(2,2) : ', model_targets%share_lfo_top_inc(2,2)
	write(uno, '(a,f12.6)')'share_lfo_top_inc(3,2) : ', model_targets%share_lfo_top_inc(3,2)
	write(uno, '(a,f12.6)')'share_lfo_top_inc(4,2) : ', model_targets%share_lfo_top_inc(4,2)
	write(uno, '(a,f12.6)')'share_lfo_top_inc(2,3) : ', model_targets%share_lfo_top_inc(2,3)
	write(uno, '(a,f12.6)')'share_lfo_top_inc(3,3) : ', model_targets%share_lfo_top_inc(3,3)
	write(uno, '(a,f12.6)')'share_lfo_top_inc(4,3) : ', model_targets%share_lfo_top_inc(4,3)
	write(uno, '(a,f12.6)')'share_lfo_quint(2,1) : ', model_targets%share_lfo_quint(2,1)
	write(uno, '(a,f12.6)')'share_lfo_quint(3,1) : ', model_targets%share_lfo_quint(3,1)
	write(uno, '(a,f12.6)')'share_lfo_quint(4,1) : ', model_targets%share_lfo_quint(4,1)
	write(uno, '(a,f12.6)')'share_lfo_quint(2,2) : ', model_targets%share_lfo_quint(2,2)
	write(uno, '(a,f12.6)')'share_lfo_quint(3,2) : ', model_targets%share_lfo_quint(3,2)
	write(uno, '(a,f12.6)')'share_lfo_quint(4,2) : ', model_targets%share_lfo_quint(4,2)
	write(uno, '(a,f12.6)')'share_lfo_quint(2,3) : ', model_targets%share_lfo_quint(2,3)
	write(uno, '(a,f12.6)')'share_lfo_quint(3,3) : ', model_targets%share_lfo_quint(3,3)
	write(uno, '(a,f12.6)')'share_lfo_quint(4,3) : ', model_targets%share_lfo_quint(4,3)
	write(uno, '(a,f12.6)')'share_lfo_quint(2,4) : ', model_targets%share_lfo_quint(2,4)
	write(uno, '(a,f12.6)')'share_lfo_quint(3,4) : ', model_targets%share_lfo_quint(3,4)
	write(uno, '(a,f12.6)')'share_lfo_quint(4,4) : ', model_targets%share_lfo_quint(4,4)
	write(uno, '(a,f12.6)')'share_lfo_quint(2,5) : ', model_targets%share_lfo_quint(2,5)
	write(uno, '(a,f12.6)')'share_lfo_quint(3,5) : ', model_targets%share_lfo_quint(3,5)
	write(uno, '(a,f12.6)')'share_lfo_quint(4,5) : ', model_targets%share_lfo_quint(4,5)
	write(uno, '(a,f12.6)')'share_lfo_quint_inc(2,1) : ', model_targets%share_lfo_quint_inc(2,1)
	write(uno, '(a,f12.6)')'share_lfo_quint_inc(3,1) : ', model_targets%share_lfo_quint_inc(3,1)
	write(uno, '(a,f12.6)')'share_lfo_quint_inc(4,1) : ', model_targets%share_lfo_quint_inc(4,1)
	write(uno, '(a,f12.6)')'share_lfo_quint_inc(2,2) : ', model_targets%share_lfo_quint_inc(2,2)
	write(uno, '(a,f12.6)')'share_lfo_quint_inc(3,2) : ', model_targets%share_lfo_quint_inc(3,2)
	write(uno, '(a,f12.6)')'share_lfo_quint_inc(4,2) : ', model_targets%share_lfo_quint_inc(4,2)
	write(uno, '(a,f12.6)')'share_lfo_quint_inc(2,3) : ', model_targets%share_lfo_quint_inc(2,3)
	write(uno, '(a,f12.6)')'share_lfo_quint_inc(3,3) : ', model_targets%share_lfo_quint_inc(3,3)
	write(uno, '(a,f12.6)')'share_lfo_quint_inc(4,3) : ', model_targets%share_lfo_quint_inc(4,3)
	write(uno, '(a,f12.6)')'share_lfo_quint_inc(2,4) : ', model_targets%share_lfo_quint_inc(2,4)
	write(uno, '(a,f12.6)')'share_lfo_quint_inc(3,4) : ', model_targets%share_lfo_quint_inc(3,4)
	write(uno, '(a,f12.6)')'share_lfo_quint_inc(4,4) : ', model_targets%share_lfo_quint_inc(4,4)
	write(uno, '(a,f12.6)')'share_lfo_quint_inc(2,5) : ', model_targets%share_lfo_quint_inc(2,5)
	write(uno, '(a,f12.6)')'share_lfo_quint_inc(3,5) : ', model_targets%share_lfo_quint_inc(3,5)
	write(uno, '(a,f12.6)')'share_lfo_quint_inc(4,5) : ', model_targets%share_lfo_quint_inc(4,5)
	! - Employment share by firm size
    write(uno, '(a)')' '
    write(uno, '(a)')'Employment share by firm size'
	write(uno, '(a,f12.6)')'share_emp(1) : ', model_targets%share_emp(1)
	write(uno, '(a,f12.6)')'share_emp(2) : ', model_targets%share_emp(2)
	write(uno, '(a,f12.6)')'share_emp(3) : ', model_targets%share_emp(3)
	write(uno, '(a,f12.6)')'share_emp(4) : ', model_targets%share_emp(4)
	! - Transitions
    write(uno, '(a)')' '
    write(uno, '(a)')'Transitions'
	write(uno, '(a,f12.6)')'trans_entre_work : ', model_targets%trans_entre_work
	write(uno, '(a,f12.6)')'trans_work_entre : ', model_targets%trans_work_entre
	write(uno, '(a,f12.6)')'trans_CP : ', model_targets%trans_CP
	write(uno, '(a,f12.6)')'trans_PC : ', model_targets%trans_PC
	write(uno, '(a,f12.6)')'trans_SC : ', model_targets%trans_SC

	write(uno, '(a,f12.6)')'trans_mat(1,1) : ', model_targets%trans_mat(1,1) 
	write(uno, '(a,f12.6)')'trans_mat(1,2) : ', model_targets%trans_mat(1,2)
	write(uno, '(a,f12.6)')'trans_mat(1,3) : ', model_targets%trans_mat(1,3)
	write(uno, '(a,f12.6)')'trans_mat(1,4) : ', model_targets%trans_mat(1,4)
	write(uno, '(a,f12.6)')'trans_mat(2,1) : ', model_targets%trans_mat(2,1) 
	write(uno, '(a,f12.6)')'trans_mat(2,2) : ', model_targets%trans_mat(2,2)
	write(uno, '(a,f12.6)')'trans_mat(2,3) : ', model_targets%trans_mat(2,3)
	write(uno, '(a,f12.6)')'trans_mat(2,4) : ', model_targets%trans_mat(2,4)
	write(uno, '(a,f12.6)')'trans_mat(3,1) : ', model_targets%trans_mat(3,1) 
	write(uno, '(a,f12.6)')'trans_mat(3,2) : ', model_targets%trans_mat(3,2)
	write(uno, '(a,f12.6)')'trans_mat(3,3) : ', model_targets%trans_mat(3,3)
	write(uno, '(a,f12.6)')'trans_mat(3,4) : ', model_targets%trans_mat(3,4)
	write(uno, '(a,f12.6)')'trans_mat(4,1) : ', model_targets%trans_mat(4,1) 
	write(uno, '(a,f12.6)')'trans_mat(4,2) : ', model_targets%trans_mat(4,2)
	write(uno, '(a,f12.6)')'trans_mat(4,3) : ', model_targets%trans_mat(4,3)
	write(uno, '(a,f12.6)')'trans_mat(4,4) : ', model_targets%trans_mat(4,4)

	! Average tax rate by income percentiles
	write(uno, '(a)')' '
    write(uno, '(a)')'ATR'
    write(uno, '(a,f12.6)')'P99.99+ : ',      model_targets%atr_vec(1)
    write(uno, '(a,f12.6)')'P99.9-P99.99 : ', model_targets%atr_vec(2)
	write(uno, '(a,f12.6)')'P99-P99.9 : ',    model_targets%atr_vec(3)
	write(uno, '(a,f12.6)')'P90-P99 : ',      model_targets%atr_vec(4)
	write(uno, '(a,f12.6)')'P50-P90 : ',      model_targets%atr_vec(5)
	write(uno, '(a,f12.6)')'Bottom 50% : ',   model_targets%atr_vec(6)

    write(uno, '(a)')' '
    write(uno, '(a)')'Checks wealth distribution'
    write(uno, '(a,f12.6)')'Q1,mass : ', mass_q(1)
    write(uno, '(a,f12.6)')'Q2,mass : ', mass_q(2)
    write(uno, '(a,f12.6)')'Q3,mass : ', mass_q(3)
    write(uno, '(a,f12.6)')'Q4,mass : ', mass_q(4)
    write(uno, '(a,f12.6)')'Q5,mass : ', mass_q(5)
    
    write(uno, '(a)')' '
    write(uno, '(a)')'Checks income distribution'
    write(uno, '(a,f12.6)')'Q1,mass : ', inc_mass_q(1)
    write(uno, '(a,f12.6)')'Q2,mass : ', inc_mass_q(2)
    write(uno, '(a,f12.6)')'Q3,mass : ', inc_mass_q(3)
    write(uno, '(a,f12.6)')'Q4,mass : ', inc_mass_q(4)
    write(uno, '(a,f12.6)')'Q5,mass : ', inc_mass_q(5)
    
    write(uno, '(a)')' '
    write(uno, '(a)')'Checks max income'
    write(uno, '(a,f12.6)')'maxinc_work : ', maxinc_work
    write(uno, '(a,f12.6)')'maxinc_ep   : ', maxinc_ep
    write(uno, '(a,f12.6)')'maxinc_es   : ', maxinc_es
    write(uno, '(a,f12.6)')'maxinc_ec   : ', maxinc_ec
    write(uno, '(a,f12.6)')'mass_a1     : ', mass_a1
    write(uno, '(a,f12.6)')'mass_a1_a3  : ', mass_a1_a3
    write(uno, '(a,f12.6)')'mass_aN     : ', mass_aN
    close(uno)
    
    open(newunit=uno, file=trim(savedir)//'inequality.txt', status='replace')
    write(uno, '(f12.6)') model_targets%gini_wealth_all, model_targets%gini_wealth_entre, model_targets%wealth_share(1), &
	model_targets%wealth_share(2), model_targets%wealth_share(3), model_targets%wealth_share(4), model_targets%gini_inc_all, &
    model_targets%gini_inc_entre, model_targets%inc_share(1), model_targets%inc_share(2), model_targets%inc_share(3), model_targets%inc_share(4)
    close(uno)
    
    open(newunit=uno, file=trim(savedir)//'occupation_quint.txt', status='replace')
    write(uno, '(f12.6)') model_targets%share_entre_top(1), model_targets%share_entre_top(2), model_targets%share_entre_top(3), &
    model_targets%share_entre_top_inc(1), model_targets%share_entre_top_inc(2), model_targets%share_entre_top_inc(3), &
    model_targets%share_entre_quint(1), model_targets%share_entre_quint(2), model_targets%share_entre_quint(3), &
    model_targets%share_entre_quint(4), model_targets%share_entre_quint(5), model_targets%share_entre_quint_inc(1), &
	model_targets%share_entre_quint_inc(2), model_targets%share_entre_quint_inc(3), model_targets%share_entre_quint_inc(4), model_targets%share_entre_quint_inc(5), &
    model_targets%share_lfo_top(2,1), model_targets%share_lfo_top(3,1), model_targets%share_lfo_top(4,1), model_targets%share_lfo_top(2,2), &
    model_targets%share_lfo_top(3,2), model_targets%share_lfo_top(4,2), model_targets%share_lfo_top(2,3), model_targets%share_lfo_top(3,3), &
	model_targets%share_lfo_top(4,3), model_targets%share_lfo_top_inc(2,1), model_targets%share_lfo_top_inc(3,1), model_targets%share_lfo_top_inc(4,1), &
	model_targets%share_lfo_top_inc(2,2), model_targets%share_lfo_top_inc(3,2), model_targets%share_lfo_top_inc(4,2), model_targets%share_lfo_top_inc(2,3), &
    model_targets%share_lfo_top_inc(3,3), model_targets%share_lfo_top_inc(4,3), model_targets%share_lfo_quint(2,1), model_targets%share_lfo_quint(3,1), &
    model_targets%share_lfo_quint(4,1), model_targets%share_lfo_quint(2,2), model_targets%share_lfo_quint(3,2), model_targets%share_lfo_quint(4,2), &
	model_targets%share_lfo_quint(2,3), model_targets%share_lfo_quint(3,3), model_targets%share_lfo_quint(4,3), model_targets%share_lfo_quint(2,4), &
	model_targets%share_lfo_quint(3,4), model_targets%share_lfo_quint(4,4), model_targets%share_lfo_quint(2,5), model_targets%share_lfo_quint(3,5), &
	model_targets%share_lfo_quint(4,5), model_targets%share_lfo_quint_inc(2,1), model_targets%share_lfo_quint_inc(3,1), model_targets%share_lfo_quint_inc(4,1), &
    model_targets%share_lfo_quint_inc(2,2), model_targets%share_lfo_quint_inc(3,2), model_targets%share_lfo_quint_inc(4,2), model_targets%share_lfo_quint_inc(2,3), &
	model_targets%share_lfo_quint_inc(3,3), model_targets%share_lfo_quint_inc(4,3), model_targets%share_lfo_quint_inc(2,4), model_targets%share_lfo_quint_inc(3,4), &
	model_targets%share_lfo_quint_inc(4,4), model_targets%share_lfo_quint_inc(2,5), model_targets%share_lfo_quint_inc(3,5), model_targets%share_lfo_quint_inc(4,5)
    close(uno)
    
    end subroutine txt_export
    !===============================================================================!
    
    subroutine mom_export(dirname,model_targets)
    ! Purpose: WRITES certain moments to external txt files in folder filenam

    !Declare inputs:
    character(LEN=*), intent(in)    :: dirname
    type(modelResults), intent(in) :: model_targets
    
    call writedim(model_targets%share_entre_quint_inc,trim(dirname)//trim('share_entre_quint_inc.txt'))
    call writedim(model_targets%share_entre_quint,trim(dirname)//trim('share_entre_quint.txt'))
    call writedim(model_targets%share_lfo_quint_inc,trim(dirname)//trim('share_lfo_quint_inc.txt'))
    call writedim(model_targets%share_lfo_quint,trim(dirname)//trim('share_lfo_quint.txt'))
    
    call writescalar(model_targets%share_entre_act,trim(dirname)//trim('share_entre_act.txt'))
    call writescalar(model_targets%share_EP_entre,trim(dirname)//trim('share_EP_entre.txt'))
    call writescalar(model_targets%share_ES_entre,trim(dirname)//trim('share_ES_entre.txt'))
    call writescalar(model_targets%share_EC_entre,trim(dirname)//trim('share_EC_entre.txt'))

    
    end subroutine mom_export
     !===============================================================================!
end module mod_targets
    
    