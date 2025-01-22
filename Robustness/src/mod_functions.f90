module mod_functions

    ! Module mod_functions contains procedures specific to this project
    ! Numerical routines are contained in mod_numerical
    
    use mod_utilities, only: readscalar
    use mod_globals  ! Make model parameters known to these procedures
    implicit none
   
contains
    
    !=======================================================================!
    subroutine fun_prices(p)
    implicit none
    type(modelPrices), intent(inout) :: p
    !Purpose:
    !Given the intitial guess for r, computes wage w and capital-labor ratio
    !in the corporate sector using the Cobb-Douglas FOCS.
    
    p%KN = (alpha/(p%r+delta))**(1.0d0/(1.0d0-alpha))
    p%w = (1.0d0-alpha)*(p%KN**alpha)
    
    ! Set y_ave to be the average labor income. 0.33 is the target for average hours worked
    p%y_ave  = p%w*0.33d0
    
    ! Compute pension benefit and social security cap as functions of avg labor income
    p%pen    = repl*p%y_ave 
    p%ss_cap = ss_cap_ratio*p%y_ave
	
	! if not benchmark - ss_cap is fixed at the benchmark level
    if ((.not. do_benchmark) .and. do_fix_ss_cap == 1) then
         call readscalar(p%ss_cap,trim(dir_bench)//trim("ss_cap.txt"))
    endif
    ! if not benchmark - pen is fixed at the benchmark level
    if ((.not. do_benchmark) .and. do_fix_pen == 1) then
        call readscalar(p%pen,trim(dir_bench)//trim("pen.txt"))
    endif
    end subroutine fun_prices
    !=======================================================================!
    subroutine fun_prices_KN(p, KN)
    implicit none
    real(8), intent(in) :: KN
    type(modelPrices), intent(out) :: p
    !Purpose:
    !Given the intitial guess for KN, computes w and r
    !in the corporate sector using the Cobb-Douglas FOCS.
    
    if (KN <= 0.0d0) then
    	write(*,*) "fun_prices_KN: input KN is negative! Setting r to 1/beta-1."
    	p%r = 1.0d0/beta - 1.0d0
    	p%w = (1.0d0-alpha)*(alpha/(p%r + delta))**(alpha/(1.0d0-alpha))
    else
    	p%r = alpha*KN**(alpha-1.0d0) - delta
    	p%w = (1.0d0-alpha)*(KN**alpha)
    endif
    
	! Set y_ave to be the average labor income. 0.33 is the target for average hours worked
	p%y_ave  = p%w*0.33d0
	p%pen    = repl*p%y_ave 
	p%ss_cap = ss_cap_ratio*p%y_ave
	
	 if ((.not. do_benchmark) .and. do_fix_ss_cap == 1) then
         call readscalar(p%ss_cap,trim(dir_bench)//trim("ss_cap.txt"))
     endif
     ! if not benchmark - pen is fixed at the benchmark level
    if ((.not. do_benchmark) .and. do_fix_pen == 1) then
        call readscalar(p%pen,trim(dir_bench)//trim("pen.txt"))
    endif
    
    end subroutine fun_prices_KN
    !=======================================================================!
    
    pure function fun_update_hsv_0(calY1,calY2,G_bench,taxes_corp,taxes_div) result(hsv_0_new)  
    ! This function computes the updated value for hsv_0 (lambda_i) to clear the
    ! government budget constraint  
    implicit none
    ! Inputs/outputs
    real(8),intent(in) :: calY1,calY2,G_bench,taxes_corp,taxes_div
    real(8) :: hsv_0_new
    
    hsv_0_new = (calY1 - (G_bench - taxes_corp - taxes_div))/calY2
    
    end function fun_update_hsv_0
    
!=======================================================================!
    
    pure function cobb_douglas(K,N) result(F)
        implicit none
        real(8), intent(in) :: K,N
        real(8) :: F
        !Purpose:
        !Cobb-Douglas prod fun used by corporate firm
    
        F = (K**alpha)*(N**(1.0d0-alpha))
    
    end function cobb_douglas
    !=======================================================================!
    
    pure function util(c,l) result(F)
        implicit none
    
        !Declare inputs:
        real(8), intent(in) :: c, l ! consumption and labor
        !Declare output:
        real(8) :: F 
        
        if (c<=0d0) then
            F = large_negative
            return
        endif
        
        if (sigma1==2d0) then
            F = -1.0d0/c - chi*l**(1.0d0+sigma2)/(1.0d0+sigma2)
        elseif (sigma1==1d0) then
            F = log(c) - chi*l**(1.0d0+sigma2)/(1.0d0+sigma2)
        else
            F = c**(1.0d0-sigma1)/(1.0d0-sigma1) - chi*l**(1.0d0+sigma2)/(1.0d0+sigma2)
        endif
        
    end function util
    !=======================================================================!
    
    pure function util_cons(c) result(F)
        implicit none
        ! PURPOSE:
        ! Utility from consumption. Compare to util above.
        !Declare inputs:
        real(8), intent(in) :: c ! consumption 
        !Declare output:
        real(8) :: F 
    
        if (c<=0d0) then
            F = large_negative
            return
        endif
    
        if (sigma1==2d0) then
            F = -1.0d0/c 
        elseif (sigma1==1d0) then
            F = log(c) 
        else
            F = c**(1.0d0-sigma1)/(1.0d0-sigma1) 
        endif
        
    end function util_cons
    !=======================================================================!
    
    pure function prodfun(k,nbar) result(F)
        implicit none
        ! F: production function f(k,nbar) (without ability "theta")
        ! where firm rents capital k and total labor nbar,
        ! Declare inputs:
        real(8), intent(in) :: k  ! capital 
        real(8), intent(in) :: nbar  ! hired labor 
        !Declare output:
        real(8) :: F

        F = (k**gamma*(nbar)**(1.0d0-gamma))**vi

    end function prodfun

    !=======================================================================!
    
    subroutine profit_max(a,theta,lambda_,r,w,k_opt,n_opt)
    ! Compute optimal capital and labor inputs
    ! Declare inputs:
    real(8), intent(in) :: a       ! Asset holdings
    real(8), intent(in) :: theta   ! Entre ability
    real(8), intent(in) :: lambda_ ! Collateral constraint
    real(8), intent(in) :: r, w    ! Prices
    ! Declare outputs:
    real(8), intent(out) :: k_opt, n_opt ! capital, labor
    ! Declare locals:
    real(8) :: k_unc
    
    ! Body of subroutine
    
    ! Compute unconstrained capital demand
    k_unc = (theta*vi)**(1d0/(1d0-vi))*(gamma/(r+delta))**(1d0+gamma*vi/(1d0-vi))*((1d0-gamma)/w)**(vi*(1d0-gamma)/(1d0-vi))
    
    ! Compute optimal capital demand
    k_opt = max(min(lambda_*a,k_unc),0d0)
    
    ! Compute optimal labor as a function of capital
    n_opt = (theta*(1.0d0-gamma)*vi/w)**(1.0d0/(1.0d0-(1.0d0-gamma)*vi))*k_opt**((gamma*vi)/(1.0d0-(1.0d0-gamma)*vi))
    
    end subroutine profit_max
    !=======================================================================!
    
    pure function fun_profit(theta_val,k,le,n,rbar,wage) result(profit)
	implicit none
	REAL(8), INTENT(in) :: theta_val,k,le,n,rbar,wage
	REAL(8) :: profit

	profit = theta_val*(k**gamma*(le+n)**(1.0d0-gamma))**vi-delta*k-rbar*k-wage*n

    end function fun_profit
    !=======================================================================!
    
    pure function cost_evasion_es(x) result(F)
        implicit none

        !Declare inputs:
        real(8), intent(in) :: x ! here x=1-\phi!!!
        !Declare outputs:
        real(8) :: F
        
        if (x>0d0  .and. cf_avoidance ==0 ) then
        ! Input share is valid and baseline model environment
            F = c0_es*x**c1_es !+10000.0d0
        else
        ! Input is invalid or counterfactual no tax avoidance environment
            F = 0d0
        endif
        
       
    end function cost_evasion_es
    !=======================================================================!
    
    pure function cost_evasion_ec(x) result(F)
       implicit none
    
       !Declare inputs:
       real(8), intent(in) :: x ! here x=\phi!!!
       !Declare outputs:
       real(8) :: F
       
       if (x>0d0 .and. cf_avoidance ==0) then 
       ! Input share is valid and baseline model environment
           F = c0_ec*x**c1_ec !+10000.0d0
       else
       ! Input is invalid or counterfactual no tax avoidance environment
           F = 0d0
       endif
       
      
    end function cost_evasion_ec
    !=======================================================================!
    
    pure function fun_tax_inc(income,y_H) result(F)
        implicit none
        ! Personal income tax
        ! Non-linear HSV tax function for progressivity
        ! plus flat marginal tax rate at the top
        
        !Declare inputs:
        real(8), intent(in) :: income ! income in levels
        real(8), intent(in) :: y_H    ! income cutoff 
        !Declare outputs:
        real(8) :: F ! total taxes paid 
        ! Declare locals:
        real(8) :: y
        
        ! Take care of negative income
        y = max(income,0d0)
        if (y<y_H) then
        	F = y - hsv_0*y**(1d0-hsv_1)
        else
        	F = y_h- hsv_0*y_h**(1d0-hsv_1) + tau_h*(y-y_h)
        endif
    
    end function fun_tax_inc 
 !=======================================================================!
    
    function fun_y_H(lambda,tau,tau_H) result(y_H)  
    ! This function returns the income threshold in the income tax function:
    ! if income<=y_H, income tax is standard HSV, if income>y_H, top tax rate applies   
    implicit none
    ! Inputs/outputs
    real(8),intent(in) :: lambda  ! Scale parameter
    real(8),intent(in) :: tau     ! Progressivity parameter
    real(8),intent(in) :: tau_H   ! Top tax rate
    ! Function result
    real(8) :: y_H
    
    y_H = (lambda*(1.0d0-tau)/(1.0d0-tau_H))**(1.0d0/tau)
    
    end function fun_y_H
    
!=======================================================================!
    
    function fun_l_H(a_val,eps_val,w,r,y_H,ss_cap) result(l_H)  
    ! Threshold for income tax function in terms of labor hours, defined as
    ! l_H such that taxable income is equal to y_H   
    implicit none
    ! Inputs/outputs
    real(8),intent(in) :: a_val,eps_val,w,r,y_H,ss_cap
    real(8) :: l_H
    ! Local variables
    real(8) :: l_ss_cap
    
    l_ss_cap = ss_cap/w*eps_val
    
    l_H = (y_H-r*a_val)/((1.0d0-tau_p)*w*eps_val)
    
    if (l_H>l_ss_cap) then
    	l_H = (y_H-r*a_val+tau_p*ss_cap)/(w*eps_val)
    endif	
    
    end function fun_l_H
!=======================================================================!   
    
    pure function fun_tax_div(div_inc) result(tax)
        implicit none
        ! Dividend income tax
        ! For an explanation and details, see memo XX (Georgi!)
        ! NOTE: We do not consider the standard deduction

        !Declare inputs:
        real(8), intent(in) :: div_inc ! income in levels 
        !Declare outputs:
        real(8) :: tax ! total taxes paid 
        !Declare locals:
        real(8) :: div_inc1
        
        ! Make sure we don't have negative income
        div_inc1 = max(div_inc,0d0)
        
        tax = tau_d*div_inc1
    
    end function fun_tax_div 
!=======================================================================!
    
    pure function fun_tax_corp(income) result(F)
        implicit none
        ! Corporate income tax
        
        !Declare inputs:
        real(8), intent(in) :: income ! income in levels 
        !Declare outputs:
        real(8) :: F ! total taxes paid 
        ! Declare locals:
        real(8) :: y 
        
        y = max(income,0d0)
        
        F = tau_c*y 
        
        
    end function fun_tax_corp
!=======================================================================!    
    pure function fun_tax_pay(income,ss_cap) result(F)
        implicit none
        ! Payroll or Social Security tax
        
        !Declare inputs:
        real(8), intent(in) :: income  ! income in levels
        real(8), intent(in) :: ss_cap  ! Social security cap in levels
        ! Globals: 
        ! tau_p
        !Declare outputs:
        real(8) :: F ! total taxes paid 
        ! Declare locals:
        real(8) :: y
        
        y = max(income,0d0)
        
        F =  tau_p*min( y, ss_cap ) 
        
    end function fun_tax_pay

!=======================================================================!
    
end module mod_functions
