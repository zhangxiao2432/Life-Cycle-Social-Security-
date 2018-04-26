!===================================================================================
!Project 1
!Zhonghui Zhang
!===================================================================================
!parameter values are chosen.
!-----------------------------------------------------------------------------------
module calibration
implicit none

integer, parameter :: jmax=65,jr=45             ! max lifespan, retirement age
integer, parameter :: dm=601                    !# of asset value

real, parameter :: beta = 1.01                  ! discount factor
real, parameter :: alpha =0.640                 ! labor share
real, parameter :: gamma = 2.000                ! intertemporal elastivity of substitution
real, parameter :: h = 0.45                     ! verage hours of work
real, parameter :: delta = 0.080                ! depreciation rate
real, parameter :: B = 1.319
real, parameter :: pi = 0.94                    ! prob of employment
real, parameter :: rho = 0.012                ! growth rate of population
real, parameter :: theta = 0.35                  ! employment replacement rate
real, parameter :: xi = 0.40                    ! unemployment replacement rate 

end module 

!--------------------------------------------------------------------------------
!Declare variables;
!--------------------------------------------------------------------------------
module global
use calibration
implicit none

integer :: i,ii,j,jj,es
real    :: w,r,benf                                  !wage/capital return/social security benefit
integer :: yloc(1:jmax,1:dm,1:2)                     !location of the optimal solution

real, dimension(1:jmax+1,1:dm)     :: d              !asset values
real, dimension(1:jmax,1:dm,1:2)   :: c_grid         !consumption of an age-j individual
real, dimension(1:jmax)            :: c              !optimal consumption for different age of agents
real, dimension(1:jmax)            :: Au             !average utility
real, dimension(1:jmax)            :: y              !optimal asset for different age of agents
real, dimension(1:jmax)            :: inc            !disposable income for different age of agents
real, dimension(1:jmax,1:2)        :: q              !disposable income of an age-j agent
real, dimension(1:jmax+1)          :: phi            !survival probability
real, dimension(1:jmax,1:dm,1:2)   :: lambda         !distribution of agents
real, dimension(1:jmax,1:dm,1:2)   :: ysol           !aggregate decision rule
real, dimension(1:jmax,1:dm,1:2)   :: csol           !optimal consumption rule
real, dimension(1:jmax+1,1:dm,1:2) :: v              !value function
real, dimension(1:jmax)            :: epsilon,epn    !efficiency profile/efficiency profile index
real, dimension(1:jmax)            :: mu,muhat       !siza/share of age group
real, dimension(2,2)               :: trans          !trans matrix
real, dimension(1:jmax,1:dm,1:2)   :: mea            !measure of age-j, asset-i group
real, dimension(1:jmax,1:dm,1:2)   :: T              !lum sum accidental bequest
real, dimension(1:jmax,1:dm,1:2)   :: usol           !utility of each group

end module 

!---------------------------------------------------------------------------------------------------------------------------
!Survival Probabilities and Efficiency Profile
!---------------------------------------------------------------------------------------------------------------------------

subroutine demographics
use calibration
use global
 
implicit none

!efficiency profile

open(unit=11,file='epsilon.txt') 
epn=0.0
do j=1,jr-1
  read(11,*) epsilon(j)
end do
close(11)  

!normalize epsilon
do j=1,jr-1
  epn(j)=epsilon(j)/sum(epsilon(1:jr-1))*(jr-1)
end do

do j = jr, jmax
  epn(j) = 0.0
end do                    

!do j = 1, jmax
!  epn(j) = 1.0
!end do    

! survival probability

open(unit=13,file='phi.txt') 
phi = 0.0
do j=1,jmax
  read(13,*) phi(j)
end do
close(13)

!do j = 1, jmax
!  phi(j) = 1.0
!end do 

muhat(1) = 1.0
do j=1,jmax-1
  muhat(j+1)=(phi(j+1)*muhat(j))/(1+rho)
end do
   
do j=1,jmax
  mu(j) = muhat(j)/sum(muhat)
end do

end subroutine 

!---------------------------------------------------------------------------------------------------------------------------
!Asset Holding
!---------------------------------------------------------------------------------------------------------------------------
subroutine Asset
use global
implicit none

real,parameter :: dstep = 0.025

d(jmax+1,:)=0
do j=1,jmax
  do i=1,dm
      d(j,i)=(i-1)*dstep
  end do
end do

end subroutine 

!--------------------------------------------------------------------------------------------------------------------------
!Value Function Iteration
!--------------------------------------------------------------------------------------------------------------------------
subroutine value_function(tau_u0, tau_s0, T_0, r_0)

use calibration
use global

implicit none
 
 
real, intent(in)                 ::  tau_u0, tau_s0, T_0, r_0
 
real, dimension(1:jmax,1:dm,1:2) :: vtemp                                ! Value function for all individual with all possible d



call Asset
call Demographics
w=alpha*B*((r_0+delta)/((1-alpha)*B))**((alpha-1)/alpha) 
benf= w*h*theta*sum(epn(1:jr-1))/(jr-1)  
    
!Disposable Income  
do j= 1,jr-1
  q(j,2) = (1-tau_u0-tau_s0)*w*h*epn(j)                                ! Disposable income for employed
  q(j,1) = xi*w*h                                                          ! Disposable income for unemployed
end do

do j=jr,jmax
  q(j,:)=benf
end do
          
!Value function iteration for age 1 to jmax
v(jmax+1,:,:)=0                                               ! Initialize value function
do j=jmax,1,-1                                               ! loop 1 begins (starts from the last period for retired people)
   do i=1,dm                                                  ! loop 2 begins (assets values for d(1) to d(601)) state variable
     do es=1,2                                                ! loop 3 begins (employment for 1 to 2) state variable
       do ii=1,dm                                             ! loop 4 begins (assets values for d(1) to d(601)) decision variable
         c_grid(j,ii,es)=(1+r_0)*d(j,i)-d(j+1,ii)+q(j,es)+T_0
           if (c_grid(j,ii,es)<=0) then
         vtemp(j,ii,es) = -1.0e16
           else
         vtemp(j,ii,es)=(c_grid(j,ii,es)**(1-gamma)-1)/(1-gamma)+beta*phi(j)*(pi*v(j+1,ii,2)+(1-pi)*v(j+1,ii,1)) 
           end if  
      end do
          v(j,i,es) = maxval(vtemp(j,:,es))
          yloc(j,i,es) = maxloc(vtemp(j,:,es),dim=1)
          ysol(j,i,es) = d(j+1,yloc(j,i,es))
          csol(j,i,es) = c_grid(j,yloc(j,i,es),es) 
    end do
  end do         
end do

end subroutine value_function

!--------------------------------------------------------------------------------------------------------------------------
!Main Program
!--------------------------------------------------------------------------------------------------------------------------
program main

use calibration
use global

implicit none


integer,parameter        :: itermax = 1000
real,parameter           :: adj = 0.4
real,parameter           :: tol = 0.005          
   
real                                 :: r_0,T_0,T1,N1,K1,r1,tau_u0,tau_s0,diffT,diffr,diffs,diffu,tau_u1,tau_s1,omiga,avinc,avcon
integer                              :: iter,ess
real, dimension(1:jmax)              :: cumphi

call asset
call demographics

!Initial Guess
   r_0 = 0.04                                
   T_0 = 0.03                                   
tau_u0 = 0.05
tau_s0 = 0.07                                  
  iter = 0

do
  call value_function (tau_u0, tau_s0, T_0, r_0)
 
!-----------------------------------------------Computing Age-Dependent Distribution---------------------------------------!
 
lambda=0.0                                                       ! Initialize the value of lambda
lambda(1,1,2)=0.94                                               ! Distribution of age 1 cohor with 0 asset and employed 
lambda(1,1,1)=0.06                                               ! Distribution of age 1 cohor with 0 asset and unemployed 
trans(1,:)=(/1-pi,pi/)
trans(2,:)=(/1-pi,pi/)                                          
  
do j=1,jr-2
  do i=1,dm                                                      ! Employment status today
    do ess=1,2                                                
      do es=1,2                                                  ! Employment status tomorrow
           jj=yloc(j,i,es)
           lambda(j+1,jj,ess)=lambda(j+1,jj,ess)+trans(es,ess)*lambda(j,i,es)
      end do
    end do
  end do
end do

do es=1,2                                                       ! Employment status today
    do i=1,dm                                                                                                 
         jj=yloc(j,i,1)
         lambda(jr,jj,1)=lambda(jr,jj,1)+lambda(jr-1,i,es)
    end do
end do

do j=jr,jmax-1
  do i=1,dm                                                   
         jj=yloc(j,i,1)
         lambda(j+1,jj,1)=lambda(j+1,jj,1)+lambda(j,i,1)
  end do
end do    
!-----------------------------------------------Aggregation-------------------------------------------!
do j=1,jmax
   do i=1,dm
      do es=1,2
         mea(j,i,es)=mu(j)*lambda(j,i,es)
      end do
   end do
end do
 
do j=1,jmax
   do i=1,dm
      do es=1,2
        T(j,i,es)=(1-phi(j))*mea(j,i,es)*ysol(j,i,es)  
      end do
   end do 
end do 

N1=0                                                
do j=1,jr-1                                   
      N1=N1+(mu(j)*sum(lambda(j,:,2))*epn(j))*h       
end do

K1=sum(mea*ysol)
T1=sum(T)
       
tau_s1 = sum(mu(jr:jmax))*benf/(sum(mea(1:jr-1,:,:)*epn(1:jr-1))*w*h)
tau_u1 = xi*sum(mu(1:jr-1))/(sum(mu(1:jr-1)*epn(1:jr-1)))
    r1 = (1-alpha)*B*(K1/N1)**(-alpha)-delta

 diffr = abs(r1-r_0)
 diffT = abs(T1-T_0)
 diffs = abs(tau_s1-tau_s0)
 diffu = abs(tau_u1-tau_u0)
 
if (diffr>tol .or. diffT>tol .or. diffs>tol .or. diffu>tol) then
      
    if (iter<itermax) then
        iter = iter+1
         r_0 = adj*r1+(1-adj)*r_0
         T_0 = adj*T1+(1-adj)*T_0
      tau_s0 = adj*tau_s1+(1-adj)*tau_s0
      tau_u0 = adj*tau_u1+(1-adj)*tau_u0

        write (*,*) '======================================'
        write (*,*) 'iteration', iter
        write (*,*) 'r error' ,diffr
        write (*,*) 't error' ,diffT
        write (*,*) 'tau_s error' ,diffs
        write (*,*) 'tau_u error' ,diffu
        write (*,*) '======================================'
    else 
        write (*,*) '======================================'
     	write (*,*) 'No convergence found'
        write (*,*) 'r error' ,diffr
        write (*,*) 't error' ,diffT
        write (*,*) 'tau_s error' ,diffs
        write (*,*) 'tau_u error' ,diffu
        write (*,*) '======================================'
	exit
    
    end if

else 
        goto 123
      end if
      end do  
     
123 continue
	write (*,*) '===============Convergence found================'
    
cumphi(1)=1
do j=1,jmax-1
  cumphi(j+1)=cumphi(j)*phi(j+1)
end do
  
do j=1,jr-1
  do i=1,dm
    do es=1,2   
       usol(j,i,es)=(csol(j,i,es)**(1-gamma)-1)/(1-gamma)
    end do
  end do
       c(j) = sum(csol(j,:,:)*lambda(j,:,:))                                 ! consumption profile
       Au(j) = beta**(j-1)*cumphi(j)*sum(usol(j,:,:)*lambda(j,:,:))
end do

do j=jr,jmax
  do i=1,dm   
       usol(j,i,1)=(csol(j,i,1)**(1-gamma)-1)/(1-gamma)
  end do
       c(j) = sum(csol(j,:,1)*lambda(j,:,1))                                 ! consumption profile
       Au(j) = beta**(j-1)*cumphi(j)*sum(usol(j,:,1)*lambda(j,:,1))
end do

omiga=sum(Au(:))
avcon=sum(c)/jmax
 
inc=0
  do j=1,jr-1
     do i=1,dm
       do es=1,2       
		     inc(j)=inc(j)+q(j,es)*lambda(j,i,es)
       end do
     end do  
  end do
  
do j=jr,jmax
  do i=1,dm      
	  inc(j)=inc(j)+q(j,1)*lambda(j,i,1)
  end do 
end do
  
avinc=(B*K1**(1-alpha)*N1**alpha-delta*K1)/sum(mu(1:44)*0.94)

!==================Print result===========================================
    
    write (*,*) '     replacement rate',theta
    write (*,*) 
    write (*,*) '             Tax rate',tau_s0
    write (*,*) 
    write (*,*) '            wage rate',w
    write (*,*) 
    write (*,*) '    return to capital',r_0
    write (*,*)
    write (*,*) '  average consumption',avcon
    write (*,*)
    write (*,*) '        capital stock',K1
    write (*,*)
    write (*,*) '       average income',avinc
    write (*,*)
    write (*,*) '      average utility',omiga
    write (*,*)  
    write (*,*) '               output',B*K1**(1-alpha)*N1**alpha
    write (*,*) 
    write (*,*) ' capital output ratio',K1/(B*K1**(1-alpha)*N1**alpha)
    write (*,*)  

    !print out
    open (unit=20,file='result.txt',action='write',status='replace')
    write (20,*) 'replacement rate',';','wage rate',';','return to capital',';',&
    'average consumption',';','capital stock',';','average income',';','average utility'
    write (20,*) theta,';',tau_s0,';',w,';',r_0,';',avcon,';',K1,';',avinc,';',omiga
    write (20,*) 'output',';','capital output ratio'
    write (20,*) B*K1**(1-alpha)*N1**alpha,';',K1/(B*K1**(1-alpha)*N1**alpha)
	write (20,*) 'Consumption'
    do j=1,jmax
    write (20,*) c(j)
    end do
    close (20)
end program