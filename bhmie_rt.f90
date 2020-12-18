!----------------------------------------------------------------------
!                           BHMIE code
!----------------------------------------------------------------------
! This is a modified version of the BHMIE code written by Bruce Draine.
!  https://www.astro.princeton.edu/~draine/scattering.html
!
! R.T. modified following points.
!- Big static arrays are now defined as dynamic array.
!- Implemented the Lentz's continued fraction method to estimate
!  the logarithmic derivative of Riccati Bessel function, 
!  which is safe for very large grains.
!- Implemented adjustable real and complex precision
!- Slightly modified xstop to match the Wiscombe (1980)'s criterion.
!
!                                              Ryo Tazaki (2020.12.16)
!----------------------------------------------------------------------
module types
implicit none
integer, parameter      :: dp = selected_real_kind(P=15)
end module types

subroutine bhmie_rt(x,refrel,nang,S1,S2,Qext,Qsca,Qback,gsca)
use types
implicit none
integer::nang,j,jj,n,nn,nstop,nmx
complex(kind=dp)::refrel,y,dini
complex(kind=dp),dimension(1:2*nang-1)::S1,S2
complex(kind=dp),allocatable,dimension(:)::d,ajy
complex(kind=dp)::an,bn,an1,bn1
complex(kind=dp)::xi1,xi,dpcx
real(kind=dp)::x,qext,qsca,gsca,qback
real(kind=dp)::en,fn,p,xstop,ymod,nu
real(kind=dp)::psi0,psi1,psi,chi0,chi1,chi
real(kind=dp)::dang,pii,theta
real(kind=dp),dimension(1:nang)::amu,pi0,pi1,pi,tau
real(kind=dp) realpart
real(kind=dp) imagpart
realpart(dpcx)=(dble(dpcx))
imagpart(dpcx)=(dimag(dpcx))

! some checks...
y     = refrel*x
ymod  = abs(y)
xstop = x + 4.05_dp * x ** (1.0/3.0) + 2.0_dp
nstop = nint(xstop)

!*** require nang.ge.1 in order to calculate scattering intensities
pii  = 4.0_dp*atan(1.0_dp)
dang = 0.0_dp
if(nang.gt.1)dang=0.5_dp*pii/dble(nang-1)
do j=1,nang
        theta =real(j-1,kind=dp)*dang
        amu(j)=cos(theta)
enddo
do j=1,nang
        pi0(j)=0.0_dp
        pi1(j)=1.0_dp
enddo
nn=2*nang-1
do j=1,nn
        S1(j)=cmplx(0.0_dp,0.0_dp)
        S2(j)=cmplx(0.0_dp,0.0_dp)
enddo

!----------------------------------------------------------------------
!       logarithmic derivative
!----------------------------------------------------------------------
! 1. Standard BHMIE algorithm ( Draine's code as well )
! Warming-up terms of 15 may not be sufficient for very large spheres.
!----------------------------------------------------------------------
!nmx = nint(max(xstop,ymod)) + 15
!nmx = nint(max(xstop,ymod)) + 100
!allocate(d(1:nmx))
!if(nmx .gt. 150000) then !1.5e5 is from BTD's code.
!        write(*,*) "Error: nmx > nmxx=",150000,"for |m|x=",ymod
!        stop
!endif
!calculate D_n(mx) by downward recurrence.
!Initial value is set as D(mx)=0+0i at n=NMX.
!d(nmx)=cmplx(0.0_dp,0.0_dp)
!do n=1,nmx-1
!        en=real(nmx-n+1,kind=dp)
!        d(nmx-n) = en/y - (1.0_dp/(d(nmx-n+1) + en/y))
!enddo

!----------------------------------------------------------------------
! 2. The Lentz's continued fraction method
!    Reference: W. J. Lentz, Appl. Opt. 15. 668. (1976)
!----------------------------------------------------------------------
! Finding initial value of the logarithmic derivative "d" at 
! the starting order of downward recurrence (n=nmx).
!----------------------------------------------------------------------
nmx = nint(max(xstop,ymod)) + 15
allocate(d(1:nmx),ajy(1:nmx))
! Equation (5) in Lentz (1976):
nu = real(nmx,kind=dp)+0.5_dp
ajy(1) = 2.0_dp * nu / y
do j=2,nmx
       ajy(j) = (-1.0_dp)**real(j-1,kind=dp)*2.0_dp*(nu+real(j-1,kind=dp))/y
enddo
dini = ajy(nmx)
do j=1,nmx-1
        dini = ajy(nmx-j) + 1.0_dp / dini
enddo
deallocate(ajy)
! Equation (3) in Lentz (1976):
d(nmx)=-real(nmx,kind=dp)/y+dini
! usual downward recurrence:
do n=1,nmx-1
        en=real(nmx-n+1,kind=dp)
        d(nmx-n) = (en/y) - (1.0_dp/(d(nmx-n+1) + en/y))
enddo
!----------------------------------------------------------------------

psi0 = cos(x)
psi1 = sin(x)
chi0 =-sin(x)
chi1 = cos(x)
xi1  = cmplx(psi1,-chi1)
!
qsca = 0.0_dp
gsca = 0.0_dp
p    =-1.0_dp  
do n=1,nstop
        en = real(n,kind=dp)
        fn = (2.0_dp*en+1.0_dp)/(en*(en+1.0_dp))

        !Calculate psi and chi via upward recurrence:
        psi = (2.0_dp*en-1.0_dp)*psi1/x - psi0
        chi = (2.0_dp*en-1.0_dp)*chi1/x - chi0
        xi = cmplx(psi,-chi)

        !*** store previous values of an and bn for use
        !    in computation of g=<cos(theta)>
        if(n.gt.1)then
                an1=an
                bn1=bn
        endif

        !*** compute an and bn:
        an = (d(n)/refrel+en/x)*psi-psi1
        an = an/((d(n)/refrel+en/x)*xi-xi1)
        bn = (refrel*d(n)+en/x)*psi-psi1
        bn = bn/((refrel*d(n)+en/x)*xi-xi1)

        !*** augment sums for qsca and g=<cos(theta)>
        qsca=qsca+(2.0_dp*en+1.0_dp)*(abs(an)**2.0+abs(bn)**2.0)
        gsca=gsca+((2.0_dp*en+1.0_dp)/(en*(en+1.0_dp)))*&
        &        (realpart(an)*realpart(bn)+imagpart(an)*imagpart(bn))
        if(n.gt.1)then
            gsca=gsca+((en-1.0_dp)*(en+1.0_dp)/en)*&
        &      (realpart(an1)*realpart(an)+imagpart(an1)*imagpart(an)+&
        &      realpart(bn1)*realpart(bn)+imagpart(bn1)*imagpart(bn))
        endif

        !c*** now calculate scattering intensity pattern
        !    first do angles from 0 to 90
        do j=1,nang
                jj=2*nang-j
                pi(j)=pi1(j)
                tau(j)=en*amu(j)*pi(j)-(en+1.0_dp)*pi0(j)
                S1(j)=S1(j)+fn*(an*pi(j)+bn*tau(j))
                S2(j)=S2(j)+fn*(an*tau(j)+bn*pi(j))
        enddo

        !*** now do angles greater than 90 using pi and tau from
        !    angles less than 90.
        !    p=1 for n=1,3,...; p=-1 for n=2,4,...
        p=-p
        do j=1,nang-1
                jj=2*nang-j
                S1(jj)=S1(jj)+fn*p*(an*pi(j)-bn*tau(j))
                S2(jj)=S2(jj)+fn*p*(bn*pi(j)-an*tau(j))
        enddo

        !Preparation for next step
        psi0 = psi1
        psi1 = psi
        chi0 = chi1
        chi1 = chi
        xi1 = cmplx(psi1,-chi1)

        !*** compute pi_n for next value of n
        !    for each angle j, compute pi_n+1
        !    from pi = pi_n , pi0 = pi_n-1
        do j=1,nang
                pi1(j)=((2.0_dp*en+1.0_dp)*amu(j)*pi(j)-(en+1.0_dp)*pi0(j))/en
                pi0(j)=pi(j)
        enddo
enddo

!*** have summed sufficient terms.
!    now compute qsca,qext,qback,and gsca
gsca=2.0_dp*gsca/qsca
qsca=(2.0_dp/(x*x))*qsca
qext=(4.0_dp/(x*x))*realpart(S1(1))
qback=(4.0_dp*(ABS(S1(2*NANG-1))/x)**2._dp) 

deallocate(d)

return
end subroutine bhmie_rt
