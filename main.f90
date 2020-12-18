!----------------------------------------------------------------------
!  Call BHMIE_RT
!
!  This calling program performes some test cases used in MIEV0.f 
!  developed by W. Wiscombe. I consider here the test case from #6 to #19.
!  The parameter set of each test case is as follows.
!
!                         n       k       x
!  --------------------------------------------
!  Test case 6  :       0.750   0.000   0.101
!  Test case 7  :       0.750   0.000   1.0e1
!  Test case 8  :       0.750   0.000   1.0e3
!  Test case 9  :       1.330   1.e-5   1.0e0
!  Test case 10 :       1.330   1.e-5   1.0e2
!  Test case 11 :       1.330   1.e-5   1.0e4
!  Test case 12 :       1.500   1.000   0.055
!  Test case 13 :       1.500   1.000   0.056
!  Test case 14 :       1.500   1.000   1.000
!  Test case 15 :       1.500   1.000   1.0e2
!  Test case 16 :       1.500   1.000   1.0e4
!  Test case 17 :       1.0e1   1.0e1   1.0e0
!  Test case 18 :       1.0e1   1.0e1   1.0e2
!  Test case 19 :       1.0e1   1.0e1   1.0e4
!
!----------------------------------------------------------------------
program call_bhmiert
use types
implicit none
integer::i,j
integer,parameter::nang=181
real(kind=dp)::qe,qs,qa,qb,gs,x
real(kind=dp)::dang,ncoefi,kcoefi
real(kind=dp),dimension(2*nang-1)::ANG,S11,S12,S33,S34
complex(kind=dp)::refrel
complex(kind=dp),dimension(2*nang-1)::S1,S2

dang=90.0_dp/dble(nang-1)
do j=1,2*nang-1
        ang(j)=dang*dble(j-1)
enddo
open(10,file="wiscombe_benchmark.out",status="unknown")
do i=6,19
        call wiscombe_testcases(i,ncoefi,kcoefi,x)
        refrel = cmplx(ncoefi,kcoefi)
        call bhmie_rt(x,refrel,nang,S1,S2,qe,qs,qb,gs)
        do j=1,2*nang-1
                S11(j)=0.5_dp*abs(S2(j))*abs(S2(j))
                S11(j)=S11(j)+0.5_dp*abs(S1(j))*abs(S1(j))
                S12(j)=0.5_dp*abs(S2(j))*abs(S2(j))
                S12(j)=S12(j)-0.5_dp*abs(S1(j))*abs(S1(j))       
                S33(j)=dble(S1(j)*conjg(S2(j)))
                S34(j)=aimag(S2(j)*conjg(S1(j)))
        enddo
        write(10,1000) " Wiscombe Test Case # = ",i
        write(10,1100) " Size parameter       = ",x
        write(10,1200) " Refractive indices m = ",ncoefi," + ",kcoefi," i"
        write(10,1100) " Qext  = ",qe
        write(10,1100) " Qsca  = ",qs
        write(10,1100) " Qback = ",qb
        write(10,1100) " g     = ",gs
        write(10,2000) "ang (deg)","S11","S12","S33","S34"
        do j=1,2*nang-1
        write(10,1300) ang(j),S11(j),S12(j),S33(j),S34(j)
        enddo
        write(10,*)
enddo

1000 format('#',A24,I5)
1100 format('#',A24,1PE15.5)
1200 format('#',A24,1PE15.5,A3,1PE15.5,A2)
1300 format(' ',1P7E15.5)
2000 format('#',7A15)
stop
end program call_bhmiert

subroutine wiscombe_testcases( itest , refre, refim, xx )
use types
integer::itest
real(kind=dp)::refre,refim,xx
if(itest .eq. 6) then
        !case 6
        refre = 0.750_dp
        refim = 0.000_dp
        xx    = 0.101_dp
elseif(itest .eq. 7) then
        !case 7
        refre = 0.750_dp
        refim = 0.000_dp
        xx    = 1.0e1_dp
elseif(itest .eq. 8) then
        !case 8
        refre = 0.750_dp
        refim = 0.000_dp
        xx    = 1.0e3_dp
elseif(itest .eq. 9) then
        !case 9
        refre = 1.330_dp
        refim = 1.0e-5_dp
        xx    = 1.0e0_dp
elseif(itest .eq. 10) then
        !case 10
        refre = 1.330_dp
        refim = 1.0e-5_dp
        xx    = 1.0e2_dp
elseif(itest .eq. 11) then
        !case 11
        refre = 1.330_dp
        !refre = 3.0_dp
        refim = 1.0e-5_dp
        xx    = 1.0e4_dp
elseif(itest .eq. 12) then
        !case 12
        refre = 1.5_dp
        refim = 1.0_dp
        xx    = 0.055_dp
elseif(itest .eq. 13) then
        !case 13
        refre = 1.5_dp
        refim = 1.0_dp
        xx    = 0.056_dp
elseif(itest .eq. 14) then
        !case 14
        refre = 1.5_dp
        refim = 1.0_dp
        xx    = 1.0_dp
elseif(itest .eq. 15) then
        !case 15
        refre = 1.5_dp
        refim = 1.0_dp
        xx    = 1.0e2_dp
elseif(itest .eq. 16) then
        !case 16
        refre = 1.5_dp
        refim = 1.0_dp
        xx    = 1.0e4_dp
elseif(itest .eq. 17) then
        !case 17
        refre = 10.0_dp
        refim = 10.0_dp
        xx    = 1.0e0_dp
elseif(itest .eq. 18) then
        !case 18
        refre = 10.0_dp
        refim = 10.0_dp
        xx    = 1.0e2_dp
elseif(itest .eq. 19) then
        !case 19
        refre = 10.0_dp
        refim = 10.0_dp
        xx    = 1.0e4_dp
else
        print *, 'error'
        stop
endif
return
end subroutine wiscombe_testcases
