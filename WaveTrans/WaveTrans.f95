program wavetrans
implicit none                                       
complex(kind=8), allocatable :: coeff(:,:)
complex(kind=16), allocatable :: cener(:)
real(kind=8), allocatable :: occ(:)
integer, allocatable :: igall(:,:)
real :: a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),vtmp(3),sumkg(3),wk(3)
real :: xnrecl,xnspin,xnprec, xnwk,xnband,ecut, Vcell, b1mag, b2mag, b3mag
real :: phi12, phi13, phi23, sinphi123, vmag, phi123
integer :: nrecl, nspin, nprec, nwk, nband, npmaxA, npmaxB, npmax, &
     npmaxC, nb1maxA, nb2maxA,nb3maxA, nb1maxB, nb2maxB, nb3maxB, nb1maxC, &
     nb2maxC, nb3maxC, nb1max, nb2max, nb3max
integer :: i, j


!----------------------------
! These invariables are for loop over k and bands
!----------------------------
integer :: irec, isp, iwk, iband, nplane, ncnt, ig3, ig2, ig1, &
           ig1p, ig2p,ig3p, iplane
real :: xnplane, gtot, etot




     
!!$*   constant 'c' below is 2m/hbar**2 in units of 1/eV Ang^2 (value is
!!$*   adjusted in final decimal places to agree with VASP value; program
!!$*   checks for discrepancy of any results between this and VASP values)

real, parameter :: c = 0.262465831d0
!!$*   data c/0.26246582250210965422d0/ 

real, parameter :: pi=3.1415926


nrecl=24
open(unit=10,file='WAVECAR',access='direct',recl=nrecl)
read(unit=10,rec=1) xnrecl,xnspin,xnprec
close(unit=10)

nrecl=nint(xnrecl)
nspin=nint(xnspin)
nprec=nint(xnprec)

if(nprec.eq.45210) then
   write(6,*) '*** error - WAVECAR_double requires complex*16'
   stop
endif


write(6,*) 
write(6,*) 'record length  =',nrecl,' spins =',nspin, ' prec flag ',nprec
open(unit=10,file='WAVECAR',access='direct',recl=nrecl, status='old')
open(unit=11,file='GCOEFF.txt')
read(unit=10,rec=2) xnwk,xnband,ecut,(a1(j),j=1,3),(a2(j),j=1,3), (a3(j),j=1,3)
nwk=nint(xnwk)
nband=nint(xnband)
allocate(occ(nband))
allocate(cener(nband))

write(6,*) 'no. k points =',nwk
write(6,*) 'no. bands =',nband
write(6,*) 'max. energy =',ecut
write(6,*) 'real space lattice vectors:'
write(6,*) 'a1 =',(a1(j),j=1,3)
write(6,*) 'a2 =',(a2(j),j=1,3)
write(6,*) 'a3 =',(a3(j),j=1,3)
write(6,*) ' '
write(11,*) nspin
write(11,*) nwk
write(11,*) nband

!-------------------------------------------
!!$*   compute reciprocal properties
!-------------------------------------------
write(6,*) ' '
call vcross(vtmp,a2,a3)
Vcell=a1(1)*vtmp(1)+a1(2)*vtmp(2)+a1(3)*vtmp(3)
write(6,*) 'volume unit cell =', Vcell
call vcross(b1,a2,a3)
call vcross(b2,a3,a1)
call vcross(b3,a1,a2)
do j=1,3
   b1(j)=2.*pi*b1(j)/Vcell
   b2(j)=2.*pi*b2(j)/Vcell
   b3(j)=2.*pi*b3(j)/Vcell
enddo
b1mag=sqrt(b1(1)**2+b1(2)**2+b1(3)**2)
b2mag=sqrt(b2(1)**2+b2(2)**2+b2(3)**2)
b3mag=sqrt(b3(1)**2+b3(2)**2+b3(3)**2)
write(6,*) 'reciprocal lattice vectors:'
write(6,*) 'b1 =',(b1(j),j=1,3)
write(6,*) 'b2 =',(b2(j),j=1,3)
write(6,*) 'b3 =',(b3(j),j=1,3)
write(6,*) 'reciprocal lattice vector magnitudes:'
write(6,*) b1mag, b2mag, b3mag
write(6,*) ' '
write(11,*) (a1(j),j=1,3)
write(11,*) (a2(j),j=1,3)
write(11,*) (a3(j),j=1,3)
write(11,*) (b1(j),j=1,3)
write(11,*) (b2(j),j=1,3)
write(11,*) (b3(j),j=1,3)



phi12=acos((b1(1)*b2(1)+b1(2)*b2(2)+b1(3)*b2(3))/(b1mag*b2mag))
call vcross(vtmp,b1,b2)
vmag=sqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b3(1)*vtmp(1)+b3(2)*vtmp(2)+b3(3)*vtmp(3))/(vmag*b3mag)
nb1maxA=(sqrt(ecut*c)/(b1mag*abs(sin(phi12))))+1
nb2maxA=(sqrt(ecut*c)/(b2mag*abs(sin(phi12))))+1
nb3maxA=(sqrt(ecut*c)/(b3mag*abs(sinphi123)))+1
npmaxA=nint(4.*pi*nb1maxA*nb2maxA*nb3maxA/3.)
      
phi13=acos((b1(1)*b3(1)+b1(2)*b3(2)+b1(3)*b3(3))/(b1mag*b3mag))
call vcross(vtmp,b1,b3)
vmag=sqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b2(1)*vtmp(1)+b2(2)*vtmp(2)+b2(3)*vtmp(3))/(vmag*b2mag)
phi123=abs(asin(sinphi123))
nb1maxB=(sqrt(ecut*c)/(b1mag*abs(sin(phi13))))+1
nb2maxB=(sqrt(ecut*c)/(b2mag*abs(sinphi123)))+1
nb3maxB=(sqrt(ecut*c)/(b3mag*abs(sin(phi13))))+1
npmaxB=nint(4.*pi*nb1maxB*nb2maxB*nb3maxB/3.)
      
phi23=acos((b2(1)*b3(1)+b2(2)*b3(2)+b2(3)*b3(3))/(b2mag*b3mag))
call vcross(vtmp,b2,b3)
vmag=sqrt(vtmp(1)**2+vtmp(2)**2+vtmp(3)**2)
sinphi123=(b1(1)*vtmp(1)+b1(2)*vtmp(2)+b1(3)*vtmp(3))/(vmag*b1mag)
phi123=abs(asin(sinphi123))
nb1maxC=(sqrt(ecut*c)/(b1mag*abs(sinphi123)))+1
nb2maxC=(sqrt(ecut*c)/(b2mag*abs(sin(phi23))))+1
nb3maxC=(sqrt(ecut*c)/(b3mag*abs(sin(phi23))))+1 
npmaxC=nint(4.*pi*nb1maxC*nb2maxC*nb3maxC/3.)

nb1max=max0(nb1maxA,nb1maxB,nb1maxC)
nb2max=max0(nb2maxA,nb2maxB,nb2maxC)
nb3max=max0(nb3maxA,nb3maxB,nb3maxC)
npmax=min0(npmaxA,npmaxB,npmaxC)

write(6,*) 'max. no. G values; 1,2,3 =',nb1max,nb2max,nb3max
write(6,*) ' '

allocate (igall(3,npmax))
write(6,*) 'estimated max. no. plane waves =',npmax
allocate (coeff(npmax,nband))

!---------------------------------------------------------------------------
!--------------------------------------------------------------------------
!!$*   Begin loops over spin, k-points and bands
irec=2
do isp=1,nspin
   write(*,*) ' '
   write(*,*) '******'
   write(*,*) 'reading spin ',isp
   do iwk=1,nwk
      irec=irec+1
      read(unit=10,rec=irec) xnplane,(wk(i),i=1,3), &
           (cener(iband),occ(iband),iband=1,nband)
      nplane=nint(xnplane)
      write(6,*) 'k point #',iwk,'  input no. of plane waves =', &
           nplane
      write(6,*) 'k value =',(wk(j),j=1,3)
      write(11,*) (wk(j),j=1,3)
      
!!$*   Calculate plane waves
      ncnt=0
      do ig3=0,2*nb3max
         ig3p=ig3
         if (ig3.gt.nb3max) ig3p=ig3-2*nb3max-1
         do ig2=0,2*nb2max
            ig2p=ig2
            if (ig2.gt.nb2max) ig2p=ig2-2*nb2max-1
            do ig1=0,2*nb1max
               ig1p=ig1
               if (ig1.gt.nb1max) ig1p=ig1-2*nb1max-1
               do j=1,3
                  sumkg(j)=(wk(1)+ig1p)*b1(j)+ &
                       (wk(2)+ig2p)*b2(j)+(wk(3)+ig3p)*b3(j)
               enddo
               gtot=sqrt(sumkg(1)**2+sumkg(2)**2+sumkg(3)**2)
               etot=gtot**2/c
               if (etot.lt.ecut) then
                  ncnt=ncnt+1
                  igall(1,ncnt)=ig1p
                  igall(2,ncnt)=ig2p
                  igall(3,ncnt)=ig3p
               end if
            enddo
         enddo
      enddo
      if (ncnt.ne.nplane) then
         write(6,*) '*** error - computed no. != input no.'
         stop
      endif
      if (ncnt.gt.npmax) then
         write(6,*) '*** error - plane wave count exceeds estimate'
         stop
      endif
      
      do iband=1,nband
         irec=irec+1
         read(unit=10,rec=irec) (coeff(iplane,iband), &
              iplane=1,nplane)
      enddo
!!$*   output G values and coefficients
      do iband=1,nband
         write(11,*) iband,nplane
         write(11,560) cener(iband),occ(iband)
560      format('( ',g14.6,' , ',g14.6,' ) ',g14.6)            
         do iplane=1,nplane
            write(11,570) (igall(j,iplane),j=1,3), &
                 coeff(iplane,iband)
570         format(3i6,'  ( ',g14.6,' , ',g14.6,' )')     
         enddo
      enddo
   enddo
enddo

stop

end program wavetrans






!!$*
!!$*   routine for computing vector cross-product
!!$*
subroutine vcross(a,b,c)
  implicit none
  real :: a(3),b(3),c(3)
  
  a(1)=b(2)*c(3)-b(3)*c(2)
  a(2)=b(3)*c(1)-b(1)*c(3)
  a(3)=b(1)*c(2)-b(2)*c(1)
  return
end subroutine vcross

