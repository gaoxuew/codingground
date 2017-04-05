program surfacestate
implicit none
!---------------------------------------------------------------------------
! Read PROCAR file, and extract the contribution from each atom
! Identify the surface state of slab
!---------------------------------------------------------------------------

INTEGER :: nkpoints, nbands, natoms
INTEGER :: kp, band, niont,  ion, k, j, nb

REAL :: pt1, pt2, pt3
REAL, allocatable :: wt(:), eig_ev(:,:), occ(:,:), oc(:,:,:,:)
!INTEGER :: ispin

REAL, PARAMETER :: e = 2.71828, pi = 3.141592


open(7,file='PROCAR') 


!    write(*,*) 'Spin polarized calculation? (no=1,yes=2):'
!    read (*,*) ispin
!    if ((ispin.ne.1).and.(ispin.ne.2)) then
!    write(*,*) ' INPUT ERROR, ispin must equal to 1 or 2 '
!    stop
!    endif
    
!    if (ispin.eq.1) then
        open(31,file='bandstructure.dat') 
!      elseif (ispin.eq.2) then
!        open(31,file='bandstructre_up.dat') 
!        open(41,file='bandstructure-dow.dat') 
!    endif
    
    read(7,*)
    read(7,104) nkpoints,nbands,natoms
104 format(16x,i3,20x,i5,19x,i4)
print *, "nkpoints, nbands, natoms", nkpoints, nbands, natoms


allocate (wt(nkpoints))
allocate (eig_ev(nkpoints, nbands), occ(nkpoints, nbands), oc(nkpoints,nbands,natoms,10))


        do k = 1,nkpoints
          read(7,*)
          read(7,105) kp,pt1,pt2,pt3,wt(k)
 105      format(10x,i3,5x,3f11.8,13x,f11.8)
!         print *, kp,pt1,pt2,pt3,wt(k)                  !testing

          read(7,*)
          do  nb = 1,nbands
            read(7,106) band,eig_ev(k,nb),occ(k,nb)
 106        format(4x,i4,9x,f14.8,7x,f12.8)
!            print *, band,eig_ev(k,nb),occ(k,nb)        ! testing

!            eig_ev(k,nb) = eig_ev(k,nb)-fermi 
!            if (eig_ev(k,nb) .gt. emax) emax = eig_ev(k,nb)
!            if (eig_ev(k,nb) .lt. emin) emin = eig_ev(k,nb)
            read(7,*) 
            read(7,*) 

            niont = natoms +1
            if (natoms .eq. 1) niont = 1
            do  ion = 1,niont
              read(7,107) (oc(k,nb,ion,j),j=1,10)
 107          format(3x,10f7.3)
            print *, (oc(k,nb,ion,j),j=1,10)
            
            enddo
            read(7,*)
          enddo
          enddo
    
    
    
    
    
    






end program surfacestate

