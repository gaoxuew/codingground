program surfacestate
implicit none
!==============================
! Read PROCAR file, and extract the contribution from each atom
! Identify the surface state of slab
!==============================

INTEGER :: nkpoints, nbands, natoms
INTEGER :: ispin

REAL, PARAMETER :: e = 2.71828, pi = 3.141592


open(7,file='PROCAR') 


    write(*,*) 'Spin polarized calculation? (no=1,yes=2):'
    read (*,*) ispin
    if ((ispin.ne.1).and.(ispin.ne.2)) then
    write(*,*) ' INPUT ERROR, ispin must equal to 1 or 2 '
    stop
    endif





print *, "Hello world!", "pi = ", pi
end program surfacestate

