! this is a hack of Shenyuan's code. It gives out an integral for looking at charge transfer. 
! Best to refer to her original source when doing new stuff

! This program is to read the charge or potential data from CHG or LOCPOT files
!    and calculate the average values in planes.
! INPUT: CHG or LOCPOT  -- read charge or local potential
!        planechg.in    -- which file to read, the orientation of the plane
! OUTPUT: planechg.dat  -- the average values of the data in each plane
! Shenyuan 10-20-2008 at LBL
! IMPORTANT: to be upgraded

PROGRAM avechgplane
IMPLICIT NONE
REAL(8) :: alatt,vec(3,3),lvec(3),vol,tt(3),vecxy(2,2),area,cutoff,shift, integral
! lvec: length of the cell
REAL(8), ALLOCATABLE :: pot(:,:,:,:),potplane(:,:),potave(:,:,:),posz(:),poszave(:,:),potsum(:)
! pot(ix,iy,iz,is),potplane(iz,is),potave(iz,iave,is),posz(iz),poszave(iz,iave)
REAL(8), ALLOCATABLE :: lave(:)
! lave(iave): average over lave along z
INTEGER :: i,j,k,ix,iy,iz,ns,is,scale,plane(3)
INTEGER :: laverage,lsub,ladd,iave,iposzave,nng,nnv,nxyz,lcutoff,ncutoff
! laverage: whether to average; iposzave: which point to average along z
! lcircle: whether average within a circle, npts: # of points in xy plane
INTEGER :: ngrid(3),atom(10),natom
! nng: # of points inplane; lnn: # of point to average along z
INTEGER, ALLOCATABLE :: lnn(:),nptsum(:,:)
! lnn(laverage),nptsum(npts,2),ix,iy, coordinates within the circle
CHARACTER(20) :: Fin,Fpot,Fout,A,output
CHARACTER(60) :: atoms

! ****************************   INITIAL   ******************************

Fin='planechg.in'
Fout='planechg.dat'
output='planechg.out'

!   ************  READ THE CONTROL PARAMETERS FROM POTPLANE.IN  **********

OPEN(10,FILE=Fin)
OPEN(11,FILE=Fout)
OPEN(20,FILE=output)
READ(10,'(A)')Fpot        ! which file to read: CHG or LOCPOT
WRITE(20,*)'Open file ',Fpot
READ(10,*)plane(3)        ! orientation of the plane: 1-x,2-y,3-z
WRITE(20,*)'The normal of the plane: ',plane(3)
READ(10,*)scale           ! scale the data by 1/vol: 1-yes,other-no
WRITE(20,*)'Whethre scale the data: ',scale
READ(10,*)shift           ! shift the data, maybe useful for LOCPOT
WRITE(20,*)'Shift the data by: ',shift
READ(10,*)ns            ! # of spin
WRITE(20,*)'number of spins: ',ns
READ(10,*)lcutoff      ! provide a cutoff for the data: 1-yes, other-no
IF(lcutoff == 1) THEN
     READ(10,*)cutoff    ! only sum over the data < cutoff
     WRITE(20,*)'cutoff of data: ',cutoff
ELSE
     cutoff=999999.0d0
END IF
READ(10,*)laverage        ! laverage over z: 0-no, other-# of average over z
IF(laverage > 0) THEN
    ALLOCATE(lave(laverage))
    ALLOCATE(lnn(laverage))
    READ(10,*)lave(1:laverage)
    WRITE(20,*)'average along z over ',lave(1:laverage), 'A'
END IF

!   ********************* READ CHG OR LOCPOT FILE ************************

OPEN(12,FILE=Fpot)
READ(12,'(A)')A
READ(12,*)alatt     ! read lattice constant
WRITE(20,*)'lattice constant ',alatt
DO i=1,3            ! read the lattice vectors
    READ(12,*)vec(i,1:3)
    lvec(i)=alatt*SQRT(vec(i,1)**2+vec(i,2)**2+vec(i,3)**2)         ! lvec(3)in A
    WRITE(20,*)vec(i,1:3),lvec(i)
END DO

READ(12,'(A)')atoms     ! read the number of atoms, max. 10 species
atoms=atoms(:41)//' 0 0 0 0 0 0 0 0 0'
READ(atoms,*)atom(1:10)
natom=SUM(atom)
READ(12,'(A)')A
DO i=1,natom            ! skip all the atoms
    READ(12,'(A)')A
END DO
READ(12,*)ngrid(1:3)        ! read the dimensions and the data
WRITE(20,*)'grid: ',ngrid(1:3)

ALLOCATE(pot(ngrid(1),ngrid(2),ngrid(3),ns))
ALLOCATE(potplane(ngrid(plane(3)),ns))
ALLOCATE(posz(ngrid(plane(3))))
!ALLOCATE(potsum(ns))

IF(laverage > 0) THEN
    ALLOCATE(potave(ngrid(plane(3)),laverage,ns))
    ALLOCATE(poszave(ngrid(plane(3)),laverage))
    ALLOCATE(potsum(ns))
END IF

READ(12,*)((((pot(ix,iy,iz,is),ix=1,ngrid(1)),iy=1,ngrid(2)),iz=1,ngrid(3)),is=1,ns)
WRITE(20,*) 'read data from ', Fpot, 'OK!'
pot=pot-shift

DO i=1,2                         ! determine the vectors in plane
    plane(i)=mod(plane(3)+i,3)
    IF(plane(i) == 0) plane(i)=3
END DO

nng=ngrid(plane(1))*ngrid(plane(2))
ALLOCATE(nptsum(2,nng))
potplane=0.0d0
DO iz=1,ngrid(plane(3))
    posz(iz)=lvec(plane(3))*DBLE(iz-0.5)/DBLE(ngrid(plane(3)))
END DO

DO is=1,ns
    DO k=1,ngrid(plane(3))      ! calculate the sum and average in planexy
        IF(lcutoff == 1) ncutoff=0
        DO j=1,ngrid(plane(2))
            do i=1,ngrid(plane(1))
                IF(plane(3) == 3) THEN
                    iz=k
                    iy=j
                    ix=i
                ELSE IF(plane(3) == 1)THEN
                    ix=k
                    iz=j
                    iy=i
                ELSE
                    iy=k
                    ix=j
                    iz=i
                END IF
                IF((lcutoff == 1) .AND. (pot(ix,iy,iz,is) > cutoff)) CYCLE
                potplane(k,is)=potplane(k,is)+pot(ix,iy,iz,is)
                IF(lcutoff == 1) ncutoff=ncutoff+1
            END DO
        END DO
        IF(lcutoff == 1) THEN
            potplane(k,is)=potplane(k,is)/DBLE(ncutoff)
        ELSE
            potplane(k,is)=potplane(k,is)/DBLE(nng) ! this is where the "average" happens
        END IF
    END DO
END DO
IF(scale == 1) THEN
    CALL crossprod(vec(plane(1),:),vec(plane(2),:),tt)
    CALL dotprod(vec(plane(3),:),tt,vol)
    vol=vol*(alatt**3)
    potplane=potplane/vol
END IF
!   potplane=potplane/(ngrid(1)*ngrid(2)*ngrid(3))
!
!   *********************** AVERAGE POTPLANE ALONG Z ********************** 

IF(laverage > 0)THEN        ! calculate the average value within lave
    DO iave=1,laverage
        lnn(iave)=NINT(lave(iave)/lvec(plane(3))*DBLE(ngrid(plane(3))))
        potsum=0.0d0
        DO i=1,lnn(iave)
            DO is=1,ns
                potsum(is)=potsum(is)+potplane(i,is)
            END DO
        END DO
        DO i=1,ngrid(plane(3))
            lsub=i
            ladd=i+lnn(iave)
            iposzave=i+(lnn(iave)-1)/2
            IF(ladd > ngrid(plane(3))) ladd=ladd-ngrid(plane(3))
            IF(iposzave > ngrid(plane(3))) iposzave=iposzave-ngrid(plane(3))
            IF(mod(lnn(iave),2) == 0) THEN
                IF(iposzave == ngrid(plane(3))) THEN
                    poszave(iposzave,iave)=posz(iposzave)+lvec(plane(3))/DBLE(ngrid(plane(3)))/2.0d0
                ELSE
                    poszave(iposzave,iave)=(posz(iposzave)+posz(iposzave+1))/2.0d0
                END IF
            ELSE
                poszave(iposzave,iave)=posz(iposzave)
            END IF
            DO is=1,ns
                potave(iposzave,iave,is)=potsum(is)/REAL(lnn(iave))
                potsum(is)=(potplane(ladd,is)-potplane(lsub,is))+potsum(is)
            END DO
        END DO
        WRITE(20,*)'Average over z ',iave,'OK!'
        WRITE(20,*)'Number of points: ',lnn(iave)
    END DO
END IF

! normalize the data
nnv = (ngrid(plane(1))*ngrid(plane(2))*ngrid(plane(3)) )  ! resuing a variable...bad practice

write(11,*) '# z, pot, potnorm, potnormint'

IF(laverage == 0) THEN
!    WRITE(11,300)ngrid(plane(1)),ngrid(plane(2)),ngrid(plane(3)),ns
    DO is=1,ns
        do i=1,ngrid(plane(3))
            integral = integral + nng*potplane(i,is)/nnv
            write(11,*) posz(i), potplane(i,is) !, nng*potplane(i,is)/nnv,  integral
        end do
        WRITE(11,'(/)')
    END DO
        WRITE(20,*)'Write data OK!'
ELSE
!    WRITE(11,300)ngrid(plane(1)),ngrid(plane(2)),ngrid(plane(3)),ns,(lnn(iave),iave=1,laverage)
    DO is =1,ns
        DO i=1,ngrid(plane(3))
            WRITE(11,100)posz(i),potplane(i,is),((poszave(i,iave),potave(i,iave,is)),iave=1,laverage)
        END DO
        WRITE(11,'(/)')
    END DO
    WRITE(20,*)'WRITE data OK!'
END IF

DEALLOCATE(pot)
DEALLOCATE(potplane)
DEALLOCATE(posz)
!DEALLOCATE(nptsum)
IF(laverage > 0)DEALLOCATE(potave,poszave,potsum)

CLOSE(10)
CLOSE(11)
CLOSE(12)
CLOSE(20)

100    FORMAT(F11.6,E20.11,6(F11.6,E20.11))
200    FORMAT(F11.6,E20.11)
300    FORMAT(20I6)
STOP
END PROGRAM avechgplane

!  ****  SUBROUTINE dotproduce  ****
SUBROUTINE dotprod(a,b,d)
IMPLICIT NONE
DOUBLE PRECISION :: a(3),b(3),d
d=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
RETURN
END SUBROUTINE dotprod

!  ****  SUBROUTINE crossprod  ****
SUBROUTINE crossprod(a,b,c)
IMPLICIT NONE
DOUBLE PRECISION :: a(3),b(3),c(3)
c(1)=a(2)*b(3)-a(3)*b(2)
c(2)=a(3)*b(1)-a(1)*b(3)
c(3)=a(1)*b(2)-a(2)*b(1)
RETURN
END SUBROUTINE crossprod
