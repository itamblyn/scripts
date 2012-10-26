!=============================================================================
!
! Utilities:
!
! (1) wfnbinasc         Originally By JLL       Last Modified 7/1/2008 (JRD)
!
!     Converts a CWFE or WFN0 file written by paratec v4.5 or later into an
!     ASCII file.
!
!==============================================================================

#include "f_defs.h"

      program wfnbinasc

      use global_m
      implicit none

      type(crystal) crys
      type(symmetry) syms
      type(kpoints) kp
      type(gspace) gvec

!---------------------------
! Local variables

      character*80 infile, outfile
      integer kx,ky,kz, &
      i,irk,j,n,nkpt,nargs, ilen, ierror

      ! iflag = 1 real version, iflag = 2 complex version
      integer iflag

      integer :: cell_symmetry
      real(DP) :: al,a(3,3),bl,b(3,3)

      integer, allocatable :: ifmin(:,:),ifmax(:,:)
      SCALAR, allocatable :: zc(:,:)

#ifndef GNU
      integer iargc
      external iargc
      external getarg
#endif

! Determine REAL version or COMPLEX version

#ifdef CPLX
      iflag=2
#else
      iflag=1
#endif

      if(iflag .eq. 1) then
        write(6,*) 'Real version is used to convert file'
      elseif(iflag .eq. 2) then
        write(6,*) 'Complex version is used to convert file'
      else
        write(0,*) 'Error in iflag! : ',iflag
      end if

!---------------------------------
! Get file names from command-line arguments

      nargs = iargc()

      if (nargs .ne. 2) then
        call die('Usage: wfnbinasc name_wfn_bin name_wfn_asc')
      endif

      call getarg(1,infile)
      call getarg(2,outfile)

! Open units

      open(unit=7,file=TRUNC(infile),form='unformatted')
      open(unit=8,file=TRUNC(outfile),form='formatted')
      write(6,*) 'Converting file ',TRUNC(infile),' into ascii format'

! Read/write data

      read(7) ((crys%bdot(i,j),i=1,3),j=1,3),bl,((b(i,j),i=1,3),j=1,3)
      write(8,*) ((crys%bdot(i,j),i=1,3),j=1,3),bl,((b(i,j),i=1,3),j=1,3)
      read(7) crys%celvol,al,((a(i,j),i=1,3),j=1,3)
      write(8,*) crys%celvol,al,((a(i,j),i=1,3),j=1,3)

      read(7) syms%ntran,cell_symmetry
      write(8,*) syms%ntran,cell_symmetry
      do n=1,syms%ntran
        read(7) ((syms%mtrx(n,i,j),i=1,3),j=1,3)
        write(8,*) ((syms%mtrx(n,i,j),i=1,3),j=1,3)
      enddo
      do n=1,syms%ntran
        read(7) (syms%tnp(n,i),i=1,3)
        write(8,*) (syms%tnp(n,i),i=1,3)
      enddo

      read(7) kp%nspin
      write(8,*) kp%nspin
      read(7) kp%nrk
      write(8,*) kp%nrk
      read(7) (kp%kgrid(i),i=1,3)
      write(8,*) (kp%kgrid(i),i=1,3)
      read(7) (kp%shift(i),i=1,3)
      write(8,*) (kp%shift(i),i=1,3)
      SAFE_ALLOCATE(kp%w, (kp%nrk))
      read(7) (kp%w(i),i=1,kp%nrk)
      write(8,*) (kp%w(i),i=1,kp%nrk)
      SAFE_DEALLOCATE(kp%w)
      SAFE_ALLOCATE(kp%rk, (3,kp%nrk))
      do j=1,kp%nrk
        read(7) (kp%rk(i,j),i=1,3)
        write(8,*) (kp%rk(i,j),i=1,3)
      enddo
      SAFE_DEALLOCATE(kp%rk)

      read(7) kp%mnband
      write(8,*) kp%mnband
      SAFE_ALLOCATE(ifmin, (kp%nrk,kp%nspin))
      SAFE_ALLOCATE(ifmax, (kp%nrk,kp%nspin))
      read(7) ((ifmin(i,j),i=1,kp%nrk),j=1,kp%nspin)
      write(8,*) ((ifmin(i,j),i=1,kp%nrk),j=1,kp%nspin)
      read(7) ((ifmax(i,j),i=1,kp%nrk),j=1,kp%nspin)
      write(8,*) ((ifmax(i,j),i=1,kp%nrk),j=1,kp%nspin)
      SAFE_DEALLOCATE(ifmin)
      SAFE_DEALLOCATE(ifmax)

      SAFE_ALLOCATE(kp%el, (kp%mnband,kp%nrk,kp%nspin))
      do n=1,kp%nspin
        do j=1,kp%nrk
          read(7) (kp%el(i,j,n),i=1,kp%mnband)
          write(8,*) (kp%el(i,j,n),i=1,kp%mnband)
        enddo
      enddo
      SAFE_DEALLOCATE(kp%el)

      read(7) (gvec%kmax(i),i=1,3)
      write(8,*) (gvec%kmax(i),i=1,3)
      read(7) gvec%ng
      write(8,*) gvec%ng
      SAFE_ALLOCATE(gvec%kx, (gvec%ng))
      SAFE_ALLOCATE(gvec%ky, (gvec%ng))
      SAFE_ALLOCATE(gvec%kz, (gvec%ng))
      do i=1,gvec%ng
        read(7) gvec%kx(i),gvec%ky(i),gvec%kz(i)
        write(8,*) gvec%kx(i),gvec%ky(i),gvec%kz(i)
      enddo
      SAFE_DEALLOCATE(gvec%kx)
      SAFE_DEALLOCATE(gvec%ky)
      SAFE_DEALLOCATE(gvec%kz)

! Output info

      write(6,*) ' crystal volume: ',crys%celvol
      write(6,*) ' number of spins: ',kp%nspin
      write(6,*) ' number of bands in file: ',kp%mnband

      do irk=1,kp%nrk
        read(7) nkpt
        write(8,*) nkpt
        do i=1,nkpt
          read(7) kx,ky,kz
          write(8,*) kx,ky,kz
        enddo
        SAFE_ALLOCATE(zc, (nkpt,kp%nspin))
        do n=1,kp%mnband
          read(7) ((zc(i,j),i=1,nkpt),j=1,kp%nspin)
          write(8,*) ((zc(i,j),i=1,nkpt),j=1,kp%nspin)
        enddo
        SAFE_DEALLOCATE(zc)
      enddo

      close(7)
      close(8)

      write(6,*) 'Done '

      end program wfnbinasc
