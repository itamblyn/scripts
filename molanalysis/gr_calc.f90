program gr_calc
  ! A program to calculate atom-atom radial distribution functions from
  ! castep simulations.
  ! 11/7/06:  Modified for non-orthorhombic boxes (LEF).
  !
  implicit none
  integer, parameter :: MAXATOM = 1000                 ! Max number of atoms.
  integer, parameter :: MAXHISTO = 5000                ! Maximum length of the histogram.
  integer, parameter :: MAXBLOCK = 5000

  double precision rx, ry, rz                          ! Displacement between atom pairs.
  double precision x(MAXATOM), y(MAXATOM), z(MAXATOM)  ! The current atom positions.
  double precision x0(MAXATOM,MAXBLOCK), &
           y0(MAXATOM,MAXBLOCK), z0(MAXATOM,MAXBLOCK)  ! The initial atom positions.
  double precision distance(MAXBLOCK)                  ! The distance traveled as a function of time.
  double precision nc_int
  integer dcount(MAXBLOCK)                             ! Number of samples for diffusion constant.
  integer nindex1                                      ! Number of atoms of type 1
  integer index1(MAXATOM)                              ! Indices of type 1 atom.
  integer nindex2                                      ! Number of atoms of type 2
  integer index2(MAXATOM)                              ! Indices of type 2 atoms.
  integer i, j, k, l, kk, ll, ia, ib, ic, it
  integer isuper                                       ! Number of super-cells used
                                                       ! in finding g(r).
                                                       ! File name.  
  character(len=80) line, file
  character(len=1) el                                  ! Element character.
  integer maxlines                                     ! Maximum number of frames in
                                                       ! trajectory file.
  integer, parameter :: in = 20                        ! Unit number of input file.
  integer natom                                        ! Number of atoms.
  integer istat                                        ! Read status.
  integer histo(MAXHISTO)                              ! Histogram of distances.
  integer histomax                                     ! Maximum occupied element of histogram file.
  double precision nc(MAXHISTO)                        ! coordination number       
  double precision time                                ! Time of current step.
  double precision ax, ay, az                          ! Components of the box vectors a, b, c.
  double precision bx, by, bz 
  double precision cx, cy, cz
  double precision dr                                  ! Step size used in the histogram.
  double precision r2, r
  double precision tmax, tmin                          ! Maximum and minimum times used in the current time block.
  double precision rl, ru, gr(MAXHISTO), dv, vfac, pi
  double precision sum
  double precision rbin(MAXHISTO)                      ! Binned radii.
  integer distance_block                               ! Averaging block for distances.
  integer nconfig                                      ! Number of configurations in the histogram.
  integer nstep                                        ! Number of trajectory steps read in.
  integer icut                                         ! Cutoff index for g(r).  Results beyond this index are not
                                                       ! reliable.
  integer iformat                                      ! Trajectory input format type.
  double precision rcut                                ! Cutoff distance for g(r).
  double precision tblock                              ! Time block size
  double precision dtime                               ! MD time step in trajectory file.
  integer nt                                           ! Number of time blocks.
  double precision alen, blen, clen                    ! length of the a, b, and c lattice vectors.
  integer itype1, itype2                               ! Atom types to calculate g(r) for.
  integer m, n, idx
  double precision ax1, ay1, az1                       ! Components of the box vectors a, b, c.
  double precision bx1, by1, bz1 
  double precision cx1, cy1, cz1
  double precision invbox(3,3)                         ! The inverse of the box matrix.
  double precision boxvol                              ! The box volume.
  double precision gasfac                              ! Ideal gas density of A-B pairs.
  double precision r1                                  ! Bond distance based on g(r) = 1.
  double precision rmin                                ! Bond distance based on minimum of g(r).
  integer k1                            
  double precision dconst                             ! Diffusion constant.
  integer offset, idxt
  integer iend

  maxlines = 100000000
  nindex1 = 0 
  ax = 0.0
  natom = 0 

  read(*,*) file
  read(*,*) iformat
  write(*,*) 'Taking input from ', file
  if ( iformat .ne. 4 ) then
     read(*,*) natom
     write(*,*) 'Number of atoms = ', natom
  end if
  read(*,*) itype1, itype2

  if ( iformat .ne. 4 ) then
     read(*,*) nindex1
     read(*,*) ( index1(j), j = 1, nindex1 )
     read(*,*) nindex2
     read(*,*) ( index2(j), j = 1, nindex2 )
     read(*,*) ax, ay, az
     read(*,*) bx, by, bz
     read(*,*) cx, cy, cz
  endif

  read(*,*) dr
  read(*,*) isuper
  read(*,*) tblock
  read(*,*) nt
  if ( iformat .eq. 2 .or. iformat .eq. 3 .or. iformat .eq. 4 ) then
     read(*,*) dtime
  endif
  read(*,*) distance_block

  if ( iformat .eq. 1 ) then
     write(*,*) 'Trajectory file is from casstep'
  else if ( iformat .eq. 2 ) then
     write(*,*) 'Trajectory file is from fireball'
  else if ( iformat .eq. 3 ) then
     write(*,*) 'Trajectory is is in xyz format'
  else if ( iformat .eq. 4 ) then
     write(*,*) 'Trajectory file is from DFTB'
  endif

  write(*,*) 'Calculating g(r) for element types', itype1, 'and', itype2
  if ( iformat .ne. 4 ) then
     write(*,*) 'Number of atoms of type 1 = ', nindex1
     write(*,*) 'Atom indices of type1 = '
     write(*,*) ( index1(j), j = 1, nindex1 )
     write(*,*) 'Number of atoms of type 2 = ', nindex2
     write(*,*) ( index2(j), j = 1, nindex2 )
     write(*,*) 'Cell vector a = ', ax, ay, az
     write(*,*) 'Cell vector b = ', bx, by, bz
     write(*,*) 'Cell vector c = ', cx, cy, cz
  else
     write(*,*) 'Atom indices and cell parameters taken from the DFTB file'
  endif

  write(*,*) 'The bin size = ', dr
  write(*,*) 'The number of periodic images in each direction = ', isuper
  write(*,*) 'The time block = ', tblock
  write(*,*) 'The number of blocks = ', nt
  write(*,*) 'The MD time step = ', dtime
  write(*,*) 'Block averaging used in diffusion constants = ', distance_block

  TIMELOOP : do it = 1, nt
     tmax = it * tblock
     tmin = (it-1) * tblock
     print*,'calculate g(r) from', tmin, 'to', tmax
     open(in, status='old', FILE=file) 
     do i = 1, MAXHISTO
        histo(i) = 0 
     end do

     do i = 1, MAXBLOCK
        distance(i) = 0.0d0
        dcount(i) = 0
     end do

     histomax = 0
     nstep = 0

     READLINES : do j = 1, maxlines
        
        if ( iformat .eq. 1 ) then
           ! castep output.
           read(in,'(A20)',iostat=istat) line
           read(in,*,iostat=istat) time

           if ( istat .ne. 0 ) then
              exit READLINES
           endif
           do k = 1 , natom 
              read(in,*,iostat=istat) x(k), y(k), z(k)
              if ( istat .ne. 0 ) then
                 exit READLINES
              endif
           end do
        else if ( iformat .eq. 2 ) then
           ! fireball output
           time = dtime * (j-1)
           read(in,80,iostat=istat) k
80         format(10x,i20/)
           if ( istat .ne. 0 ) then
              exit READLINES
           endif
           if ( k .ne. natom ) then
              write(*,*) 'Error: atom number check failed', k, natom
              exit READLINES
           endif
           do k = 1 , natom 
              read(in,90,iostat=istat) x(k), y(k), z(k)
90            format(9X,f14.6,f14.6,f14.6)
              if ( istat .ne. 0 ) then
                 exit READLINES
              endif
           end do
        else if ( iformat .eq. 3 ) then
           ! xyz format.
           time = dtime * (j-1)
           read(in,91,iostat=istat) k
91         format(i5/)
           if ( istat .ne. 0 ) then
              exit READLINES
           endif
           if ( k .ne. natom ) then
              write(*,*) 'Error: atom number check failed', k, natom
           endif
           do k = 1 , natom 
              read(in,*,iostat=istat) el, x(k), y(k), z(k)
              if ( istat .ne. 0 ) then
                 exit READLINES
              endif
           end do
        else if ( iformat .eq. 4 ) then
           ! DFTB Gen file.
           time = dtime * (j-1)
           read(in,91,iostat=istat) k
           if ( istat .ne. 0 ) then
              exit READLINES
           endif
           if ( natom .eq. 0 ) then
              natom = k 
           else if ( k .ne. natom ) then
              write(*,*) 'Error: atom number check failed', k, natom
           endif
           l = 1
           m = 1
           do k = 1 , natom 
              read(in,*,iostat=istat) n, idx, x(k), y(k), z(k)
              if ( istat .ne. 0 ) then
                 exit READLINES
              endif
              if ( k .ne. n ) then
                 write(*,*) 'Failed index consistency check'
                 stop
              end if
              if ( idx .eq. itype1 ) then
                 index1(l) = n
                 l = l + 1
              endif
              if ( idx .eq. itype2 ) then
                 index2(m) = n 
                 m = m + 1
              endif
              if ( mod(j-1,distance_block) .eq. 0 ) then
                 ll = (j-1)/distance_block + 1
                 if ( ll .gt. MAXBLOCK ) then
                    stop 'Recompile with a larger MAXBLOCK or use a larger diffusion time block input parameter'
                 endif
                 x0(k,ll) = x(k) 
                 y0(k,ll) = y(k)
                 z0(k,ll) = z(k)
              endif
           end do
           if ( nindex1 .gt. 0 ) then
              if ( l - 1 .ne. nindex1 ) then
                 write(*,*) 'Match on number of index 1 failed', l, nindex1
                 stop
              endif
              if ( m - 1 .ne. nindex2 ) then
                 write(*,*) 'Match on number of index 2 failed'
                 stop
              endif
           else
              nindex1 = l - 1
              nindex2 = m - 1 

              write(*,*) 'Number of atoms of type1 = ', nindex1
              write(*,*) 'Atom indices of type1 = '
              do k = 0, nindex1, 10
                 if ( k + 10 > nindex1 ) then
                    iend = nindex1 - k
                 else
                    iend = 10
                 endif
                 write(*,'(10I6)') ( index1(l+k), l = 1, iend )
              end do
              write(*,*) 'Number of atoms of type 2 = ', nindex2
              write(*,*) 'Atom indices of type2 = '
              do k = 0, nindex2, 10
                 if ( k + 10 > nindex2 ) then
                    iend = nindex2 - k
                 else
                    iend = 10
                 endif
                 write(*,'(10I6)') ( index2(l+k), l = 1, iend )
              end do
           endif

           !             Dummy numbers read first.
           !          Read in the box parameters.           
           read(in,*,iostat=istat) ax1, ay1, az1
           read(in,*,iostat=istat) ax1, ay1, az1
           read(in,*,iostat=istat) bx1, by1, bz1
           read(in,*,iostat=istat) cx1, cy1, cz1
!           if ( ax .gt. 0.0 ) then
!              if ( abs(ax1-ax) > 1.0e-10 .or. abs(ay1-ay) > 1.0e-10 .or. &
!                   abs(az1-az) > 1.0e-10 ) then
!                 write(*,*) 'Box vector A match test failed'
!                 stop
!              endif
!              if ( abs(bx1-bx) > 1.0e-10 .or. abs(by1-by) > 1.0e-10 .or. &
!                   abs(bz1-bz) > 1.0e-10 ) then
!                 write(*,*) 'Box vector B match test failed'
!                 stop
!              endif
!              if ( abs(cx1-cx) > 1.0e-10 .or. abs(cy1-cy) > 1.0e-10 .or. &
!                   abs(cz1-cz) > 1.0e-10 ) then
!                 write(*,*) 'Box vector C match test failed'
!                 stop
!              endif
!           else
              ax = ax1
              ay = ay1
              az = az1
              bx = bx1
              by = by1
              bz = bz1
              cx = cx1
              cy = cy1
              cz = cz1
              !write(*,1011) ax, ay, az
1011          format('The box vector a =',f11.4, f11.4, f11.4)

              !write(*,1012) bx, by, bz
1012          format('The box vector b =',f11.4, f11.4, f11.4)

              !write(*,1013) cx, cy, cz
1013          format('The box vector c =',f11.4, f11.4, f11.4)

!           endif
           if ( istat .ne. 0 ) then
              exit READLINES
           endif
        end if
        ! Calculate diffusion before wrapping atoms into unit cell.
        k = 1 + j / distance_block
        do m = 1, k
           ! Loop over initial conditions.
           offset = distance_block * (m-1)
           if ( offset .lt. j ) then
              idxt = (j-offset)/distance_block
              do l = 1, nindex1 
                 ll = index1(l)
                 distance(idxt) = distance(idxt) + &
                      (x0(ll,m) - x(ll))**2 + (y0(ll,m)-y(ll))**2 &
                      + (z0(ll,m)-z(ll))**2
                 dcount(idxt) = dcount(idxt) + 1
              end do
           end if
        end do
        call inversebox(ax, ay, az, bx, by, bz, cx, cy, cz, invbox, boxvol)
        do k = 1, natom
           ! Wrap atoms into the primitive unit cell.
           call wrap_in_box(ax, ay, az, bx, by, bz, cx, cy, cz, invbox, &
                x(k), y(k), z(k))
        end do
        if ( time .lt. tmin ) then
           cycle
        else if ( time .gt. tmax ) then
           exit READLINES
        endif
        do l = 1, nindex1
           ll = index1(l) 
           do k = 1, nindex2
              kk = index2(k)
              do ia = -isuper, isuper, 1
                 do ib = -isuper, isuper, 1
                    do ic = -isuper, isuper, 1
                       if ( kk .eq. ll .and. ia .eq. 0 .and. ib .eq. 0 &
                            .and. ic .eq. 0 ) then
                          ! avoid identical atoms.
                          cycle
                       endif

                       rx = x(ll) - x(kk)
                       ry = y(ll) - y(kk)
                       rz = z(ll) - z(kk)

                       rx = rx + ia * ax + ib * bx + ic * cx
                       ry = ry + ia * ay + ib * by + ic * cy
                       rz = rz + ia * az + ib * bz + ic * cz 


                       r2 = rx * rx + ry * ry + rz * rz

                       r  = sqrt(r2)

                       i = int(r/dr) + 1

                       if ( i .lt. MAXHISTO ) then
                          histo(i) = histo(i) + 1
                          if ( i .gt. histomax ) then
                             histomax = i 
                          endif
                       else
                          write(*,*) 'The histo array needs to be made bigger'
                          stop
                       end if

                    end do
                 end do
              end do

           end do
        end do
        nstep = nstep + 1
     end do READLINES
     nc(:) = real(histo(:))/nstep

     write(*,* ) 'Number of steps read in = ', nstep

     write(*,*) 'Distance^2 for element 1' 
     write(*,*) ' '
     write(*,*) 'Time (ps)  Distance^2 (Ang.^2)'

     do i = 1, nstep/distance_block
        distance(i) = distance(i) / dcount(i)
        write(*,89) dtime * i * distance_block, distance(i) 
89      format(e13.5, ' ' , e13.5)
     end do
     write(*,*) ' '
     write(*,*) 'Diffusion constant estimates:'
     write(*,*) 'Time Diffusion const. (cm^2 / s)'
     do i = 2, 5
        k = (i * nstep) / (5 * distance_block) 
        k1 = ((i-1) * nstep) / (5 * distance_block)
        dconst = distance(k) - distance(k1)
        dconst = dconst / (dtime * distance_block * (k - k1) )
        dconst = dconst * 1.0e-04 / 6.0
        write(*,89), dtime * distance_block * k, dconst
     end do
     
     alen = sqrt(ax * ax + ay * ay + az * az) 
     blen = sqrt(bx * bx + by * by + bz * bz) 
     clen = sqrt(cx * cx + cy * cy + cz * cz) 

     rcut = alen
     if ( rcut .gt. blen ) then
        rcut = blen 
     end if
     if ( rcut .gt. clen ) then
        rcut = clen 
     end if

     nconfig = 0 
     do i = 1, histomax
        ru = i * dr
        if ( ru .gt. rcut ) then
           icut = i - 1
           exit
        endif
        nconfig = nconfig + histo(i)
     end do

     if ( nconfig .gt. 0 ) then
        pi = 4.0 * atan(1.0)
        ! debug ! put in the correct factor.
        vfac = boxvol

        do i = 1, histomax
           ru = i * dr
           rl = (i-1) * dr
           dv = 4.0 * pi * (ru * ru * ru - rl * rl * rl) / 3.0
           !          gasfac is the number of A-B pairs expected in dv if 
           !          the system was an ideal gas.
           if ( itype1 .ne. itype2 ) then
              gasfac = nindex1 * nindex2 * dv / boxvol
           else 
              gasfac = nindex1 * (nindex1 - 1) * dv / boxvol
           endif
           gr(i) = histo(i) / (nstep * gasfac)
           rbin(i) = (ru + rl) / 2.0d0
        end do

        write(*,*) ' '
        write(*,*) 'r     G(r)'
        open (file = 'gr.dat', status = 'replace', unit = 10)
        do i = 1, icut
           write(10,99) rbin(i), gr(i)
99         format(e13.5, ',' , e13.5)
        end do
        close (10)
        open (file = 'nc.dat', status = 'replace', unit = 10)
        nc_int = 0d0
        do i = 1, icut
           nc_int = nc_int + nc(i)/nindex1
           write(10,99) rbin(i), nc_int
        end do
        write(*,101) 
101     format(/)
        close (10)

        write(*,*) 'Bond definition second value where: g(r) = 1'
        do i = 1, histomax - 1
           if ( gr(i) .gt. 1.0 .and. gr(i+1) .lt. 1.0 ) then
              r1 = rbin(i) + (1.0d0-gr(i)) * (dr / (gr(i+1)-gr(i)))
              if ( r1 .lt. 2.2 ) then
                 write(*,103) r1
              else
                 write(*,*) 'No covalent bond found.'
              endif
              exit
           endif
        end do

        write(*,*) 'Bond definition first minimum of g(r) where g(r) < 1 '
        do i = 1, histomax - 2
           if ( gr(i) > gr(i+1) .and. gr(i+1) < gr(i+2) .and. gr(i+1) < 1.0 ) then
              rmin = rbin(i+1)
              if ( rmin .lt. 2.2 ) then
                 write(*,103) rmin
103              format('Bond distance =', f7.3)
              else 
                 write(*,*) 'No covalent bond found'
              endif
              exit
           endif
        end do

     endif
     close(in)
  end do TIMELOOP

  stop
end program gr_calc

  
subroutine inversebox(ax, ay, az, bx, by, bz, cx, cy, cz, invbox, boxvol)
  ! Calculate the inverse box maxtrix, used for finding
  !  crystal coordinates.
  implicit none
  double precision ax, ay, az ! Box vector a.
  double precision bx, by, bz ! Box vector b.
  double precision cx, cy, cz ! Box vector c.
  double precision invbox(3,3)
  double precision xhlp
  double precision boxvol

  xhlp=-az*by*cx+ ay*bz*cx 
  xhlp=xhlp+az*bx*cy 
  xhlp=xhlp-ax*bz*cy 
  xhlp=xhlp-ay*bx*cz 
  xhlp=xhlp+ax*by*cz 

  invbox(1,1) =(-bz*cy+by*cz)/xhlp 
  invbox(2,1) =(bz*cx-bx*cz)/xhlp 
  invbox(3,1) =(-by*cx+bx*cy)/xhlp 

  invbox(1,2)=(az*cy-ay*cz)/xhlp 
  invbox(2,2)=(-az*cx+ax*cz)/xhlp 
  invbox(3,2)=(ay*cx-ax*cy)/xhlp 

  invbox(1,3)=(-az*by+ay*bz)/xhlp 
  invbox(2,3)=(az*bx-ax*bz)/xhlp 
  invbox(3,3)=(-ay*bx+ax*by)/xhlp 

  boxvol = dabs(xhlp)
  return
end subroutine inversebox




subroutine wrap_in_box(ax, ay, az, bx, by, bz, cx, cy, cz, invbox, x, y, z)
  ! This function wraps the given point x, y, z into the primitive simulation
  ! box. 
  implicit none
  double precision ax, ay, az, bx, by, bz, cx, cy, cz
  double precision invbox(3,3)
  double precision x, y, z
  double precision ca, cb, cc 

!  write(*,*) 'Original x, y, z: ', x, y, z
  ca = invbox(1,1) * x + invbox(2,1) * y + invbox(3,1) * z 
  cb = invbox(1,2) * x + invbox(2,2) * y + invbox(3,2) * z 
  cc = invbox(1,3) * x + invbox(2,3) * y + invbox(3,3) * z 

  ca = ca - nint(ca) 
  cb = cb - nint(cb) 
  cc = cc - nint(cc) 

  x = ca * ax + cb * bx + cc * cx 
  y = ca * ay + cb * by + cc * cy 
  z = ca * az + cb * bz + cc * cz 

!  write(*,*) 'New x, y, z: ', x, y, z
  return
end subroutine wrap_in_box
