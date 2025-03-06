module Utils
   use globals
   implicit none

   public :: printChar, trimQuotes, ranNum, seedRandom, ranVec

contains

   subroutine printChar(c, len)
      integer, intent(in) :: len
      character(len=*), intent(in) :: c
      integer :: i
      do i = 1, len
         write (*, '(A)', advance='no') c
      end do
      print *
   end subroutine printChar

   subroutine trimQuotes(s)
      character(len=*), intent(inout) :: s
      if (s(1:1) == '"' .and. s(len_trim(s):len_trim(s)) == '"') then
         s = s(2:len_trim(s) - 1)
      end if
   end subroutine trimQuotes

   function ranNum()
      real(8) :: ranNum
      call random_number(ranNum)
   end function ranNum

   function rangeRanNum(min, max)
      real(8), intent(in) :: min, max
      real(8) :: rangeRanNum
      rangeRanNum = ranNum()
      rangeRanNum = rangeRanNum*(max - min) + min
   end function rangeRanNum

   subroutine seedRandom(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: seed_array(:)

      call random_seed(size=n)
      allocate (seed_array(n))
      ! Fill the seed array with the seed value
      do i = 1, n
         seed_array(i) = seed
      end do
      call random_seed(put=seed_array)
   end subroutine seedRandom

   subroutine createOutputDir(dirName)
      character(len=*), intent(in) :: dirName
      character(len=50) :: cmd
      logical :: dirExists, dirIsDir

      inquire (file=dirName, exist=dirExists)
      if (.not. dirExists) then
         call system('mkdir '//trim(adjustl(dirName)))
         call system('mkdir '//trim(adjustl(dirName))//'/coords')
         call system('mkdir '//trim(adjustl(dirName))//'/props')
      end if
      ! cmd = 'mkdir output/'//trim(adjustl(dirName))
      ! call system(cmd)
   end subroutine createOutputDir

   ! du(k) = du(k)*rangeRanNum(-p%lambda, p%lambda)
   subroutine ranVec(vec)
      real(8) :: ran,norm,sum
      real(8), dimension(3), intent(out) :: vec
      integer :: i

      norm=2
      do while (norm > 1.0)
         sum=0
         do i = 1,3
            ran = ranNum()
            vec(i) = 2*(ran - 0.5)
            sum=sum+vec(i)**2
         end do
         norm=sum
      end do
      norm=sqrt(norm)
      do i = 1,3
         vec(i) = vec(i)/norm
      end do

   end subroutine ranVec

   subroutine saveVTU(p, cfg, step)
      type(Particles), intent(in) :: p
      type(ConfigFile), intent(inout) :: cfg
      integer, intent(in) :: step
      character(len=50) :: stepChar, nChar, lChar, etaChar, fileName
      character(len=100) :: fullPath
      integer :: ioStatus, i

      ! Convert integer i to character string
      write (stepChar, '(I0)') step
      write (nChar, '(I0)') p%nParticles
      write (lChar, '(F0.2)') p%l
      write (etaChar, '(F0.2)') p%eta

      fileName = trim(cfg%fileName)//'_'//trim(stepChar)//'.vtu'
      fullPath = trim(cfg%dirName)//'/coords/'//trim(fileName)

      open (unit=10, file=fullPath, status='replace', action='write', iostat=ioStatus)
      if (ioStatus /= 0) then
         print *, 'Error opening file:', fullPath
         stop
      end if

      write (10, '(A)') '<?xml version="1.0"?>'
      write (10, '(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
      write (10, '(A)') '  <UnstructuredGrid>'
      write (10, '(A)') '    <Piece NumberOfPoints="'//trim(adjustl(nChar))//'" NumberOfCells="0">'
      write (10, '(A)') '      <Points>'
      write (10, '(A)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
      do i = 1, p%nParticles
         write (10, '(3F10.5, 1X)') p%r(i, 1), p%r(i, 2), p%r(i, 3)
      end do
      write (10, '(A)') '        </DataArray>'
      write (10, '(A)') '      </Points>'
      write (10, '(A)') '      <PointData Vectors="orientations">'
      write (10, '(A)') '        <DataArray type="Float64" Name="orientations" NumberOfComponents="3" format="ascii">'
      do i = 1, p%nParticles
         write (10, '(3F10.5)') p%u(i, 1), p%u(i, 2), p%u(i, 3)
      end do
      write (10, '(A)') '        </DataArray>'
      write (10, '(A)') '      </PointData>'
      ! write (10, '(A)') '      <FieldData>'
      ! write (10, '(A)') '        <DataArray type="Float64" Name="lBox" NumberOfComponents="1" format="ascii">'
      ! write (10, '(F10.5)') p%lBox
      ! write (10, '(A)') '        </DataArray>'
      ! write (10, '(A)') '        <DataArray type="Float64" Name="eta" NumberOfComponents="1" format="ascii">'
      ! write (10, '(F10.5)') p%eta
      ! write (10, '(A)') '        </DataArray>'
      ! write (10, '(A)') '      </FieldData>'
      write (10, '(A)') '    </Piece>'
      write (10, '(A)') '  </UnstructuredGrid>'
      write (10, '(A)') '</VTKFile>'

      close (10)
   end subroutine saveVTU


   subroutine readVTU(dir,fileName,p)
      type(Particles), intent(inout) :: p
      character(len=*), intent(in) :: fileName
      character(len=*), intent(in) :: dir
      character(len=50) :: file
      integer :: i, ioStatus, nParticles,pos1,pos2
      character(len=100) :: line, text

      file = trim(dir)//'/coords/'//trim(fileName)//'.vtu'
      print *, 'Reading props file: ', file
      open (unit=10, file=file, status='old', action='read', iostat=ioStatus)
      if (ioStatus /= 0) then
         print *, 'Error opening file:', file
         stop
      end if

      ! Skip the header lines
      do i = 1, 3
         read (10, '(A)') line
      end do

      read (10, '(A)') line

      pos1 = index(line, "NumberOfPoints=")
      pos1 = pos1+1 + len("NumberOfPoints=")
      text = line(pos1:len(line))
      pos2 = index(text, ' ')
      text = text(1:pos2-2)

      read (text, *) nParticles

      if (nParticles /= p%nParticles) then
         print *, 'Number of particles in file does not match the number of particles specified in the config, exiting.'
         stop
      end if

      allocate (p%r(nParticles, 3))
      allocate (p%u(nParticles, 3))

      ! skip to the start of point data
      do i = 1, 2
         read (10, '(A)') line
      end do

      ! read the particle positions
      do i = 1, nParticles
         read(10, '(A)') line
         line = trim(line)
         read (line, *) p%r(i, :)
      end do

      ! skip to the start of orientation data
      do i = 1, 4
         read (10, '(A)') line
      end do

      ! read the particle orientations
      do i = 1, nParticles
         read(10, '(A)') line
         line = trim(line)
         read (line, *) p%u(i, :)
      end do

      ! do i = 1,4
      !    read(10, '(A)') line
      ! end do

      ! read(10,'(A)') line
      ! read(line,*) p%lBox

      ! do i = 1,2
      !    read(10,'(A)') line
      ! end do

      ! read(10,'(A)') line
      ! read(line,*) p%eta

      close (10)

   end subroutine readVTU


   subroutine copyFile(src, dst)
      character(len=*), intent(in) :: src, dst
      character(len=100) :: cmd
      cmd = 'cp '//trim(adjustl(src))//' '//trim(adjustl(dst))
      call system(cmd)
   end subroutine copyFile

   subroutine saveXYZ(p,cfg,step)
      type(Particles), intent(in) :: p
      type(ConfigFile), intent(in) :: cfg
      integer, intent(in) :: step

      character(len=50) :: stepChar, nChar, fileName
      character(len=100) :: fullPath
      integer :: ioStatus, i

      ! Convert integer i to character string
      write (stepChar, '(I0)') step
      write (nChar, '(I0)') p%nParticles

      fileName = trim(cfg%fileName)//'_'//trim(stepChar)//'.xyz'
      fullPath = trim(cfg%dirName)//'/coords/'//trim(fileName)

      print *, 'Saving XYZ file: ', fileName
      print *, 'Output Directory: ', trim(cfg%dirName)
      open (unit=10, file=fullPath, status='replace', action='write', iostat=ioStatus)
      if (ioStatus /= 0) then
         print *, 'Error opening file:', fullPath
         stop
      end if

      write (10, '(A)') trim(adjustl(nChar))
      write (10, '(A)') 'Spherocylinder Position and Orientation'



      do i=1,p%nParticles
         write (10, '(A, 3F10.5, 3F10.5)') 'SC ', p%r(i, 1), p%r(i, 2), p%r(i, 3), p%u(i, 1), p%u(i, 2), p%u(i, 3)
      end do

      close(10)
   end subroutine saveXYZ


   subroutine readXYZ(dir,fileName,p)
      type(Particles), intent(inout) :: p
      character(len=*), intent(in) :: fileName,dir
      integer :: i, ioStatus, nParticles
      character(len=100) :: line
      character(len=50) :: var,file

      file = trim(dir)//'/coords/'//trim(fileName)//'.xyz'
      open (unit=10, file=file, status='old', action='read', iostat=ioStatus)
      if (ioStatus /= 0) then
         print *, 'Error opening file:', file
         stop
      end if

      read (10, *) nParticles

      if (nParticles /= p%nParticles) then
         print *, 'Number of particles in file does not match the number of particles specified in the config, exiting.'
         stop
      end if

      allocate (p%r(nParticles, 3))
      allocate (p%u(nParticles, 3))

      ! skip the comment line
      read(10,*) line
      print *, line

      do i = 1, nParticles
         read(10,'(A)') line
         read(line,*) var, p%r(i, 1), p%r(i, 2), p%r(i, 3), p%u(i, 1), p%u(i, 2), p%u(i, 3)
      end do

      close (10)

   end subroutine readXYZ

   subroutine saveProps(p,cfg,step)
      type(Particles), intent(in) :: p
      type(ConfigFile), intent(in) :: cfg
      integer, intent(in) :: step

      character(len=50) :: stepChar, nChar,etaChar,lBox ,fileName,drmaxChar,dvmaxChar,lambdaChar
      character(len=100) :: fullPath
      integer :: ioStatus

      ! Convert integer i to character string
      write (stepChar, '(I0)') step
      write (nChar, '(I0)') p%nParticles
      write(etaChar,'(F10.5)') p%eta
      write(lBox,'(F10.5)') p%lBox
      write(drmaxChar,'(F10.5)') p%drMax
      write(dvmaxChar,'(F10.5)') p%dvMax
      write(lambdaChar,'(F10.5)') p%lambda


      fileName = trim(cfg%fileName)//'_'//trim(stepChar)//'.dat'
      fullPath = trim(cfg%dirName)//'/props/'//trim(fileName)

      open (unit=10, file=fullPath, status='replace', action='write', iostat=ioStatus)
      if (ioStatus /= 0) then
         print *, 'Error opening file:', fullPath
         stop
      end if

      write (10,"(A,A)") "step=",trim(adjustl(stepChar))
      write (10, '(A, A)') 'lBox=', trim(adjustl(lBox))
      write (10, '(A, A)') 'eta=', trim(adjustl(etaChar))
      write (10, '(A, A)') 'nParticles=', trim(adjustl(nChar))
      write (10, '(A, A)') 'drMax=', trim(adjustl(drmaxChar))
      write (10, '(A, A)') 'dvMax=', trim(adjustl(dvmaxChar))
      write (10, '(A, A)') 'lambda=', trim(adjustl(lambdaChar))

   end subroutine saveProps


   subroutine readProps(dir,fileName,p)
      type(Particles), intent(inout) :: p
      character(len=*), intent(in) :: fileName,dir
      integer :: ioStatus,pos,i
      character(len=100) :: line
      character(len=50) :: var,value,file

      file = trim(dir)//'/props/'//trim(fileName)//'.dat'
      print *, 'Reading props file: ', file
      open (unit=10, file=file, status='old', action='read', iostat=ioStatus)
      if (ioStatus /= 0) then
         print *, 'Error opening file:', file
         stop
      end if

      do while (line /= "")
         read(10,'(A)',iostat=ioStatus) line
         if (ioStatus /= 0) exit
         if (line == '') exit ! check if at the end of the file
         pos = index(line, '=')
         if (pos > 0) then
            var = trim(adjustl(line(1:pos - 1)))
            value = trim(adjustl(line(pos + 1:)))
            print *, trim(var), trim(value)

            select case (var)
             case ('lBox')
               read(value,*) p%lBox
             case ('eta')
               read(value,*) p%eta
             case ('nParticles')
               read(value,*) p%nParticles
             case ('drMax')
               read(value,*) p%drMax
             case ('dvMax')
               read(value,*) p%dvMax
             case ('lambda')
               read(value,*) p%lambda
            end select
         end if
      end do
      ! do i = 1, 7
      !    read(10,'(A)') line
      !    if (line == '') exit ! check if at the end of the file
      !    pos = index(line, '=')
      !    if (pos > 0) then
      !       var = trim(adjustl(line(1:pos - 1)))
      !       value = trim(adjustl(line(pos + 1:)))
      !       print *, trim(var), trim(value)

      !       select case (var)
      !        case ('lBox')
      !          read(value,*) p%lBox
      !        case ('eta')
      !          read(value,*) p%eta
      !        case ('nParticles')
      !          read(value,*) p%nParticles
      !        case ('drMax')
      !          read(value,*) p%drMax
      !        case ('dvMax')
      !          read(value,*) p%dvMax
      !        case ('lambda')
      !          read(value,*) p%lambda
      !       end select
      !    end if
      ! end do
   end subroutine readProps


   subroutine saveState(p,cfg,step)
      type(Particles), intent(in) :: p
      type(ConfigFile), intent(inout) :: cfg
      integer, intent(in) :: step
      call saveVTU(p,cfg,step)
      call saveProps(p,cfg,step)
   end subroutine saveState


   subroutine readState(dir,fileName,p)
      type(Particles), intent(inout) :: p
      character(len=*), intent(in) :: dir,fileName

      call readVTU(dir,fileName,p)
      call readProps(dir,fileName,p)

   end subroutine readState

   subroutine saveCSV(p, rhoAvg, errorRho, cfg, step)
    type(Particles), intent(in) :: p
    type(ConfigFile), intent(in) :: cfg
    integer, intent(in) :: step
    real(8), intent(in) :: rhoAvg, errorRho
    real(8) :: errorEta,etaAvg
    character(len=50) :: fileName
    character(len=100) :: fullPath
    integer :: ioStatus, unit, fileSize

    etaAvg = rhoAvg * p%v0 
    errorEta = errorRho * p%v0

    fileName = 'output.csv'
    fullPath = trim(cfg%dirName)//'/'//trim(fileName)

    ! Open the file in append mode
    open (unit=10, file=fullPath, status='unknown', action='write', position='append', iostat=ioStatus)
    if (ioStatus /= 0) then
        print *, 'Error opening file:', fullPath
        stop
    end if

    ! Check if the file is empty and write the header if it is
    inquire(unit=10, size=fileSize)
    if (fileSize == 0) then
        write (10, '(A)') 'step,etaAvg,errorEta,rhoAvg,errorRho,lBox,drMax,dvMax,lambda'
    end if

    ! Write the data
    write (10, '(I0)', advance='no') step
    write (10, '(A)', advance='no') ','
    write (10, '(F10.5)', advance='no') etaAvg
    write (10, '(A)', advance='no') ','
    write (10, '(F10.5)', advance='no') errorEta
    write (10, '(A)', advance='no') ','
    write (10, '(F10.5)', advance='no') rhoAvg
    write (10, '(A)', advance='no') ','
    write (10, '(F10.5)', advance='no') errorRho
    write (10, '(A)', advance='no') ','
    write (10, '(F10.5)', advance='no') p%lBox
    write (10, '(A)', advance='no') ','
    write (10, '(F10.5)', advance='no') p%drMax
    write (10, '(A)', advance='no') ','
    write (10, '(F10.5)', advance='no') p%dvMax
    write (10, '(A)', advance='no') ','
    write (10, '(F10.5)', advance='no') p%lambda
    write (10, '(A)') ''

    close (10)
end subroutine saveCSV



   function calcEigen(p) result(eig)
      type(Particles), intent(in) :: p
      type(Eigen):: eig
      real(8) :: sq, sqSq
      real(8) :: qSum(3, 3), qRun(3, 3),frac
      integer :: ii, jj, kk, nrot



      qSum = 0
      do ii = 1, p%nParticles
         do jj = 1, 3
            do kk = 1, 3
               if (jj == kk) then
                  frac = 0.5
               else
                  frac = 0
               end if
               qSum(jj, kk) = qSum(jj, kk) + 1.5*p%u(ii, jj)*p%u(ii, kk) - frac
            end do
         end do
      end do
      qRun = qRun + qSum
      qSum = qSum/p%nParticles

      call jacobi(qSum, 3, 3, eig%seval, eig%sevec, nrot)
      call eig2srt(eig%seval, eig%sevec, 3, 3)

      sq = sq - 2*eig%seval(2)
      sqSq = sqSq + 4*eig%seval(2)**2
   end function calcEigen


!********************************************************************
!  ** numerical recipes routine to diagonalise a matrix            **
!********************************************************************
   SUBROUTINE jacobi(a, n, np, d, v, nrot)
      IMPLICIT NONE
      INTEGER n, np, nrot, NMAX
      real*8 a(np, np), d(np), v(np, np)
      PARAMETER(NMAX=500)
      INTEGER i, ip, iq, j
      real*8 c, g, h, s, sm, t, tau, theta, tresh, b(NMAX), z(NMAX)
      do 12 ip = 1, n
         do 11 iq = 1, n
            v(ip, iq) = 0.
11       continue
         v(ip, ip) = 1.
12    continue
      do 13 ip = 1, n
         b(ip) = a(ip, ip)
         d(ip) = b(ip)
         z(ip) = 0.
13    continue
      nrot = 0
      do 24 i = 1, 50
         sm = 0.
         do 15 ip = 1, n - 1
            do 14 iq = ip + 1, n
               sm = sm + abs(a(ip, iq))
14          continue
15       continue
         if (sm .eq. 0.) return
         if (i .lt. 4) then
            tresh = 0.2*sm/n**2
         else
            tresh = 0.
         end if
         do 22 ip = 1, n - 1
            do 21 iq = ip + 1, n
               g = 100.*abs(a(ip, iq))
               if ((i .gt. 4) .and. (abs(d(ip)) &
                  *g .eq. abs(d(ip))) .and. (abs(d(iq)) + g .eq. abs(d(iq)))) then
                  a(ip, iq) = 0.
               else if (abs(a(ip, iq)) .gt. tresh) then
                  h = d(iq) - d(ip)
                  if (abs(h) + g .eq. abs(h)) then
                     t = a(ip, iq)/h
                  else
                     theta = 0.5*h/a(ip, iq)
                     t = 1./(abs(theta) + sqrt(1.+theta**2))
                     if (theta .lt. 0.) t = -t
                  end if
                  c = 1./sqrt(1 + t**2)
                  s = t*c
                  tau = s/(1.+c)
                  h = t*a(ip, iq)
                  z(ip) = z(ip) - h
                  z(iq) = z(iq) + h
                  d(ip) = d(ip) - h
                  d(iq) = d(iq) + h
                  a(ip, iq) = 0.
                  do 16 j = 1, ip - 1
                     g = a(j, ip)
                     h = a(j, iq)
                     a(j, ip) = g - s*(h + g*tau)
                     a(j, iq) = h + s*(g - h*tau)
16                continue
                  do 17 j = ip + 1, iq - 1
                     g = a(ip, j)
                     h = a(j, iq)
                     a(ip, j) = g - s*(h + g*tau)
                     a(j, iq) = h + s*(g - h*tau)
17                continue
                  do 18 j = iq + 1, n
                     g = a(ip, j)
                     h = a(iq, j)
                     a(ip, j) = g - s*(h + g*tau)
                     a(iq, j) = h + s*(g - h*tau)
18                continue
                  do 19 j = 1, n
                     g = v(j, ip)
                     h = v(j, iq)
                     v(j, ip) = g - s*(h + g*tau)
                     v(j, iq) = h + s*(g - h*tau)
19                continue
                  nrot = nrot + 1
               end if
21          continue
22       continue
         do 23 ip = 1, n
            b(ip) = b(ip) + z(ip)
            d(ip) = b(ip)
            z(ip) = 0.
23       continue
24    continue
      stop 'too many iterations in jacobi'
!      pause 'too many iterations in jacobi'
!      return
   END subroutine jacobi
!  (C) Copr. 1986-92 Numerical Recipes Software Y.#?u25.){2p49.
!********************************************************************
!  ** numerical recipes routine to sort eigenvalues (smallest first)
!********************************************************************
   SUBROUTINE eig2srt(d, v, n, np)
      IMPLICIT NONE
!     as eigsrt but opposite order
      INTEGER n, np
      real*8 d(np), v(np, np)
      INTEGER i, j, k
      real*8 p
      do 13 i = 1, n - 1
         k = i
         p = d(i)
         do 11 j = i + 1, n
            if (d(j) .lt. p) then
               k = j
               p = d(j)
            end if
11       continue
         if (k .ne. i) then
            d(k) = d(i)
            d(i) = p
            do 12 j = 1, n
               p = v(j, i)
               v(j, i) = v(j, k)
               v(j, k) = p
12          continue
         end if
13    continue
      return
   END subroutine eig2srt
end module Utils
