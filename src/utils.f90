module Utils
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

      inquire (file="output", exist=dirExists)
      if (.not. dirExists) then
         call system('mkdir output')
      end if
      cmd = 'mkdir output/'//trim(adjustl(dirName))
      call system(cmd)
   end subroutine createOutputDir

   subroutine ranVec(vec)
      real(8) :: ran
      real(8), dimension(3), intent(out) :: vec
      integer :: i

      do i = 1, 3
         ran = ranNum()
         vec(i) = ran
      end do
      vec = vec/sqrt(sum(vec**2))
   end subroutine ranVec

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
11          continue
            v(ip, ip) = 1.
12          continue
            do 13 ip = 1, n
               b(ip) = a(ip, ip)
               d(ip) = b(ip)
               z(ip) = 0.
13             continue
               nrot = 0
               do 24 i = 1, 50
                  sm = 0.
                  do 15 ip = 1, n - 1
                     do 14 iq = ip + 1, n
                        sm = sm + abs(a(ip, iq))
14                      continue
15                      continue
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
16                                  continue
                                    do 17 j = ip + 1, iq - 1
                                       g = a(ip, j)
                                       h = a(j, iq)
                                       a(ip, j) = g - s*(h + g*tau)
                                       a(j, iq) = h + s*(g - h*tau)
17                                     continue
                                       do 18 j = iq + 1, n
                                          g = a(ip, j)
                                          h = a(iq, j)
                                          a(ip, j) = g - s*(h + g*tau)
                                          a(iq, j) = h + s*(g - h*tau)
18                                        continue
                                          do 19 j = 1, n
                                             g = v(j, ip)
                                             h = v(j, iq)
                                             v(j, ip) = g - s*(h + g*tau)
                                             v(j, iq) = h + s*(g - h*tau)
19                                           continue
                                             nrot = nrot + 1
                                             end if
21                                           continue
22                                           continue
                                             do 23 ip = 1, n
                                                b(ip) = b(ip) + z(ip)
                                                d(ip) = b(ip)
                                                z(ip) = 0.
23                                              continue
24                                              continue
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
11                                                       continue
                                                         if (k .ne. i) then
                                                            d(k) = d(i)
                                                            d(i) = p
                                                            do 12 j = 1, n
                                                               p = v(j, i)
                                                               v(j, i) = v(j, k)
                                                               v(j, k) = p
12                                                             continue
                                                               end if
13                                                             continue
                                                               return
                                                               END subroutine eig2srt
                                                               end module Utils
