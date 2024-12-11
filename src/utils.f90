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
      real :: ranNum
      call random_number(ranNum)
   end function ranNum

   function rangeRanNum(min, max)
      real, intent(in) :: min, max
      real :: rangeRanNum
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
      real :: ran
      real, dimension(3), intent(out) :: vec
      integer :: i

      do i = 1, 3
         ran = ranNum()
         vec(i) = ran
      end do
      vec = vec/sqrt(sum(vec**2))
   end subroutine ranVec

end module Utils
