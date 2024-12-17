module Config
   use Utils
   implicit none
   public :: ConfigFile, readConfig

   type :: ConfigFile
      integer :: nParticles, seed, nCycles, nDump, nBlock, nOrder, nAdjust
      real(8) :: l, pressure, eta, w0, kappa, rCut, beta, drMax, dvMax, lambda
      character(len=50) :: setup, dirName, fileName
      real(8), allocatable :: q(:)
      logical :: seeded
   end type ConfigFile

contains

   subroutine readConfig(filename, cfg)
      character(len=*), intent(in) :: filename
      type(ConfigFile), intent(out) :: cfg
      integer :: unit, ios, pos, numElements, i
      character(len=100) :: line
      character(len=50) :: var
      character(len=50) :: value

      ! open config file
      open (newunit=unit, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         print *, "Error opening file: ", filename
         stop
      end if

      !   loop through the file
      do
         read (unit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         pos = index(line, '=')
         if (pos > 0) then

            var = trim(adjustl(line(1:pos - 1)))
            value = trim(adjustl(line(pos + 1:)))

            select case (var)
            case ("nParticles")
               read (value, *) cfg%nParticles
            case ("seed")
               read (value, *) cfg%seed
            case ("nCycles")
               read (value, *) cfg%nCycles
            case ("nDump")
               read (value, *) cfg%nDump
            case ("nBlock")
               read (value, *) cfg%nBlock
            case ("nOrder")
               read (value, *) cfg%nOrder
            case ("l")
               read (value, *) cfg%l
            case ("pressure")
               read (value, *) cfg%pressure
            case ("eta")
               read (value, *) cfg%eta
            case ("w0")
               read (value, *) cfg%w0
            case ("kappa")
               read (value, *) cfg%kappa
            case ("rCut")
               read (value, *) cfg%rCut
            case ("beta")
               read (value, *) cfg%beta
            case ("drMax")
               read (value, *) cfg%drMax
            case ("dvMax")
               read (value, *) cfg%dvMax
            case ("lambda")
               read (value, *) cfg%lambda
            case ("setup")
               call trimQuotes(value)
               cfg%setup = value
            case ("dirName")
               call trimQuotes(value)
               cfg%dirName = trim(value)
            case ("fileName")
               call trimQuotes(value)
               cfg%fileName = value
            case ("seeded")
               read (value, *) cfg%seeded
            case ("nAdjust")
               read (value, *) cfg%nAdjust
            case ("q")
               value = trim(adjustl(value))
               numElements = count([(value(i:i) == ',', i=1, len(value))]) + 1
               allocate (cfg%q(numElements))
               read (value, *) cfg%q

               if (sum(cfg%q) /= 0) then
                  print *, "Error: q does not sum to 1"
                  stop
               end if
            end select

         end if

      end do

      close (unit)
   end subroutine readConfig

end module Config
