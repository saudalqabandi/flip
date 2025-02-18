module Particle
   use Config
   use Utils
   use globals

   implicit none
contains

   function initParticles(cfg) result(p)
      type(ConfigFile), intent(in) :: cfg
      type(Particles) :: p

      real(8) :: pi, third, lBoxMin

      pi = 4.0*atan(1.0)
      third = 1.0/3.0
      lBoxMin = 2*(cfg%l + 1)

      ! assign values from config struct
      p%nParticles = cfg%nParticles
      p%l = cfg%l
      p%pressure = cfg%pressure
      p%eta = cfg%eta
      p%w0 = cfg%w0
      p%kappa = cfg%kappa
      p%rCut = cfg%rCut
      p%beta = cfg%beta
      p%drMax = cfg%drMax
      p%lambda = cfg%lambda
      p%setup = cfg%setup
      p%q = cfg%q

      ! calculate derived values
      p%v0 = (pi/6.0) + 0.25*pi*p%l
      p%pressure = p%pressure/p%v0
      p%rho = p%eta/p%v0
      p%vOld = p%nparticles/p%rho
      p%lBox = (p%vOld)**third

      p%dvMax = cfg%dvMax * p%vOld

      p%over = .false.

      ! setup the particles
      print *, 'Setting up particles...'
      call setupParticles(p,cfg)

   end function initParticles

   subroutine setupParticles(p,cfg)
      type(Particles), intent(inout) :: p
      type(ConfigFile), intent(in) :: cfg
      if (p%setup == 'random') then
         print *, 'Setting up random lattice'
         call setupRandom(p)
      else if (p%setup == 'cubic') then
         print *, 'Setting up cubic lattice'
         call setupCubic(p)

      else if (p%setup =='restart') then
         print *, "Setting up restart simulation"
         call setupRestart(p,cfg)
      else
         print *, 'Invalid setup type: ', p%setup
         stop
      end if
   end subroutine setupParticles

   subroutine setupRandom(p)
      type(Particles), intent(inout) :: p

      integer :: i, j
      real(8), dimension(3) :: vec

      ! allocate arrays
      allocate (p%r(p%nParticles, 3))
      allocate (p%u(p%nParticles, 3))

      do i = 1, p%nParticles
         do j = 1, 3
            p%r(i, j) = p%lBox*ranNum()
         end do
         call ranVec(vec)
         p%u(i, :) = vec
      end do
   end subroutine setupRandom

   subroutine setupCubic(p)
      real(8) :: third, step, halfBox, x, y, z
      integer :: i, j, k, n, index
      real(8), dimension(3) :: vec

      type(Particles), intent(inout) :: p

      third = 1.0/3.0
      n = anint(p%nParticles**third)

      if (n**3 /= p%nParticles) then
         print *, "Number of particles is not a perfect cube, adjusting from ", p%nParticles, " to ", n**3
         p%nParticles = n**3
         p%vOld = p%nParticles/p%rho
         p%lBox = (p%vOld)**third
      end if

      ! allocate arrays
      allocate (p%r(p%nParticles, 3))
      allocate (p%u(p%nParticles, 3))

      step = p%lBox/n
      ! step = p%l + 1.1
      halfBox = p%lBox/2.0

      index = 1
      do i = 1, n
         do j = 1, n
            do k = 1, n
               x = -halfBox + (i - 0.5)*step
               y = -halfBox + (j - 0.5)*step
               z = -halfBox + (k - 0.5)*step

               p%r(index, :) = (/x, y, z/)
               index = index + 1
            end do
         end do
      end do

      do i = 1, p%nParticles
         call ranVec(vec)
         p%u(i, :) = vec
      end do
   end subroutine setupCubic

   subroutine setupRestart(p,cfg)
      type(Particles), intent(inout) :: p
      type(ConfigFile), intent(in) :: cfg
      character(len=50) :: file,dir


      ! file = trim(cfg%dirName)//'/coords/'//trim(cfg%fileName)//'_'//trim(adjustl(cfg%restartStep))//'.xyz'
      dir = trim(cfg%restartDir)
      file = trim(cfg%fileName)//'_'//trim(adjustl(cfg%restartStep))
      print *, 'Reading file: ', file
      ! call readVTU(file,p)

      ! print *, 'Reading state from: ', dir, file
      call readState(dir,file,p)

      p%vOld = p%lBox**3
      p%rho = p%eta/p%v0

      ! recalculate parameters
      ! p%rho = p%eta/p%v0
      ! p%vOld = p%nparticles/p%rho
      ! p%lBox = (p%vOld)**(1.0/3.0)
   end subroutine setupRestart
end module Particle
