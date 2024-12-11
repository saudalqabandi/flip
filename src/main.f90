program flip
   use Config
   use Utils
   use Particle
   use Potential
   use Accumalators
   use Moves
   implicit none

   real :: start, finish, elapsed
   character(len=8) :: date
   character(len=10) :: time
   real :: pot, potNew, potOld, ran, delta, vNew, rhoNew, lBoxNew, scale, third, dH
   real, allocatable :: rNew(:, :), uNew(:, :), riNew(:), uiNew(:)
   real :: rij(3)
   integer :: k, n, i, d
   logical :: accept
   character(len=25) :: moveType

   type(ConfigFile) :: cfg
   type(Particles) :: p
   type(Accumalator) :: acc

   start = 0.0
   finish = 0.0
   elapsed = 0.0
   date = ''
   time = ''
   third = 1.0/3.0

   call cpu_time(start)
   call readConfig('config.cfg', cfg)

   ! seed the random number generator if requested
   if (cfg%seeded) then
      call seedRandom(cfg%seed)
   end if

   ! add date and time to output directory name
   call date_and_time(date, time)
   cfg%dirName = trim(cfg%dirName)//'_'//trim(date)//'_'//trim(time(1:6))
   call createOutputDir(cfg%dirName)

   ! initialize the particles
   p = initParticles(cfg)
   call saveVTK(p, cfg, 0)
   call checkOverlap(p%r, p%u, p)

   if (p%over) then
      print *, 'Overlap detected in initial configuraiton. Exiting.'
      stop
   end if

   p%potential = calcPotential(p, p%r, p%u)

   call printChar('=', 40)
   print *, 'Initial potential: ', p%potential
   call printChar('=', 40)

   call initAccumalators(acc)
   do n = 1, cfg%nCycles

      do k = 1, p%nParticles
         i = int(ranNum()*p%nParticles) + 1

         ran = ranNum()
         rNew = p%r
         uNew = p%u
         riNew = rNew(i, :)
         uiNew = uNew(i, :)

         if (ran < 1.0/3.0) then
            call translateMove(p, i, riNew)
            moveType = 'translate'
         else if (ran < 2.0/3.0) then
            call rotateMove(p, i, uiNew)
            moveType = 'rotate'
         else
            call flipMove(p, i, uiNew)
            moveType = 'flip'
         end if

         rNew(i, :) = riNew
         uNew(i, :) = uiNew

         call singleParticleOverlap(p, rNew, uNew, i)

         if (p%over) then
            accept = .false.
            call updateAccumalators(acc, moveType, accept)
            cycle
         end if

         potOld = singleParticlePotential(p, p%r, p%u, i)
         potNew = singleParticlePotential(p, rNew, uNew, i)

         if (p%over) then
            accept = .false.
            call updateAccumalators(acc, moveType, accept)
            cycle
         end if

         delta = potNew - potOld
         call metropolis(delta, p, accept)
         call updateAccumalators(acc, moveType, accept)

         if (accept) then
            p%r = rNew
            p%u = uNew
            p%potential = p%potential + delta
         end if
      end do

      vNew = p%vOld + p%dvMax*(2.0*ranNum() - 1.0)
      rhoNew = p%nParticles/vNew
      lBoxNew = vNew**third
      scale = lBoxNew/p%lBox
      rNew = p%r*scale
      uNew = p%u

      call checkOverlap(rNew, uNew, p)
      ! skip if overlap detected
      if (p%over) then
         accept = .false.
         call updateAccumalators(acc, moveType, accept)
         cycle
      end if

      potOld = p%potential
      potNew = calcPotential(p, rNew, uNew)
      delta = potNew - potOld

      dH = delta + p%pressure*(vNew - p%vOld) - p%nParticles*log(vNew/p%vOld)

      ! print *, "potNew: ", potNew

      ! print *, "potential: ", p%potential

      moveType = 'vol'
      call metropolis(dH, p, accept)
      call updateAccumalators(acc, moveType, accept)

      ! print *, 'Cycle: ', n, ' Potential: ', p%potential, ' potNew: ', potNew

      if (accept) then
         ! print *, "Accepted volume change"
         p%r = rNew
         p%u = uNew
         p%potential = potNew
         p%vOld = vNew
         p%rho = rhoNew
         p%lBox = lBoxNew
      end if

      if (mod(n, cfg%nDump) == 0) then
         print '(A, I10,A,F12.2)', 'Cycle: ', n, ' Potential: ', p%potential
         call printAccumalators(acc)
         call saveVTK(p, cfg, n)

         call printChar('-', 40)
         print *, 'drMax: ', p%drMax
         print *, 'dvMax: ', p%dvMax
         print *, 'lambda: ', p%lambda
         print *, "lBox: ", p%lBox
         print *, "vOld: ", p%vOld
         print *, "rho: ", p%rho
         call printChar('-', 40)
      end if

      if (mod(n, cfg%nAdjust) == 0) then

         if (acc%ratioTrans > 0.5) then
            p%drMax = p%drMax*1.1
         else
            p%drMax = p%drMax*0.9
         end if

         if (acc%ratioRot > 0.5) then
            p%lambda = p%lambda*1.1
         else
            p%lambda = p%lambda*0.9
         end if

         if (acc%ratioVol > 0.5) then
            p%dvMax = p%dvMax*1.1
         else
            p%dvMax = p%dvMax*0.9
         end if

      end if

   end do

   call cpu_time(finish)
   elapsed = finish - start
   print *, 'Elapsed time: ', elapsed, ' seconds'

end program flip
