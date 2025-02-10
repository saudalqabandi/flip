program flip
   use Config
   use Utils
   use Particle
   use Potential
   use Overlap
   use Accumalators
   use Moves
   use globals
   implicit none

   real(8) :: start, finish, elapsed
   character(len=8) :: date
   character(len=10) :: time
   real(8) :: pot, potNew, potOld, ran, delta, vNew, rhoNew, lBoxNew, scale, third, dH, blockEnergy, blockVolume, runEnergy, runVolume, tmpEnergy, tmpVolume, blockDensity, totalEnergy, totalDensity, totalDensitySq,rhoAvg, rhoAvgSq, errorRho
   real(8), allocatable :: rNew(:, :), uNew(:, :), riNew(:), uiNew(:)
   real(8) :: rij(3)
   integer :: k, n, i, blockCount, j
   logical :: accept,dirExists
   character(len=25) :: moveType
   integer :: ii, jj, kk, count, nrot
   real(8) :: frac, seval(3), sevec(3, 3), qSum(3, 3), qRun(3, 3), sq, sqSq

   type(ConfigFile) :: cfg
   type(Particles) :: p
   type(Accumalator) :: acc
   type(Eigen) :: eig

   start = 0.0
   finish = 0.0
   elapsed = 0.0
   date = ''
   time = ''
   third = 1.0/3.0

   totalDensity = 0
   totalDensitySq = 0
   blockCount = 0
   blockEnergy = 0
   blockVolume = 0
   runEnergy = 0
   runVolume = 0
   totalEnergy = 0


   call cpu_time(start)
   call readConfig('config.cfg', cfg)

   ! seed the random number generator if requested
   if (cfg%seeded) then
      call seedRandom(cfg%seed)
   end if

   ! initialize the particles
   p = initParticles(cfg)

   if (cfg%setup == 'restart') then
      inquire(file=trim(cfg%dirName)//'/config.txt', exist=dirExists)
      if (dirExists) then
         cfg%dirName = trim(cfg%dirName)//"_restart"
         print *, 'Restart directory: ', cfg%dirName
      end if
   end if

   call createOutputDir(cfg%dirName)
   ! call saveVTU(p, cfg, 0)
   call saveState(p, cfg, 0)
   call carlosCheckOverlap(p%r, p%u, p, p%lBox)
   ! call checkOverlap(p%r,p%u,p)

   if (p%over) then
      print *, 'Overlap detected in initial configuraiton. Exiting.'
      stop
   end if

   ! p%potential = calcPotential(p, p%r, p%u)

   call printChar('=', 40)
   print *, 'Initial potential: ', p%potential
   call printChar('=', 40)

   call zeroAccumalators(acc)

   do n = 1, cfg%nCycles
      do k = 1, p%nParticles
         i = int(ranNum()*p%nParticles) + 1

         ran = ranNum()
         rNew = p%r
         uNew = p%u
         riNew = rNew(i, :)
         uiNew = uNew(i, :)

         if (ran < 1.0/2.0) then
            call translateMove(p, i, riNew)
            moveType = 'translate'
         else
            call rotateMove(p, i, uiNew)
            moveType = 'rotate'
            ! else
            !    call flipMove(p, i, uiNew)
            !    moveType = 'flip'
         end if

         rNew(i, :) = riNew
         uNew(i, :) = uiNew

         ! call singleParticleOverlap(rNew, uNew, p, i)
         call carlosSingleParticleOverlap(rNew, uNew, p, i,p%lBox)

         if (p%over) then
            accept = .false.
            call updateAccumalators(acc, moveType, accept)
            cycle
         else
            accept = .true.
            call updateAccumalators(acc, moveType, accept)
         end if

         if (accept) then
            p%r = rNew
            p%u = uNew
            p%potential = p%potential + delta
         end if
      end do

      ! vNew = p%vOld+(rangeRanNum(-p%dvMax, p%dvMax)*p%vOld)
      vNew = p%vOld+(rangeRanNum(-p%dvMax, p%dvMax))
      rhoNew = p%nParticles/vNew
      lBoxNew = vNew**third
      scale = lBoxNew/p%lBox
      rNew = p%r*scale
      uNew = p%u
      moveType = 'vol'

      call carlosCheckOverlap(rNew, uNew, p,lBoxNew)
      ! skip if overlap detected
      if (p%over) then
         accept = .false.
         call updateAccumalators(acc, moveType, accept)
      else
         delta = 0;
         dH = delta + p%pressure*(vNew - p%vOld) - p%nParticles*log(vNew/p%vOld)
         call metropolis(dH, p, accept)
         call updateAccumalators(acc, moveType, accept)
      end if

      if (accept) then
         p%r = rNew
         p%u = uNew
         p%potential = potNew
         p%vOld = vNew
         p%rho = rhoNew
         p%lBox = lBoxNew
         p%eta = rhoNew * p%v0
      end if

      if (mod(n, cfg%nDump) == 0) then
         print '(A, I10,A,F12.2)', 'Cycle: ', n, ' Potential: ', p%potential
         call printAccumalators(acc)

         call saveState(p, cfg, n)

         call printChar('-', 40)
         print *, 'drMax: ', p%drMax
         print *, 'dvMax: ', p%dvMax
         print *, 'lambda: ', p%lambda
         print *, "lBox: ", p%lBox
         print *, "vOld: ", p%vOld
         print *, "rho: ", p%rho
         print *, "eta: ", p%eta
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

      blockVolume = blockVolume + p%vold
      runVolume = runVolume + p%vold

      if (mod(n, cfg%nBlock) == 0) then
         blockCount = blockCount + 1

         blockVolume = blockVolume/(cfg%nBlock)
         tmpVolume = runVolume/n
         blockDensity = p%nParticles/blockVolume
         totalDensity = totalDensity + blockDensity
         totalDensitySq = totalDensitySq + blockDensity**2

         write (20, '(i7, 1x, 2(f12.3, 1x))') blockCount, blockEnergy, tmpEnergy
         write (*, *) 'Block average volume = ', blockVolume
         write (*, *) 'Block average density = ', blockDensity
         write (*, *) 'Run average volume = ', tmpVolume
         write (*, *) 'Run average density = ', p%nParticles/tmpVolume
         write (20, '(i7, 1x, 2(f12.3, 1x))') blockCount, blockVolume, tmpVolume
         call zeroAccumalators(acc)
      end if

      if (mod(n, cfg%nOrder) == 0) then
         count = count + 1

         eig = calcEigen(p)

         write (*, *) ' SORTED OVERALL EIGENVECTORS FOR Q'
         do i = 1, 3
            write (*, '(1X,A,I3,A,F12.6)') 'EIGENVALUE', i, '=', eig%seval(i)
            write (*, *) ' EIGENVECTOR'
            write (*, '(10X,5F12.6)') (eig%sevec(j, i), j=1, 3)
         end do
         print *, 'Nematic order parameter configuration', -2*eig%seval(2)
      end if

   end do

   rhoAvg = totalDensity/blockCount
   rhoAvgSq = totalDensitySq/blockCount
   errorRho = sqrt((rhoAvgSq - rhoAvg**2)/real(blockCount - 1))

   call printChar('=', 40)
   print *, 'Final potential: ', p%potential
   call printChar('=', 40)

   call saveState(p, cfg, cfg%nCycles)

   call printChar('-', 40)
   print *, 'drMax: ', p%drMax
   print *, 'dvMax: ', p%dvMax
   print *, 'lambda: ', p%lambda

   call printChar('-', 80)
   print *, 'Average density: ', rhoAvg, ' +/- ', errorRho
   print *, "Average eta: ", rhoAvg*p%v0 , ' +/- ', errorRho*p%v0
   call printChar('-', 80)

   print *, "Final volume: ", p%vOld
   print *, "Compressibility, Z: ", p%pressure/rhoAvg


   call printChar('-', 80)
   print *, "Final Eigenvalues"
   eig = calcEigen(p)

   write (*, *) ' SORTED OVERALL EIGENVECTORS FOR Q'
   do i = 1, 3
      write (*, '(1X,A,I3,A,F12.6)') 'EIGENVALUE', i, '=', eig%seval(i)
      write (*, *) ' EIGENVECTOR'
      write (*, '(10X,5F12.6)') (eig%sevec(j, i), j=1, 3)
   end do
   print *, 'Nematic order parameter configuration', -2*eig%seval(2)

   call printChar('-',80)

   ! call copyFile('log.out', trim(cfg%dirName)//'/log.out')
   ! call copyFile('config.cfg',trim(cfg%dirName)//'/config.txt')

   call cpu_time(finish)
   elapsed = finish - start
   print *, 'Elapsed time: ', elapsed, ' seconds'

end program flip
