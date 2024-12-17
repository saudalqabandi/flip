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
   real :: pot, potNew, potOld, ran, delta, vNew, rhoNew, lBoxNew, scale, third, dH, blockEnergy, blockVolume, runEnergy, runVolume, tmpEnergy, tmpVolume, blockDensity, totalEnergy, totalDensity, totalDensitySq,rhoAvg, rhoAvgSq, errorRho, rhoStar,rhoCp
   real, allocatable :: rNew(:, :), uNew(:, :), riNew(:), uiNew(:)
   real :: rij(3)
   integer :: k, n, i, blockCount, j
   logical :: accept
   character(len=25) :: moveType
   integer :: ii, jj, kk, count, nrot
   real(8) :: frac, seval(3), sevec(3, 3), qSum(3, 3), qRun(3, 3), sq, sqSq

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
   ! call checkOverlap(p%r, p%u, p)
   call carlosCheckOverlap(p,p%r,p%u)

   if (p%over) then
      print *, 'Overlap detected in initial configuraiton. Exiting.'
      stop
   end if

   ! p%potential = calcPotential(p, p%r, p%u)

   call printChar('=', 40)
   print *, 'Initial potential: ', p%potential
   call printChar('=', 40)

   call initAccumalators(acc)

   rhoCp = 2*(sqrt(2.0) + p%l*sqrt(3.0))

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

         ! callcheckOverlap(rNew, uNew, p)
         ! call singleParticleOverlap(p, rNew, uNew, i)
         call carlosSingleOverlap(p, rNew, uNew, i)

         ! if (p%over) then
         !    accept = .false.
         !    call updateAccumalators(acc, moveType, accept)
         !    cycle
         ! end if

         ! potOld = singleParticlePotential(p, p%r, p%u, i)
         ! potNew = singleParticlePotential(p, rNew, uNew, i)

         if (p%over) then
            accept = .false.
            call updateAccumalators(acc, moveType, accept)
            cycle
         else
            accept = .true.
            call updateAccumalators(acc, moveType, accept)
         end if

         !    delta = potNew - potOld
         !    call metropolis(delta, p, accept)
         ! call updateAccumalators(acc, moveType, accept)

         if (accept) then
            p%r = rNew
            p%u = uNew
            p%potential = p%potential + delta
         end if
      end do

      if (n <= 500) then
         vNew = p%vOld + p%dvMax*(2.0*ranNum() - 1.0)
      else
         vNew = p%vOld - (p%dvMax*ranNum())
      end if

      rhoNew = p%nParticles/vNew
      lBoxNew = vNew**third
      scale = lBoxNew/p%lBox
      rNew = p%r*scale
      uNew = p%u
      moveType = 'vol'

      ! call checkOverlap(rNew, uNew, p)
      call carlosCheckOverlap(p, rNew, uNew)
      ! skip if overlap detected
      if (p%over) then
         accept = .false.
         call updateAccumalators(acc, moveType, accept)
         cycle
      else
         accept = .true.
         call updateAccumalators(acc, moveType, accept)
      end if

      ! potOld = p%potential
      ! potNew = calcPotential(p, rNew, uNew)
      ! delta = potNew - potOld

      ! dH = delta + p%pressure*(vNew - p%vOld) - p%nParticles*log(vNew/p%vOld)

      ! print *, "potNew: ", potNew

      ! print *, "potential: ", p%potential

      ! call metropolis(dH, p, accept)
      ! call updateAccumalators(acc, moveType, accept)

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

         rhoStar = p%rho/rhoCp

         call printChar('-', 40)
         print *, 'drMax: ', p%drMax
         print *, 'dvMax: ', p%dvMax
         print *, 'lambda: ', p%lambda
         print *, "lBox: ", p%lBox
         print *, "vOld: ", p%vOld
         print *, "rho: ", p%rho
         print *, 'rhoStar: ', rhoStar

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

      ! blockEnergy = blockEnergy + p%potential
      blockVolume = blockVolume + p%vold

      ! runEnergy = runEnergy + p%potential
      runVolume = runVolume + p%vold

      if (mod(n, cfg%nBlock) == 0) then
         blockCount = blockCount + 1
         ! blockEnergy = blockEnergy/(cfg%nBlock*p%nParticles)
         ! totalEnergy = totalEnergy + blockEnergy
         ! tmpEnergy = runEnergy/(n*p%nParticles)

         blockVolume = blockVolume/(cfg%nBlock)
         tmpVolume = runVolume/n
         blockDensity = p%nParticles/blockVolume
         totalDensity = totalDensity + blockDensity
         totalDensitySq = totalDensitySq + blockDensity**2

         write (*, *) 'Cycle number = ', n
         ! write (*, *) 'Block average beta*energy/n = ', blockEnergy
         ! write (*, *) 'Run average beta*energy/n = ', tmpEnergy
         write (20, '(i7, 1x, 2(f12.3, 1x))') blockCount, blockEnergy, tmpEnergy
         write (*, *) 'Block average volume = ', blockVolume
         write (*, *) 'Block average density = ', blockDensity
         write (*, *) 'Run average volume = ', tmpVolume
         write (*, *) 'Run average density = ', p%nParticles/tmpVolume
         write (20, '(i7, 1x, 2(f12.3, 1x))') blockCount, blockVolume, tmpVolume

      end if

      if (mod(n, cfg%nOrder) == 0) then
         count = count + 1

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

         call jacobi(qSum, 3, 3, seval, sevec, nrot)
         call eig2srt(seval, sevec, 3, 3)

         sq = sq - 2*seval(2)
         sqSq = sqSq + 4*seval(2)**2
         write (*, *) ' SORTED OVERALL EIGENVECTORS FOR Q'
         do i = 1, 3
            write (*, '(1X,A,I3,A,F12.6)') 'EIGENVALUE', i, '=', SEVAL(I)
            write (*, *) ' EIGENVECTOR'
            write (*, '(10X,5F12.6)') (SEVEC(j, i), j=1, 3)
         end do
         print *, 'Nematic order parameter configuration', -2*seval(2)
      end if

   end do

   rhoAvg = totalDensity/blockCount
   rhoAvgSq = totalDensitySq/blockCount
   errorRho = sqrt((rhoAvgSq - rhoAvg**2)/real(blockCount - 1))

   call printChar('=', 40)
   print *, 'Final potential: ', p%potential
   call printChar('=', 40)

   call printAccumalators(acc)
   call saveVTK(p, cfg, cfg%nCycles)

   print *, 'rhoAvg: ', rhoAvg
   print *, 'rhoAvgSq: ', rhoAvgSq
   print *, 'errorRho: ', errorRho

   print *, 'Average density: ', rhoAvg, ' +/- ', errorRho
   call cpu_time(finish)
   elapsed = finish - start
   print *, 'Elapsed time: ', elapsed, ' seconds'

end program flip
