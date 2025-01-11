module globals

   type :: Particles
      integer :: nParticles, nCells
      real(8):: lBox, l, pressure, eta, w0, kappa, rCut, beta, drMax, dvMax, lambda, v0, rho, vOld, potential
      logical :: over
      real(8), allocatable :: r(:, :), u(:, :), q(:)
      character(len=50) :: setup
   end type Particles


   type :: ConfigFile
      integer :: nParticles, seed, nCycles, nDump, nBlock, nOrder, nAdjust
      real(8) :: l, pressure, eta, w0, kappa, rCut, beta, drMax, dvMax, lambda
      character(len=50) :: setup, dirName, fileName,restartStep
      real(8), allocatable :: q(:)
      logical :: seeded
   end type ConfigFile

   type :: Accumalator
      integer :: nMoves, nTrans, nRot, nFlip, nVol, nAccept, nTransAccept, nRotAccept, nFlipAccept, nVolAccept
      real(8) :: ratioTrans, ratioRot, ratioFlip, ratioVol, ratio
   end type Accumalator

end module globals
