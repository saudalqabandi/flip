module Moves
   use globals
   use Utils
   implicit none

contains

   subroutine translateMove(p, i, riNew)
      type(Particles), intent(inout) :: p
      integer, intent(in) :: i
      real(8), intent(out) :: riNew(3)
      integer :: d

      do d = 1, 3
         riNew(d) = p%r(i, d) + rangeRanNum(-p%drMax, p%drMax)
         riNew(d) = riNew(d) - anint(riNew(d)/p%lBox)*p%lBox
      end do
   end subroutine translateMove

   subroutine rotateMove(p, i, uiNew)
      type(Particles), intent(inout) :: p
      integer, intent(in) :: i
      real(8), intent(inout) :: uiNew(3)
      real(8) :: du(3)
      integer :: k

      call ranVec(du)

      do k = 1, 3
         ! du(k) = du(k)*rangeRanNum(-p%lambda, p%lambda)
         du(k) = du(k)*p%lambda
      end do

      uiNew = p%u(i, :) + du
      uiNew = uiNew/sqrt(sum(uiNew**2))
   end subroutine rotateMove

   subroutine flipMove(p, i, uiNew)
      type(Particles), intent(inout) :: p
      integer, intent(in) :: i
      real(8), intent(inout) :: uiNew(3)
      uiNew = -p%u(i, :)
   end subroutine flipMove


   subroutine volumeMove(p,rNew)
      type(Particles), intent(in) :: p
      real(8), intent(inout) :: rNew(3,p%nParticles)
      real(8) :: lBoxNew, vNew, dV,third,rhoNew, scale

      third = 1.0/3.0
      dV = rangeRanNum(-p%dvMax, p%dvMax) * p%vOld
      vNew = p%vOld + dV
      rhoNew = p%nParticles/vNew
      lBoxNew = vNew**third
      scale = lBoxNew/p%lBox
      rNew = rNew*scale
      ! rNew = p%r*scale
   end subroutine volumeMove

end module Moves
