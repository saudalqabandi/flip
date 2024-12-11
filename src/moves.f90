module Moves
   use Particle
   use Utils

   implicit none

contains

   subroutine translateMove(p, i, riNew)
      type(Particles), intent(inout) :: p
      integer, intent(in) :: i
      real, intent(out) :: riNew(3)
      real :: dr(3)
      dr = p%drMax*(2.0*ranNum() - 1.0)
      riNew= p%r(i, :) + dr
   end subroutine translateMove

   subroutine rotateMove(p, i, uiNew)
      type(Particles), intent(inout) :: p
      integer, intent(in) :: i
      real, intent(inout) :: uiNew(3)
      real :: du(3)
      integer :: k

      call ranVec(du)

      do k = 1, 3
         du(k) = du(k)*rangeRanNum(-p%lambda, p%lambda)
      end do

      uiNew = p%u(i, :) + du
      uiNew = uiNew/sqrt(sum(uiNew**2))
   end subroutine rotateMove

   subroutine flipMove(p, i, uiNew)
      type(Particles), intent(inout) :: p
      integer, intent(in) :: i
      real, intent(inout) :: uiNew(3)
      uiNew = -p%u(i, :)
   end subroutine flipMove

end module Moves
