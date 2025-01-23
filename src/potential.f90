module Potential
   use globals
   use Utils
   implicit none

contains
   function nematicPotential(p, u) result(pot)
      type(Particles), intent(in) :: p
      real(8), intent(in) :: u(p%nParticles, 3)
      integer :: i
      real(8) :: pot
      pot = 0.0

      do i = 1, p%nParticles
         pot = pot + p%w0*(1.5*u(i, 2)*u(i, 2) - 0.5)
      end do
   end function nematicPotential

   function pairPotential(p, r, u, i, j) result(pot)
      type(Particles), intent(inout) :: p
      integer, intent(in) :: i, j
      real(8), intent(in) :: r(p%nParticles, 3), u(p%nParticles, 3)
      real:: rq(2, size(p%q), 3), rij(3), rijPb(3)
      real(8) :: rr
      integer:: k, qSign, qSize, q1, q2
      real(8) :: pot
      pot = 0.0

      qSize = size(p%q)
      do k = 1, qSize
         qSign = dsign(1.0d0, p%q(k))
         rq(1, k, :) = r(i, :) + qSign*(1.0/qSize)*p%l*u(i, :)
         rq(2, k, :) = r(j, :) + qSign*(1.0/qSize)*p%l*u(j, :)
      end do

      do q1 = 1, qSize
         do q2 = 1, qSize
            rij = rq(1, q1, :) - rq(2, q2, :)
            rijPb = rij - p%lBox*anint(rij/p%lBox)
            rr = sqrt(sum(rijPb**2))
            if (rr < 1.0) then
               p%over = .true.
            end if
            pot = pot + p%q(q1)*p%q(q2)*exp(-p%kappa*rr)/rr
         end do
      end do

   end function pairPotential

   function calcPotential(p, r, u) result(pot)
      type(Particles), intent(inout) :: p
      real(8), intent(in) :: r(p%nParticles, 3), u(p%nParticles, 3)
      integer :: i, j
      real(8) :: pot

      pot = 0.0
      do i = 1, p%nParticles
         do j = 1, p%nParticles
            if (i == j) cycle
            pot = pot + pairPotential(p, r, u, i, j)
         end do
      end do

      pot = pot + nematicPotential(p, u)
   end function calcPotential

   subroutine metropolis(delta, p, accept)
      real(8), intent(in) :: delta
      logical, intent(out) :: accept
      type(Particles), intent(inout) :: p

      real(8) :: betaDelta, ran

      betaDelta = p%beta*delta
      ran = ranNum()

      if (betaDelta < 0.0 .or. exp(-betaDelta) > ran) then
         accept = .true.
      else
         accept = .false.
      end if
   end subroutine metropolis

   function singleParticlePotential(p, r, u, i) result(pot)
      type(Particles), intent(inout) :: p
      real(8), intent(in) :: r(p%nParticles, 3), u(p%nParticles, 3)
      integer, intent(in) :: i
      integer :: j
      real(8) :: pot
      pot = 0.0

      do j = 1, p%nParticles
         if (i == j) cycle
         pot = pot + pairPotential(p, r, u, i, j)
      end do

      pot = pot + p%w0*(1.5*u(i, 2)*u(i, 2) - 0.5)

   end function singleParticlePotential

end module Potential

