module Potential
   use Particle
   implicit none

   public :: checkOverlap, pairOverlap, pairPotential
contains
   function distSq(rijSq, rei, rej, eij, l) result(sijSq)
      implicit none
      real             :: sijSq ! Returns squared distance between line segments
      real, intent(in) :: rijSq ! Squared centre-centre distance
      real, intent(in) :: rei    ! Scalar product rij.ei where ei is unit vector along i
      real, intent(in) :: rej    ! Scalar product rij.ej where ej is unit vector along j
      real, intent(in) :: eij    ! Scalar product ei.ej
      real, intent(in) :: l    ! Line segment length

      real            :: sinSq, ci, cj, ai, aj, di, dj, halfL
      real, parameter :: tol = 1.e-6

      sinSq = 1.0 - eij**2 ! Squared sine of angle between line segments
      halfL = l/2.0    ! Half the line segment length

      IF (sinSq < tol) THEN ! Guard against nearly-parallel lines
         ci = -rei
         cj = rej
      ELSE
         ci = (-rei + eij*rej)/sinSq
         cj = (rej - eij*rei)/sinSq
      END IF

      ai = abs(ci)
      aj = abs(cj)
      IF (ai > halfL) ci = sign(halfL, ci)
      IF (aj > halfL) cj = sign(halfL, cj)

      IF (ai > aj) THEN
         cj = rej + ci*eij
      ELSE
         ci = -rei + cj*eij
      END IF

      ai = abs(ci)
      aj = abs(cj)
      IF (ai > halfL) ci = sign(halfL, ci)
      IF (aj > halfL) cj = sign(halfL, cj)

      di = 2.0*rei + ci - cj*eij
      dj = -2.0*rej + cj - ci*eij

      sijSq = rijSq + ci*di + cj*dj
   end function distSq

   subroutine checkOverlap(r, u, p)
      type(Particles), intent(inout) :: p
      real, intent(in) :: r(p%nParticles, 3), u(p%nParticles, 3)
      integer :: i, j

      do i = 1, p%nParticles
         do j = i + 1, p%nParticles
            if (pairOverlap(r, u, p, i, j)) then
               print *, 'Overlap detected between particles ', i, ' and ', j
               p%over = .true.
               return
            end if
         end do
      end do
      p%over = .false.
   end subroutine checkOverlap

   function pairOverlap(r, u, p, i, j) result(overlap)
      type(Particles), intent(in) :: p
      real, intent(in) :: r(p%nParticles, 3), u(p%nParticles, 3)
      integer, intent(in) :: i, j

      real :: range, rangeBoxSq, lBoxSq, rijSq, rei, rej, eij, sijSq
      real :: rij(3), rijPb(3), rijPbBox(3)
      logical :: overlap

      overlap = .false.

      range = 1 + p%l
      rangeBoxSq = (range/p%lBox)**2
      lBoxSq = p%lBox**2

      rij = r(i, :) - r(j, :)
      rijPb = rij - p%lBox*anint(rij/p%lBox)
      rijPbBox = rijPb*p%lBox
      rijSq = sum(rijPb**2)*lBoxSq

      if (rijSq > rangeBoxSq*lBoxSq) then
         return
      end if

      rei = dot_product(rijPbBox, u(i, :))
      rej = dot_product(rijPbBox, u(j, :))
      eij = dot_product(u(i, :), u(j, :))

      sijSq = distSq(rijSq, rei, rej, eij, p%lBox)

      if (sijSq < 1) then
         overlap = .true.
         return
      end if
   end function pairOverlap

   function nematicPotential(p, u) result(pot)
      type(Particles), intent(in) :: p
      real, intent(in) :: u(p%nParticles, 3)
      integer :: i
      real :: pot
      pot = 0.0

      do i = 1, p%nParticles
         pot = pot + p%w0*(1.5*u(i, 2)*u(i, 2) - 0.5)
      end do
   end function nematicPotential

   function pairPotential(p, r, u, i, j) result(pot)
      type(Particles), intent(inout) :: p
      integer, intent(in) :: i, j
      real, intent(in) :: r(p%nParticles, 3), u(p%nParticles, 3)
      real:: rq(2, size(p%q), 3), rij(3), rijPb(3)
      real :: rr
      integer:: k, qSign, qSize, q1, q2
      real :: pot
      pot = 0.0

      qSize = size(p%q)
      do k = 1, qSize
         qSign = sign(1.0, p%q(k))
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
      real, intent(in) :: r(p%nParticles, 3), u(p%nParticles, 3)
      integer :: i, j
      real :: pot

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
      real, intent(in) :: delta
      logical, intent(out) :: accept
      type(Particles), intent(inout) :: p

      real :: betaDelta, ran

      betaDelta = p%beta*delta
      ran = ranNum()

      if (betaDelta < 0.0 .or. exp(-betaDelta) > ran) then
         accept = .true.
      else
         accept = .false.
      end if

   end subroutine metropolis

   subroutine singleParticleOverlap(p, r, u, i)
      type(Particles), intent(inout) :: p
      real, intent(in) :: r(p%nParticles, 3), u(p%nParticles, 3)
      integer, intent(in) :: i
      integer :: j

      do j = 1, p%nParticles
         if (i == j) cycle
         if (pairOverlap(r, u, p, i, j)) then
            ! print *, 'Overlap detected between particles ', i, ' and ', j
            p%over = .true.
            return
         end if
      end do
      p%over = .false.

   end subroutine singleParticleOverlap

   function singleParticlePotential(p, r, u, i) result(pot)
      type(Particles), intent(inout) :: p
      real, intent(in) :: r(p%nParticles, 3), u(p%nParticles, 3)
      integer, intent(in) :: i
      integer :: j
      real :: pot
      pot = 0.0

      do j = 1, p%nParticles
         if (i == j) cycle
         pot = pot + pairPotential(p, r, u, i, j)
      end do

      pot = pot + p%w0*(1.5*u(i, 2)*u(i, 2) - 0.5)

   end function singleParticlePotential

end module Potential
