module Potential
   use Particle
   implicit none

   ! public :: checkOverlap, pairOverlap, pairPotential, carlosCheckOverlap, carlosSingleOverlap
contains
   function distSq(rijSq, rei, rej, eij, l) result(sijSq)
      implicit none
      real(8)             :: sijSq ! Returns squared distance between line segments
      real(8), intent(in) :: rijSq ! Squared centre-centre distance
      real(8), intent(in) :: rei    ! Scalar product rij.ei where ei is unit vector along i
      real(8), intent(in) :: rej    ! Scalar product rij.ej where ej is unit vector along j
      real(8), intent(in) :: eij    ! Scalar product ei.ej
      real(8), intent(in) :: l    ! Line segment length

      real            :: sinSq, ci, cj, ai, aj, di, dj, halfL
      real(8), parameter :: tol = 1.e-6

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
      real(8), intent(in) :: r(p%nParticles, 3), u(p%nParticles, 3)
      integer :: i, j

      do i = 1, p%nParticles
         do j = i + 1, p%nParticles
            if (pairOverlap(r, u, p, i, j)) then
               ! print *, 'Overlap detected between particles ', i, ' and ', j
               p%over = .true.
               return
            end if
         end do
      end do
      p%over = .false.
   end subroutine checkOverlap

   function pairOverlap(r, u, p, i, j) result(overlap)
      type(Particles), intent(in) :: p
      real(8), intent(in) :: r(p%nParticles, 3), u(p%nParticles, 3)
      integer, intent(in) :: i, j

      real(8) :: range, rangeBoxSq, lBoxSq, rijSq, rei, rej, eij, sijSq
      real(8) :: rij(3), rijPb(3), rijPbBox(3)
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

   subroutine singleParticleOverlap(p, r, u, i)
      type(Particles), intent(inout) :: p
      real(8), intent(in) :: r(p%nParticles, 3), u(p%nParticles, 3)
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

   real function minDistance(rij, rr, ei, ej, l)
      implicit none
      real(8) :: rui, ruj, uij
      real(8) :: sinsq, lambda, beta
      real(8), intent(in) ::  rr, l
      real(8), intent(in) ::  rij(3), ei(3), ej(3)

      minDistance = 0.0d0

      rui = dot_product(rij, ei)
      ruj = dot_product(rij, ej)
      uij = dot_product(ei, ej)

      sinsq = 1.0d0 - uij*uij

      if (sinsq > 0.0d0) then
         lambda = (ruj*uij - rui)/sinsq
         beta = (ruj - rui*uij)/sinsq
      else
         lambda = 0.0d0
         beta = 0.0d0
      end if

      if (abs(lambda) <= l .and. abs(beta) <= l) then
         MinDistance = rr + lambda*lambda + beta*beta + 2.0d0*lambda*rui &
                       - 2.0d0*beta*ruj - 2.0d0*lambda*beta*uij
         return
      end if

      if (abs(lambda) >= abs(beta)) then
         lambda = sign(l, lambda)
         beta = ruj + lambda*uij
      else
         beta = sign(l, beta)
         lambda = beta*uij - rui
      end if

      if (abs(beta) > l) beta = sign(l, beta)

      if (abs(lambda) > l) lambda = sign(l, lambda)

      minDistance = rr + lambda*lambda + beta*beta + 2.0d0*lambda*rui &
                    - 2.0d0*beta*ruj - 2.0d0*lambda*beta*uij

      return

   end function minDistance

   function carlosPairOverlap(p, r, u, i, j) result(overlap)
      type(Particles), intent(inout) :: p
      real(8), intent(in) :: r(p%nParticles, 3), u(p%nParticles, 3)
      integer :: i, j
      real(8) :: dd, rr, rij(3), ei(3), ej(3)
      logical :: overlap

      dd = 0.0

      print *, 'Checking overlap between particles ', i, ' and ', j
      print *, 'r(i, :) = ', r(i, :)
      print *, 'r(j, :) = ', r(j, :)

      rij = r(i, :) - r(j, :)
      rij = rij - p%lBox*anint(rij/p%lBox)
      rr = dot_product(rij, rij)
      ei = u(i, :)
      ej = u(j, :)

      print *, 'rij = ', rij

      print *, 'rr = ', rr
      dd = minDistance(rij, rr, ei, ej, p%l)

      print *, 'dd = ', dd

      if (dd < 1.0) then
         overlap = .true.
      else
         overlap = .false.
      end if

      return
   end function carlosPairOverlap

   subroutine carlosCheckOverlap(p, r, u)
      type(Particles), intent(inout) :: p
      real(8), intent(in) :: r(p%nParticles, 3), u(p%nParticles, 3)
      integer :: i, j

      do i = 1, p%nParticles
         do j = 1, p%nParticles
            if (i == j) cycle
            if (carlosPairOverlap(p, r, u, i, j)) then
               print *, 'Overlap detected between particles ', i, ' and ', j
               p%over = .true.
               return
            end if
         end do
      end do
      p%over = .false.
   end subroutine carlosCheckOverlap

   subroutine carlosSingleOverlap(p, r, u, i)
      type(Particles), intent(inout) :: p
      real(8), intent(in) :: r(p%nParticles, 3), u(p%nParticles, 3)
      integer, intent(in) :: i
      integer :: j

      do j = 1, p%nParticles
         if (j /= i) then
            if (carlosPairOverlap(p, r, u, i, j)) then
               p%over = .true.
               return
            end if
         end if
      end do
      p%over = .false.
   end subroutine carlosSingleOverlap

end module Potential

