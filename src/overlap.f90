module Overlap
   use globals
   implicit none

contains
   function dist_sq(rij_sq, rei, rej, eij, ell) result(sij_sq)
      implicit none
      real(8) :: sij_sq ! Returns squared distance between line segments
      real(8), intent(in) :: rij_sq ! Squared centre-centre distance
      real(8), intent(in) :: rei    ! Scalar product rij.ei where ei is unit vector along i
      real(8), intent(in) :: rej    ! Scalar product rij.ej where ej is unit vector along j
      real(8), intent(in) :: eij    ! Scalar product ei.ej
      real(8), intent(in) :: ell    ! Line segment length

      real(8) :: sin_sq, ci, cj, ai, aj, di, dj, ell2
      real(8), parameter :: tol = 1.0e-6

      sin_sq = 1.0 - eij**2 ! Squared sine of angle between line segments
      ell2 = ell / 2.0      ! Half the line segment length

      if (sin_sq < tol) then ! Guard against nearly-parallel lines
         ci = -rei
         cj = rej
      else
         ci = (-rei + eij * rej) / sin_sq
         cj = (rej - eij * rei) / sin_sq
      end if

      ai = abs(ci)
      aj = abs(cj)
      if (ai > ell2) ci = sign(ell2, ci)
      if (aj > ell2) cj = sign(ell2, cj)

      if (ai > aj) then
         cj = rej + ci * eij
      else
         ci = -rei + cj * eij
      end if

      ai = abs(ci)
      aj = abs(cj)
      if (ai > ell2) ci = sign(ell2, ci)
      if (aj > ell2) cj = sign(ell2, cj)

      di = 2.0 * rei + ci - cj * eij
      dj = -2.0 * rej + cj - ci * eij

      sij_sq = rij_sq + ci * di + cj * dj

   end function dist_sq


   function pairOverlap(r, u, p, i, j) result(overlap)
      type(Particles), intent(inout) :: p
      real(8), intent(in) :: r(p%nParticles, 3), u(p%nParticles, 3)
      integer, intent(in) :: i, j
      logical :: overlap

      real(8) :: rij(3), rij_sq, rei, rej, eij, ell, sij_sq,range,range_box_sq

      range = 1+p%l
      range_box_sq = (range/p%lBox)**2

      rij(:) = r(j, :) - r(i, :)
      rij(:) = rij(:) - anint(rij(:)) !PBC in box=1 units
      rij_sq = dot_product(rij, rij)

      if (rij_sq > range_box_sq) then
         overlap = .false.
         return
      end if

      rij_sq = rij_sq *p%lBox**2
      rij = rij * p%lBox
      rei = dot_product(rij, u(i, :))
      rej = dot_product(rij, u(j, :))
      eij = dot_product(u(i, :), u(j, :))
      ell = p%l

      sij_sq = dist_sq(rij_sq, rei, rej, eij, ell)

      overlap = sij_sq < 1

   end function pairOverlap

   subroutine checkOverlap(r, u, p)
      type(Particles), intent(inout) :: p
      real(8), intent(in) :: r(p%nParticles, 3), u(p%nParticles, 3)
      integer :: i, j

      do i = 1, p%nParticles
         do j = i + 1, p%nParticles
            if (pairOverlap(r, u, p, i, j)) then
               p%over = .true.
               return
            end if
         end do
      end do
      p%over = .false.
   end subroutine checkOverlap


   subroutine singleParticleOverlap(r, u, p, i)
      type(Particles), intent(inout) :: p
      real(8), intent(in) :: r(p%nParticles, 3), u(p%nParticles, 3)
      integer, intent(in) :: i

      integer :: j

      do j = 1, p%nParticles
         if (i /= j) then
            if (pairOverlap(r, u, p, i, j)) then
               ! print *, 'Overlap detected between particles ', i, ' and ', j
               p%over = .true.
               return
            end if
         end if
      end do
      p%over = .false.
   end subroutine singleParticleOverlap


end module Overlap
