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
               print *, "Overlap between particles in Allen and Tildesly routine: ",i,j
               print *, "r(i): ",r(i,:), "u(i): ",u(i,:)
               print *, "r(j): ",r(j,:), "u(j): ",u(j,:)
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

   subroutine carlosSingleParticleOverlap(r,u,p,i,lBox)
      type(Particles), intent(inout) :: p
      real(8), intent(in) :: r(p%nParticles, 3), u(p%nParticles, 3),lBox
      integer, intent(in) :: i

      integer :: j

      do j = 1, p%nParticles
         if (i /= j) then
            if (carlosPairOverlap(r,u,p,i,j,lBox)) then
               ! print *, 'Overlap detected between particles ', i, ' and ', j
               p%over = .true.
               return
            end if
         end if
      end do

      p%over = .false.
   end subroutine carlosSingleParticleOverlap

   subroutine carlosCheckOverlap(r,u,p,lBox)
      type(Particles), intent(inout) :: p
      real(8), intent(in) :: r(p%nParticles, 3), u(p%nParticles, 3),lBox
      integer :: i, j

      do i = 1, p%nParticles
         do j = i + 1, p%nParticles
            if (carlosPairOverlap(r,u,p,i,j,lBox)) then
               p%over = .true.
               return
            end if
         end do
      end do
      p%over = .false.
   end subroutine carlosCheckOverlap


   function carlosPairOverlap(r,u,p,i,j,lBox) result(overlap)
      implicit none
      type(Particles), intent(inout) :: p
      real(8), intent(in) :: r(p%nParticles, 3), u(p%nParticles, 3),lBox
      integer, intent(in) :: i, j
      logical :: overlap
      real(8) :: rij(3), rr, dd,dd2, range,ui(3),uj(3)


      range = 1.0d0 + p%l
      overlap = .false.

      rij = r(j,:) - r(i,:)
      rij = rij - anint(rij/lBox)*lBox
      rr = dot_product(rij,rij)


      ui = u(i,:)
      uj = u(j,:)

      if (rr <= range**2) then
         dd = minDistance(rij,rr,ui,uj,p%l*0.5)
         dd2 = shortestDistance(rij,ui,uj,p%l*0.5)

         if (dd < 1.0) then
            overlap = .true.
         end if

      end if

      return
   end function carlosPairOverlap

   real(8) function minDistance( rij, rr, eParti, ePartj, ld)
      implicit none
      real(8) :: rui, ruj, uij
      real(8) :: sinsq, lambda, beta
      real(8), intent(in) ::  rr, ld
      real(8), intent(in) ::  rij(3), ePartI(3), ePartJ(3)


      minDistance = 0.0d0

      ! ... spherocylinder-spherocylinder criterion

      rui = dot_product(rij,ePartI)
      ruj = dot_product(rij,ePartJ)
      uij = dot_product(ePartI,ePartJ)


      sinsq = 1.0d0 - uij*uij

      if( sinsq > 0.0d0 ) then
         lambda = ( ruj*uij - rui )/sinsq
         beta = ( ruj - rui*uij )/sinsq
      else
         lambda = 0.0d0
         beta = 0.0d0
      end if

      if( dabs( lambda ) <= ld .and. dabs( beta ) <= ld ) then
         MinDistance = rr + lambda*lambda + beta*beta + 2.0d0*lambda*rui &
            - 2.0d0*beta*ruj - 2.0d0*lambda*beta*uij
         return
      end if

      if( dabs( lambda ) >=  dabs( beta ) ) then
         lambda = dsign( ld, lambda )
         beta = ruj + lambda*uij
      else
         beta = dsign( ld, beta )
         lambda = beta*uij - rui
      end if

      if( dabs( beta ) > ld ) beta = dsign( ld, beta )

      if( dabs( lambda ) > ld ) lambda = dsign( ld, lambda )


      minDistance = rr + lambda*lambda + beta*beta + 2.0d0*lambda*rui &
         - 2.0d0*beta*ruj - 2.0d0*lambda*beta*uij

      return


   end function minDistance


   function shortestDistance(rij,ui,uj,halfL) result (distance)
      real(8), intent(in) :: rij(3), ui(3), uj(3), halfL
      real(8) :: rr, rui, ruj, uij
      real(8) :: xLanda, xMu,aux1,aux2,cc,distance
      real(8) :: tol = 1.0d-6

      rr = dot_product(rij,rij)
      rui = dot_product(rij,ui)
      ruj = dot_product(rij,uj)
      uij = dot_product(ui,uj)

      cc = 1.0d0 - uij**2

   end function shortestDistance


end module Overlap
