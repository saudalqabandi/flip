module Accumalators
   use globals
   implicit none
contains

   subroutine zeroAccumalators(acc)
      type(Accumalator), intent(out) :: acc
      acc%nMoves = 0
      acc%nTrans = 0
      acc%nRot = 0
      acc%nFlip = 0
      acc%nVol = 0
      acc%nAccept = 0
      acc%nTransAccept = 0
      acc%nRotAccept = 0
      acc%nFlipAccept = 0
      acc%nVolAccept = 0
      acc%ratioTrans = 0.0
      acc%ratioRot = 0.0
      acc%ratioFlip = 0.0
   end subroutine zeroAccumalators

   subroutine updateAccumalators(acc, moveType, accept)
      type(Accumalator), intent(inout) :: acc
      character(len=*), intent(in) :: moveType
      logical, intent(in) :: accept

      select case (moveType)

       case ("translate")
         acc%nTrans = acc%nTrans + 1
         if (accept) acc%nTransAccept = acc%nTransAccept + 1
         acc%ratioTrans = acc%nTransAccept/real(acc%nTrans)

       case ("rotate")
         acc%nRot = acc%nRot + 1
         if (accept) acc%nRotAccept = acc%nRotAccept + 1
         acc%ratioRot = acc%nRotAccept/real(acc%nRot)

       case ("flip")
         acc%nFlip = acc%nFlip + 1
         if (accept) acc%nFlipAccept = acc%nFlipAccept + 1
         acc%ratioFlip = acc%nFlipAccept/real(acc%nFlip)

       case ("vol")
         acc%nVol = acc%nVol + 1
         if (accept) acc%nVolAccept = acc%nVolAccept + 1
         acc%ratioVol = acc%nVolAccept/real(acc%nVol)

      end select

      acc%nMoves = acc%nMoves + 1
      if (accept) acc%nAccept = acc%nAccept + 1
      acc%ratio = acc%nAccept/real(acc%nMoves)

   end subroutine updateAccumalators

   subroutine printAccumalators(acc)
      type(Accumalator), intent(in) :: acc
      print '(A, I10, A, I10, A, F6.2)', 'Moves:        ', acc%nMoves, ' Accepted: ', acc%nAccept, ' Ratio: ', acc%ratio
      print '(A, I10, A, I10, A, F6.4)', 'Translations: ', acc%nTrans, ' Accepted: ', acc%nTransAccept, ' Ratio: ', acc%ratioTrans
      print '(A, I10, A, I10, A, F6.4)', 'Rotations:    ', acc%nRot, ' Accepted: ', acc%nRotAccept, ' Ratio: ', acc%ratioRot
      print '(A, I10, A, I10, A, F6.4)', 'Flips:        ', acc%nFlip, ' Accepted: ', acc%nFlipAccept, ' Ratio: ', acc%ratioFlip
      print '(A, I10, A, I10, A, F6.4)', 'Volume changes:', acc%nVol, ' Accepted: ', acc%nVolAccept, ' Ratio: ', acc%ratioVol
   end subroutine printAccumalators

end module Accumalators
