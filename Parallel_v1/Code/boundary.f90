! Author: Jaume Garcia

module boundary
! Mòdul que conté les dues subrutines amb les condicions periòdiques de contorn 
implicit none
contains

subroutine PBC1(L,x)
! Subrutina que calcula les boundary conditions aplicant la minimum image convention
! i si una partícula surt de la caixa, entra per l'altre costat.
! Argument d'entrada: tamany de la caixa L
! Argument de sortida: posició de la partícula, x.
! Adreçada al càlcul de la força LJ.
implicit none
double precision x,L
if (x > L/2.d0) then
   x = x - L  
endif  
if (x < -L/2.d0) then
   x = x + L
endif    
return 

end subroutine 
end module boundary

