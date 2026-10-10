c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine derivs(x,j,djds,beta,dbetadr,gamma)
	USE SET_KIND_MODULE
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Subroutine to calculate the value of the derivitive
c djds at the given j point.  Called in subroutines charactistics and
c runge_kutta
c
c Altered 2/21/96 DLM Updated to F90 standard
c
c--------------------------------------------------------------------
c
      implicit none
c
      REAL(KIND=LDP), dimension(2) :: j,djds
      REAL(KIND=LDP) x,r,mu,beta,dbetadr,gamma
c
      r=j(1)
      mu=j(2)
c
      djds(1)=gamma*(mu+beta)
      djds(2)=gamma*(1.0_LDP-mu*mu)*((1.0_LDP+beta*mu)/r-
     *     gamma*gamma*(mu+beta)*dbetadr)
c
      return
      end
c
