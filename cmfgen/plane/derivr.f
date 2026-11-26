c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine derivr(x,j,djdr,beta,dbetadr,gamma)
	USE SET_KIND_MODULE
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Subroutine to calculate the value of the derivitive
c djdr at the given j point.  Called in subroutines characteristics and
c runge_kutta
c
c Altered 2/21/96 DLM Updated to F90 standard
c
c--------------------------------------------------------------------
c
      implicit none
c
      REAL(KIND=LDP), dimension(2) :: j,djdr
      REAL(KIND=LDP) x,r,mu,beta,dbetadr,gamma
c
      r=x
      mu=j(2)
c
      djdr(1)=1.0_LDP/gamma/(mu+beta)
      djdr(2)=(1.0_LDP-mu*mu)/(mu+beta)*((1.0_LDP+beta*mu)/r-
     *     gamma*gamma*(mu+beta)*dbetadr)
c
      return
      end
c
