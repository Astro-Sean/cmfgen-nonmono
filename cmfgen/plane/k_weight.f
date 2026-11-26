c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine k_weight(nw,angle,temp)
	USE SET_KIND_MODULE
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c Calculate integration weights for the radiation pressure, K.
c (see description in routine weights and blue derivation book)
c
c Written    8-95 DLM
c
c Altered    8-95 DLM Made more general.  Checks now depend on end points
c                     instead of 0->1.
c
c Altered 2/21/96 DLM Updated to F90 standard
c
c--------------------------------------------------------------------
c
      implicit none
c
c Number of angles: passed from calling routine
c
      integer :: nw
c
c Local variables
c
      integer :: j
c
      REAL(KIND=LDP) :: a1,a2,sum,check
      REAL(KIND=LDP), dimension(nw) :: temp(nw)
c
c Angles
c
      REAL(KIND=LDP), dimension(nw) :: angle
c
c--------------------------------------------------------------------
c
      temp(:)=0.0_LDP
c
      do j=1,nw-1
        a1=(angle(j)*angle(j)+
     *      angle(j)*angle(j+1)+
     *      angle(j+1)*angle(j+1))/6.0_LDP
        a2=(angle(j)*angle(j)*angle(j)+
     *      angle(j)*angle(j)*angle(j+1)+
     *      angle(j)*angle(j+1)*angle(j+1)+
     *      angle(j+1)*angle(j+1)*angle(j+1))/8.0_LDP
        temp(j)=temp(j)-angle(j+1)*a1+a2
        temp(j+1)=temp(j+1)+angle(j)*a1-a2
      enddo
c
c Check that 1/2*int[mu^2*dmu]=[mu_max^3-mu+min^3]/6
c
      sum=0.0_LDP
      do j=1,nw
        sum=sum+temp(j)
      enddo
c
      check=(angle(1)**3-angle(nw)**3)/6.0_LDP
c
      if(abs(check-sum).gt.1.0E-12_LDP)then
        print*,' K integration weights need to be normalized (1)'
        print*,' error=',abs(check-sum)
      endif
c
c Check that 1/2*int[mu^3*dmu]=[mu_max^4-mu_min^4]/8
c
      sum=0.0_LDP
      do j=1,nw
        sum=sum+temp(j)*angle(j)
      enddo
c
      check=(angle(1)**4-angle(nw)**4)/8.0_LDP
c
      if(abs(check-sum).gt.1.0E-12_LDP)then
        print*,' K integration weights need to be normalized (2)'
        print*,' error=',abs(check-sum)
      endif
c
      return
      end
