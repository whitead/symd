
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Stochastic velocity rescale, as described in
! Bussi, Donadio and Parrinello, J. Chem. Phys. 126, 014101 (2007)
!
! This subroutine implements Eq.(A7) and returns the new value for the kinetic energy,
! which can be used to rescale the velocities.
! The procedure can be applied to all atoms or to smaller groups.
! If it is applied to intersecting groups in sequence, the kinetic energy
! that is given as an input (kk) has to be up-to-date with respect to the previous rescalings.
!
! When applied to the entire system, and when performing standard molecular dynamics (fixed c.o.m. (center of mass))
! the degrees of freedom of the c.o.m. have to be discarded in the calculation of ndeg,
! and the c.o.m. momentum HAS TO BE SET TO ZERO.
! When applied to subgroups, one can chose to:
! (a) calculate the subgroup kinetic energy in the usual reference frame, and count the c.o.m. in ndeg
! (b) calculate the subgroup kinetic energy with respect to its c.o.m. motion, discard the c.o.m. in ndeg
!     and apply the rescale factor with respect to the subgroup c.o.m. velocity.
! They should be almost equivalent.
! If the subgroups are expected to move one respect to the other, the choice (b) should be better.
!
! If a null relaxation time is required (taut=0.0), the procedure reduces to an istantaneous
! randomization of the kinetic energy, as described in paragraph IIA.
!
! HOW TO CALCULATE THE EFFECTIVE-ENERGY DRIFT
! The effective-energy (htilde) drift can be used to check the integrator against discretization errors.
! The easiest recipe is:
! htilde = h + conint
! where h is the total energy (kinetic + potential)
! and conint is a quantity accumulated along the trajectory as minus the sum of all the increments of kinetic
! energy due to the thermostat.
!
function resamplekin(kk,sigma,ndeg,taut)
  implicit none
  real*8               :: resamplekin
  real*8,  intent(in)  :: kk    ! present value of the kinetic energy of the atoms to be thermalized (in arbitrary units)
  real*8,  intent(in)  :: sigma ! target average value of the kinetic energy (ndeg k_b T/2)  (in the same units as kk)
  integer, intent(in)  :: ndeg  ! number of degrees of freedom of the atoms to be thermalized
  real*8,  intent(in)  :: taut  ! relaxation time of the thermostat, in units of 'how often this routine is called'
  real*8 :: factor,rr
  real*8, external :: gasdev
  if(taut>0.1) then
    factor=exp(-1.0/taut)
  else
    factor=0.0
  end if
  rr = gasdev()
  resamplekin = kk + (1.0-factor)* (sigma*(sumnoises(ndeg-1)+rr**2)/ndeg-kk) &
               + 2.0*rr*sqrt(kk*sigma/ndeg*(1.0-factor)*factor)

contains 

double precision function sumnoises(nn)
  implicit none
  integer, intent(in) :: nn
! returns the sum of n independent gaussian noises squared
! (i.e. equivalent to summing the square of the return values of nn calls to gasdev)
  real*8, external :: gamdev,gasdev
  if(nn==0) then
    sumnoises=0.0
  else if(nn==1) then
    sumnoises=gasdev()**2
  else if(modulo(nn,2)==0) then
    sumnoises=2.0*gamdev(nn/2)
  else
    sumnoises=2.0*gamdev((nn-1)/2) + gasdev()**2
  end if
end function sumnoises

end function resamplekin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THE FOLLOWING ROUTINES ARE TRANSCRIBED FROM NUMERICAL RECIPES

double precision function gamdev(ia)
! gamma-distributed random number, implemented as described in numerical recipes

implicit none
integer, intent(in) :: ia
integer j
real*8 am,e,s,v1,v2,x,y
real*8, external :: ran1
if(ia.lt.1)pause 'bad argument in gamdev'
if(ia.lt.6)then
  x=1.
  do 11 j=1,ia
    x=x*ran1()
11  continue
  x=-log(x)
else
1 v1=2.*ran1()-1.
    v2=2.*ran1()-1.
  if(v1**2+v2**2.gt.1.)goto 1
    y=v2/v1
    am=ia-1
    s=sqrt(2.*am+1.)
    x=s*y+am
  if(x.le.0.)goto 1
    e=(1.+y**2)*exp(am*log(x/am)-s*y)
  if(ran1().gt.e)goto 1
endif
gamdev=x
end

double precision function gasdev()
! gaussian-distributed random number, implemented as described in numerical recipes

implicit none
integer, save :: iset = 0
real*8, save :: gset
real*8, external :: ran1
real*8 fac,rsq,v1,v2
if(iset==0) then
1       v1=2.*ran1()-1.0d0
  v2=2.*ran1()-1.0d0
  rsq=v1**2+v2**2
  if(rsq.ge.1..or.rsq.eq.0.)goto 1
  fac=sqrt(-2.*log(rsq)/rsq)
  gset=v1*fac
  gasdev=v2*fac
  iset=1
else
  gasdev=gset
  iset=0
end if
end

FUNCTION ran1()
! random number generator
INTEGER IA,IM,IQ,IR,NTAB,NDIV
REAL ran1,AM,EPS,RNMX
PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
  NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
INTEGER j,k,iv(NTAB),iy
SAVE iv,iy
DATA iv /NTAB*0/, iy /0/
INTEGER, SAVE :: idum=0 !! ATTENTION: THE SEED IS HARDCODED
if (idum.le.0.or.iy.eq.0) then
  idum=max(-idum,1)
  do 11 j=NTAB+8,1,-1
    k=idum/IQ
    idum=IA*(idum-k*IQ)-IR*k
    if (idum.lt.0) idum=idum+IM
    if (j.le.NTAB) iv(j)=idum
11  continue
  iy=iv(1)
endif
k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if (idum.lt.0) idum=idum+IM
j=1+iy/NDIV
iy=iv(j)
iv(j)=idum
ran1=min(AM*iy,RNMX)
return
      END

