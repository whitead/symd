program rdf
implicit none
integer           :: iostat
logical           :: first=.true.
integer           :: natoms,natoms_check
real              :: cell(3),cell_check(3)
integer           :: iatom,jatom
real, allocatable :: positions(:,:)
real, allocatable :: histogram(:)
real              :: distance(3),distance_pbc(3),d
character(10)     :: atomname
integer           :: istep,ibigstep
integer           :: jnow,norma,ilag,jold
real              :: rmax,deltar
integer           :: nbin,ibin
real, parameter   :: pi=3.14159265358979323844
real              :: nideal

nbin=200
allocate(histogram(nbin))
histogram=0
rmax=3.0
deltar=rmax/nbin
istep=-1
do
  read(*,*,iostat=iostat) natoms
  if(iostat/=0) exit
  read(*,*) cell
  if(first) then
    natoms_check=natoms
    cell_check=cell
    allocate(positions(3,natoms))
    first=.false.
  end if
  if(natoms/=natoms_check) stop
  if(any(cell/=cell_check)) stop
  do iatom=1,natoms
    read(*,*) atomname,positions(:,iatom)
  end do
  do iatom=1,natoms-1
    do jatom=iatom+1,natoms
      distance=positions(:,iatom)-positions(:,jatom)
      call pbc(cell,distance,distance_pbc)
      d=sqrt(sum(distance_pbc**2))
      ibin=int(d/deltar)+1
      if(ibin<=nbin .and. ibin>0) then
        histogram(ibin)=histogram(ibin)+2.0
      end if
    end do
  end do
  istep=istep+1
end do

do ibin=1,nbin
  nideal=4.0*pi*natoms/product(cell)/3.0*deltar**3*(ibin**3-(ibin-1)**3)
  write(*,"(2g16.5)") (ibin-0.5)*deltar,histogram(ibin)/nideal/istep/natoms
end do


end program rdf

subroutine pbc(cell,vin,vout)
! apply periodic boundary condition to a vector
  implicit none
  real, intent(in) :: cell(3)
  real, intent(in) :: vin(3)
  real, intent(out) :: vout(3)
  integer :: i
  do i=1,3
    vout(i)=vin(i)-nint(vin(i)/cell(i))*cell(i)
  end do
end subroutine pbc
