program diffusion
implicit none
integer           :: iostat
logical           :: first=.true.
integer           :: natoms,natoms_check
real              :: cell(3),cell_check(3)
integer           :: iatom
integer           :: i
real, allocatable :: positions(:,:),positions0(:,:)
real              :: displacement
character(10)     :: atomname
integer           :: istep

displacement=0.0
istep=-1
do
  read(*,*,iostat=iostat) natoms
  if(iostat/=0) exit
  read(*,*) cell
  if(first) then
    natoms_check=natoms
    cell_check=cell
    allocate(positions(3,natoms))
    allocate(positions0(3,natoms))
  end if
  if(natoms/=natoms_check) stop
  if(any(cell/=cell_check)) stop
  do iatom=1,natoms
    read(*,*) atomname,positions(:,iatom)
  end do
  do i=1,3
    positions(i,:)=positions(i,:)-sum(positions(i,:))/natoms
  end do
  if(first) then
    positions0=positions
    first=.false.
  end if
  istep=istep+1
  write(*,*) istep,sum((positions-positions0)**2)/size(positions)
end do

end program diffusion
