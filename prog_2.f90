program prandtl
implicit none
integer, parameter:: i0 = 12 ! output unit
integer i,j,s,smax,ni,nj,n
real(8) H, L,u0,nu,e,eps,eps1,G_in,G
real(8) dx,dy
real(8),allocatable :: a(:),b(:),c(:),d(:),per(:)
real(8),allocatable :: u(:,:),uold(:,:),v(:,:),vold(:,:),p(:,:),pold(:,:),c_f(:,:),r_n(:,:),x(:,:),y(:,:)
smax=300 
ni=51
nj=21
H=0.1
L=1
u0=1.0 
nu=0.01
e=1
s=1
eps=0.00001
eps1=0.00001
G_in=u0*H
allocate(a(nj),b(nj),c(nj),d(nj),per(nj))
allocate(u(ni,nj),uold(ni,nj),v(ni,nj),vold(ni,nj),p(ni,nj),pold(ni,nj),r_n(ni,nj),c_f(ni,nj),x(ni,nj),y(ni,nj))
dx=L/(ni-1)
dy=H/(nj-1)
do i=1,ni
do j=1,nj
x(i,j)=dx*(i-1)
y(i,j)=dy*(j-1)
end do
end do
call boundary(u,v,p)
do i=2,ni
G=0.0
p(i,1:nj)=p(i-1,1:nj)
vold(i,1:nj)=v(i-1,1:nj)
uold(i,1:nj)=u(i-1,1:nj)
n=1 
do while(((abs(G_in-G)/G_in)>eps1).AND.(n<smax))
s=1
e=1
do while((e>eps).AND.(s<smax))
a(1)=0.0
b(1)=1.0
c(1)=0.0
d(1)=0.0
a(nj)=-1.0
b(nj)=1.0
c(nj)=0.0
d(nj)=0.0
do j=2,nj-1
a(j)=(-vold(i,j-1)/(2*dy))-nu/(dy**2)
b(j)=(uold(i,j)/dx + 2.0*nu/(dy**2))
c(j)=(vold(i,j+1)/(2.0*dy)-nu/(dy**2))
d(j)=(uold(i-1,j)**2)/dx-(p(i,j)-p(i-1,j))/dx
end do
call progonka (per,nj,a,b,c,d)
e=0.0
do j=1,nj
u(i,j)=per(j)
e=max(e,(abs(u(i,j)-uold(i,j)))/uold(i,j))
end do
do j=2,nj
v(i,j)= v(i,j-1)-dy/2.0*((u(i,j)-u(i-1,j))/dx+(u(i,j-1)-u(i-1,j-1))/dx)
e=max(e,(abs(v(i,j)-vold(i,j)))/vold(i,j))
end do
v(i,1)=0
do j=1,nj
uold(i,j)= u(i,j)
vold(i,j)=v(i,j)
end do
s=s+1
end do
G=0.0
do j=1,nj
G=G+u(i,j)*dy
end do
do j=1,nj
p(i,j)=p(i,j)-((G_in*(G_in-G)*0.1)/(H**2))
end do
n=n+1
end do
print*, x(1,j),G
end do

do i=1,ni
r_n(i,:)=1000
!print*, p(i,:)
end do

do i=1,ni
do j=1,nj
pold(i,j)=abs(p(ni-i+1,nj-j+1))
end do
end do

call C_f_calculation(nu,ni,nj,dy,u0,u,c_f)

write(*,*) 'Output data' 
Open(i0,FILE='Results_incompressible.plt')
Call Output_Fields(i0,ni,nj,x,y,u,v,p,r_n,c_f)

 !print*,x(1:ni,1)
 !print*,y(1,1:nj)

Close(i0)
deallocate(x,y,u,v,p,c_f)

contains

subroutine progonka (z,n,a,b,c,d)
implicit none
integer n,i
real(8),dimension(n) :: z,a,b,c,d
real(8),dimension(n) :: alp,bet
alp(1)=-c(1)/b(1)
bet(1)=d(1)/b(1)
do i=2,n-1
alp(i)=-c(i)/(b(i)+a(i)*alp(i-1))
bet(i)=(d(i)-a(i)*bet(i-1))/(b(i)+a(i)*alp(i-1))
end do
z(n)=(d(n)-a(n)*bet(n-1))/(b(n)+a(n)*alp(n-1))
do i=n-1,1,-1
z(i)=alp(i)*z(i+1)+bet(i)
end do
end subroutine

subroutine boundary (u,v,p)
real(8),dimension(ni,nj):: u,v,p
u(1,1:nj)=u0
u(1:ni,1)=0
v(1:ni,1)=0
u(1:ni,nj)=u0
v(1,1:nj)=0
p(1,1:nj)=0
end subroutine

subroutine C_f_calculation(nu,ni,nj,dy,u0,u,c_f)
implicit none
integer ni, nj, i
real(8) :: nu, dy, u0
real(8) :: u(ni,nj), c_f(ni, nj)
do i = 1, ni
c_f(i, :) = 2*nu/(u0**2.0)*((u(i,2)-u(i,1))/dy)    
end do
end subroutine

subroutine Output_Fields(i0,ni,nj,x,y,u,v,p,r_n,c_f)
implicit none
integer ni,nj,i0, j
real(8),dimension(ni,nj):: u,v,p,c_f,r_n,x,y
Write(i0,*) 'VARIABLES = "X", "Y", "U", "V", "P", "RO", "C_f" ' 
Write(i0,*) 'ZONE i=',ni,', j=',nj, ', DATAPACKING=BLOCK'
Write(i0,'(100E25.16)') x(1:ni,1:nj) 
Write(i0,'(100E25.16)') y(1:ni,1:nj)
Write(i0,'(100E25.16)') u(1:ni,1:nj)
Write(i0,'(100E25.16)') v(1:ni,1:nj)
Write(i0,'(100E25.16)') p(1:ni,1:nj)
Write(i0,'(100E25.16)') r_n(1:ni,1:nj)
Write(i0,'(100E25.16)') c_f(1:ni,1:nj)
end subroutine
end program


   
      