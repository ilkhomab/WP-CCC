program protonH 

use precision_type
use mpi
use data_base,only:v,vmat_dir,vmat_ex,vmat_obk,en_t,n_t,l_t,m_t,num_t,zint,nzint,vmat_dir_ar,vmat_ex_ar,vmat_obk_ar,mppic,b,key_c,key_num,wig_d,nstep,step,cu,omega,nomega,vmat_dir_om,vmat_ex_om,vmat_ex0_om,homega,zero,projectile_charge,probfile,wffile,ejenfile,nb
use omp_lib
use openacc
implicit none
include 'mpif.h'
integer i,ist_f, ndp, nchx,l,m,n,f,j,ii,ix,iz,ixmax,izmax,istep,ia,info
parameter (ndp = 30)
real(kind=id)::t1,t2,t3,t4,t5,exp_integ,dir_init,plgndr,alf,xx,rho,zz,xlam,rylm,z,x,dx,dz,y,xmin,xmax,zmin,zmax,radwf,tt,temp,ttt,uu,bigr,smallr,rmr,xrmr,yrmr,temp2,temp3,cosrmr,cosang,cosr,sinr,rj,ry,rjp,ryp
real(kind=id):: mpreal,ss,direct_num,direct_amp,exch_amp,exch0_amp,rcond
complex(kind=id)::znum,zsum,c_phase  
complex(kind=id),allocatable,dimension(:)::zam
complex(kind=id),allocatable,dimension(:,:)::vmatL,vmatR,af
complex(kind=id),allocatable,dimension(:)::work,xxx
real(kind=id),allocatable,dimension(:)::r,c,rwork,swork
integer,allocatable,dimension(:)::ipiv
real(kind=8),dimension(1) :: ferr, berr 
character*1 chx(ndp+200)
integer,parameter::nmax_f=100,nmax_i=16
real(kind=id),dimension(1:nmax_f,1:nmax_i):: dir_el
!real(kind=id),dimension(1:nmax_f):: x,w 
real(kind=id),dimension(0:nmax_f):: sj,dj 
real(kind=id),dimension(1:5,1:5):: resr,resi 
character (len=132) :: skip
logical::ex
character*1::fact,trans,equed

myid=-1
ntasks=1
! MPI initialization and environment characteristics
! note: make these your first (or nearly so) executable statements
call MPI_INIT( ierr )
call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
call MPI_COMM_SIZE( MPI_COMM_WORLD, ntasks, ierr )
nomp = min(OMP_GET_NUM_PROCS(),OMP_GET_MAX_THREADS()) 
nodes = ntasks / nomp
if (mod(myid,nomp).eq.0) then
   nodeid = myid/nomp + 1
if(myid==0) then
   print*,'nomp',nomp
   print*,'nodes',nodes
   print*,'ntasks',ntasks
endif
   print*,'myid',myid
   print*,'nodeid',nodeid

!projectile_charge=-1.e0_id
!inquire(file='prob_2s00',exist=probfile)
inquire(file='prob',exist=probfile)
if(myid==0)then
  print*,'probfile=',probfile
endif
  
if(.not.probfile)then
   call acc_set_device_num(myid,acc_device_nvidia)
endif


!inquire(file='wffile',exist=wffile)

!inquire(file='ejenfile',exist=ejenfile)

print*,'myid111111=',myid

call readin
!call get_c8


call get_cross_sections

endif

call mpi_finalize(ierr)

stop
 30   format(i5,20000g14.6)
stop
end 
