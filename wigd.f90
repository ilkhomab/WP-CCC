module precision_type

integer, parameter :: p4 = selected_real_kind(4)
integer, parameter :: p8 = selected_real_kind(8)
integer, parameter :: p16 = selected_real_kind(8)
!integer, parameter :: id = selected_real_kind(16)
integer, parameter :: id = selected_real_kind(8)
integer, parameter :: mpoud=2*id
!integer, parameter :: mpoud_p16=16
integer, parameter :: mpoud_p16=32
real(kind=id) :: mppic
real(kind=p16) :: mppic16
real(kind=id),parameter :: hr=27.2116_id
real(kind=id),parameter ::pi=3.1415926535897932384626433832795028841971693993751_id
real(kind=id),parameter :: pifc=0.797884560802865355879892119869_id

end module precision_type

module mpi
implicit none 
integer::ngpu,myid,ierr,ntasks,nodes,nomp,nodeid,min_rho,max_rho 
end module mpi

!program test
!  
!  use mpmodule
!  implicit none
!  integer::l,m1,m2,ndp
!  type (mp_real)::wigd,t,ww

!  ndp=100
!  call mpinit (ndp)
!  mpoud=ndp
!
!  l=5;m1=-4;m2=-3
!
!  t=(1_id)
!
!  ww=wigd(l,m1,m2,t)
!
!  call mpwrite(6,ww)
!
!end program test
module acc_wig
  use precision_type
  implicit none

  contains

function wigd(l,m1,m2,t)
  
!$acc routine
  
  implicit none
  integer::l,m1,m2
  real(kind=id)::wigd,t

if(l==0.and.m1==0.and.m2==0) then

  wigd=1._id

endif

if(l==1.and.m1==-1.and.m2==-1) then

  wigd=Cos(t*0.5_id)**2

endif

if(l==1.and.m1==-1.and.m2==0) then

  wigd=Cos(t*0.5_id)*((-1.4142135623730950488016887242096980785696718753769_id))*Sin(t*0.5_id)

endif

if(l==1.and.m1==-1.and.m2==1) then

  wigd=Sin(t*0.5_id)**2

endif

if(l==1.and.m1==0.and.m2==-1) then

  wigd=Cos(t*0.5_id)*1.41421356237309504880168872420969807856967187537695_id*Sin(t*0.5_id)

endif

if(l==1.and.m1==0.and.m2==0) then

  wigd=Cos(t)

endif

if(l==1.and.m1==0.and.m2==1) then

  wigd=Cos(t*0.5_id)*((-1.4142135623730950488016887242096980785696718753769_id))*Sin(t*0.5_id)

endif

if(l==1.and.m1==1.and.m2==-1) then

  wigd=Sin(t*0.5_id)**2

endif

if(l==1.and.m1==1.and.m2==0) then

  wigd=Cos(t*0.5_id)*1.41421356237309504880168872420969807856967187537695_id*Sin(t*0.5_id)

endif

if(l==1.and.m1==1.and.m2==1) then

  wigd=Cos(t*0.5_id)**2

endif

if(l==2.and.m1==-2.and.m2==-2) then

  wigd=Cos(t*0.5_id)**4

endif

if(l==2.and.m1==-2.and.m2==-1) then

  wigd=Cos(t*0.5_id)**3*(-2._id)*Sin(t*0.5_id)

endif

if(l==2.and.m1==-2.and.m2==0) then

  wigd=Cos(t*0.5_id)**2*2.4494897427831780981972840747058913919659474806567_id*Sin(t*0.5_id)**2

endif

if(l==2.and.m1==-2.and.m2==1) then

  wigd=Cos(t*0.5_id)*(-2._id)*Sin(t*0.5_id)**3

endif

if(l==2.and.m1==-2.and.m2==2) then

  wigd=Sin(t*0.5_id)**4

endif

if(l==2.and.m1==-1.and.m2==-2) then

  wigd=Cos(t*0.5_id)**3*2._id*Sin(t*0.5_id)

endif

if(l==2.and.m1==-1.and.m2==-1) then

  wigd=Cos(t*0.5_id)**2*0.5_id*(-2._id + Cos(t)*4._id)

endif

if(l==2.and.m1==-1.and.m2==0) then

  wigd=Cos(t)*Cos(t*0.5_id)*(-2.4494897427831780981972840747058913919659474806567_id)*Sin(t*0.5_id)

endif

if(l==2.and.m1==-1.and.m2==1) then

  wigd=0.5_id*(2._id + Cos(t)*4._id)*Sin(t*0.5_id)**2

endif

if(l==2.and.m1==-1.and.m2==2) then

  wigd=Cos(t*0.5_id)*(-2._id)*Sin(t*0.5_id)**3

endif

if(l==2.and.m1==0.and.m2==-2) then

  wigd=Cos(t*0.5_id)**2*2.4494897427831780981972840747058913919659474806567_id*Sin(t*0.5_id)**2

endif

if(l==2.and.m1==0.and.m2==-1) then

  wigd=Cos(t)*Cos(t*0.5_id)*2.4494897427831780981972840747058913919659474806567_id*Sin(t*0.5_id)

endif

if(l==2.and.m1==0.and.m2==0) then

  wigd=0.5_id*(-1._id + Cos(t)**2*3._id)

endif

if(l==2.and.m1==0.and.m2==1) then

  wigd=Cos(t)*Cos(t*0.5_id)*(-2.4494897427831780981972840747058913919659474806567_id)*Sin(t*0.5_id)

endif

if(l==2.and.m1==0.and.m2==2) then

  wigd=Cos(t*0.5_id)**2*2.4494897427831780981972840747058913919659474806567_id*Sin(t*0.5_id)**2

endif

if(l==2.and.m1==1.and.m2==-2) then

  wigd=Cos(t*0.5_id)*2._id*Sin(t*0.5_id)**3

endif

if(l==2.and.m1==1.and.m2==-1) then

  wigd=0.5_id*(2._id + Cos(t)*4._id)*Sin(t*0.5_id)**2

endif

if(l==2.and.m1==1.and.m2==0) then

  wigd=Cos(t)*Cos(t*0.5_id)*2.4494897427831780981972840747058913919659474806567_id*Sin(t*0.5_id)

endif

if(l==2.and.m1==1.and.m2==1) then

  wigd=Cos(t*0.5_id)**2*0.5_id*(-2._id + Cos(t)*4._id)

endif

if(l==2.and.m1==1.and.m2==2) then

  wigd=Cos(t*0.5_id)**3*(-2._id)*Sin(t*0.5_id)

endif

if(l==2.and.m1==2.and.m2==-2) then

  wigd=Sin(t*0.5_id)**4

endif

if(l==2.and.m1==2.and.m2==-1) then

  wigd=Cos(t*0.5_id)*2._id*Sin(t*0.5_id)**3

endif

if(l==2.and.m1==2.and.m2==0) then

  wigd=Cos(t*0.5_id)**2*2.4494897427831780981972840747058913919659474806567_id*Sin(t*0.5_id)**2

endif

if(l==2.and.m1==2.and.m2==1) then

  wigd=Cos(t*0.5_id)**3*2._id*Sin(t*0.5_id)

endif

if(l==2.and.m1==2.and.m2==2) then

  wigd=Cos(t*0.5_id)**4

endif

if(l==3.and.m1==-3.and.m2==-3) then

  wigd=Cos(t*0.5_id)**6

endif

if(l==3.and.m1==-3.and.m2==-2) then

  wigd=Cos(t*0.5_id)**5*(-2.4494897427831780981972840747058913919659474806567_id)*Sin(t*0.5_id)

endif

if(l==3.and.m1==-3.and.m2==-1) then

  wigd=Cos(t*0.5_id)**4*3.8729833462074168851792653997823996108329217052916_id*Sin(t*0.5_id)**2

endif

if(l==3.and.m1==-3.and.m2==0) then

  wigd=Cos(t*0.5_id)**3*(-4.4721359549995793928183473374625524708812367192231_id)*Sin(t*0.5_id)**3

endif

if(l==3.and.m1==-3.and.m2==1) then

  wigd=Cos(t*0.5_id)**2*3.8729833462074168851792653997823996108329217052916_id*Sin(t*0.5_id)**4

endif

if(l==3.and.m1==-3.and.m2==2) then

  wigd=Cos(t*0.5_id)*(-2.4494897427831780981972840747058913919659474806567_id)*Sin(t*0.5_id)**5

endif

if(l==3.and.m1==-3.and.m2==3) then

  wigd=Sin(t*0.5_id)**6

endif

if(l==3.and.m1==-2.and.m2==-3) then

  wigd=Cos(t*0.5_id)**5*2.4494897427831780981972840747058913919659474806567_id*Sin(t*0.5_id)

endif

if(l==3.and.m1==-2.and.m2==-2) then

  wigd=Cos(t*0.5_id)**4*0.5_id*(-4._id + Cos(t)*6._id)

endif

if(l==3.and.m1==-2.and.m2==-1) then

  wigd=Cos(t*0.5_id)**3*(-0.7905694150420948329997233861081796334298887848313_id)*(-2._id + Cos(t)*6._id)*Sin(t*0.5_id)

endif

if(l==3.and.m1==-2.and.m2==0) then

  wigd=Cos(t)*Cos(t*0.5_id)**2*5.4772255750516611345696978280080213395274469499798_id*Sin(t*0.5_id)**2

endif

if(l==3.and.m1==-2.and.m2==1) then

  wigd=Cos(t*0.5_id)*(-0.7905694150420948329997233861081796334298887848313_id)*(2._id + Cos(t)*6._id)*Sin(t*0.5_id)**3

endif

if(l==3.and.m1==-2.and.m2==2) then

  wigd=0.5_id*(4._id + Cos(t)*6._id)*Sin(t*0.5_id)**4

endif

if(l==3.and.m1==-2.and.m2==3) then

  wigd=Cos(t*0.5_id)*(-2.4494897427831780981972840747058913919659474806567_id)*Sin(t*0.5_id)**5

endif

if(l==3.and.m1==-1.and.m2==-3) then

  wigd=Cos(t*0.5_id)**4*3.8729833462074168851792653997823996108329217052916_id*Sin(t*0.5_id)**2

endif

if(l==3.and.m1==-1.and.m2==-2) then

  wigd=Cos(t*0.5_id)**3*0.7905694150420948329997233861081796334298887848313_id*(-2._id + Cos(t)*6._id)*Sin(t*0.5_id)

endif

if(l==3.and.m1==-1.and.m2==-1) then

  wigd=Cos(t*0.5_id)**2*(1._id + (Cos(t) -1._id)**2*3.75_id + (Cos(t) -1._id)*5._id)

endif

if(l==3.and.m1==-1.and.m2==0) then

  wigd=Cos(t*0.5_id)*(-1.1547005383792515290182975610039149112952035025403_id)*(3._id + (Cos(t) -1._id)**2*3.75_id + (Cos(t) -1._id)*7.5_id)*Sin(t*0.5_id)

endif

if(l==3.and.m1==-1.and.m2==1) then

  wigd=((Cos(t) -1._id)**2*3.75_id + 6._id + (Cos(t) -1._id)*10._id)*Sin(t*0.5_id)**2

endif

if(l==3.and.m1==-1.and.m2==2) then

  wigd=Cos(t*0.5_id)*(-0.7905694150420948329997233861081796334298887848313_id)*(2._id + Cos(t)*6._id)*Sin(t*0.5_id)**3

endif

if(l==3.and.m1==-1.and.m2==3) then

  wigd=Cos(t*0.5_id)**2*3.8729833462074168851792653997823996108329217052916_id*Sin(t*0.5_id)**4

endif

if(l==3.and.m1==0.and.m2==-3) then

  wigd=Cos(t*0.5_id)**3*4.4721359549995793928183473374625524708812367192231_id*Sin(t*0.5_id)**3

endif

if(l==3.and.m1==0.and.m2==-2) then

  wigd=Cos(t)*Cos(t*0.5_id)**2*5.4772255750516611345696978280080213395274469499798_id*Sin(t*0.5_id)**2

endif

if(l==3.and.m1==0.and.m2==-1) then

  wigd=Cos(t*0.5_id)*1.1547005383792515290182975610039149112952035025403_id*(3._id + (Cos(t) -1._id)**2*3.75_id + (Cos(t) -1._id)*7.5_id)*Sin(t*0.5_id)

endif

if(l==3.and.m1==0.and.m2==0) then

  wigd=0.5_id*(Cos(t)*(-3._id) + Cos(t)**3*5._id)

endif

if(l==3.and.m1==0.and.m2==1) then

  wigd=Cos(t*0.5_id)*(-1.1547005383792515290182975610039149112952035025403_id)*(3._id + (Cos(t) -1._id)**2*3.75_id + (Cos(t) -1._id)*7.5_id)*Sin(t*0.5_id)

endif

if(l==3.and.m1==0.and.m2==2) then

  wigd=Cos(t)*Cos(t*0.5_id)**2*5.4772255750516611345696978280080213395274469499798_id*Sin(t*0.5_id)**2

endif

if(l==3.and.m1==0.and.m2==3) then

  wigd=Cos(t*0.5_id)**3*(-4.4721359549995793928183473374625524708812367192231_id)*Sin(t*0.5_id)**3

endif

if(l==3.and.m1==1.and.m2==-3) then

  wigd=Cos(t*0.5_id)**2*3.8729833462074168851792653997823996108329217052916_id*Sin(t*0.5_id)**4

endif

if(l==3.and.m1==1.and.m2==-2) then

  wigd=Cos(t*0.5_id)*0.7905694150420948329997233861081796334298887848313_id*(2._id + Cos(t)*6._id)*Sin(t*0.5_id)**3

endif

if(l==3.and.m1==1.and.m2==-1) then

  wigd=((Cos(t) -1._id)**2*3.75_id + 6._id + (Cos(t) -1._id)*10._id)*Sin(t*0.5_id)**2

endif

if(l==3.and.m1==1.and.m2==0) then

  wigd=Cos(t*0.5_id)*1.1547005383792515290182975610039149112952035025403_id*(3._id + (Cos(t) -1._id)**2*3.75_id + (Cos(t) -1._id)*7.5_id)*Sin(t*0.5_id)

endif

if(l==3.and.m1==1.and.m2==1) then

  wigd=Cos(t*0.5_id)**2*(1._id + (Cos(t) -1._id)**2*3.75_id + (Cos(t) -1._id)*5._id)

endif

if(l==3.and.m1==1.and.m2==2) then

  wigd=Cos(t*0.5_id)**3*(-0.7905694150420948329997233861081796334298887848313_id)*(-2._id + Cos(t)*6._id)*Sin(t*0.5_id)

endif

if(l==3.and.m1==1.and.m2==3) then

  wigd=Cos(t*0.5_id)**4*3.8729833462074168851792653997823996108329217052916_id*Sin(t*0.5_id)**2

endif

if(l==3.and.m1==2.and.m2==-3) then

  wigd=Cos(t*0.5_id)*2.4494897427831780981972840747058913919659474806567_id*Sin(t*0.5_id)**5

endif

if(l==3.and.m1==2.and.m2==-2) then

  wigd=0.5_id*(4._id + Cos(t)*6._id)*Sin(t*0.5_id)**4

endif

if(l==3.and.m1==2.and.m2==-1) then

  wigd=Cos(t*0.5_id)*0.7905694150420948329997233861081796334298887848313_id*(2._id + Cos(t)*6._id)*Sin(t*0.5_id)**3

endif

if(l==3.and.m1==2.and.m2==0) then

  wigd=Cos(t)*Cos(t*0.5_id)**2*5.4772255750516611345696978280080213395274469499798_id*Sin(t*0.5_id)**2

endif

if(l==3.and.m1==2.and.m2==1) then

  wigd=Cos(t*0.5_id)**3*0.7905694150420948329997233861081796334298887848313_id*(-2._id + Cos(t)*6._id)*Sin(t*0.5_id)

endif

if(l==3.and.m1==2.and.m2==2) then

  wigd=Cos(t*0.5_id)**4*0.5_id*(-4._id + Cos(t)*6._id)

endif

if(l==3.and.m1==2.and.m2==3) then

  wigd=Cos(t*0.5_id)**5*(-2.4494897427831780981972840747058913919659474806567_id)*Sin(t*0.5_id)

endif

if(l==3.and.m1==3.and.m2==-3) then

  wigd=Sin(t*0.5_id)**6

endif

if(l==3.and.m1==3.and.m2==-2) then

  wigd=Cos(t*0.5_id)*2.4494897427831780981972840747058913919659474806567_id*Sin(t*0.5_id)**5

endif

if(l==3.and.m1==3.and.m2==-1) then

  wigd=Cos(t*0.5_id)**2*3.8729833462074168851792653997823996108329217052916_id*Sin(t*0.5_id)**4

endif

if(l==3.and.m1==3.and.m2==0) then

  wigd=Cos(t*0.5_id)**3*4.4721359549995793928183473374625524708812367192231_id*Sin(t*0.5_id)**3

endif

if(l==3.and.m1==3.and.m2==1) then

  wigd=Cos(t*0.5_id)**4*3.8729833462074168851792653997823996108329217052916_id*Sin(t*0.5_id)**2

endif

if(l==3.and.m1==3.and.m2==2) then

  wigd=Cos(t*0.5_id)**5*2.4494897427831780981972840747058913919659474806567_id*Sin(t*0.5_id)

endif

if(l==3.and.m1==3.and.m2==3) then

  wigd=Cos(t*0.5_id)**6

endif

if(l==4.and.m1==-4.and.m2==-4) then

  wigd=Cos(t*0.5_id)**8

endif

if(l==4.and.m1==-4.and.m2==-3) then

  wigd=Cos(t*0.5_id)**7*((-2.8284271247461900976033774484193961571393437507539_id))*Sin(t*0.5_id)

endif

if(l==4.and.m1==-4.and.m2==-2) then

  wigd=Cos(t*0.5_id)**6*5.2915026221291811810032315072785208514205183661649_id*Sin(t*0.5_id)**2

endif

if(l==4.and.m1==-4.and.m2==-1) then

  wigd=Cos(t*0.5_id)**5*((-7.483314773547882771167497464633098603512039615557_id))*Sin(t*0.5_id)**3

endif

if(l==4.and.m1==-4.and.m2==0) then

  wigd=Cos(t*0.5_id)**4*8.3666002653407554797817202578518748939281536929867_id*Sin(t*0.5_id)**4

endif

if(l==4.and.m1==-4.and.m2==1) then

  wigd=Cos(t*0.5_id)**3*((-7.483314773547882771167497464633098603512039615557_id))*Sin(t*0.5_id)**5

endif

if(l==4.and.m1==-4.and.m2==2) then

  wigd=Cos(t*0.5_id)**2*5.2915026221291811810032315072785208514205183661649_id*Sin(t*0.5_id)**6

endif

if(l==4.and.m1==-4.and.m2==3) then

  wigd=Cos(t*0.5_id)*((-2.8284271247461900976033774484193961571393437507539_id))*Sin(t*0.5_id)**7

endif

if(l==4.and.m1==-4.and.m2==4) then

  wigd=Sin(t*0.5_id)**8

endif

if(l==4.and.m1==-3.and.m2==-4) then

  wigd=Cos(t*0.5_id)**7*2.8284271247461900976033774484193961571393437507539_id*Sin(t*0.5_id)

endif

if(l==4.and.m1==-3.and.m2==-3) then

  wigd=Cos(t*0.5_id)**6*0.5_id*(-6._id + Cos(t)*8._id)

endif

if(l==4.and.m1==-3.and.m2==-2) then

  wigd=Cos(t*0.5_id)**5*((-0.9354143466934853463959371830791373254390049519447_id))*(-4._id + Cos(t)*8._id)*Sin(t*0.5_id)

endif

if(l==4.and.m1==-3.and.m2==-1) then

  wigd=Cos(t*0.5_id)**4*1.3228756555322952952508078768196302128551295915412_id*(-2._id + Cos(t)*8._id)*Sin(t*0.5_id)**2

endif

if(l==4.and.m1==-3.and.m2==0) then

  wigd=Cos(t)*Cos(t*0.5_id)**3*((-11.832159566199232085134656583123234096831002461589_id))*Sin(t*0.5_id)**3

endif

if(l==4.and.m1==-3.and.m2==1) then

  wigd=Cos(t*0.5_id)**2*1.3228756555322952952508078768196302128551295915412_id*(2._id + Cos(t)*8._id)*Sin(t*0.5_id)**4

endif

if(l==4.and.m1==-3.and.m2==2) then

  wigd=Cos(t*0.5_id)*(-0.9354143466934853463959371830791373254390049519447_id)*(4._id + Cos(t)*8._id)*Sin(t*0.5_id)**5

endif

if(l==4.and.m1==-3.and.m2==3) then

  wigd=0.5_id*(6._id + Cos(t)*8._id)*Sin(t*0.5_id)**6

endif

if(l==4.and.m1==-3.and.m2==4) then

  wigd=Cos(t*0.5_id)*((-2.8284271247461900976033774484193961571393437507539_id))*Sin(t*0.5_id)**7

endif

if(l==4.and.m1==-2.and.m2==-4) then

  wigd=Cos(t*0.5_id)**6*5.2915026221291811810032315072785208514205183661649_id*Sin(t*0.5_id)**2

endif

if(l==4.and.m1==-2.and.m2==-3) then

  wigd=Cos(t*0.5_id)**5*0.9354143466934853463959371830791373254390049519447_id*(-4._id + Cos(t)*8._id)*Sin(t*0.5_id)

endif

if(l==4.and.m1==-2.and.m2==-2) then

  wigd=Cos(t*0.5_id)**4*(1._id + (Cos(t) -1._id)*7._id + (Cos(t) -1._id)**2*7._id)

endif

if(l==4.and.m1==-2.and.m2==-1) then

  wigd=Cos(t*0.5_id)**3*((-1.4142135623730950488016887242096980785696718753769_id))*(3._id + (Cos(t) -1._id)**2*7._id + (Cos(t) -1._id)*10.5_id)*Sin(t*0.5_id)

endif

if(l==4.and.m1==-2.and.m2==0) then

  wigd=Cos(t*0.5_id)**2*1.58113883008418966599944677221635926685977756966261_id*(6._id + (Cos(t) -1._id)**2*7._id + (Cos(t) -1._id)*14._id)*Sin(t*0.5_id)**2

endif

if(l==4.and.m1==-2.and.m2==1) then

  wigd=Cos(t*0.5_id)*((-1.4142135623730950488016887242096980785696718753769_id))*((Cos(t) -1._id)**2*7._id + 10._id + (Cos(t) -1._id)*17.5_id)*Sin(t*0.5_id)**3

endif

if(l==4.and.m1==-2.and.m2==2) then

  wigd=((Cos(t) -1._id)**2*7._id + 15._id + (Cos(t) -1._id)*21._id)*Sin(t*0.5_id)**4

endif

if(l==4.and.m1==-2.and.m2==3) then

  wigd=Cos(t*0.5_id)*((-0.9354143466934853463959371830791373254390049519447_id))*(4._id + Cos(t)*8._id)*Sin(t*0.5_id)**5

endif

if(l==4.and.m1==-2.and.m2==4) then

  wigd=Cos(t*0.5_id)**2*5.2915026221291811810032315072785208514205183661649_id*Sin(t*0.5_id)**6

endif

if(l==4.and.m1==-1.and.m2==-4) then

  wigd=Cos(t*0.5_id)**5*7.483314773547882771167497464633098603512039615557_id*Sin(t*0.5_id)**3

endif

if(l==4.and.m1==-1.and.m2==-3) then

  wigd=Cos(t*0.5_id)**4*1.3228756555322952952508078768196302128551295915412_id*(-2._id + Cos(t)*8._id)*Sin(t*0.5_id)**2

endif

if(l==4.and.m1==-1.and.m2==-2) then

  wigd=Cos(t*0.5_id)**3*1.41421356237309504880168872420969807856967187537695_id*(3._id + (Cos(t) -1._id)**2*7._id + (Cos(t) -1._id)*10.5_id)*Sin(t*0.5_id)

endif

if(l==4.and.m1==-1.and.m2==-1) then

  wigd=Cos(t*0.5_id)**2*(1._id + (Cos(t) -1._id)**3*7._id + (Cos(t) -1._id)*9._id + (Cos(t) -1._id)**2*15.75_id)

endif

if(l==4.and.m1==-1.and.m2==0) then

  wigd=Cos(t*0.5_id)*(-1.1180339887498948482045868343656381177203091798058_id)*(4._id + (Cos(t) -1._id)**3*7._id + (Cos(t) -1._id)*18._id + (Cos(t) -1._id)**2*21._id)*Sin(t*0.5_id)

endif

if(l==4.and.m1==-1.and.m2==1) then

  wigd=((Cos(t) -1._id)**3*7._id + 10._id + (Cos(t) -1._id)**2*26.25_id + (Cos(t) -1._id)*30._id)*Sin(t*0.5_id)**2

endif

if(l==4.and.m1==-1.and.m2==2) then

  wigd=Cos(t*0.5_id)*(-1.4142135623730950488016887242096980785696718753769_id)*((Cos(t) -1._id)**2*7._id + 10._id + (Cos(t) -1._id)*17.5_id)*Sin(t*0.5_id)**3

endif

if(l==4.and.m1==-1.and.m2==3) then

  wigd=Cos(t*0.5_id)**2*1.3228756555322952952508078768196302128551295915412_id*(2._id + Cos(t)*8._id)*Sin(t*0.5_id)**4

endif

if(l==4.and.m1==-1.and.m2==4) then

  wigd=Cos(t*0.5_id)**3*(-7.483314773547882771167497464633098603512039615557_id)*Sin(t*0.5_id)**5

endif

if(l==4.and.m1==0.and.m2==-4) then

  wigd=Cos(t*0.5_id)**4*8.3666002653407554797817202578518748939281536929867_id*Sin(t*0.5_id)**4

endif

if(l==4.and.m1==0.and.m2==-3) then

  wigd=Cos(t)*Cos(t*0.5_id)**3*11.832159566199232085134656583123234096831002461589_id*Sin(t*0.5_id)**3

endif

if(l==4.and.m1==0.and.m2==-2) then

  wigd=Cos(t*0.5_id)**2*1.58113883008418966599944677221635926685977756966261_id*(6._id + (Cos(t) -1._id)**2*7._id + (Cos(t) -1._id)*14._id)*Sin(t*0.5_id)**2

endif

if(l==4.and.m1==0.and.m2==-1) then

  wigd=Cos(t*0.5_id)*1.1180339887498948482045868343656381177203091798058_id*(4._id + (Cos(t) -1._id)**3*7._id + (Cos(t) -1._id)*18._id + (Cos(t) -1._id)**2*21._id)*Sin(t*0.5_id)

endif

if(l==4.and.m1==0.and.m2==0) then

  wigd=0.125_id*(Cos(t)**2*(-30._id) + 3._id + Cos(t)**4*35._id)

endif

if(l==4.and.m1==0.and.m2==1) then

  wigd=Cos(t*0.5_id)*(-1.1180339887498948482045868343656381177203091798058_id)*(4._id + (Cos(t) -1._id)**3*7._id + (Cos(t) -1._id)*18._id + (Cos(t) -1._id)**2*21._id)*Sin(t*0.5_id)

endif

if(l==4.and.m1==0.and.m2==2) then

  wigd=Cos(t*0.5_id)**2*1.58113883008418966599944677221635926685977756966261_id*(6._id + (Cos(t) -1._id)**2*7._id + (Cos(t) -1._id)*14._id)*Sin(t*0.5_id)**2

endif

if(l==4.and.m1==0.and.m2==3) then

  wigd=Cos(t)*Cos(t*0.5_id)**3*(-11.832159566199232085134656583123234096831002461589_id)*Sin(t*0.5_id)**3

endif

if(l==4.and.m1==0.and.m2==4) then

  wigd=Cos(t*0.5_id)**4*8.3666002653407554797817202578518748939281536929867_id*Sin(t*0.5_id)**4

endif

if(l==4.and.m1==1.and.m2==-4) then

  wigd=Cos(t*0.5_id)**3*7.483314773547882771167497464633098603512039615557_id*Sin(t*0.5_id)**5

endif

if(l==4.and.m1==1.and.m2==-3) then

  wigd=Cos(t*0.5_id)**2*1.3228756555322952952508078768196302128551295915412_id*(2._id + Cos(t)*8._id)*Sin(t*0.5_id)**4

endif

if(l==4.and.m1==1.and.m2==-2) then

  wigd=Cos(t*0.5_id)*1.41421356237309504880168872420969807856967187537695_id*((Cos(t) -1._id)**2*7._id + 10._id + (Cos(t) -1._id)*17.5_id)*Sin(t*0.5_id)**3

endif

if(l==4.and.m1==1.and.m2==-1) then

  wigd=((Cos(t) -1._id)**3*7._id + 10._id + (Cos(t) -1._id)**2*26.25_id + (Cos(t) -1._id)*30._id)*Sin(t*0.5_id)**2

endif

if(l==4.and.m1==1.and.m2==0) then

  wigd=Cos(t*0.5_id)*1.1180339887498948482045868343656381177203091798058_id*(4._id + (Cos(t) -1._id)**3*7._id + (Cos(t) -1._id)*18._id + (Cos(t) -1._id)**2*21._id)*Sin(t*0.5_id)

endif

if(l==4.and.m1==1.and.m2==1) then

  wigd=Cos(t*0.5_id)**2*(1._id + (Cos(t) -1._id)**3*7._id + (Cos(t) -1._id)*9._id + (Cos(t) -1._id)**2*15.75_id)

endif

if(l==4.and.m1==1.and.m2==2) then

  wigd=Cos(t*0.5_id)**3*(-1.4142135623730950488016887242096980785696718753769_id)*(3._id + (Cos(t) -1._id)**2*7._id + (Cos(t) -1._id)*10.5_id)*Sin(t*0.5_id)

endif

if(l==4.and.m1==1.and.m2==3) then

  wigd=Cos(t*0.5_id)**4*1.3228756555322952952508078768196302128551295915412_id*(-2._id + Cos(t)*8._id)*Sin(t*0.5_id)**2

endif

if(l==4.and.m1==1.and.m2==4) then

  wigd=Cos(t*0.5_id)**5*(-7.483314773547882771167497464633098603512039615557_id)*Sin(t*0.5_id)**3

endif

if(l==4.and.m1==2.and.m2==-4) then

  wigd=Cos(t*0.5_id)**2*5.2915026221291811810032315072785208514205183661649_id*Sin(t*0.5_id)**6

endif

if(l==4.and.m1==2.and.m2==-3) then

  wigd=Cos(t*0.5_id)*0.9354143466934853463959371830791373254390049519447_id*(4._id + Cos(t)*8._id)*Sin(t*0.5_id)**5

endif

if(l==4.and.m1==2.and.m2==-2) then

  wigd=((Cos(t) -1._id)**2*7._id + 15._id + (Cos(t) -1._id)*21._id)*Sin(t*0.5_id)**4

endif

if(l==4.and.m1==2.and.m2==-1) then

  wigd=Cos(t*0.5_id)*1.41421356237309504880168872420969807856967187537695_id*((Cos(t) -1._id)**2*7._id + 10._id + (Cos(t) -1._id)*17.5_id)*Sin(t*0.5_id)**3

endif

if(l==4.and.m1==2.and.m2==0) then

  wigd=Cos(t*0.5_id)**2*1.58113883008418966599944677221635926685977756966261_id*(6._id + (Cos(t) -1._id)**2*7._id + (Cos(t) -1._id)*14._id)*Sin(t*0.5_id)**2

endif

if(l==4.and.m1==2.and.m2==1) then

  wigd=Cos(t*0.5_id)**3*1.41421356237309504880168872420969807856967187537695_id*(3._id + (Cos(t) -1._id)**2*7._id + (Cos(t) -1._id)*10.5_id)*Sin(t*0.5_id)

endif

if(l==4.and.m1==2.and.m2==2) then

  wigd=Cos(t*0.5_id)**4*(1._id + (Cos(t) -1._id)*7._id + (Cos(t) -1._id)**2*7._id)

endif

if(l==4.and.m1==2.and.m2==3) then

  wigd=Cos(t*0.5_id)**5*(-0.9354143466934853463959371830791373254390049519447_id)*(-4._id + Cos(t)*8._id)*Sin(t*0.5_id)

endif

if(l==4.and.m1==2.and.m2==4) then

  wigd=Cos(t*0.5_id)**6*5.2915026221291811810032315072785208514205183661649_id*Sin(t*0.5_id)**2

endif

if(l==4.and.m1==3.and.m2==-4) then

  wigd=Cos(t*0.5_id)*2.8284271247461900976033774484193961571393437507539_id*Sin(t*0.5_id)**7

endif

if(l==4.and.m1==3.and.m2==-3) then

  wigd=0.5_id*(6._id + Cos(t)*8._id)*Sin(t*0.5_id)**6

endif

if(l==4.and.m1==3.and.m2==-2) then

  wigd=Cos(t*0.5_id)*0.9354143466934853463959371830791373254390049519447_id*(4._id + Cos(t)*8._id)*Sin(t*0.5_id)**5

endif

if(l==4.and.m1==3.and.m2==-1) then

  wigd=Cos(t*0.5_id)**2*1.3228756555322952952508078768196302128551295915412_id*(2._id + Cos(t)*8._id)*Sin(t*0.5_id)**4

endif

if(l==4.and.m1==3.and.m2==0) then

  wigd=Cos(t)*Cos(t*0.5_id)**3*11.832159566199232085134656583123234096831002461589_id*Sin(t*0.5_id)**3

endif

if(l==4.and.m1==3.and.m2==1) then

  wigd=Cos(t*0.5_id)**4*1.3228756555322952952508078768196302128551295915412_id*(-2._id + Cos(t)*8._id)*Sin(t*0.5_id)**2

endif

if(l==4.and.m1==3.and.m2==2) then

  wigd=Cos(t*0.5_id)**5*0.9354143466934853463959371830791373254390049519447_id*(-4._id + Cos(t)*8._id)*Sin(t*0.5_id)

endif

if(l==4.and.m1==3.and.m2==3) then

  wigd=Cos(t*0.5_id)**6*0.5_id*(-6._id + Cos(t)*8._id)

endif

if(l==4.and.m1==3.and.m2==4) then

  wigd=Cos(t*0.5_id)**7*(-2.8284271247461900976033774484193961571393437507539_id)*Sin(t*0.5_id)

endif

if(l==4.and.m1==4.and.m2==-4) then

  wigd=Sin(t*0.5_id)**8

endif

if(l==4.and.m1==4.and.m2==-3) then

  wigd=Cos(t*0.5_id)*2.8284271247461900976033774484193961571393437507539_id*Sin(t*0.5_id)**7

endif

if(l==4.and.m1==4.and.m2==-2) then

  wigd=Cos(t*0.5_id)**2*5.2915026221291811810032315072785208514205183661649_id*Sin(t*0.5_id)**6

endif

if(l==4.and.m1==4.and.m2==-1) then

  wigd=Cos(t*0.5_id)**3*7.483314773547882771167497464633098603512039615557_id*Sin(t*0.5_id)**5

endif

if(l==4.and.m1==4.and.m2==0) then

  wigd=Cos(t*0.5_id)**4*8.3666002653407554797817202578518748939281536929867_id*Sin(t*0.5_id)**4

endif

if(l==4.and.m1==4.and.m2==1) then

  wigd=Cos(t*0.5_id)**5*7.483314773547882771167497464633098603512039615557_id*Sin(t*0.5_id)**3

endif

if(l==4.and.m1==4.and.m2==2) then

  wigd=Cos(t*0.5_id)**6*5.2915026221291811810032315072785208514205183661649_id*Sin(t*0.5_id)**2

endif

if(l==4.and.m1==4.and.m2==3) then

  wigd=Cos(t*0.5_id)**7*2.8284271247461900976033774484193961571393437507539_id*Sin(t*0.5_id)

endif

if(l==4.and.m1==4.and.m2==4) then

  wigd=Cos(t*0.5_id)**8

endif

if(l==5.and.m1==-5.and.m2==-5) then

  wigd=Cos(t*0.5_id)**10

endif

if(l==5.and.m1==-5.and.m2==-4) then

  wigd=Cos(t*0.5_id)**9*(-3.1622776601683793319988935444327185337195551393252_id)*Sin(t*0.5_id)

endif

if(l==5.and.m1==-5.and.m2==-3) then

  wigd=Cos(t*0.5_id)**8*6.708203932499369089227521006193828706321855078835_id*Sin(t*0.5_id)**2

endif

if(l==5.and.m1==-5.and.m2==-2) then

  wigd=Cos(t*0.5_id)**7*(-10.95445115010332226913939565601604267905489389996_id)*Sin(t*0.5_id)**3

endif

if(l==5.and.m1==-5.and.m2==-1) then

  wigd=Cos(t*0.5_id)**6*14.491376746189438573718664157169771723140132874759_id*Sin(t*0.5_id)**4

endif

if(l==5.and.m1==-5.and.m2==0) then

  wigd=Cos(t*0.5_id)**5*(-15.874507866387543543009694521835562554261555098495_id)*Sin(t*0.5_id)**5

endif

if(l==5.and.m1==-5.and.m2==1) then

  wigd=Cos(t*0.5_id)**4*14.491376746189438573718664157169771723140132874759_id*Sin(t*0.5_id)**6

endif

if(l==5.and.m1==-5.and.m2==2) then

  wigd=Cos(t*0.5_id)**3*(-10.95445115010332226913939565601604267905489389996_id)*Sin(t*0.5_id)**7

endif

if(l==5.and.m1==-5.and.m2==3) then

  wigd=Cos(t*0.5_id)**2*6.708203932499369089227521006193828706321855078835_id*Sin(t*0.5_id)**8

endif

if(l==5.and.m1==-5.and.m2==4) then

  wigd=Cos(t*0.5_id)*(-3.1622776601683793319988935444327185337195551393252_id)*Sin(t*0.5_id)**9

endif

if(l==5.and.m1==-5.and.m2==5) then

  wigd=Sin(t*0.5_id)**10

endif

if(l==5.and.m1==-4.and.m2==-5) then

  wigd=Cos(t*0.5_id)**9*3.1622776601683793319988935444327185337195551393252_id*Sin(t*0.5_id)

endif

if(l==5.and.m1==-4.and.m2==-4) then

  wigd=Cos(t*0.5_id)**8*0.5_id*(-8._id + Cos(t)*10._id)

endif

if(l==5.and.m1==-4.and.m2==-3) then

  wigd=Cos(t*0.5_id)**7*(-1.0606601717798212866012665431572735589272539065327_id)*(-6._id + Cos(t)*10._id)*Sin(t*0.5_id)

endif

if(l==5.and.m1==-4.and.m2==-2) then

  wigd=Cos(t*0.5_id)**6*1.73205080756887729352744634150587236694280525381038_id*(-4._id + Cos(t)*10._id)*Sin(t*0.5_id)**2

endif

if(l==5.and.m1==-4.and.m2==-1) then

  wigd=Cos(t*0.5_id)**5*(-2.291287847477920003294023596864004244492228288384_id)*(-2._id + Cos(t)*10._id)*Sin(t*0.5_id)**3

endif

if(l==5.and.m1==-4.and.m2==0) then

  wigd=Cos(t)*Cos(t*0.5_id)**4*25.09980079602226643934516077355562468178446107896_id*Sin(t*0.5_id)**4

endif

if(l==5.and.m1==-4.and.m2==1) then

  wigd=Cos(t*0.5_id)**3*(-2.291287847477920003294023596864004244492228288384_id)*(2._id + Cos(t)*10._id)*Sin(t*0.5_id)**5

endif

if(l==5.and.m1==-4.and.m2==2) then

  wigd=Cos(t*0.5_id)**2*1.73205080756887729352744634150587236694280525381038_id*(4._id + Cos(t)*10._id)*Sin(t*0.5_id)**6

endif

if(l==5.and.m1==-4.and.m2==3) then

  wigd=Cos(t*0.5_id)*(-1.0606601717798212866012665431572735589272539065327_id)*(6._id + Cos(t)*10._id)*Sin(t*0.5_id)**7

endif

if(l==5.and.m1==-4.and.m2==4) then

  wigd=0.5_id*(8._id + Cos(t)*10._id)*Sin(t*0.5_id)**8

endif

if(l==5.and.m1==-4.and.m2==5) then

  wigd=Cos(t*0.5_id)*(-3.1622776601683793319988935444327185337195551393252_id)*Sin(t*0.5_id)**9

endif

if(l==5.and.m1==-3.and.m2==-5) then

  wigd=Cos(t*0.5_id)**8*6.708203932499369089227521006193828706321855078835_id*Sin(t*0.5_id)**2

endif

if(l==5.and.m1==-3.and.m2==-4) then

  wigd=Cos(t*0.5_id)**7*1.0606601717798212866012665431572735589272539065327_id*(-6._id + Cos(t)*10._id)*Sin(t*0.5_id)

endif

if(l==5.and.m1==-3.and.m2==-3) then

  wigd=Cos(t*0.5_id)**6*(1._id + (Cos(t) -1._id)*9._id + (Cos(t) -1._id)**2*11.25_id)

endif

if(l==5.and.m1==-3.and.m2==-2) then

  wigd=Cos(t*0.5_id)**5*(-1.6329931618554520654648560498039275946439649871044_id)*(3._id + (Cos(t) -1._id)**2*11.25_id + (Cos(t) -1._id)*13.5_id)*Sin(t*0.5_id)

endif

if(l==5.and.m1==-3.and.m2==-1) then

  wigd=Cos(t*0.5_id)**4*2.1602468994692867436553224786959988859017347690194_id*(6._id + (Cos(t) -1._id)**2*11.25_id + (Cos(t) -1._id)*18._id)*Sin(t*0.5_id)**2

endif

if(l==5.and.m1==-3.and.m2==0) then

  wigd=Cos(t*0.5_id)**3*(-2.3664319132398464170269313166246468193662004923177_id)*(10._id + (Cos(t) -1._id)**2*11.25_id + (Cos(t) -1._id)*22.5_id)*Sin(t*0.5_id)**3

endif

if(l==5.and.m1==-3.and.m2==1) then

  wigd=Cos(t*0.5_id)**2*2.1602468994692867436553224786959988859017347690194_id*((Cos(t) -1._id)**2*11.25_id + 15._id + (Cos(t) -1._id)*27._id)*Sin(t*0.5_id)**4

endif

if(l==5.and.m1==-3.and.m2==2) then

  wigd=Cos(t*0.5_id)*(-1.6329931618554520654648560498039275946439649871044_id)*((Cos(t) -1._id)**2*11.25_id + 21._id + (Cos(t) -1._id)*31.5_id)*Sin(t*0.5_id)**5

endif

if(l==5.and.m1==-3.and.m2==3) then

  wigd=((Cos(t) -1._id)**2*11.25_id + 28._id + (Cos(t) -1._id)*36._id)*Sin(t*0.5_id)**6

endif

if(l==5.and.m1==-3.and.m2==4) then

  wigd=Cos(t*0.5_id)*(-1.0606601717798212866012665431572735589272539065327_id)*(6._id + Cos(t)*10._id)*Sin(t*0.5_id)**7

endif

if(l==5.and.m1==-3.and.m2==5) then

  wigd=Cos(t*0.5_id)**2*6.708203932499369089227521006193828706321855078835_id*Sin(t*0.5_id)**8

endif

if(l==5.and.m1==-2.and.m2==-5) then

  wigd=Cos(t*0.5_id)**7*10.95445115010332226913939565601604267905489389996_id*Sin(t*0.5_id)**3

endif

if(l==5.and.m1==-2.and.m2==-4) then

  wigd=Cos(t*0.5_id)**6*1.73205080756887729352744634150587236694280525381038_id*(-4._id + Cos(t)*10._id)*Sin(t*0.5_id)**2

endif

if(l==5.and.m1==-2.and.m2==-3) then

  wigd=Cos(t*0.5_id)**5*1.6329931618554520654648560498039275946439649871044_id*(3._id + (Cos(t) -1._id)**2*11.25_id + (Cos(t) -1._id)*13.5_id)*Sin(t*0.5_id)

endif

if(l==5.and.m1==-2.and.m2==-2) then

  wigd=Cos(t*0.5_id)**4*(1._id + (Cos(t) -1._id)*12._id + (Cos(t) -1._id)**3*15._id + (Cos(t) -1._id)**2*27._id)

endif

if(l==5.and.m1==-2.and.m2==-1) then

  wigd=Cos(t*0.5_id)**3*(-1.3228756555322952952508078768196302128551295915412_id)*(4._id + (Cos(t) -1._id)**3*15._id + (Cos(t) -1._id)*24._id + (Cos(t) -1._id)**2*36._id)*Sin(t*0.5_id)

endif

if(l==5.and.m1==-2.and.m2==0) then

  wigd=Cos(t*0.5_id)**2*1.4491376746189438573718664157169771723140132874759_id*(10._id + (Cos(t) -1._id)**3*15._id + (Cos(t) -1._id)*40._id + (Cos(t) -1._id)**2*45._id)*Sin(t*0.5_id)**2

endif

if(l==5.and.m1==-2.and.m2==1) then

  wigd=Cos(t*0.5_id)*(-1.3228756555322952952508078768196302128551295915412_id)*((Cos(t) -1._id)**3*15._id + 20._id + (Cos(t) -1._id)**2*54._id + (Cos(t) -1._id)*60._id)*Sin(t*0.5_id)**3

endif

if(l==5.and.m1==-2.and.m2==2) then

  wigd=((Cos(t) -1._id)**3*15._id + 35._id + (Cos(t) -1._id)**2*63._id + (Cos(t) -1._id)*84._id)*Sin(t*0.5_id)**4

endif

if(l==5.and.m1==-2.and.m2==3) then

  wigd=Cos(t*0.5_id)*(-1.6329931618554520654648560498039275946439649871044_id)*((Cos(t) -1._id)**2*11.25_id + 21._id + (Cos(t) -1._id)*31.5_id)*Sin(t*0.5_id)**5

endif

if(l==5.and.m1==-2.and.m2==4) then

  wigd=Cos(t*0.5_id)**2*1.73205080756887729352744634150587236694280525381038_id*(4._id + Cos(t)*10._id)*Sin(t*0.5_id)**6

endif

if(l==5.and.m1==-2.and.m2==5) then

  wigd=Cos(t*0.5_id)**3*(-10.95445115010332226913939565601604267905489389996_id)*Sin(t*0.5_id)**7

endif

if(l==5.and.m1==-1.and.m2==-5) then

  wigd=Cos(t*0.5_id)**6*14.491376746189438573718664157169771723140132874759_id*Sin(t*0.5_id)**4

endif

if(l==5.and.m1==-1.and.m2==-4) then

  wigd=Cos(t*0.5_id)**5*2.291287847477920003294023596864004244492228288384_id*(-2._id + Cos(t)*10._id)*Sin(t*0.5_id)**3

endif

if(l==5.and.m1==-1.and.m2==-3) then

  wigd=Cos(t*0.5_id)**4*2.1602468994692867436553224786959988859017347690194_id*(6._id + (Cos(t) -1._id)**2*11.25_id + (Cos(t) -1._id)*18._id)*Sin(t*0.5_id)**2

endif

if(l==5.and.m1==-1.and.m2==-2) then

  wigd=Cos(t*0.5_id)**3*1.3228756555322952952508078768196302128551295915412_id*(4._id + (Cos(t) -1._id)**3*15._id + (Cos(t) -1._id)*24._id + (Cos(t) -1._id)**2*36._id)*Sin(t*0.5_id)

endif

if(l==5.and.m1==-1.and.m2==-1) then

  wigd=Cos(t*0.5_id)**2*(1._id + (Cos(t) -1._id)**4*13.125_id + (Cos(t) -1._id)*14._id + (Cos(t) -1._id)**2*42._id + (Cos(t) -1._id)**3*42._id)

endif

if(l==5.and.m1==-1.and.m2==0) then

  wigd=Cos(t*0.5_id)*(-1.095445115010332226913939565601604267905489389996_id)*(5._id + (Cos(t) -1._id)**4*13.125_id + (Cos(t) -1._id)*35._id + (Cos(t) -1._id)**3*52.5_id + (Cos(t) -1._id)**2*70._id)*Sin(t*0.5_id)

endif

if(l==5.and.m1==-1.and.m2==1) then

  wigd=((Cos(t) -1._id)**4*13.125_id + 15._id + (Cos(t) -1._id)**3*63._id + (Cos(t) -1._id)*70._id + (Cos(t) -1._id)**2*105._id)*Sin(t*0.5_id)**2

endif

if(l==5.and.m1==-1.and.m2==2) then

  wigd=Cos(t*0.5_id)*(-1.3228756555322952952508078768196302128551295915412_id)*((Cos(t) -1._id)**3*15._id + 20._id + (Cos(t) -1._id)**2*54._id + (Cos(t) -1._id)*60._id)*Sin(t*0.5_id)**3

endif

if(l==5.and.m1==-1.and.m2==3) then

  wigd=Cos(t*0.5_id)**2*2.1602468994692867436553224786959988859017347690194_id*((Cos(t) -1._id)**2*11.25_id + 15._id + (Cos(t) -1._id)*27._id)*Sin(t*0.5_id)**4

endif

if(l==5.and.m1==-1.and.m2==4) then

  wigd=Cos(t*0.5_id)**3*(-2.291287847477920003294023596864004244492228288384_id)*(2._id + Cos(t)*10._id)*Sin(t*0.5_id)**5

endif

if(l==5.and.m1==-1.and.m2==5) then

  wigd=Cos(t*0.5_id)**4*14.491376746189438573718664157169771723140132874759_id*Sin(t*0.5_id)**6

endif

if(l==5.and.m1==0.and.m2==-5) then

  wigd=Cos(t*0.5_id)**5*15.874507866387543543009694521835562554261555098495_id*Sin(t*0.5_id)**5

endif

if(l==5.and.m1==0.and.m2==-4) then

  wigd=Cos(t)*Cos(t*0.5_id)**4*25.09980079602226643934516077355562468178446107896_id*Sin(t*0.5_id)**4

endif

if(l==5.and.m1==0.and.m2==-3) then

  wigd=Cos(t*0.5_id)**3*2.3664319132398464170269313166246468193662004923177_id*(10._id + (Cos(t) -1._id)**2*11.25_id + (Cos(t) -1._id)*22.5_id)*Sin(t*0.5_id)**3

endif

if(l==5.and.m1==0.and.m2==-2) then

  wigd=Cos(t*0.5_id)**2*1.4491376746189438573718664157169771723140132874759_id*(10._id + (Cos(t) -1._id)**3*15._id + (Cos(t) -1._id)*40._id + (Cos(t) -1._id)**2*45._id)*Sin(t*0.5_id)**2

endif

if(l==5.and.m1==0.and.m2==-1) then

  wigd=Cos(t*0.5_id)*1.09544511501033222691393956560160426790548938999597_id*(5._id + (Cos(t) -1._id)**4*13.125_id + (Cos(t) -1._id)*35._id + (Cos(t) -1._id)**3*52.5_id + (Cos(t) -1._id)**2*70._id)*Sin(t*0.5_id)

endif

if(l==5.and.m1==0.and.m2==0) then

  wigd=0.125_id*(Cos(t)**3*(-70._id) + Cos(t)*15._id + Cos(t)**5*63._id)

endif

if(l==5.and.m1==0.and.m2==1) then

  wigd=Cos(t*0.5_id)*(-1.095445115010332226913939565601604267905489389996_id)*(5._id + (Cos(t) -1._id)**4*13.125_id + (Cos(t) -1._id)*35._id + (Cos(t) -1._id)**3*52.5_id + (Cos(t) -1._id)**2*70._id)*Sin(t*0.5_id)

endif

if(l==5.and.m1==0.and.m2==2) then

  wigd=Cos(t*0.5_id)**2*1.4491376746189438573718664157169771723140132874759_id*(10._id + (Cos(t) -1._id)**3*15._id + (Cos(t) -1._id)*40._id + (Cos(t) -1._id)**2*45._id)*Sin(t*0.5_id)**2

endif

if(l==5.and.m1==0.and.m2==3) then

  wigd=Cos(t*0.5_id)**3*(-2.3664319132398464170269313166246468193662004923177_id)*(10._id + (Cos(t) -1._id)**2*11.25_id + (Cos(t) -1._id)*22.5_id)*Sin(t*0.5_id)**3

endif

if(l==5.and.m1==0.and.m2==4) then

  wigd=Cos(t)*Cos(t*0.5_id)**4*25.09980079602226643934516077355562468178446107896_id*Sin(t*0.5_id)**4

endif

if(l==5.and.m1==0.and.m2==5) then

  wigd=Cos(t*0.5_id)**5*(-15.874507866387543543009694521835562554261555098495_id)*Sin(t*0.5_id)**5

endif

if(l==5.and.m1==1.and.m2==-5) then

  wigd=Cos(t*0.5_id)**4*14.491376746189438573718664157169771723140132874759_id*Sin(t*0.5_id)**6

endif

if(l==5.and.m1==1.and.m2==-4) then

  wigd=Cos(t*0.5_id)**3*2.291287847477920003294023596864004244492228288384_id*(2._id + Cos(t)*10._id)*Sin(t*0.5_id)**5

endif

if(l==5.and.m1==1.and.m2==-3) then

  wigd=Cos(t*0.5_id)**2*2.1602468994692867436553224786959988859017347690194_id*((Cos(t) -1._id)**2*11.25_id + 15._id + (Cos(t) -1._id)*27._id)*Sin(t*0.5_id)**4

endif

if(l==5.and.m1==1.and.m2==-2) then

  wigd=Cos(t*0.5_id)*1.3228756555322952952508078768196302128551295915412_id*((Cos(t) -1._id)**3*15._id + 20._id + (Cos(t) -1._id)**2*54._id + (Cos(t) -1._id)*60._id)*Sin(t*0.5_id)**3

endif

if(l==5.and.m1==1.and.m2==-1) then

  wigd=((Cos(t) -1._id)**4*13.125_id + 15._id + (Cos(t) -1._id)**3*63._id + (Cos(t) -1._id)*70._id + (Cos(t) -1._id)**2*105._id)*Sin(t*0.5_id)**2

endif

if(l==5.and.m1==1.and.m2==0) then

  wigd=Cos(t*0.5_id)*1.09544511501033222691393956560160426790548938999597_id*(5._id + (Cos(t) -1._id)**4*13.125_id + (Cos(t) -1._id)*35._id + (Cos(t) -1._id)**3*52.5_id + (Cos(t) -1._id)**2*70._id)*Sin(t*0.5_id)

endif

if(l==5.and.m1==1.and.m2==1) then

  wigd=Cos(t*0.5_id)**2*(1._id + (Cos(t) -1._id)**4*13.125_id + (Cos(t) -1._id)*14._id + (Cos(t) -1._id)**2*42._id + (Cos(t) -1._id)**3*42._id)

endif

if(l==5.and.m1==1.and.m2==2) then

  wigd=Cos(t*0.5_id)**3*(-1.3228756555322952952508078768196302128551295915412_id)*(4._id + (Cos(t) -1._id)**3*15._id + (Cos(t) -1._id)*24._id + (Cos(t) -1._id)**2*36._id)*Sin(t*0.5_id)

endif

if(l==5.and.m1==1.and.m2==3) then

  wigd=Cos(t*0.5_id)**4*2.1602468994692867436553224786959988859017347690194_id*(6._id + (Cos(t) -1._id)**2*11.25_id + (Cos(t) -1._id)*18._id)*Sin(t*0.5_id)**2

endif

if(l==5.and.m1==1.and.m2==4) then

  wigd=Cos(t*0.5_id)**5*(-2.291287847477920003294023596864004244492228288384_id)*(-2._id + Cos(t)*10._id)*Sin(t*0.5_id)**3

endif

if(l==5.and.m1==1.and.m2==5) then

  wigd=Cos(t*0.5_id)**6*14.491376746189438573718664157169771723140132874759_id*Sin(t*0.5_id)**4

endif

if(l==5.and.m1==2.and.m2==-5) then

  wigd=Cos(t*0.5_id)**3*10.95445115010332226913939565601604267905489389996_id*Sin(t*0.5_id)**7

endif

if(l==5.and.m1==2.and.m2==-4) then

  wigd=Cos(t*0.5_id)**2*1.73205080756887729352744634150587236694280525381038_id*(4._id + Cos(t)*10._id)*Sin(t*0.5_id)**6

endif

if(l==5.and.m1==2.and.m2==-3) then

  wigd=Cos(t*0.5_id)*1.6329931618554520654648560498039275946439649871044_id*((Cos(t) -1._id)**2*11.25_id + 21._id + (Cos(t) -1._id)*31.5_id)*Sin(t*0.5_id)**5

endif

if(l==5.and.m1==2.and.m2==-2) then

  wigd=((Cos(t) -1._id)**3*15._id + 35._id + (Cos(t) -1._id)**2*63._id + (Cos(t) -1._id)*84._id)*Sin(t*0.5_id)**4

endif

if(l==5.and.m1==2.and.m2==-1) then

  wigd=Cos(t*0.5_id)*1.3228756555322952952508078768196302128551295915412_id*((Cos(t) -1._id)**3*15._id + 20._id + (Cos(t) -1._id)**2*54._id + (Cos(t) -1._id)*60._id)*Sin(t*0.5_id)**3

endif

if(l==5.and.m1==2.and.m2==0) then

  wigd=Cos(t*0.5_id)**2*1.4491376746189438573718664157169771723140132874759_id*(10._id + (Cos(t) -1._id)**3*15._id + (Cos(t) -1._id)*40._id + (Cos(t) -1._id)**2*45._id)*Sin(t*0.5_id)**2

endif

if(l==5.and.m1==2.and.m2==1) then

  wigd=Cos(t*0.5_id)**3*1.3228756555322952952508078768196302128551295915412_id*(4._id + (Cos(t) -1._id)**3*15._id + (Cos(t) -1._id)*24._id + (Cos(t) -1._id)**2*36._id)*Sin(t*0.5_id)

endif

if(l==5.and.m1==2.and.m2==2) then

  wigd=Cos(t*0.5_id)**4*(1._id + (Cos(t) -1._id)*12._id + (Cos(t) -1._id)**3*15._id + (Cos(t) -1._id)**2*27._id)

endif

if(l==5.and.m1==2.and.m2==3) then

  wigd=Cos(t*0.5_id)**5*(-1.6329931618554520654648560498039275946439649871044_id)*(3._id + (Cos(t) -1._id)**2*11.25_id + (Cos(t) -1._id)*13.5_id)*Sin(t*0.5_id)

endif

if(l==5.and.m1==2.and.m2==4) then

  wigd=Cos(t*0.5_id)**6*1.73205080756887729352744634150587236694280525381038_id*(-4._id + Cos(t)*10._id)*Sin(t*0.5_id)**2

endif

if(l==5.and.m1==2.and.m2==5) then

  wigd=Cos(t*0.5_id)**7*(-10.95445115010332226913939565601604267905489389996_id)*Sin(t*0.5_id)**3

endif

if(l==5.and.m1==3.and.m2==-5) then

  wigd=Cos(t*0.5_id)**2*6.708203932499369089227521006193828706321855078835_id*Sin(t*0.5_id)**8

endif

if(l==5.and.m1==3.and.m2==-4) then

  wigd=Cos(t*0.5_id)*1.0606601717798212866012665431572735589272539065327_id*(6._id + Cos(t)*10._id)*Sin(t*0.5_id)**7

endif

if(l==5.and.m1==3.and.m2==-3) then

  wigd=((Cos(t) -1._id)**2*11.25_id + 28._id + (Cos(t) -1._id)*36._id)*Sin(t*0.5_id)**6

endif

if(l==5.and.m1==3.and.m2==-2) then

  wigd=Cos(t*0.5_id)*1.6329931618554520654648560498039275946439649871044_id*((Cos(t) -1._id)**2*11.25_id + 21._id + (Cos(t) -1._id)*31.5_id)*Sin(t*0.5_id)**5

endif

if(l==5.and.m1==3.and.m2==-1) then

  wigd=Cos(t*0.5_id)**2*2.1602468994692867436553224786959988859017347690194_id*((Cos(t) -1._id)**2*11.25_id + 15._id + (Cos(t) -1._id)*27._id)*Sin(t*0.5_id)**4

endif

if(l==5.and.m1==3.and.m2==0) then

  wigd=Cos(t*0.5_id)**3*2.3664319132398464170269313166246468193662004923177_id*(10._id + (Cos(t) -1._id)**2*11.25_id + (Cos(t) -1._id)*22.5_id)*Sin(t*0.5_id)**3

endif

if(l==5.and.m1==3.and.m2==1) then

  wigd=Cos(t*0.5_id)**4*2.1602468994692867436553224786959988859017347690194_id*(6._id + (Cos(t) -1._id)**2*11.25_id + (Cos(t) -1._id)*18._id)*Sin(t*0.5_id)**2

endif

if(l==5.and.m1==3.and.m2==2) then

  wigd=Cos(t*0.5_id)**5*1.6329931618554520654648560498039275946439649871044_id*(3._id + (Cos(t) -1._id)**2*11.25_id + (Cos(t) -1._id)*13.5_id)*Sin(t*0.5_id)

endif

if(l==5.and.m1==3.and.m2==3) then

  wigd=Cos(t*0.5_id)**6*(1._id + (Cos(t) -1._id)*9._id + (Cos(t) -1._id)**2*11.25_id)

endif

if(l==5.and.m1==3.and.m2==4) then

  wigd=Cos(t*0.5_id)**7*(-1.0606601717798212866012665431572735589272539065327_id)*(-6._id + Cos(t)*10._id)*Sin(t*0.5_id)

endif

if(l==5.and.m1==3.and.m2==5) then

  wigd=Cos(t*0.5_id)**8*6.708203932499369089227521006193828706321855078835_id*Sin(t*0.5_id)**2

endif

if(l==5.and.m1==4.and.m2==-5) then

  wigd=Cos(t*0.5_id)*3.1622776601683793319988935444327185337195551393252_id*Sin(t*0.5_id)**9

endif

if(l==5.and.m1==4.and.m2==-4) then

  wigd=0.5_id*(8._id + Cos(t)*10._id)*Sin(t*0.5_id)**8

endif

if(l==5.and.m1==4.and.m2==-3) then

  wigd=Cos(t*0.5_id)*1.0606601717798212866012665431572735589272539065327_id*(6._id + Cos(t)*10._id)*Sin(t*0.5_id)**7

endif

if(l==5.and.m1==4.and.m2==-2) then

  wigd=Cos(t*0.5_id)**2*1.73205080756887729352744634150587236694280525381038_id*(4._id + Cos(t)*10._id)*Sin(t*0.5_id)**6

endif

if(l==5.and.m1==4.and.m2==-1) then

  wigd=Cos(t*0.5_id)**3*2.291287847477920003294023596864004244492228288384_id*(2._id + Cos(t)*10._id)*Sin(t*0.5_id)**5

endif

if(l==5.and.m1==4.and.m2==0) then

  wigd=Cos(t)*Cos(t*0.5_id)**4*25.09980079602226643934516077355562468178446107896_id*Sin(t*0.5_id)**4

endif

if(l==5.and.m1==4.and.m2==1) then

  wigd=Cos(t*0.5_id)**5*2.291287847477920003294023596864004244492228288384_id*(-2._id + Cos(t)*10._id)*Sin(t*0.5_id)**3

endif

if(l==5.and.m1==4.and.m2==2) then

  wigd=Cos(t*0.5_id)**6*1.73205080756887729352744634150587236694280525381038_id*(-4._id + Cos(t)*10._id)*Sin(t*0.5_id)**2

endif

if(l==5.and.m1==4.and.m2==3) then

  wigd=Cos(t*0.5_id)**7*1.0606601717798212866012665431572735589272539065327_id*(-6._id + Cos(t)*10._id)*Sin(t*0.5_id)

endif

if(l==5.and.m1==4.and.m2==4) then

  wigd=Cos(t*0.5_id)**8*0.5_id*(-8._id + Cos(t)*10._id)

endif

if(l==5.and.m1==4.and.m2==5) then

  wigd=Cos(t*0.5_id)**9*(-3.1622776601683793319988935444327185337195551393252_id)*Sin(t*0.5_id)

endif

if(l==5.and.m1==5.and.m2==-5) then

  wigd=Sin(t*0.5_id)**10

endif

if(l==5.and.m1==5.and.m2==-4) then

  wigd=Cos(t*0.5_id)*3.1622776601683793319988935444327185337195551393252_id*Sin(t*0.5_id)**9

endif

if(l==5.and.m1==5.and.m2==-3) then

  wigd=Cos(t*0.5_id)**2*6.708203932499369089227521006193828706321855078835_id*Sin(t*0.5_id)**8

endif

if(l==5.and.m1==5.and.m2==-2) then

  wigd=Cos(t*0.5_id)**3*10.95445115010332226913939565601604267905489389996_id*Sin(t*0.5_id)**7

endif

if(l==5.and.m1==5.and.m2==-1) then

  wigd=Cos(t*0.5_id)**4*14.491376746189438573718664157169771723140132874759_id*Sin(t*0.5_id)**6

endif

if(l==5.and.m1==5.and.m2==0) then

  wigd=Cos(t*0.5_id)**5*15.874507866387543543009694521835562554261555098495_id*Sin(t*0.5_id)**5

endif

if(l==5.and.m1==5.and.m2==1) then

  wigd=Cos(t*0.5_id)**6*14.491376746189438573718664157169771723140132874759_id*Sin(t*0.5_id)**4

endif

if(l==5.and.m1==5.and.m2==2) then

  wigd=Cos(t*0.5_id)**7*10.95445115010332226913939565601604267905489389996_id*Sin(t*0.5_id)**3

endif

if(l==5.and.m1==5.and.m2==3) then

  wigd=Cos(t*0.5_id)**8*6.708203932499369089227521006193828706321855078835_id*Sin(t*0.5_id)**2

endif

if(l==5.and.m1==5.and.m2==4) then

  wigd=Cos(t*0.5_id)**9*3.1622776601683793319988935444327185337195551393252_id*Sin(t*0.5_id)

endif

if(l==5.and.m1==5.and.m2==5) then

  wigd=Cos(t*0.5_id)**10

endif

end function wigd

end module acc_wig
