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

if(l==6.and.m1==-6.and.m2==-6) then
   wigd=Cos(t*0.5_id)**12
endif

if(l==6.and.m1==-6.and.m2==-5) then
   wigd=Cos(t*0.5_id)**11*-3.4641016151377545870548926830117447338856105076208_id*Sin(t*0.5_id)
endif

if(l==6.and.m1==-6.and.m2==-4) then
   wigd=Cos(t*0.5_id)**10*8.1240384046359603604598835682660403485042040867253_id*Sin(t*0.5_id)**2
endif

if(l==6.and.m1==-6.and.m2==-3) then
   wigd=Cos(t*0.5_id)**9*-14.832396974191325897422794881601426121959808638195_id*Sin(t*0.5_id)**3
endif

if(l==6.and.m1==-6.and.m2==-2) then
   wigd=Cos(t*0.5_id)**8*22.248595461286988846134192322402139182939712957293_id*Sin(t*0.5_id)**4
endif

if(l==6.and.m1==-6.and.m2==-1) then
   wigd=Cos(t*0.5_id)**7*-28.14249455894057732739378068126679768352937012047_id*Sin(t*0.5_id)**5
endif

if(l==6.and.m1==-6.and.m2==0) then
   wigd=Cos(t*0.5_id)**6*30.397368307141327263337488412079290075164652923404_id*Sin(t*0.5_id)**6
endif

if(l==6.and.m1==-6.and.m2==1) then
   wigd=Cos(t*0.5_id)**5*-28.14249455894057732739378068126679768352937012047_id*Sin(t*0.5_id)**7
endif

if(l==6.and.m1==-6.and.m2==2) then
   wigd=Cos(t*0.5_id)**4*22.248595461286988846134192322402139182939712957293_id*Sin(t*0.5_id)**8
endif

if(l==6.and.m1==-6.and.m2==3) then
   wigd=Cos(t*0.5_id)**3*-14.832396974191325897422794881601426121959808638195_id*Sin(t*0.5_id)**9
endif

if(l==6.and.m1==-6.and.m2==4) then
   wigd=Cos(t*0.5_id)**2*8.1240384046359603604598835682660403485042040867253_id*Sin(t*0.5_id)**10
endif

if(l==6.and.m1==-6.and.m2==5) then
   wigd=Cos(t*0.5_id)*-3.4641016151377545870548926830117447338856105076208_id*Sin(t*0.5_id)**11
endif

if(l==6.and.m1==-6.and.m2==6) then
   wigd=Sin(t*0.5_id)**12
endif

if(l==6.and.m1==-5.and.m2==-6) then
   wigd=Cos(t*0.5_id)**11*3.4641016151377545870548926830117447338856105076208_id*Sin(t*0.5_id)
endif

if(l==6.and.m1==-5.and.m2==-5) then
   wigd=Cos(t*0.5_id)**10*0.5_id*(-10._id + Cos(t)*12._id)
endif

if(l==6.and.m1==-5.and.m2==-4) then
   wigd=Cos(t*0.5_id)**9*-1.1726039399558573886414075283861165701470570883529_id*(-8._id + Cos(t)*12._id)*Sin(t*0.5_id)
endif

if(l==6.and.m1==-5.and.m2==-3) then
   wigd=Cos(t*0.5_id)**8*2.14087209644418817001948307960665673657789208526_id*(-6._id + Cos(t)*12._id)*Sin(t*0.5_id)**2
endif

if(l==6.and.m1==-5.and.m2==-2) then
   wigd=Cos(t*0.5_id)**7*-3.2113081446662822550292246194099851048668381278901_id*(-4._id + Cos(t)*12._id)*Sin(t*0.5_id)**3
endif

if(l==6.and.m1==-5.and.m2==-1) then
   wigd=Cos(t*0.5_id)**6*4.0620192023179801802299417841330201742521020433627_id*(-2._id + Cos(t)*12._id)*Sin(t*0.5_id)**4
endif

if(l==6.and.m1==-5.and.m2==0) then
   wigd=Cos(t)*Cos(t*0.5_id)**5*-52.649786324352732362438329844497857365255260965327_id*Sin(t*0.5_id)**5
endif

if(l==6.and.m1==-5.and.m2==1) then
   wigd=Cos(t*0.5_id)**4*4.0620192023179801802299417841330201742521020433627_id*(2._id + Cos(t)*12._id)*Sin(t*0.5_id)**6
endif

if(l==6.and.m1==-5.and.m2==2) then
   wigd=Cos(t*0.5_id)**3*-3.2113081446662822550292246194099851048668381278901_id*(4._id + Cos(t)*12._id)*Sin(t*0.5_id)**7
endif

if(l==6.and.m1==-5.and.m2==3) then
   wigd=Cos(t*0.5_id)**2*2.14087209644418817001948307960665673657789208526_id*(6._id + Cos(t)*12._id)*Sin(t*0.5_id)**8
endif

if(l==6.and.m1==-5.and.m2==4) then
   wigd=Cos(t*0.5_id)*-1.1726039399558573886414075283861165701470570883529_id*(8._id + Cos(t)*12._id)*Sin(t*0.5_id)**9
endif

if(l==6.and.m1==-5.and.m2==5) then
   wigd=0.5_id*(10._id + Cos(t)*12._id)*Sin(t*0.5_id)**10
endif

if(l==6.and.m1==-5.and.m2==6) then
   wigd=Cos(t*0.5_id)*-3.4641016151377545870548926830117447338856105076208_id*Sin(t*0.5_id)**11
endif

if(l==6.and.m1==-4.and.m2==-6) then
   wigd=Cos(t*0.5_id)**10*8.1240384046359603604598835682660403485042040867253_id*Sin(t*0.5_id)**2
endif

if(l==6.and.m1==-4.and.m2==-5) then
   wigd=Cos(t*0.5_id)**9*1.1726039399558573886414075283861165701470570883529_id*(-8._id + Cos(t)*12._id)*Sin(t*0.5_id)
endif

if(l==6.and.m1==-4.and.m2==-4) then
   wigd=Cos(t*0.5_id)**8*(1._id + (Cos(t) + -1._id)*11._id + (Cos(t) + -1._id)**2*16.5_id)
endif

if(l==6.and.m1==-4.and.m2==-3) then
   wigd=Cos(t*0.5_id)**7*-1.8257418583505537115232326093360071131758156499933_id*(3._id + (Cos(t) + -1._id)*16.5_id + (Cos(t) + -1._id)**2*16.5_id)*Sin(t*0.5_id)
endif

if(l==6.and.m1==-4.and.m2==-2) then
   wigd=Cos(t*0.5_id)**6*2.7386127875258305672848489140040106697637234749899_id*(6._id + (Cos(t) + -1._id)**2*16.5_id + (Cos(t) + -1._id)*22._id)*Sin(t*0.5_id)**2
endif

if(l==6.and.m1==-4.and.m2==-1) then
   wigd=Cos(t*0.5_id)**5*-3.4641016151377545870548926830117447338856105076208_id*(10._id + (Cos(t) + -1._id)**2*16.5_id + (Cos(t) + -1._id)*27.5_id)*Sin(t*0.5_id)**3
endif

if(l==6.and.m1==-4.and.m2==0) then
   wigd=Cos(t*0.5_id)**4*3.7416573867739413855837487323165493017560198077787_id*(15._id + (Cos(t) + -1._id)**2*16.5_id + (Cos(t) + -1._id)*33._id)*Sin(t*0.5_id)**4
endif

if(l==6.and.m1==-4.and.m2==1) then
   wigd=Cos(t*0.5_id)**3*-3.4641016151377545870548926830117447338856105076208_id*((Cos(t) + -1._id)**2*16.5_id + 21._id + (Cos(t) + -1._id)*38.5_id)*Sin(t*0.5_id)**5
endif

if(l==6.and.m1==-4.and.m2==2) then
   wigd=Cos(t*0.5_id)**2*2.7386127875258305672848489140040106697637234749899_id*((Cos(t) + -1._id)**2*16.5_id + 28._id + (Cos(t) + -1._id)*44._id)*Sin(t*0.5_id)**6
endif

if(l==6.and.m1==-4.and.m2==3) then
   wigd=Cos(t*0.5_id)*-1.8257418583505537115232326093360071131758156499933_id*((Cos(t) + -1._id)**2*16.5_id + 36._id + (Cos(t) + -1._id)*49.5_id)*Sin(t*0.5_id)**7
endif

if(l==6.and.m1==-4.and.m2==4) then
   wigd=((Cos(t) + -1._id)**2*16.5_id + 45._id + (Cos(t) + -1._id)*55._id)*Sin(t*0.5_id)**8
endif

if(l==6.and.m1==-4.and.m2==5) then
   wigd=Cos(t*0.5_id)*-1.1726039399558573886414075283861165701470570883529_id*(8._id + Cos(t)*12._id)*Sin(t*0.5_id)**9
endif

if(l==6.and.m1==-4.and.m2==6) then
   wigd=Cos(t*0.5_id)**2*8.1240384046359603604598835682660403485042040867253_id*Sin(t*0.5_id)**10
endif

if(l==6.and.m1==-3.and.m2==-6) then
   wigd=Cos(t*0.5_id)**9*14.832396974191325897422794881601426121959808638195_id*Sin(t*0.5_id)**3
endif

if(l==6.and.m1==-3.and.m2==-5) then
   wigd=Cos(t*0.5_id)**8*2.14087209644418817001948307960665673657789208526_id*(-6._id + Cos(t)*12._id)*Sin(t*0.5_id)**2
endif

if(l==6.and.m1==-3.and.m2==-4) then
   wigd=Cos(t*0.5_id)**7*1.82574185835055371152323260933600711317581564999328_id*(3._id + (Cos(t) + -1._id)*16.5_id + (Cos(t) + -1._id)**2*16.5_id)*Sin(t*0.5_id)
endif

if(l==6.and.m1==-3.and.m2==-3) then
   wigd=Cos(t*0.5_id)**6*(1._id + (Cos(t) + -1._id)*15._id + (Cos(t) + -1._id)**3*27.5_id + (Cos(t) + -1._id)**2*41.25_id)
endif

if(l==6.and.m1==-3.and.m2==-2) then
   wigd=Cos(t*0.5_id)**5*-1.5_id*(4._id + (Cos(t) + -1._id)**3*27.5_id + (Cos(t) + -1._id)*30._id + (Cos(t) + -1._id)**2*55._id)*Sin(t*0.5_id)
endif

if(l==6.and.m1==-3.and.m2==-1) then
   wigd=Cos(t*0.5_id)**4*1.8973665961010275991993361266596311202317330835951_id*(10._id + (Cos(t) + -1._id)**3*27.5_id + (Cos(t) + -1._id)*50._id + (Cos(t) + -1._id)**2*68.75_id)*Sin(t*0.5_id)**2
endif

if(l==6.and.m1==-3.and.m2==0) then
   wigd=Cos(t*0.5_id)**3*-2.049390153191919676644207736104210398147006532691_id*(20._id + (Cos(t) + -1._id)**3*27.5_id +&
& (Cos(t) + -1._id)*75._id + (Cos(t) + -1._id)**2*82.5_id)*Sin(t*0.5_id)**3
endif

if(l==6.and.m1==-3.and.m2==1) then
   wigd=Cos(t*0.5_id)**2*1.8973665961010275991993361266596311202317330835951_id*((Cos(t) + -1._id)**3*27.5_id + 35._id +&
& (Cos(t) + -1._id)**2*96.25_id + (Cos(t) + -1._id)*105._id)*Sin(t*0.5_id)**4
endif

if(l==6.and.m1==-3.and.m2==2) then
   wigd=Cos(t*0.5_id)*-1.5_id*((Cos(t) + -1._id)**3*27.5_id + 56._id + (Cos(t) + -1._id)**2*110._id + (Cos(t) + -1._id)*&
&140._id)*Sin(t*0.5_id)**5
endif

if(l==6.and.m1==-3.and.m2==3) then
   wigd=((Cos(t) + -1._id)**3*27.5_id + 84._id + (Cos(t) + -1._id)**2*123.75_id + (Cos(t) + -1._id)*180._id)*Sin(t*0.5_i&
&d)**6
endif

if(l==6.and.m1==-3.and.m2==4) then
   wigd=Cos(t*0.5_id)*-1.8257418583505537115232326093360071131758156499933_id*((Cos(t) + -1._id)**2*16.5_id + 36._id + (&
&Cos(t) + -1._id)*49.5_id)*Sin(t*0.5_id)**7
endif

if(l==6.and.m1==-3.and.m2==5) then
   wigd=Cos(t*0.5_id)**2*2.14087209644418817001948307960665673657789208526_id*(6._id + Cos(t)*12._id)*Sin(t*0.5_id)**8
endif

if(l==6.and.m1==-3.and.m2==6) then
   wigd=Cos(t*0.5_id)**3*-14.832396974191325897422794881601426121959808638195_id*Sin(t*0.5_id)**9
endif

if(l==6.and.m1==-2.and.m2==-6) then
   wigd=Cos(t*0.5_id)**8*22.248595461286988846134192322402139182939712957293_id*Sin(t*0.5_id)**4
endif

if(l==6.and.m1==-2.and.m2==-5) then
   wigd=Cos(t*0.5_id)**7*3.2113081446662822550292246194099851048668381278901_id*(-4._id + Cos(t)*12._id)*Sin(t*0.5_id)**&
&3
endif

if(l==6.and.m1==-2.and.m2==-4) then
   wigd=Cos(t*0.5_id)**6*2.7386127875258305672848489140040106697637234749899_id*(6._id + (Cos(t) + -1._id)**2*16.5_id + &
&(Cos(t) + -1._id)*22._id)*Sin(t*0.5_id)**2
endif

if(l==6.and.m1==-2.and.m2==-3) then
   wigd=Cos(t*0.5_id)**5*1.5_id*(4._id + (Cos(t) + -1._id)**3*27.5_id + (Cos(t) + -1._id)*30._id + (Cos(t) + -1._id)**2*&
&55._id)*Sin(t*0.5_id)
endif

if(l==6.and.m1==-2.and.m2==-2) then
   wigd=Cos(t*0.5_id)**4*(1._id + (Cos(t) + -1._id)*18._id + (Cos(t) + -1._id)**4*30.9375_id + (Cos(t) + -1._id)**2*67.5&
&_id + (Cos(t) + -1._id)**3*82.5_id)
endif

if(l==6.and.m1==-2.and.m2==-1) then
   wigd=Cos(t*0.5_id)**3*-1.2649110640673517327995574177730874134878220557301_id*(5._id + (Cos(t) + -1._id)**4*30.9375_i&
&d + (Cos(t) + -1._id)*45._id + (Cos(t) + -1._id)**3*103.125_id + (Cos(t) + -1._id)**2*112.5_id)*Sin(t*0.5_id)
endif

if(l==6.and.m1==-2.and.m2==0) then
   wigd=Cos(t*0.5_id)**2*1.3662601021279464510961384907361402654313376884606_id*(15._id + (Cos(t) + -1._id)**4*30.9375_i&
&d + (Cos(t) + -1._id)*90._id + (Cos(t) + -1._id)**3*123.75_id + (Cos(t) + -1._id)**2*168.75_id)*Sin(t*0.5_id)**2
endif

if(l==6.and.m1==-2.and.m2==1) then
   wigd=Cos(t*0.5_id)*-1.2649110640673517327995574177730874134878220557301_id*((Cos(t) + -1._id)**4*30.9375_id + 35._id &
&+ (Cos(t) + -1._id)**3*144.375_id + (Cos(t) + -1._id)*157.5_id + (Cos(t) + -1._id)**2*236.25_id)*Sin(t*0.5_id)**3
endif

if(l==6.and.m1==-2.and.m2==2) then
   wigd=((Cos(t) + -1._id)**4*30.9375_id + 70._id + (Cos(t) + -1._id)**3*165._id + (Cos(t) + -1._id)*252._id + (Cos(t) +&
& -1._id)**2*315._id)*Sin(t*0.5_id)**4
endif

if(l==6.and.m1==-2.and.m2==3) then
   wigd=Cos(t*0.5_id)*-1.5_id*((Cos(t) + -1._id)**3*27.5_id + 56._id + (Cos(t) + -1._id)**2*110._id + (Cos(t) + -1._id)*&
&140._id)*Sin(t*0.5_id)**5
endif

if(l==6.and.m1==-2.and.m2==4) then
   wigd=Cos(t*0.5_id)**2*2.7386127875258305672848489140040106697637234749899_id*((Cos(t) + -1._id)**2*16.5_id + 28._id +&
& (Cos(t) + -1._id)*44._id)*Sin(t*0.5_id)**6
endif

if(l==6.and.m1==-2.and.m2==5) then
   wigd=Cos(t*0.5_id)**3*-3.2113081446662822550292246194099851048668381278901_id*(4._id + Cos(t)*12._id)*Sin(t*0.5_id)**&
&7
endif

if(l==6.and.m1==-2.and.m2==6) then
   wigd=Cos(t*0.5_id)**4*22.248595461286988846134192322402139182939712957293_id*Sin(t*0.5_id)**8
endif

if(l==6.and.m1==-1.and.m2==-6) then
   wigd=Cos(t*0.5_id)**7*28.14249455894057732739378068126679768352937012047_id*Sin(t*0.5_id)**5
endif

if(l==6.and.m1==-1.and.m2==-5) then
   wigd=Cos(t*0.5_id)**6*4.0620192023179801802299417841330201742521020433627_id*(-2._id + Cos(t)*12._id)*Sin(t*0.5_id)**&
&4
endif

if(l==6.and.m1==-1.and.m2==-4) then
   wigd=Cos(t*0.5_id)**5*3.4641016151377545870548926830117447338856105076208_id*(10._id + (Cos(t) + -1._id)**2*16.5_id +&
& (Cos(t) + -1._id)*27.5_id)*Sin(t*0.5_id)**3
endif

if(l==6.and.m1==-1.and.m2==-3) then
   wigd=Cos(t*0.5_id)**4*1.8973665961010275991993361266596311202317330835951_id*(10._id + (Cos(t) + -1._id)**3*27.5_id +&
& (Cos(t) + -1._id)*50._id + (Cos(t) + -1._id)**2*68.75_id)*Sin(t*0.5_id)**2
endif

if(l==6.and.m1==-1.and.m2==-2) then
   wigd=Cos(t*0.5_id)**3*1.2649110640673517327995574177730874134878220557301_id*(5._id + (Cos(t) + -1._id)**4*30.9375_id&
& + (Cos(t) + -1._id)*45._id + (Cos(t) + -1._id)**3*103.125_id + (Cos(t) + -1._id)**2*112.5_id)*Sin(t*0.5_id)
endif

if(l==6.and.m1==-1.and.m2==-1) then
   wigd=Cos(t*0.5_id)**2*(1._id + (Cos(t) + -1._id)*20._id + (Cos(t) + -1._id)**5*24.75_id + (Cos(t) + -1._id)**2*90._id&
& + (Cos(t) + -1._id)**4*103.125_id + (Cos(t) + -1._id)**3*150._id)
endif

if(l==6.and.m1==-1.and.m2==0) then
   wigd=Cos(t*0.5_id)*-1.0801234497346433718276612393479994429508673845097_id*(6._id + (Cos(t) + -1._id)**5*24.75_id + (&
&Cos(t) + -1._id)*60._id + (Cos(t) + -1._id)**4*123.75_id + (Cos(t) + -1._id)**2*180._id + (Cos(t) + -1._id)**3*225._id)*Sin(t*0.5_id)
endif

if(l==6.and.m1==-1.and.m2==1) then
   wigd=(21._id + (Cos(t) + -1._id)**5*24.75_id + (Cos(t) + -1._id)*140._id + (Cos(t) + -1._id)**4*144.375_id + (Cos(t) &
&+ -1._id)**2*315._id + (Cos(t) + -1._id)**3*315._id)*Sin(t*0.5_id)**2
endif

if(l==6.and.m1==-1.and.m2==2) then
   wigd=Cos(t*0.5_id)*-1.2649110640673517327995574177730874134878220557301_id*((Cos(t) + -1._id)**4*30.9375_id + 35._id &
&+ (Cos(t) + -1._id)**3*144.375_id + (Cos(t) + -1._id)*157.5_id + (Cos(t) + -1._id)**2*236.25_id)*Sin(t*0.5_id)**3
endif

if(l==6.and.m1==-1.and.m2==3) then
   wigd=Cos(t*0.5_id)**2*1.8973665961010275991993361266596311202317330835951_id*((Cos(t) + -1._id)**3*27.5_id + 35._id +&
& (Cos(t) + -1._id)**2*96.25_id + (Cos(t) + -1._id)*105._id)*Sin(t*0.5_id)**4
endif

if(l==6.and.m1==-1.and.m2==4) then
   wigd=Cos(t*0.5_id)**3*-3.4641016151377545870548926830117447338856105076208_id*((Cos(t) + -1._id)**2*16.5_id + 21._id &
&+ (Cos(t) + -1._id)*38.5_id)*Sin(t*0.5_id)**5
endif

if(l==6.and.m1==-1.and.m2==5) then
   wigd=Cos(t*0.5_id)**4*4.0620192023179801802299417841330201742521020433627_id*(2._id + Cos(t)*12._id)*Sin(t*0.5_id)**6&
&
endif

if(l==6.and.m1==-1.and.m2==6) then
   wigd=Cos(t*0.5_id)**5*-28.14249455894057732739378068126679768352937012047_id*Sin(t*0.5_id)**7
endif

if(l==6.and.m1==0.and.m2==-6) then
   wigd=Cos(t*0.5_id)**6*30.397368307141327263337488412079290075164652923404_id*Sin(t*0.5_id)**6
endif

if(l==6.and.m1==0.and.m2==-5) then
   wigd=Cos(t)*Cos(t*0.5_id)**5*52.649786324352732362438329844497857365255260965327_id*Sin(t*0.5_id)**5
endif

if(l==6.and.m1==0.and.m2==-4) then
   wigd=Cos(t*0.5_id)**4*3.7416573867739413855837487323165493017560198077787_id*(15._id + (Cos(t) + -1._id)**2*16.5_id +&
& (Cos(t) + -1._id)*33._id)*Sin(t*0.5_id)**4
endif

if(l==6.and.m1==0.and.m2==-3) then
   wigd=Cos(t*0.5_id)**3*2.049390153191919676644207736104210398147006532691_id*(20._id + (Cos(t) + -1._id)**3*27.5_id + &
&(Cos(t) + -1._id)*75._id + (Cos(t) + -1._id)**2*82.5_id)*Sin(t*0.5_id)**3
endif

if(l==6.and.m1==0.and.m2==-2) then
   wigd=Cos(t*0.5_id)**2*1.3662601021279464510961384907361402654313376884606_id*(15._id + (Cos(t) + -1._id)**4*30.9375_i&
&d + (Cos(t) + -1._id)*90._id + (Cos(t) + -1._id)**3*123.75_id + (Cos(t) + -1._id)**2*168.75_id)*Sin(t*0.5_id)**2
endif

if(l==6.and.m1==0.and.m2==-1) then
   wigd=Cos(t*0.5_id)*1.08012344973464337182766123934799944295086738450972_id*(6._id + (Cos(t) + -1._id)**5*24.75_id + (&
&Cos(t) + -1._id)*60._id + (Cos(t) + -1._id)**4*123.75_id + (Cos(t) + -1._id)**2*180._id + (Cos(t) + -1._id)**3*225._id)*Sin(t*0.5_id)
endif

if(l==6.and.m1==0.and.m2==0) then
   wigd=0.0625_id*(Cos(t)**4*-315._id + -5._id + Cos(t)**2*105._id + Cos(t)**6*231._id)
endif

if(l==6.and.m1==0.and.m2==1) then
   wigd=Cos(t*0.5_id)*-1.0801234497346433718276612393479994429508673845097_id*(6._id + (Cos(t) + -1._id)**5*24.75_id + (&
&Cos(t) + -1._id)*60._id + (Cos(t) + -1._id)**4*123.75_id + (Cos(t) + -1._id)**2*180._id + (Cos(t) + -1._id)**3*225._id)*Sin(t*0.5_id)
endif

if(l==6.and.m1==0.and.m2==2) then
   wigd=Cos(t*0.5_id)**2*1.3662601021279464510961384907361402654313376884606_id*(15._id + (Cos(t) + -1._id)**4*30.9375_i&
&d + (Cos(t) + -1._id)*90._id + (Cos(t) + -1._id)**3*123.75_id + (Cos(t) + -1._id)**2*168.75_id)*Sin(t*0.5_id)**2
endif

if(l==6.and.m1==0.and.m2==3) then
   wigd=Cos(t*0.5_id)**3*-2.049390153191919676644207736104210398147006532691_id*(20._id + (Cos(t) + -1._id)**3*27.5_id +&
& (Cos(t) + -1._id)*75._id + (Cos(t) + -1._id)**2*82.5_id)*Sin(t*0.5_id)**3
endif

if(l==6.and.m1==0.and.m2==4) then
   wigd=Cos(t*0.5_id)**4*3.7416573867739413855837487323165493017560198077787_id*(15._id + (Cos(t) + -1._id)**2*16.5_id +&
& (Cos(t) + -1._id)*33._id)*Sin(t*0.5_id)**4
endif

if(l==6.and.m1==0.and.m2==5) then
   wigd=Cos(t)*Cos(t*0.5_id)**5*-52.649786324352732362438329844497857365255260965327_id*Sin(t*0.5_id)**5
endif

if(l==6.and.m1==0.and.m2==6) then
   wigd=Cos(t*0.5_id)**6*30.397368307141327263337488412079290075164652923404_id*Sin(t*0.5_id)**6
endif

if(l==6.and.m1==1.and.m2==-6) then
   wigd=Cos(t*0.5_id)**5*28.14249455894057732739378068126679768352937012047_id*Sin(t*0.5_id)**7
endif

if(l==6.and.m1==1.and.m2==-5) then
   wigd=Cos(t*0.5_id)**4*4.0620192023179801802299417841330201742521020433627_id*(2._id + Cos(t)*12._id)*Sin(t*0.5_id)**6&
&
endif

if(l==6.and.m1==1.and.m2==-4) then
   wigd=Cos(t*0.5_id)**3*3.4641016151377545870548926830117447338856105076208_id*((Cos(t) + -1._id)**2*16.5_id + 21._id +&
& (Cos(t) + -1._id)*38.5_id)*Sin(t*0.5_id)**5
endif

if(l==6.and.m1==1.and.m2==-3) then
   wigd=Cos(t*0.5_id)**2*1.8973665961010275991993361266596311202317330835951_id*((Cos(t) + -1._id)**3*27.5_id + 35._id +&
& (Cos(t) + -1._id)**2*96.25_id + (Cos(t) + -1._id)*105._id)*Sin(t*0.5_id)**4
endif

if(l==6.and.m1==1.and.m2==-2) then
   wigd=Cos(t*0.5_id)*1.2649110640673517327995574177730874134878220557301_id*((Cos(t) + -1._id)**4*30.9375_id + 35._id +&
& (Cos(t) + -1._id)**3*144.375_id + (Cos(t) + -1._id)*157.5_id + (Cos(t) + -1._id)**2*236.25_id)*Sin(t*0.5_id)**3
endif

if(l==6.and.m1==1.and.m2==-1) then
   wigd=(21._id + (Cos(t) + -1._id)**5*24.75_id + (Cos(t) + -1._id)*140._id + (Cos(t) + -1._id)**4*144.375_id + (Cos(t) &
&+ -1._id)**2*315._id + (Cos(t) + -1._id)**3*315._id)*Sin(t*0.5_id)**2
endif

if(l==6.and.m1==1.and.m2==0) then
   wigd=Cos(t*0.5_id)*1.08012344973464337182766123934799944295086738450972_id*(6._id + (Cos(t) + -1._id)**5*24.75_id + (&
&Cos(t) + -1._id)*60._id + (Cos(t) + -1._id)**4*123.75_id + (Cos(t) + -1._id)**2*180._id + (Cos(t) + -1._id)**3*225._id)*Sin(t*0.5_id)
endif

if(l==6.and.m1==1.and.m2==1) then
   wigd=Cos(t*0.5_id)**2*(1._id + (Cos(t) + -1._id)*20._id + (Cos(t) + -1._id)**5*24.75_id + (Cos(t) + -1._id)**2*90._id&
& + (Cos(t) + -1._id)**4*103.125_id + (Cos(t) + -1._id)**3*150._id)
endif

if(l==6.and.m1==1.and.m2==2) then
   wigd=Cos(t*0.5_id)**3*-1.2649110640673517327995574177730874134878220557301_id*(5._id + (Cos(t) + -1._id)**4*30.9375_i&
&d + (Cos(t) + -1._id)*45._id + (Cos(t) + -1._id)**3*103.125_id + (Cos(t) + -1._id)**2*112.5_id)*Sin(t*0.5_id)
endif

if(l==6.and.m1==1.and.m2==3) then
   wigd=Cos(t*0.5_id)**4*1.8973665961010275991993361266596311202317330835951_id*(10._id + (Cos(t) + -1._id)**3*27.5_id +&
& (Cos(t) + -1._id)*50._id + (Cos(t) + -1._id)**2*68.75_id)*Sin(t*0.5_id)**2
endif

if(l==6.and.m1==1.and.m2==4) then
   wigd=Cos(t*0.5_id)**5*-3.4641016151377545870548926830117447338856105076208_id*(10._id + (Cos(t) + -1._id)**2*16.5_id &
&+ (Cos(t) + -1._id)*27.5_id)*Sin(t*0.5_id)**3
endif

if(l==6.and.m1==1.and.m2==5) then
   wigd=Cos(t*0.5_id)**6*4.0620192023179801802299417841330201742521020433627_id*(-2._id + Cos(t)*12._id)*Sin(t*0.5_id)**&
&4
endif

if(l==6.and.m1==1.and.m2==6) then
   wigd=Cos(t*0.5_id)**7*-28.14249455894057732739378068126679768352937012047_id*Sin(t*0.5_id)**5
endif

if(l==6.and.m1==2.and.m2==-6) then
   wigd=Cos(t*0.5_id)**4*22.248595461286988846134192322402139182939712957293_id*Sin(t*0.5_id)**8
endif

if(l==6.and.m1==2.and.m2==-5) then
   wigd=Cos(t*0.5_id)**3*3.2113081446662822550292246194099851048668381278901_id*(4._id + Cos(t)*12._id)*Sin(t*0.5_id)**7&
&
endif

if(l==6.and.m1==2.and.m2==-4) then
   wigd=Cos(t*0.5_id)**2*2.7386127875258305672848489140040106697637234749899_id*((Cos(t) + -1._id)**2*16.5_id + 28._id +&
& (Cos(t) + -1._id)*44._id)*Sin(t*0.5_id)**6
endif

if(l==6.and.m1==2.and.m2==-3) then
   wigd=Cos(t*0.5_id)*1.5_id*((Cos(t) + -1._id)**3*27.5_id + 56._id + (Cos(t) + -1._id)**2*110._id + (Cos(t) + -1._id)*1&
&40._id)*Sin(t*0.5_id)**5
endif

if(l==6.and.m1==2.and.m2==-2) then
   wigd=((Cos(t) + -1._id)**4*30.9375_id + 70._id + (Cos(t) + -1._id)**3*165._id + (Cos(t) + -1._id)*252._id + (Cos(t) +&
& -1._id)**2*315._id)*Sin(t*0.5_id)**4
endif

if(l==6.and.m1==2.and.m2==-1) then
   wigd=Cos(t*0.5_id)*1.2649110640673517327995574177730874134878220557301_id*((Cos(t) + -1._id)**4*30.9375_id + 35._id +&
& (Cos(t) + -1._id)**3*144.375_id + (Cos(t) + -1._id)*157.5_id + (Cos(t) + -1._id)**2*236.25_id)*Sin(t*0.5_id)**3
endif

if(l==6.and.m1==2.and.m2==0) then
   wigd=Cos(t*0.5_id)**2*1.3662601021279464510961384907361402654313376884606_id*(15._id + (Cos(t) + -1._id)**4*30.9375_i&
&d + (Cos(t) + -1._id)*90._id + (Cos(t) + -1._id)**3*123.75_id + (Cos(t) + -1._id)**2*168.75_id)*Sin(t*0.5_id)**2
endif

if(l==6.and.m1==2.and.m2==1) then
   wigd=Cos(t*0.5_id)**3*1.2649110640673517327995574177730874134878220557301_id*(5._id + (Cos(t) + -1._id)**4*30.9375_id&
& + (Cos(t) + -1._id)*45._id + (Cos(t) + -1._id)**3*103.125_id + (Cos(t) + -1._id)**2*112.5_id)*Sin(t*0.5_id)
endif

if(l==6.and.m1==2.and.m2==2) then
   wigd=Cos(t*0.5_id)**4*(1._id + (Cos(t) + -1._id)*18._id + (Cos(t) + -1._id)**4*30.9375_id + (Cos(t) + -1._id)**2*67.5&
&_id + (Cos(t) + -1._id)**3*82.5_id)
endif

if(l==6.and.m1==2.and.m2==3) then
   wigd=Cos(t*0.5_id)**5*-1.5_id*(4._id + (Cos(t) + -1._id)**3*27.5_id + (Cos(t) + -1._id)*30._id + (Cos(t) + -1._id)**2&
&*55._id)*Sin(t*0.5_id)
endif

if(l==6.and.m1==2.and.m2==4) then
   wigd=Cos(t*0.5_id)**6*2.7386127875258305672848489140040106697637234749899_id*(6._id + (Cos(t) + -1._id)**2*16.5_id + &
&(Cos(t) + -1._id)*22._id)*Sin(t*0.5_id)**2
endif

if(l==6.and.m1==2.and.m2==5) then
   wigd=Cos(t*0.5_id)**7*-3.2113081446662822550292246194099851048668381278901_id*(-4._id + Cos(t)*12._id)*Sin(t*0.5_id)*&
&*3
endif

if(l==6.and.m1==2.and.m2==6) then
   wigd=Cos(t*0.5_id)**8*22.248595461286988846134192322402139182939712957293_id*Sin(t*0.5_id)**4
endif

if(l==6.and.m1==3.and.m2==-6) then
   wigd=Cos(t*0.5_id)**3*14.832396974191325897422794881601426121959808638195_id*Sin(t*0.5_id)**9
endif

if(l==6.and.m1==3.and.m2==-5) then
   wigd=Cos(t*0.5_id)**2*2.14087209644418817001948307960665673657789208526_id*(6._id + Cos(t)*12._id)*Sin(t*0.5_id)**8
endif

if(l==6.and.m1==3.and.m2==-4) then
   wigd=Cos(t*0.5_id)*1.82574185835055371152323260933600711317581564999328_id*((Cos(t) + -1._id)**2*16.5_id + 36._id + (&
&Cos(t) + -1._id)*49.5_id)*Sin(t*0.5_id)**7
endif

if(l==6.and.m1==3.and.m2==-3) then
   wigd=((Cos(t) + -1._id)**3*27.5_id + 84._id + (Cos(t) + -1._id)**2*123.75_id + (Cos(t) + -1._id)*180._id)*Sin(t*0.5_i&
&d)**6
endif

if(l==6.and.m1==3.and.m2==-2) then
   wigd=Cos(t*0.5_id)*1.5_id*((Cos(t) + -1._id)**3*27.5_id + 56._id + (Cos(t) + -1._id)**2*110._id + (Cos(t) + -1._id)*1&
&40._id)*Sin(t*0.5_id)**5
endif

if(l==6.and.m1==3.and.m2==-1) then
   wigd=Cos(t*0.5_id)**2*1.8973665961010275991993361266596311202317330835951_id*((Cos(t) + -1._id)**3*27.5_id + 35._id +&
& (Cos(t) + -1._id)**2*96.25_id + (Cos(t) + -1._id)*105._id)*Sin(t*0.5_id)**4
endif

if(l==6.and.m1==3.and.m2==0) then
   wigd=Cos(t*0.5_id)**3*2.049390153191919676644207736104210398147006532691_id*(20._id + (Cos(t) + -1._id)**3*27.5_id + &
&(Cos(t) + -1._id)*75._id + (Cos(t) + -1._id)**2*82.5_id)*Sin(t*0.5_id)**3
endif

if(l==6.and.m1==3.and.m2==1) then
   wigd=Cos(t*0.5_id)**4*1.8973665961010275991993361266596311202317330835951_id*(10._id + (Cos(t) + -1._id)**3*27.5_id +&
& (Cos(t) + -1._id)*50._id + (Cos(t) + -1._id)**2*68.75_id)*Sin(t*0.5_id)**2
endif

if(l==6.and.m1==3.and.m2==2) then
   wigd=Cos(t*0.5_id)**5*1.5_id*(4._id + (Cos(t) + -1._id)**3*27.5_id + (Cos(t) + -1._id)*30._id + (Cos(t) + -1._id)**2*&
&55._id)*Sin(t*0.5_id)
endif

if(l==6.and.m1==3.and.m2==3) then
   wigd=Cos(t*0.5_id)**6*(1._id + (Cos(t) + -1._id)*15._id + (Cos(t) + -1._id)**3*27.5_id + (Cos(t) + -1._id)**2*41.25_i&
&d)
endif

if(l==6.and.m1==3.and.m2==4) then
   wigd=Cos(t*0.5_id)**7*-1.8257418583505537115232326093360071131758156499933_id*(3._id + (Cos(t) + -1._id)*16.5_id + (C&
&os(t) + -1._id)**2*16.5_id)*Sin(t*0.5_id)
endif

if(l==6.and.m1==3.and.m2==5) then
   wigd=Cos(t*0.5_id)**8*2.14087209644418817001948307960665673657789208526_id*(-6._id + Cos(t)*12._id)*Sin(t*0.5_id)**2
endif

if(l==6.and.m1==3.and.m2==6) then
   wigd=Cos(t*0.5_id)**9*-14.832396974191325897422794881601426121959808638195_id*Sin(t*0.5_id)**3
endif

if(l==6.and.m1==4.and.m2==-6) then
   wigd=Cos(t*0.5_id)**2*8.1240384046359603604598835682660403485042040867253_id*Sin(t*0.5_id)**10
endif

if(l==6.and.m1==4.and.m2==-5) then
   wigd=Cos(t*0.5_id)*1.1726039399558573886414075283861165701470570883529_id*(8._id + Cos(t)*12._id)*Sin(t*0.5_id)**9
endif

if(l==6.and.m1==4.and.m2==-4) then
   wigd=((Cos(t) + -1._id)**2*16.5_id + 45._id + (Cos(t) + -1._id)*55._id)*Sin(t*0.5_id)**8
endif

if(l==6.and.m1==4.and.m2==-3) then
   wigd=Cos(t*0.5_id)*1.82574185835055371152323260933600711317581564999328_id*((Cos(t) + -1._id)**2*16.5_id + 36._id + (&
&Cos(t) + -1._id)*49.5_id)*Sin(t*0.5_id)**7
endif

if(l==6.and.m1==4.and.m2==-2) then
   wigd=Cos(t*0.5_id)**2*2.7386127875258305672848489140040106697637234749899_id*((Cos(t) + -1._id)**2*16.5_id + 28._id +&
& (Cos(t) + -1._id)*44._id)*Sin(t*0.5_id)**6
endif

if(l==6.and.m1==4.and.m2==-1) then
   wigd=Cos(t*0.5_id)**3*3.4641016151377545870548926830117447338856105076208_id*((Cos(t) + -1._id)**2*16.5_id + 21._id +&
& (Cos(t) + -1._id)*38.5_id)*Sin(t*0.5_id)**5
endif

if(l==6.and.m1==4.and.m2==0) then
   wigd=Cos(t*0.5_id)**4*3.7416573867739413855837487323165493017560198077787_id*(15._id + (Cos(t) + -1._id)**2*16.5_id +&
& (Cos(t) + -1._id)*33._id)*Sin(t*0.5_id)**4
endif

if(l==6.and.m1==4.and.m2==1) then
   wigd=Cos(t*0.5_id)**5*3.4641016151377545870548926830117447338856105076208_id*(10._id + (Cos(t) + -1._id)**2*16.5_id +&
& (Cos(t) + -1._id)*27.5_id)*Sin(t*0.5_id)**3
endif

if(l==6.and.m1==4.and.m2==2) then
   wigd=Cos(t*0.5_id)**6*2.7386127875258305672848489140040106697637234749899_id*(6._id + (Cos(t) + -1._id)**2*16.5_id + &
&(Cos(t) + -1._id)*22._id)*Sin(t*0.5_id)**2
endif

if(l==6.and.m1==4.and.m2==3) then
   wigd=Cos(t*0.5_id)**7*1.82574185835055371152323260933600711317581564999328_id*(3._id + (Cos(t) + -1._id)*16.5_id + (C&
&os(t) + -1._id)**2*16.5_id)*Sin(t*0.5_id)
endif

if(l==6.and.m1==4.and.m2==4) then
   wigd=Cos(t*0.5_id)**8*(1._id + (Cos(t) + -1._id)*11._id + (Cos(t) + -1._id)**2*16.5_id)
endif

if(l==6.and.m1==4.and.m2==5) then
   wigd=Cos(t*0.5_id)**9*-1.1726039399558573886414075283861165701470570883529_id*(-8._id + Cos(t)*12._id)*Sin(t*0.5_id)
endif

if(l==6.and.m1==4.and.m2==6) then
   wigd=Cos(t*0.5_id)**10*8.1240384046359603604598835682660403485042040867253_id*Sin(t*0.5_id)**2
endif

if(l==6.and.m1==5.and.m2==-6) then
   wigd=Cos(t*0.5_id)*3.4641016151377545870548926830117447338856105076208_id*Sin(t*0.5_id)**11
endif

if(l==6.and.m1==5.and.m2==-5) then
   wigd=0.5_id*(10._id + Cos(t)*12._id)*Sin(t*0.5_id)**10
endif

if(l==6.and.m1==5.and.m2==-4) then
   wigd=Cos(t*0.5_id)*1.1726039399558573886414075283861165701470570883529_id*(8._id + Cos(t)*12._id)*Sin(t*0.5_id)**9
endif

if(l==6.and.m1==5.and.m2==-3) then
   wigd=Cos(t*0.5_id)**2*2.14087209644418817001948307960665673657789208526_id*(6._id + Cos(t)*12._id)*Sin(t*0.5_id)**8
endif

if(l==6.and.m1==5.and.m2==-2) then
   wigd=Cos(t*0.5_id)**3*3.2113081446662822550292246194099851048668381278901_id*(4._id + Cos(t)*12._id)*Sin(t*0.5_id)**7&
&
endif

if(l==6.and.m1==5.and.m2==-1) then
   wigd=Cos(t*0.5_id)**4*4.0620192023179801802299417841330201742521020433627_id*(2._id + Cos(t)*12._id)*Sin(t*0.5_id)**6&
&
endif

if(l==6.and.m1==5.and.m2==0) then
   wigd=Cos(t)*Cos(t*0.5_id)**5*52.649786324352732362438329844497857365255260965327_id*Sin(t*0.5_id)**5
endif

if(l==6.and.m1==5.and.m2==1) then
   wigd=Cos(t*0.5_id)**6*4.0620192023179801802299417841330201742521020433627_id*(-2._id + Cos(t)*12._id)*Sin(t*0.5_id)**&
&4
endif

if(l==6.and.m1==5.and.m2==2) then
   wigd=Cos(t*0.5_id)**7*3.2113081446662822550292246194099851048668381278901_id*(-4._id + Cos(t)*12._id)*Sin(t*0.5_id)**&
&3
endif

if(l==6.and.m1==5.and.m2==3) then
   wigd=Cos(t*0.5_id)**8*2.14087209644418817001948307960665673657789208526_id*(-6._id + Cos(t)*12._id)*Sin(t*0.5_id)**2
endif

if(l==6.and.m1==5.and.m2==4) then
   wigd=Cos(t*0.5_id)**9*1.1726039399558573886414075283861165701470570883529_id*(-8._id + Cos(t)*12._id)*Sin(t*0.5_id)
endif

if(l==6.and.m1==5.and.m2==5) then
   wigd=Cos(t*0.5_id)**10*0.5_id*(-10._id + Cos(t)*12._id)
endif

if(l==6.and.m1==5.and.m2==6) then
   wigd=Cos(t*0.5_id)**11*-3.4641016151377545870548926830117447338856105076208_id*Sin(t*0.5_id)
endif

if(l==6.and.m1==6.and.m2==-6) then
   wigd=Sin(t*0.5_id)**12
endif

if(l==6.and.m1==6.and.m2==-5) then
   wigd=Cos(t*0.5_id)*3.4641016151377545870548926830117447338856105076208_id*Sin(t*0.5_id)**11
endif

if(l==6.and.m1==6.and.m2==-4) then
   wigd=Cos(t*0.5_id)**2*8.1240384046359603604598835682660403485042040867253_id*Sin(t*0.5_id)**10
endif

if(l==6.and.m1==6.and.m2==-3) then
   wigd=Cos(t*0.5_id)**3*14.832396974191325897422794881601426121959808638195_id*Sin(t*0.5_id)**9
endif

if(l==6.and.m1==6.and.m2==-2) then
   wigd=Cos(t*0.5_id)**4*22.248595461286988846134192322402139182939712957293_id*Sin(t*0.5_id)**8
endif

if(l==6.and.m1==6.and.m2==-1) then
   wigd=Cos(t*0.5_id)**5*28.14249455894057732739378068126679768352937012047_id*Sin(t*0.5_id)**7
endif

if(l==6.and.m1==6.and.m2==0) then
   wigd=Cos(t*0.5_id)**6*30.397368307141327263337488412079290075164652923404_id*Sin(t*0.5_id)**6
endif

if(l==6.and.m1==6.and.m2==1) then
   wigd=Cos(t*0.5_id)**7*28.14249455894057732739378068126679768352937012047_id*Sin(t*0.5_id)**5
endif

if(l==6.and.m1==6.and.m2==2) then
   wigd=Cos(t*0.5_id)**8*22.248595461286988846134192322402139182939712957293_id*Sin(t*0.5_id)**4
endif

if(l==6.and.m1==6.and.m2==3) then
   wigd=Cos(t*0.5_id)**9*14.832396974191325897422794881601426121959808638195_id*Sin(t*0.5_id)**3
endif

if(l==6.and.m1==6.and.m2==4) then
   wigd=Cos(t*0.5_id)**10*8.1240384046359603604598835682660403485042040867253_id*Sin(t*0.5_id)**2
endif

if(l==6.and.m1==6.and.m2==5) then
   wigd=Cos(t*0.5_id)**11*3.4641016151377545870548926830117447338856105076208_id*Sin(t*0.5_id)
endif

if(l==6.and.m1==6.and.m2==6) then
   wigd=Cos(t*0.5_id)**12
endif

if(l==7.and.m1==-7.and.m2==-7) then
   wigd=Cos(t*0.5_id)**14
endif

if(l==7.and.m1==-7.and.m2==-6) then
   wigd=Cos(t*0.5_id)**13*-3.7416573867739413855837487323165493017560198077787_id*Sin(t*0.5_id)
endif

if(l==7.and.m1==-7.and.m2==-5) then
   wigd=Cos(t*0.5_id)**12*9.5393920141694564915262158602322654025462342525055_id*Sin(t*0.5_id)**2
endif

if(l==7.and.m1==-7.and.m2==-4) then
   wigd=Cos(t*0.5_id)**11*-19.078784028338912983052431720464530805092468505011_id*Sin(t*0.5_id)**3
endif

if(l==7.and.m1==-7.and.m2==-3) then
   wigd=Cos(t*0.5_id)**10*31.63858403911274914310629158480098308708005351898_id*Sin(t*0.5_id)**4
endif

if(l==7.and.m1==-7.and.m2==-2) then
   wigd=Cos(t*0.5_id)**9*-44.743714642394187341373897459828840314950565642978_id*Sin(t*0.5_id)**5
endif

if(l==7.and.m1==-7.and.m2==-1) then
   wigd=Cos(t*0.5_id)**8*54.799635035281028776516933082406334274967746248937_id*Sin(t*0.5_id)**6
endif

if(l==7.and.m1==-7.and.m2==0) then
   wigd=Cos(t*0.5_id)**7*-58.583274063507239250292076223037325707919457357862_id*Sin(t*0.5_id)**7
endif

if(l==7.and.m1==-7.and.m2==1) then
   wigd=Cos(t*0.5_id)**6*54.799635035281028776516933082406334274967746248937_id*Sin(t*0.5_id)**8
endif

if(l==7.and.m1==-7.and.m2==2) then
   wigd=Cos(t*0.5_id)**5*-44.743714642394187341373897459828840314950565642978_id*Sin(t*0.5_id)**9
endif

if(l==7.and.m1==-7.and.m2==3) then
   wigd=Cos(t*0.5_id)**4*31.63858403911274914310629158480098308708005351898_id*Sin(t*0.5_id)**10
endif

if(l==7.and.m1==-7.and.m2==4) then
   wigd=Cos(t*0.5_id)**3*-19.078784028338912983052431720464530805092468505011_id*Sin(t*0.5_id)**11
endif

if(l==7.and.m1==-7.and.m2==5) then
   wigd=Cos(t*0.5_id)**2*9.5393920141694564915262158602322654025462342525055_id*Sin(t*0.5_id)**12
endif

if(l==7.and.m1==-7.and.m2==6) then
   wigd=Cos(t*0.5_id)*-3.7416573867739413855837487323165493017560198077787_id*Sin(t*0.5_id)**13
endif

if(l==7.and.m1==-7.and.m2==7) then
   wigd=Sin(t*0.5_id)**14
endif

if(l==7.and.m1==-6.and.m2==-7) then
   wigd=Cos(t*0.5_id)**13*3.7416573867739413855837487323165493017560198077787_id*Sin(t*0.5_id)
endif

if(l==7.and.m1==-6.and.m2==-6) then
   wigd=Cos(t*0.5_id)**12*0.5_id*(-12._id + Cos(t)*14._id)
endif

if(l==7.and.m1==-6.and.m2==-5) then
   wigd=Cos(t*0.5_id)**11*-1.2747548783981962075070560272556954973909427365249_id*(-10._id + Cos(t)*14._id)*Sin(t*0.5_id&
&)
endif

if(l==7.and.m1==-6.and.m2==-4) then
   wigd=Cos(t*0.5_id)**10*2.5495097567963924150141120545113909947818854730498_id*(-8._id + Cos(t)*14._id)*Sin(t*0.5_id)*&
&*2
endif

if(l==7.and.m1==-6.and.m2==-3) then
   wigd=Cos(t*0.5_id)**9*-4.2278836313219407245431294806941488436657699679026_id*(-6._id + Cos(t)*14._id)*Sin(t*0.5_id)*&
&*3
endif

if(l==7.and.m1==-6.and.m2==-2) then
   wigd=Cos(t*0.5_id)**8*5.9791303715506990105649203780978083069955431425686_id*(-4._id + Cos(t)*14._id)*Sin(t*0.5_id)**&
&4
endif

if(l==7.and.m1==-6.and.m2==-1) then
   wigd=Cos(t*0.5_id)**7*-7.322909257938404906286509527879665713489932169733_id*(-2._id + Cos(t)*14._id)*Sin(t*0.5_id)**&
&5
endif

if(l==7.and.m1==-6.and.m2==0) then
   wigd=Cos(t)*Cos(t*0.5_id)**6*109.59927007056205755303386616481266854993549249787_id*Sin(t*0.5_id)**6
endif

if(l==7.and.m1==-6.and.m2==1) then
   wigd=Cos(t*0.5_id)**5*-7.322909257938404906286509527879665713489932169733_id*(2._id + Cos(t)*14._id)*Sin(t*0.5_id)**7&
&
endif

if(l==7.and.m1==-6.and.m2==2) then
   wigd=Cos(t*0.5_id)**4*5.9791303715506990105649203780978083069955431425686_id*(4._id + Cos(t)*14._id)*Sin(t*0.5_id)**8&
&
endif

if(l==7.and.m1==-6.and.m2==3) then
   wigd=Cos(t*0.5_id)**3*-4.2278836313219407245431294806941488436657699679026_id*(6._id + Cos(t)*14._id)*Sin(t*0.5_id)**&
&9
endif

if(l==7.and.m1==-6.and.m2==4) then
   wigd=Cos(t*0.5_id)**2*2.5495097567963924150141120545113909947818854730498_id*(8._id + Cos(t)*14._id)*Sin(t*0.5_id)**1&
&0
endif

if(l==7.and.m1==-6.and.m2==5) then
   wigd=Cos(t*0.5_id)*-1.2747548783981962075070560272556954973909427365249_id*(10._id + Cos(t)*14._id)*Sin(t*0.5_id)**11&
&
endif

if(l==7.and.m1==-6.and.m2==6) then
   wigd=0.5_id*(12._id + Cos(t)*14._id)*Sin(t*0.5_id)**12
endif

if(l==7.and.m1==-6.and.m2==7) then
   wigd=Cos(t*0.5_id)*-3.7416573867739413855837487323165493017560198077787_id*Sin(t*0.5_id)**13
endif

if(l==7.and.m1==-5.and.m2==-7) then
   wigd=Cos(t*0.5_id)**12*9.5393920141694564915262158602322654025462342525055_id*Sin(t*0.5_id)**2
endif

if(l==7.and.m1==-5.and.m2==-6) then
   wigd=Cos(t*0.5_id)**11*1.2747548783981962075070560272556954973909427365249_id*(-10._id + Cos(t)*14._id)*Sin(t*0.5_id)&
&
endif

if(l==7.and.m1==-5.and.m2==-5) then
   wigd=Cos(t*0.5_id)**10*(1._id + (Cos(t) + -1._id)*13._id + (Cos(t) + -1._id)**2*22.75_id)
endif

if(l==7.and.m1==-5.and.m2==-4) then
   wigd=Cos(t*0.5_id)**9*-2._id*(3._id + (Cos(t) + -1._id)*19.5_id + (Cos(t) + -1._id)**2*22.75_id)*Sin(t*0.5_id)
endif

if(l==7.and.m1==-5.and.m2==-3) then
   wigd=Cos(t*0.5_id)**8*3.3166247903553998491149327366706866839270885455894_id*(6._id + (Cos(t) + -1._id)**2*22.75_id +&
& (Cos(t) + -1._id)*26._id)*Sin(t*0.5_id)**2
endif

if(l==7.and.m1==-5.and.m2==-2) then
   wigd=Cos(t*0.5_id)**7*-4.6904157598234295545656301135444662805882283534117_id*(10._id + (Cos(t) + -1._id)**2*22.75_id&
& + (Cos(t) + -1._id)*32.5_id)*Sin(t*0.5_id)**3
endif

if(l==7.and.m1==-5.and.m2==-1) then
   wigd=Cos(t*0.5_id)**6*5.7445626465380286598506114682189293182202644579828_id*(15._id + (Cos(t) + -1._id)**2*22.75_id &
&+ (Cos(t) + -1._id)*39._id)*Sin(t*0.5_id)**4
endif

if(l==7.and.m1==-5.and.m2==0) then
   wigd=Cos(t*0.5_id)**5*-6.1411957886299077254982887358519179988468052765271_id*(21._id + (Cos(t) + -1._id)**2*22.75_id&
& + (Cos(t) + -1._id)*45.5_id)*Sin(t*0.5_id)**5
endif

if(l==7.and.m1==-5.and.m2==1) then
   wigd=Cos(t*0.5_id)**4*5.7445626465380286598506114682189293182202644579828_id*((Cos(t) + -1._id)**2*22.75_id + 28._id &
&+ (Cos(t) + -1._id)*52._id)*Sin(t*0.5_id)**6
endif

if(l==7.and.m1==-5.and.m2==2) then
   wigd=Cos(t*0.5_id)**3*-4.6904157598234295545656301135444662805882283534117_id*((Cos(t) + -1._id)**2*22.75_id + 36._id&
& + (Cos(t) + -1._id)*58.5_id)*Sin(t*0.5_id)**7
endif

if(l==7.and.m1==-5.and.m2==3) then
   wigd=Cos(t*0.5_id)**2*3.3166247903553998491149327366706866839270885455894_id*((Cos(t) + -1._id)**2*22.75_id + 45._id &
&+ (Cos(t) + -1._id)*65._id)*Sin(t*0.5_id)**8
endif

if(l==7.and.m1==-5.and.m2==4) then
   wigd=Cos(t*0.5_id)*-2._id*((Cos(t) + -1._id)**2*22.75_id + 55._id + (Cos(t) + -1._id)*71.5_id)*Sin(t*0.5_id)**9
endif

if(l==7.and.m1==-5.and.m2==5) then
   wigd=((Cos(t) + -1._id)**2*22.75_id + 66._id + (Cos(t) + -1._id)*78._id)*Sin(t*0.5_id)**10
endif

if(l==7.and.m1==-5.and.m2==6) then
   wigd=Cos(t*0.5_id)*-1.2747548783981962075070560272556954973909427365249_id*(10._id + Cos(t)*14._id)*Sin(t*0.5_id)**11&
&
endif

if(l==7.and.m1==-5.and.m2==7) then
   wigd=Cos(t*0.5_id)**2*9.5393920141694564915262158602322654025462342525055_id*Sin(t*0.5_id)**12
endif

if(l==7.and.m1==-4.and.m2==-7) then
   wigd=Cos(t*0.5_id)**11*19.078784028338912983052431720464530805092468505011_id*Sin(t*0.5_id)**3
endif

if(l==7.and.m1==-4.and.m2==-6) then
   wigd=Cos(t*0.5_id)**10*2.5495097567963924150141120545113909947818854730498_id*(-8._id + Cos(t)*14._id)*Sin(t*0.5_id)*&
&*2
endif

if(l==7.and.m1==-4.and.m2==-5) then
   wigd=Cos(t*0.5_id)**9*2._id*(3._id + (Cos(t) + -1._id)*19.5_id + (Cos(t) + -1._id)**2*22.75_id)*Sin(t*0.5_id)
endif

if(l==7.and.m1==-4.and.m2==-4) then
   wigd=Cos(t*0.5_id)**8*(1._id + (Cos(t) + -1._id)*18._id + (Cos(t) + -1._id)**3*45.5_id + (Cos(t) + -1._id)**2*58.5_id&
&)
endif

if(l==7.and.m1==-4.and.m2==-3) then
   wigd=Cos(t*0.5_id)**7*-1.6583123951776999245574663683353433419635442727947_id*(4._id + (Cos(t) + -1._id)*36._id + (Co&
&s(t) + -1._id)**3*45.5_id + (Cos(t) + -1._id)**2*78._id)*Sin(t*0.5_id)
endif

if(l==7.and.m1==-4.and.m2==-2) then
   wigd=Cos(t*0.5_id)**6*2.3452078799117147772828150567722331402941141767059_id*(10._id + (Cos(t) + -1._id)**3*45.5_id +&
& (Cos(t) + -1._id)*60._id + (Cos(t) + -1._id)**2*97.5_id)*Sin(t*0.5_id)**2
endif

if(l==7.and.m1==-4.and.m2==-1) then
   wigd=Cos(t*0.5_id)**5*-2.8722813232690143299253057341094646591101322289914_id*(20._id + (Cos(t) + -1._id)**3*45.5_id &
&+ (Cos(t) + -1._id)*90._id + (Cos(t) + -1._id)**2*117._id)*Sin(t*0.5_id)**3
endif

if(l==7.and.m1==-4.and.m2==0) then
   wigd=Cos(t*0.5_id)**4*3.0705978943149538627491443679259589994234026382635_id*(35._id + (Cos(t) + -1._id)**3*45.5_id +&
& (Cos(t) + -1._id)*126._id + (Cos(t) + -1._id)**2*136.5_id)*Sin(t*0.5_id)**4
endif

if(l==7.and.m1==-4.and.m2==1) then
   wigd=Cos(t*0.5_id)**3*-2.8722813232690143299253057341094646591101322289914_id*((Cos(t) + -1._id)**3*45.5_id + 56._id &
&+ (Cos(t) + -1._id)**2*156._id + (Cos(t) + -1._id)*168._id)*Sin(t*0.5_id)**5
endif

if(l==7.and.m1==-4.and.m2==2) then
   wigd=Cos(t*0.5_id)**2*2.3452078799117147772828150567722331402941141767059_id*((Cos(t) + -1._id)**3*45.5_id + 84._id +&
& (Cos(t) + -1._id)**2*175.5_id + (Cos(t) + -1._id)*216._id)*Sin(t*0.5_id)**6
endif

if(l==7.and.m1==-4.and.m2==3) then
   wigd=Cos(t*0.5_id)*-1.6583123951776999245574663683353433419635442727947_id*((Cos(t) + -1._id)**3*45.5_id + 120._id + &
&(Cos(t) + -1._id)**2*195._id + (Cos(t) + -1._id)*270._id)*Sin(t*0.5_id)**7
endif

if(l==7.and.m1==-4.and.m2==4) then
   wigd=((Cos(t) + -1._id)**3*45.5_id + 165._id + (Cos(t) + -1._id)**2*214.5_id + (Cos(t) + -1._id)*330._id)*Sin(t*0.5_i&
&d)**8
endif

if(l==7.and.m1==-4.and.m2==5) then
   wigd=Cos(t*0.5_id)*-2._id*((Cos(t) + -1._id)**2*22.75_id + 55._id + (Cos(t) + -1._id)*71.5_id)*Sin(t*0.5_id)**9
endif

if(l==7.and.m1==-4.and.m2==6) then
   wigd=Cos(t*0.5_id)**2*2.5495097567963924150141120545113909947818854730498_id*(8._id + Cos(t)*14._id)*Sin(t*0.5_id)**1&
&0
endif

if(l==7.and.m1==-4.and.m2==7) then
   wigd=Cos(t*0.5_id)**3*-19.078784028338912983052431720464530805092468505011_id*Sin(t*0.5_id)**11
endif

if(l==7.and.m1==-3.and.m2==-7) then
   wigd=Cos(t*0.5_id)**10*31.63858403911274914310629158480098308708005351898_id*Sin(t*0.5_id)**4
endif

if(l==7.and.m1==-3.and.m2==-6) then
   wigd=Cos(t*0.5_id)**9*4.2278836313219407245431294806941488436657699679026_id*(-6._id + Cos(t)*14._id)*Sin(t*0.5_id)**&
&3
endif

if(l==7.and.m1==-3.and.m2==-5) then
   wigd=Cos(t*0.5_id)**8*3.3166247903553998491149327366706866839270885455894_id*(6._id + (Cos(t) + -1._id)**2*22.75_id +&
& (Cos(t) + -1._id)*26._id)*Sin(t*0.5_id)**2
endif

if(l==7.and.m1==-3.and.m2==-4) then
   wigd=Cos(t*0.5_id)**7*1.6583123951776999245574663683353433419635442727947_id*(4._id + (Cos(t) + -1._id)*36._id + (Cos&
&(t) + -1._id)**3*45.5_id + (Cos(t) + -1._id)**2*78._id)*Sin(t*0.5_id)
endif

if(l==7.and.m1==-3.and.m2==-3) then
   wigd=Cos(t*0.5_id)**6*(1._id + (Cos(t) + -1._id)*22._id + (Cos(t) + -1._id)**4*62.5625_id + (Cos(t) + -1._id)**2*99._&
&id + (Cos(t) + -1._id)**3*143._id)
endif

if(l==7.and.m1==-3.and.m2==-2) then
   wigd=Cos(t*0.5_id)**5*-1.4142135623730950488016887242096980785696718753769_id*(5._id + (Cos(t) + -1._id)*55._id + (Co&
&s(t) + -1._id)**4*62.5625_id + (Cos(t) + -1._id)**2*165._id + (Cos(t) + -1._id)**3*178.75_id)*Sin(t*0.5_id)
endif

if(l==7.and.m1==-3.and.m2==-1) then
   wigd=Cos(t*0.5_id)**4*1.73205080756887729352744634150587236694280525381038_id*(15._id + (Cos(t) + -1._id)**4*62.5625_&
&id + (Cos(t) + -1._id)*110._id + (Cos(t) + -1._id)**3*214.5_id + (Cos(t) + -1._id)**2*247.5_id)*Sin(t*0.5_id)**2
endif

if(l==7.and.m1==-3.and.m2==0) then
   wigd=Cos(t*0.5_id)**3*-1.8516401995451029231331335531679990450586298020167_id*(35._id + (Cos(t) + -1._id)**4*62.5625_&
&id + (Cos(t) + -1._id)*192.5_id + (Cos(t) + -1._id)**3*250.25_id + (Cos(t) + -1._id)**2*346.5_id)*Sin(t*0.5_id)**3
endif

if(l==7.and.m1==-3.and.m2==1) then
   wigd=Cos(t*0.5_id)**2*1.73205080756887729352744634150587236694280525381038_id*((Cos(t) + -1._id)**4*62.5625_id + 70._&
&id + (Cos(t) + -1._id)**3*286._id + (Cos(t) + -1._id)*308._id + (Cos(t) + -1._id)**2*462._id)*Sin(t*0.5_id)**4
endif

if(l==7.and.m1==-3.and.m2==2) then
   wigd=Cos(t*0.5_id)*-1.4142135623730950488016887242096980785696718753769_id*((Cos(t) + -1._id)**4*62.5625_id + 126._id&
& + (Cos(t) + -1._id)**3*321.75_id + (Cos(t) + -1._id)*462._id + (Cos(t) + -1._id)**2*594._id)*Sin(t*0.5_id)**5
endif

if(l==7.and.m1==-3.and.m2==3) then
   wigd=((Cos(t) + -1._id)**4*62.5625_id + 210._id + (Cos(t) + -1._id)**3*357.5_id + (Cos(t) + -1._id)*660._id + (Cos(t)&
& + -1._id)**2*742.5_id)*Sin(t*0.5_id)**6
endif

if(l==7.and.m1==-3.and.m2==4) then
   wigd=Cos(t*0.5_id)*-1.6583123951776999245574663683353433419635442727947_id*((Cos(t) + -1._id)**3*45.5_id + 120._id + &
&(Cos(t) + -1._id)**2*195._id + (Cos(t) + -1._id)*270._id)*Sin(t*0.5_id)**7
endif

if(l==7.and.m1==-3.and.m2==5) then
   wigd=Cos(t*0.5_id)**2*3.3166247903553998491149327366706866839270885455894_id*((Cos(t) + -1._id)**2*22.75_id + 45._id &
&+ (Cos(t) + -1._id)*65._id)*Sin(t*0.5_id)**8
endif

if(l==7.and.m1==-3.and.m2==6) then
   wigd=Cos(t*0.5_id)**3*-4.2278836313219407245431294806941488436657699679026_id*(6._id + Cos(t)*14._id)*Sin(t*0.5_id)**&
&9
endif

if(l==7.and.m1==-3.and.m2==7) then
   wigd=Cos(t*0.5_id)**4*31.63858403911274914310629158480098308708005351898_id*Sin(t*0.5_id)**10
endif

if(l==7.and.m1==-2.and.m2==-7) then
   wigd=Cos(t*0.5_id)**9*44.743714642394187341373897459828840314950565642978_id*Sin(t*0.5_id)**5
endif

if(l==7.and.m1==-2.and.m2==-6) then
   wigd=Cos(t*0.5_id)**8*5.9791303715506990105649203780978083069955431425686_id*(-4._id + Cos(t)*14._id)*Sin(t*0.5_id)**&
&4
endif

if(l==7.and.m1==-2.and.m2==-5) then
   wigd=Cos(t*0.5_id)**7*4.6904157598234295545656301135444662805882283534117_id*(10._id + (Cos(t) + -1._id)**2*22.75_id &
&+ (Cos(t) + -1._id)*32.5_id)*Sin(t*0.5_id)**3
endif

if(l==7.and.m1==-2.and.m2==-4) then
   wigd=Cos(t*0.5_id)**6*2.3452078799117147772828150567722331402941141767059_id*(10._id + (Cos(t) + -1._id)**3*45.5_id +&
& (Cos(t) + -1._id)*60._id + (Cos(t) + -1._id)**2*97.5_id)*Sin(t*0.5_id)**2
endif

if(l==7.and.m1==-2.and.m2==-3) then
   wigd=Cos(t*0.5_id)**5*1.41421356237309504880168872420969807856967187537695_id*(5._id + (Cos(t) + -1._id)*55._id + (Co&
&s(t) + -1._id)**4*62.5625_id + (Cos(t) + -1._id)**2*165._id + (Cos(t) + -1._id)**3*178.75_id)*Sin(t*0.5_id)
endif

if(l==7.and.m1==-2.and.m2==-2) then
   wigd=Cos(t*0.5_id)**4*(1._id + (Cos(t) + -1._id)*25._id + (Cos(t) + -1._id)**5*62.5625_id + (Cos(t) + -1._id)**2*137.&
&5_id + (Cos(t) + -1._id)**4*223.4375_id + (Cos(t) + -1._id)**3*275._id)
endif

if(l==7.and.m1==-2.and.m2==-1) then
   wigd=Cos(t*0.5_id)**3*-1.2247448713915890490986420373529456959829737403283_id*(6._id + (Cos(t) + -1._id)**5*62.5625_i&
&d + (Cos(t) + -1._id)*75._id + (Cos(t) + -1._id)**4*268.125_id + (Cos(t) + -1._id)**2*275._id + (Cos(t) + -1._id)**3*412.5_id)*Sin(t*0.5_id)
endif

if(l==7.and.m1==-2.and.m2==0) then
   wigd=Cos(t*0.5_id)**2*1.3093073414159542875965849124937167111384161647908_id*(21._id + (Cos(t) + -1._id)**5*62.5625_i&
&d + (Cos(t) + -1._id)*175._id + (Cos(t) + -1._id)**4*312.8125_id + (Cos(t) + -1._id)**2*481.25_id + (Cos(t) + -1._id)**3*577.5_id)*Sin(t*0.5_id)**2
endif

if(l==7.and.m1==-2.and.m2==1) then
   wigd=Cos(t*0.5_id)*-1.2247448713915890490986420373529456959829737403283_id*(56._id + (Cos(t) + -1._id)**5*62.5625_id &
&+ (Cos(t) + -1._id)*350._id + (Cos(t) + -1._id)**4*357.5_id + (Cos(t) + -1._id)**2*770._id + (Cos(t) + -1._id)**3*770._id)*Sin(t*0.5_id)**3
endif

if(l==7.and.m1==-2.and.m2==2) then
   wigd=((Cos(t) + -1._id)**5*62.5625_id + 126._id + (Cos(t) + -1._id)**4*402.1875_id + (Cos(t) + -1._id)*630._id + (Cos&
&(t) + -1._id)**3*990._id + (Cos(t) + -1._id)**2*1155._id)*Sin(t*0.5_id)**4
endif

if(l==7.and.m1==-2.and.m2==3) then
   wigd=Cos(t*0.5_id)*-1.4142135623730950488016887242096980785696718753769_id*((Cos(t) + -1._id)**4*62.5625_id + 126._id&
& + (Cos(t) + -1._id)**3*321.75_id + (Cos(t) + -1._id)*462._id + (Cos(t) + -1._id)**2*594._id)*Sin(t*0.5_id)**5
endif

if(l==7.and.m1==-2.and.m2==4) then
   wigd=Cos(t*0.5_id)**2*2.3452078799117147772828150567722331402941141767059_id*((Cos(t) + -1._id)**3*45.5_id + 84._id +&
& (Cos(t) + -1._id)**2*175.5_id + (Cos(t) + -1._id)*216._id)*Sin(t*0.5_id)**6
endif

if(l==7.and.m1==-2.and.m2==5) then
   wigd=Cos(t*0.5_id)**3*-4.6904157598234295545656301135444662805882283534117_id*((Cos(t) + -1._id)**2*22.75_id + 36._id&
& + (Cos(t) + -1._id)*58.5_id)*Sin(t*0.5_id)**7
endif

if(l==7.and.m1==-2.and.m2==6) then
   wigd=Cos(t*0.5_id)**4*5.9791303715506990105649203780978083069955431425686_id*(4._id + Cos(t)*14._id)*Sin(t*0.5_id)**8&
&
endif

if(l==7.and.m1==-2.and.m2==7) then
   wigd=Cos(t*0.5_id)**5*-44.743714642394187341373897459828840314950565642978_id*Sin(t*0.5_id)**9
endif

if(l==7.and.m1==-1.and.m2==-7) then
   wigd=Cos(t*0.5_id)**8*54.799635035281028776516933082406334274967746248937_id*Sin(t*0.5_id)**6
endif

if(l==7.and.m1==-1.and.m2==-6) then
   wigd=Cos(t*0.5_id)**7*7.322909257938404906286509527879665713489932169733_id*(-2._id + Cos(t)*14._id)*Sin(t*0.5_id)**5&
&
endif

if(l==7.and.m1==-1.and.m2==-5) then
   wigd=Cos(t*0.5_id)**6*5.7445626465380286598506114682189293182202644579828_id*(15._id + (Cos(t) + -1._id)**2*22.75_id &
&+ (Cos(t) + -1._id)*39._id)*Sin(t*0.5_id)**4
endif

if(l==7.and.m1==-1.and.m2==-4) then
   wigd=Cos(t*0.5_id)**5*2.8722813232690143299253057341094646591101322289914_id*(20._id + (Cos(t) + -1._id)**3*45.5_id +&
& (Cos(t) + -1._id)*90._id + (Cos(t) + -1._id)**2*117._id)*Sin(t*0.5_id)**3
endif

if(l==7.and.m1==-1.and.m2==-3) then
   wigd=Cos(t*0.5_id)**4*1.73205080756887729352744634150587236694280525381038_id*(15._id + (Cos(t) + -1._id)**4*62.5625_&
&id + (Cos(t) + -1._id)*110._id + (Cos(t) + -1._id)**3*214.5_id + (Cos(t) + -1._id)**2*247.5_id)*Sin(t*0.5_id)**2
endif

if(l==7.and.m1==-1.and.m2==-2) then
   wigd=Cos(t*0.5_id)**3*1.22474487139158904909864203735294569598297374032834_id*(6._id + (Cos(t) + -1._id)**5*62.5625_i&
&d + (Cos(t) + -1._id)*75._id + (Cos(t) + -1._id)**4*268.125_id + (Cos(t) + -1._id)**2*275._id + (Cos(t) + -1._id)**3*412.5_id)*Sin(t*0.5_id)
endif

if(l==7.and.m1==-1.and.m2==-1) then
   wigd=Cos(t*0.5_id)**2*(1._id + (Cos(t) + -1._id)*27._id + (Cos(t) + -1._id)**6*46.921875_id + (Cos(t) + -1._id)**2*16&
&8.75_id + (Cos(t) + -1._id)**5*241.3125_id + (Cos(t) + -1._id)**3*412.5_id + (Cos(t) + -1._id)**4*464.0625_id)
endif

if(l==7.and.m1==-1.and.m2==0) then
   wigd=Cos(t*0.5_id)*-1.0690449676496975387382139235190140862160056593654_id*(7._id + (Cos(t) + -1._id)**6*46.921875_id&
& + (Cos(t) + -1._id)*94.5_id + (Cos(t) + -1._id)**5*281.53125_id + (Cos(t) + -1._id)**2*393.75_id + (Cos(t) + -1._id)**4*649.6875_id + (Cos(t) + -1._id)**3*721.875_id)*Sin(t*0.5_id)
endif

if(l==7.and.m1==-1.and.m2==1) then
   wigd=(28._id + (Cos(t) + -1._id)**6*46.921875_id + (Cos(t) + -1._id)*252._id + (Cos(t) + -1._id)**5*321.75_id + (Cos(&
&t) + -1._id)**2*787.5_id + (Cos(t) + -1._id)**4*866.25_id + (Cos(t) + -1._id)**3*1155._id)*Sin(t*0.5_id)**2
endif

if(l==7.and.m1==-1.and.m2==2) then
   wigd=Cos(t*0.5_id)*-1.2247448713915890490986420373529456959829737403283_id*(56._id + (Cos(t) + -1._id)**5*62.5625_id &
&+ (Cos(t) + -1._id)*350._id + (Cos(t) + -1._id)**4*357.5_id + (Cos(t) + -1._id)**2*770._id + (Cos(t) + -1._id)**3*770._id)*Sin(t*0.5_id)**3
endif

if(l==7.and.m1==-1.and.m2==3) then
   wigd=Cos(t*0.5_id)**2*1.73205080756887729352744634150587236694280525381038_id*((Cos(t) + -1._id)**4*62.5625_id + 70._&
&id + (Cos(t) + -1._id)**3*286._id + (Cos(t) + -1._id)*308._id + (Cos(t) + -1._id)**2*462._id)*Sin(t*0.5_id)**4
endif

if(l==7.and.m1==-1.and.m2==4) then
   wigd=Cos(t*0.5_id)**3*-2.8722813232690143299253057341094646591101322289914_id*((Cos(t) + -1._id)**3*45.5_id + 56._id &
&+ (Cos(t) + -1._id)**2*156._id + (Cos(t) + -1._id)*168._id)*Sin(t*0.5_id)**5
endif

if(l==7.and.m1==-1.and.m2==5) then
   wigd=Cos(t*0.5_id)**4*5.7445626465380286598506114682189293182202644579828_id*((Cos(t) + -1._id)**2*22.75_id + 28._id &
&+ (Cos(t) + -1._id)*52._id)*Sin(t*0.5_id)**6
endif

if(l==7.and.m1==-1.and.m2==6) then
   wigd=Cos(t*0.5_id)**5*-7.322909257938404906286509527879665713489932169733_id*(2._id + Cos(t)*14._id)*Sin(t*0.5_id)**7&
&
endif

if(l==7.and.m1==-1.and.m2==7) then
   wigd=Cos(t*0.5_id)**6*54.799635035281028776516933082406334274967746248937_id*Sin(t*0.5_id)**8
endif

if(l==7.and.m1==0.and.m2==-7) then
   wigd=Cos(t*0.5_id)**7*58.583274063507239250292076223037325707919457357862_id*Sin(t*0.5_id)**7
endif

if(l==7.and.m1==0.and.m2==-6) then
   wigd=Cos(t)*Cos(t*0.5_id)**6*109.59927007056205755303386616481266854993549249787_id*Sin(t*0.5_id)**6
endif

if(l==7.and.m1==0.and.m2==-5) then
   wigd=Cos(t*0.5_id)**5*6.1411957886299077254982887358519179988468052765271_id*(21._id + (Cos(t) + -1._id)**2*22.75_id &
&+ (Cos(t) + -1._id)*45.5_id)*Sin(t*0.5_id)**5
endif

if(l==7.and.m1==0.and.m2==-4) then
   wigd=Cos(t*0.5_id)**4*3.0705978943149538627491443679259589994234026382635_id*(35._id + (Cos(t) + -1._id)**3*45.5_id +&
& (Cos(t) + -1._id)*126._id + (Cos(t) + -1._id)**2*136.5_id)*Sin(t*0.5_id)**4
endif

if(l==7.and.m1==0.and.m2==-3) then
   wigd=Cos(t*0.5_id)**3*1.8516401995451029231331335531679990450586298020167_id*(35._id + (Cos(t) + -1._id)**4*62.5625_i&
&d + (Cos(t) + -1._id)*192.5_id + (Cos(t) + -1._id)**3*250.25_id + (Cos(t) + -1._id)**2*346.5_id)*Sin(t*0.5_id)**3
endif

if(l==7.and.m1==0.and.m2==-2) then
   wigd=Cos(t*0.5_id)**2*1.3093073414159542875965849124937167111384161647908_id*(21._id + (Cos(t) + -1._id)**5*62.5625_i&
&d + (Cos(t) + -1._id)*175._id + (Cos(t) + -1._id)**4*312.8125_id + (Cos(t) + -1._id)**2*481.25_id + (Cos(t) + -1._id)**3*577.5_id)*Sin(t*0.5_id)**2
endif

if(l==7.and.m1==0.and.m2==-1) then
   wigd=Cos(t*0.5_id)*1.0690449676496975387382139235190140862160056593654_id*(7._id + (Cos(t) + -1._id)**6*46.921875_id &
&+ (Cos(t) + -1._id)*94.5_id + (Cos(t) + -1._id)**5*281.53125_id + (Cos(t) + -1._id)**2*393.75_id + (Cos(t) + -1._id)**4*649.6875_id + (Cos(t) + -1._id)**3*721.875_id)*Sin(t*0.5_id)
endif

if(l==7.and.m1==0.and.m2==0) then
   wigd=0.0625_id*(Cos(t)**5*-693._id + Cos(t)*-35._id + Cos(t)**3*315._id + Cos(t)**7*429._id)
endif

if(l==7.and.m1==0.and.m2==1) then
   wigd=Cos(t*0.5_id)*-1.0690449676496975387382139235190140862160056593654_id*(7._id + (Cos(t) + -1._id)**6*46.921875_id&
& + (Cos(t) + -1._id)*94.5_id + (Cos(t) + -1._id)**5*281.53125_id + (Cos(t) + -1._id)**2*393.75_id + (Cos(t) + -1._id)**4*649.6875_id + (Cos(t) + -1._id)**3*721.875_id)*Sin(t*0.5_id)
endif

if(l==7.and.m1==0.and.m2==2) then
   wigd=Cos(t*0.5_id)**2*1.3093073414159542875965849124937167111384161647908_id*(21._id + (Cos(t) + -1._id)**5*62.5625_i&
&d + (Cos(t) + -1._id)*175._id + (Cos(t) + -1._id)**4*312.8125_id + (Cos(t) + -1._id)**2*481.25_id + (Cos(t) + -1._id)**3*577.5_id)*Sin(t*0.5_id)**2
endif

if(l==7.and.m1==0.and.m2==3) then
   wigd=Cos(t*0.5_id)**3*-1.8516401995451029231331335531679990450586298020167_id*(35._id + (Cos(t) + -1._id)**4*62.5625_&
&id + (Cos(t) + -1._id)*192.5_id + (Cos(t) + -1._id)**3*250.25_id + (Cos(t) + -1._id)**2*346.5_id)*Sin(t*0.5_id)**3
endif

if(l==7.and.m1==0.and.m2==4) then
   wigd=Cos(t*0.5_id)**4*3.0705978943149538627491443679259589994234026382635_id*(35._id + (Cos(t) + -1._id)**3*45.5_id +&
& (Cos(t) + -1._id)*126._id + (Cos(t) + -1._id)**2*136.5_id)*Sin(t*0.5_id)**4
endif

if(l==7.and.m1==0.and.m2==5) then
   wigd=Cos(t*0.5_id)**5*-6.1411957886299077254982887358519179988468052765271_id*(21._id + (Cos(t) + -1._id)**2*22.75_id&
& + (Cos(t) + -1._id)*45.5_id)*Sin(t*0.5_id)**5
endif

if(l==7.and.m1==0.and.m2==6) then
   wigd=Cos(t)*Cos(t*0.5_id)**6*109.59927007056205755303386616481266854993549249787_id*Sin(t*0.5_id)**6
endif

if(l==7.and.m1==0.and.m2==7) then
   wigd=Cos(t*0.5_id)**7*-58.583274063507239250292076223037325707919457357862_id*Sin(t*0.5_id)**7
endif

if(l==7.and.m1==1.and.m2==-7) then
   wigd=Cos(t*0.5_id)**6*54.799635035281028776516933082406334274967746248937_id*Sin(t*0.5_id)**8
endif

if(l==7.and.m1==1.and.m2==-6) then
   wigd=Cos(t*0.5_id)**5*7.322909257938404906286509527879665713489932169733_id*(2._id + Cos(t)*14._id)*Sin(t*0.5_id)**7
endif

if(l==7.and.m1==1.and.m2==-5) then
   wigd=Cos(t*0.5_id)**4*5.7445626465380286598506114682189293182202644579828_id*((Cos(t) + -1._id)**2*22.75_id + 28._id &
&+ (Cos(t) + -1._id)*52._id)*Sin(t*0.5_id)**6
endif

if(l==7.and.m1==1.and.m2==-4) then
   wigd=Cos(t*0.5_id)**3*2.8722813232690143299253057341094646591101322289914_id*((Cos(t) + -1._id)**3*45.5_id + 56._id +&
& (Cos(t) + -1._id)**2*156._id + (Cos(t) + -1._id)*168._id)*Sin(t*0.5_id)**5
endif

if(l==7.and.m1==1.and.m2==-3) then
   wigd=Cos(t*0.5_id)**2*1.73205080756887729352744634150587236694280525381038_id*((Cos(t) + -1._id)**4*62.5625_id + 70._&
&id + (Cos(t) + -1._id)**3*286._id + (Cos(t) + -1._id)*308._id + (Cos(t) + -1._id)**2*462._id)*Sin(t*0.5_id)**4
endif

if(l==7.and.m1==1.and.m2==-2) then
   wigd=Cos(t*0.5_id)*1.22474487139158904909864203735294569598297374032834_id*(56._id + (Cos(t) + -1._id)**5*62.5625_id &
&+ (Cos(t) + -1._id)*350._id + (Cos(t) + -1._id)**4*357.5_id + (Cos(t) + -1._id)**2*770._id + (Cos(t) + -1._id)**3*770._id)*Sin(t*0.5_id)**3
endif

if(l==7.and.m1==1.and.m2==-1) then
   wigd=(28._id + (Cos(t) + -1._id)**6*46.921875_id + (Cos(t) + -1._id)*252._id + (Cos(t) + -1._id)**5*321.75_id + (Cos(&
&t) + -1._id)**2*787.5_id + (Cos(t) + -1._id)**4*866.25_id + (Cos(t) + -1._id)**3*1155._id)*Sin(t*0.5_id)**2
endif

if(l==7.and.m1==1.and.m2==0) then
   wigd=Cos(t*0.5_id)*1.0690449676496975387382139235190140862160056593654_id*(7._id + (Cos(t) + -1._id)**6*46.921875_id &
&+ (Cos(t) + -1._id)*94.5_id + (Cos(t) + -1._id)**5*281.53125_id + (Cos(t) + -1._id)**2*393.75_id + (Cos(t) + -1._id)**4*649.6875_id + (Cos(t) + -1._id)**3*721.875_id)*Sin(t*0.5_id)
endif

if(l==7.and.m1==1.and.m2==1) then
   wigd=Cos(t*0.5_id)**2*(1._id + (Cos(t) + -1._id)*27._id + (Cos(t) + -1._id)**6*46.921875_id + (Cos(t) + -1._id)**2*16&
&8.75_id + (Cos(t) + -1._id)**5*241.3125_id + (Cos(t) + -1._id)**3*412.5_id + (Cos(t) + -1._id)**4*464.0625_id)
endif

if(l==7.and.m1==1.and.m2==2) then
   wigd=Cos(t*0.5_id)**3*-1.2247448713915890490986420373529456959829737403283_id*(6._id + (Cos(t) + -1._id)**5*62.5625_i&
&d + (Cos(t) + -1._id)*75._id + (Cos(t) + -1._id)**4*268.125_id + (Cos(t) + -1._id)**2*275._id + (Cos(t) + -1._id)**3*412.5_id)*Sin(t*0.5_id)
endif

if(l==7.and.m1==1.and.m2==3) then
   wigd=Cos(t*0.5_id)**4*1.73205080756887729352744634150587236694280525381038_id*(15._id + (Cos(t) + -1._id)**4*62.5625_&
&id + (Cos(t) + -1._id)*110._id + (Cos(t) + -1._id)**3*214.5_id + (Cos(t) + -1._id)**2*247.5_id)*Sin(t*0.5_id)**2
endif

if(l==7.and.m1==1.and.m2==4) then
   wigd=Cos(t*0.5_id)**5*-2.8722813232690143299253057341094646591101322289914_id*(20._id + (Cos(t) + -1._id)**3*45.5_id &
&+ (Cos(t) + -1._id)*90._id + (Cos(t) + -1._id)**2*117._id)*Sin(t*0.5_id)**3
endif

if(l==7.and.m1==1.and.m2==5) then
   wigd=Cos(t*0.5_id)**6*5.7445626465380286598506114682189293182202644579828_id*(15._id + (Cos(t) + -1._id)**2*22.75_id &
&+ (Cos(t) + -1._id)*39._id)*Sin(t*0.5_id)**4
endif

if(l==7.and.m1==1.and.m2==6) then
   wigd=Cos(t*0.5_id)**7*-7.322909257938404906286509527879665713489932169733_id*(-2._id + Cos(t)*14._id)*Sin(t*0.5_id)**&
&5
endif

if(l==7.and.m1==1.and.m2==7) then
   wigd=Cos(t*0.5_id)**8*54.799635035281028776516933082406334274967746248937_id*Sin(t*0.5_id)**6
endif

if(l==7.and.m1==2.and.m2==-7) then
   wigd=Cos(t*0.5_id)**5*44.743714642394187341373897459828840314950565642978_id*Sin(t*0.5_id)**9
endif

if(l==7.and.m1==2.and.m2==-6) then
   wigd=Cos(t*0.5_id)**4*5.9791303715506990105649203780978083069955431425686_id*(4._id + Cos(t)*14._id)*Sin(t*0.5_id)**8&
&
endif

if(l==7.and.m1==2.and.m2==-5) then
   wigd=Cos(t*0.5_id)**3*4.6904157598234295545656301135444662805882283534117_id*((Cos(t) + -1._id)**2*22.75_id + 36._id &
&+ (Cos(t) + -1._id)*58.5_id)*Sin(t*0.5_id)**7
endif

if(l==7.and.m1==2.and.m2==-4) then
   wigd=Cos(t*0.5_id)**2*2.3452078799117147772828150567722331402941141767059_id*((Cos(t) + -1._id)**3*45.5_id + 84._id +&
& (Cos(t) + -1._id)**2*175.5_id + (Cos(t) + -1._id)*216._id)*Sin(t*0.5_id)**6
endif

if(l==7.and.m1==2.and.m2==-3) then
   wigd=Cos(t*0.5_id)*1.41421356237309504880168872420969807856967187537695_id*((Cos(t) + -1._id)**4*62.5625_id + 126._id&
& + (Cos(t) + -1._id)**3*321.75_id + (Cos(t) + -1._id)*462._id + (Cos(t) + -1._id)**2*594._id)*Sin(t*0.5_id)**5
endif

if(l==7.and.m1==2.and.m2==-2) then
   wigd=((Cos(t) + -1._id)**5*62.5625_id + 126._id + (Cos(t) + -1._id)**4*402.1875_id + (Cos(t) + -1._id)*630._id + (Cos&
&(t) + -1._id)**3*990._id + (Cos(t) + -1._id)**2*1155._id)*Sin(t*0.5_id)**4
endif

if(l==7.and.m1==2.and.m2==-1) then
   wigd=Cos(t*0.5_id)*1.22474487139158904909864203735294569598297374032834_id*(56._id + (Cos(t) + -1._id)**5*62.5625_id &
&+ (Cos(t) + -1._id)*350._id + (Cos(t) + -1._id)**4*357.5_id + (Cos(t) + -1._id)**2*770._id + (Cos(t) + -1._id)**3*770._id)*Sin(t*0.5_id)**3
endif

if(l==7.and.m1==2.and.m2==0) then
   wigd=Cos(t*0.5_id)**2*1.3093073414159542875965849124937167111384161647908_id*(21._id + (Cos(t) + -1._id)**5*62.5625_i&
&d + (Cos(t) + -1._id)*175._id + (Cos(t) + -1._id)**4*312.8125_id + (Cos(t) + -1._id)**2*481.25_id + (Cos(t) + -1._id)**3*577.5_id)*Sin(t*0.5_id)**2
endif

if(l==7.and.m1==2.and.m2==1) then
   wigd=Cos(t*0.5_id)**3*1.22474487139158904909864203735294569598297374032834_id*(6._id + (Cos(t) + -1._id)**5*62.5625_i&
&d + (Cos(t) + -1._id)*75._id + (Cos(t) + -1._id)**4*268.125_id + (Cos(t) + -1._id)**2*275._id + (Cos(t) + -1._id)**3*412.5_id)*Sin(t*0.5_id)
endif

if(l==7.and.m1==2.and.m2==2) then
   wigd=Cos(t*0.5_id)**4*(1._id + (Cos(t) + -1._id)*25._id + (Cos(t) + -1._id)**5*62.5625_id + (Cos(t) + -1._id)**2*137.&
&5_id + (Cos(t) + -1._id)**4*223.4375_id + (Cos(t) + -1._id)**3*275._id)
endif

if(l==7.and.m1==2.and.m2==3) then
   wigd=Cos(t*0.5_id)**5*-1.4142135623730950488016887242096980785696718753769_id*(5._id + (Cos(t) + -1._id)*55._id + (Co&
&s(t) + -1._id)**4*62.5625_id + (Cos(t) + -1._id)**2*165._id + (Cos(t) + -1._id)**3*178.75_id)*Sin(t*0.5_id)
endif

if(l==7.and.m1==2.and.m2==4) then
   wigd=Cos(t*0.5_id)**6*2.3452078799117147772828150567722331402941141767059_id*(10._id + (Cos(t) + -1._id)**3*45.5_id +&
& (Cos(t) + -1._id)*60._id + (Cos(t) + -1._id)**2*97.5_id)*Sin(t*0.5_id)**2
endif

if(l==7.and.m1==2.and.m2==5) then
   wigd=Cos(t*0.5_id)**7*-4.6904157598234295545656301135444662805882283534117_id*(10._id + (Cos(t) + -1._id)**2*22.75_id&
& + (Cos(t) + -1._id)*32.5_id)*Sin(t*0.5_id)**3
endif

if(l==7.and.m1==2.and.m2==6) then
   wigd=Cos(t*0.5_id)**8*5.9791303715506990105649203780978083069955431425686_id*(-4._id + Cos(t)*14._id)*Sin(t*0.5_id)**&
&4
endif

if(l==7.and.m1==2.and.m2==7) then
   wigd=Cos(t*0.5_id)**9*-44.743714642394187341373897459828840314950565642978_id*Sin(t*0.5_id)**5
endif

if(l==7.and.m1==3.and.m2==-7) then
   wigd=Cos(t*0.5_id)**4*31.63858403911274914310629158480098308708005351898_id*Sin(t*0.5_id)**10
endif

if(l==7.and.m1==3.and.m2==-6) then
   wigd=Cos(t*0.5_id)**3*4.2278836313219407245431294806941488436657699679026_id*(6._id + Cos(t)*14._id)*Sin(t*0.5_id)**9&
&
endif

if(l==7.and.m1==3.and.m2==-5) then
   wigd=Cos(t*0.5_id)**2*3.3166247903553998491149327366706866839270885455894_id*((Cos(t) + -1._id)**2*22.75_id + 45._id &
&+ (Cos(t) + -1._id)*65._id)*Sin(t*0.5_id)**8
endif

if(l==7.and.m1==3.and.m2==-4) then
   wigd=Cos(t*0.5_id)*1.6583123951776999245574663683353433419635442727947_id*((Cos(t) + -1._id)**3*45.5_id + 120._id + (&
&Cos(t) + -1._id)**2*195._id + (Cos(t) + -1._id)*270._id)*Sin(t*0.5_id)**7
endif

if(l==7.and.m1==3.and.m2==-3) then
   wigd=((Cos(t) + -1._id)**4*62.5625_id + 210._id + (Cos(t) + -1._id)**3*357.5_id + (Cos(t) + -1._id)*660._id + (Cos(t)&
& + -1._id)**2*742.5_id)*Sin(t*0.5_id)**6
endif

if(l==7.and.m1==3.and.m2==-2) then
   wigd=Cos(t*0.5_id)*1.41421356237309504880168872420969807856967187537695_id*((Cos(t) + -1._id)**4*62.5625_id + 126._id&
& + (Cos(t) + -1._id)**3*321.75_id + (Cos(t) + -1._id)*462._id + (Cos(t) + -1._id)**2*594._id)*Sin(t*0.5_id)**5
endif

if(l==7.and.m1==3.and.m2==-1) then
   wigd=Cos(t*0.5_id)**2*1.73205080756887729352744634150587236694280525381038_id*((Cos(t) + -1._id)**4*62.5625_id + 70._&
&id + (Cos(t) + -1._id)**3*286._id + (Cos(t) + -1._id)*308._id + (Cos(t) + -1._id)**2*462._id)*Sin(t*0.5_id)**4
endif

if(l==7.and.m1==3.and.m2==0) then
   wigd=Cos(t*0.5_id)**3*1.8516401995451029231331335531679990450586298020167_id*(35._id + (Cos(t) + -1._id)**4*62.5625_i&
&d + (Cos(t) + -1._id)*192.5_id + (Cos(t) + -1._id)**3*250.25_id + (Cos(t) + -1._id)**2*346.5_id)*Sin(t*0.5_id)**3
endif

if(l==7.and.m1==3.and.m2==1) then
   wigd=Cos(t*0.5_id)**4*1.73205080756887729352744634150587236694280525381038_id*(15._id + (Cos(t) + -1._id)**4*62.5625_&
&id + (Cos(t) + -1._id)*110._id + (Cos(t) + -1._id)**3*214.5_id + (Cos(t) + -1._id)**2*247.5_id)*Sin(t*0.5_id)**2
endif

if(l==7.and.m1==3.and.m2==2) then
   wigd=Cos(t*0.5_id)**5*1.41421356237309504880168872420969807856967187537695_id*(5._id + (Cos(t) + -1._id)*55._id + (Co&
&s(t) + -1._id)**4*62.5625_id + (Cos(t) + -1._id)**2*165._id + (Cos(t) + -1._id)**3*178.75_id)*Sin(t*0.5_id)
endif

if(l==7.and.m1==3.and.m2==3) then
   wigd=Cos(t*0.5_id)**6*(1._id + (Cos(t) + -1._id)*22._id + (Cos(t) + -1._id)**4*62.5625_id + (Cos(t) + -1._id)**2*99._&
&id + (Cos(t) + -1._id)**3*143._id)
endif

if(l==7.and.m1==3.and.m2==4) then
   wigd=Cos(t*0.5_id)**7*-1.6583123951776999245574663683353433419635442727947_id*(4._id + (Cos(t) + -1._id)*36._id + (Co&
&s(t) + -1._id)**3*45.5_id + (Cos(t) + -1._id)**2*78._id)*Sin(t*0.5_id)
endif

if(l==7.and.m1==3.and.m2==5) then
   wigd=Cos(t*0.5_id)**8*3.3166247903553998491149327366706866839270885455894_id*(6._id + (Cos(t) + -1._id)**2*22.75_id +&
& (Cos(t) + -1._id)*26._id)*Sin(t*0.5_id)**2
endif

if(l==7.and.m1==3.and.m2==6) then
   wigd=Cos(t*0.5_id)**9*-4.2278836313219407245431294806941488436657699679026_id*(-6._id + Cos(t)*14._id)*Sin(t*0.5_id)*&
&*3
endif

if(l==7.and.m1==3.and.m2==7) then
   wigd=Cos(t*0.5_id)**10*31.63858403911274914310629158480098308708005351898_id*Sin(t*0.5_id)**4
endif

if(l==7.and.m1==4.and.m2==-7) then
   wigd=Cos(t*0.5_id)**3*19.078784028338912983052431720464530805092468505011_id*Sin(t*0.5_id)**11
endif

if(l==7.and.m1==4.and.m2==-6) then
   wigd=Cos(t*0.5_id)**2*2.5495097567963924150141120545113909947818854730498_id*(8._id + Cos(t)*14._id)*Sin(t*0.5_id)**1&
&0
endif

if(l==7.and.m1==4.and.m2==-5) then
   wigd=Cos(t*0.5_id)*2._id*((Cos(t) + -1._id)**2*22.75_id + 55._id + (Cos(t) + -1._id)*71.5_id)*Sin(t*0.5_id)**9
endif

if(l==7.and.m1==4.and.m2==-4) then
   wigd=((Cos(t) + -1._id)**3*45.5_id + 165._id + (Cos(t) + -1._id)**2*214.5_id + (Cos(t) + -1._id)*330._id)*Sin(t*0.5_i&
&d)**8
endif

if(l==7.and.m1==4.and.m2==-3) then
   wigd=Cos(t*0.5_id)*1.6583123951776999245574663683353433419635442727947_id*((Cos(t) + -1._id)**3*45.5_id + 120._id + (&
&Cos(t) + -1._id)**2*195._id + (Cos(t) + -1._id)*270._id)*Sin(t*0.5_id)**7
endif

if(l==7.and.m1==4.and.m2==-2) then
   wigd=Cos(t*0.5_id)**2*2.3452078799117147772828150567722331402941141767059_id*((Cos(t) + -1._id)**3*45.5_id + 84._id +&
& (Cos(t) + -1._id)**2*175.5_id + (Cos(t) + -1._id)*216._id)*Sin(t*0.5_id)**6
endif

if(l==7.and.m1==4.and.m2==-1) then
   wigd=Cos(t*0.5_id)**3*2.8722813232690143299253057341094646591101322289914_id*((Cos(t) + -1._id)**3*45.5_id + 56._id +&
& (Cos(t) + -1._id)**2*156._id + (Cos(t) + -1._id)*168._id)*Sin(t*0.5_id)**5
endif

if(l==7.and.m1==4.and.m2==0) then
   wigd=Cos(t*0.5_id)**4*3.0705978943149538627491443679259589994234026382635_id*(35._id + (Cos(t) + -1._id)**3*45.5_id +&
& (Cos(t) + -1._id)*126._id + (Cos(t) + -1._id)**2*136.5_id)*Sin(t*0.5_id)**4
endif

if(l==7.and.m1==4.and.m2==1) then
   wigd=Cos(t*0.5_id)**5*2.8722813232690143299253057341094646591101322289914_id*(20._id + (Cos(t) + -1._id)**3*45.5_id +&
& (Cos(t) + -1._id)*90._id + (Cos(t) + -1._id)**2*117._id)*Sin(t*0.5_id)**3
endif

if(l==7.and.m1==4.and.m2==2) then
   wigd=Cos(t*0.5_id)**6*2.3452078799117147772828150567722331402941141767059_id*(10._id + (Cos(t) + -1._id)**3*45.5_id +&
& (Cos(t) + -1._id)*60._id + (Cos(t) + -1._id)**2*97.5_id)*Sin(t*0.5_id)**2
endif

if(l==7.and.m1==4.and.m2==3) then
   wigd=Cos(t*0.5_id)**7*1.6583123951776999245574663683353433419635442727947_id*(4._id + (Cos(t) + -1._id)*36._id + (Cos&
&(t) + -1._id)**3*45.5_id + (Cos(t) + -1._id)**2*78._id)*Sin(t*0.5_id)
endif

if(l==7.and.m1==4.and.m2==4) then
   wigd=Cos(t*0.5_id)**8*(1._id + (Cos(t) + -1._id)*18._id + (Cos(t) + -1._id)**3*45.5_id + (Cos(t) + -1._id)**2*58.5_id&
&)
endif

if(l==7.and.m1==4.and.m2==5) then
   wigd=Cos(t*0.5_id)**9*-2._id*(3._id + (Cos(t) + -1._id)*19.5_id + (Cos(t) + -1._id)**2*22.75_id)*Sin(t*0.5_id)
endif

if(l==7.and.m1==4.and.m2==6) then
   wigd=Cos(t*0.5_id)**10*2.5495097567963924150141120545113909947818854730498_id*(-8._id + Cos(t)*14._id)*Sin(t*0.5_id)*&
&*2
endif

if(l==7.and.m1==4.and.m2==7) then
   wigd=Cos(t*0.5_id)**11*-19.078784028338912983052431720464530805092468505011_id*Sin(t*0.5_id)**3
endif

if(l==7.and.m1==5.and.m2==-7) then
   wigd=Cos(t*0.5_id)**2*9.5393920141694564915262158602322654025462342525055_id*Sin(t*0.5_id)**12
endif

if(l==7.and.m1==5.and.m2==-6) then
   wigd=Cos(t*0.5_id)*1.2747548783981962075070560272556954973909427365249_id*(10._id + Cos(t)*14._id)*Sin(t*0.5_id)**11
endif

if(l==7.and.m1==5.and.m2==-5) then
   wigd=((Cos(t) + -1._id)**2*22.75_id + 66._id + (Cos(t) + -1._id)*78._id)*Sin(t*0.5_id)**10
endif

if(l==7.and.m1==5.and.m2==-4) then
   wigd=Cos(t*0.5_id)*2._id*((Cos(t) + -1._id)**2*22.75_id + 55._id + (Cos(t) + -1._id)*71.5_id)*Sin(t*0.5_id)**9
endif

if(l==7.and.m1==5.and.m2==-3) then
   wigd=Cos(t*0.5_id)**2*3.3166247903553998491149327366706866839270885455894_id*((Cos(t) + -1._id)**2*22.75_id + 45._id &
&+ (Cos(t) + -1._id)*65._id)*Sin(t*0.5_id)**8
endif

if(l==7.and.m1==5.and.m2==-2) then
   wigd=Cos(t*0.5_id)**3*4.6904157598234295545656301135444662805882283534117_id*((Cos(t) + -1._id)**2*22.75_id + 36._id &
&+ (Cos(t) + -1._id)*58.5_id)*Sin(t*0.5_id)**7
endif

if(l==7.and.m1==5.and.m2==-1) then
   wigd=Cos(t*0.5_id)**4*5.7445626465380286598506114682189293182202644579828_id*((Cos(t) + -1._id)**2*22.75_id + 28._id &
&+ (Cos(t) + -1._id)*52._id)*Sin(t*0.5_id)**6
endif

if(l==7.and.m1==5.and.m2==0) then
   wigd=Cos(t*0.5_id)**5*6.1411957886299077254982887358519179988468052765271_id*(21._id + (Cos(t) + -1._id)**2*22.75_id &
&+ (Cos(t) + -1._id)*45.5_id)*Sin(t*0.5_id)**5
endif

if(l==7.and.m1==5.and.m2==1) then
   wigd=Cos(t*0.5_id)**6*5.7445626465380286598506114682189293182202644579828_id*(15._id + (Cos(t) + -1._id)**2*22.75_id &
&+ (Cos(t) + -1._id)*39._id)*Sin(t*0.5_id)**4
endif

if(l==7.and.m1==5.and.m2==2) then
   wigd=Cos(t*0.5_id)**7*4.6904157598234295545656301135444662805882283534117_id*(10._id + (Cos(t) + -1._id)**2*22.75_id &
&+ (Cos(t) + -1._id)*32.5_id)*Sin(t*0.5_id)**3
endif

if(l==7.and.m1==5.and.m2==3) then
   wigd=Cos(t*0.5_id)**8*3.3166247903553998491149327366706866839270885455894_id*(6._id + (Cos(t) + -1._id)**2*22.75_id +&
& (Cos(t) + -1._id)*26._id)*Sin(t*0.5_id)**2
endif

if(l==7.and.m1==5.and.m2==4) then
   wigd=Cos(t*0.5_id)**9*2._id*(3._id + (Cos(t) + -1._id)*19.5_id + (Cos(t) + -1._id)**2*22.75_id)*Sin(t*0.5_id)
endif

if(l==7.and.m1==5.and.m2==5) then
   wigd=Cos(t*0.5_id)**10*(1._id + (Cos(t) + -1._id)*13._id + (Cos(t) + -1._id)**2*22.75_id)
endif

if(l==7.and.m1==5.and.m2==6) then
   wigd=Cos(t*0.5_id)**11*-1.2747548783981962075070560272556954973909427365249_id*(-10._id + Cos(t)*14._id)*Sin(t*0.5_id&
&)
endif

if(l==7.and.m1==5.and.m2==7) then
   wigd=Cos(t*0.5_id)**12*9.5393920141694564915262158602322654025462342525055_id*Sin(t*0.5_id)**2
endif

if(l==7.and.m1==6.and.m2==-7) then
   wigd=Cos(t*0.5_id)*3.7416573867739413855837487323165493017560198077787_id*Sin(t*0.5_id)**13
endif

if(l==7.and.m1==6.and.m2==-6) then
   wigd=0.5_id*(12._id + Cos(t)*14._id)*Sin(t*0.5_id)**12
endif

if(l==7.and.m1==6.and.m2==-5) then
   wigd=Cos(t*0.5_id)*1.2747548783981962075070560272556954973909427365249_id*(10._id + Cos(t)*14._id)*Sin(t*0.5_id)**11
endif

if(l==7.and.m1==6.and.m2==-4) then
   wigd=Cos(t*0.5_id)**2*2.5495097567963924150141120545113909947818854730498_id*(8._id + Cos(t)*14._id)*Sin(t*0.5_id)**1&
&0
endif

if(l==7.and.m1==6.and.m2==-3) then
   wigd=Cos(t*0.5_id)**3*4.2278836313219407245431294806941488436657699679026_id*(6._id + Cos(t)*14._id)*Sin(t*0.5_id)**9&
&
endif

if(l==7.and.m1==6.and.m2==-2) then
   wigd=Cos(t*0.5_id)**4*5.9791303715506990105649203780978083069955431425686_id*(4._id + Cos(t)*14._id)*Sin(t*0.5_id)**8&
&
endif

if(l==7.and.m1==6.and.m2==-1) then
   wigd=Cos(t*0.5_id)**5*7.322909257938404906286509527879665713489932169733_id*(2._id + Cos(t)*14._id)*Sin(t*0.5_id)**7
endif

if(l==7.and.m1==6.and.m2==0) then
   wigd=Cos(t)*Cos(t*0.5_id)**6*109.59927007056205755303386616481266854993549249787_id*Sin(t*0.5_id)**6
endif

if(l==7.and.m1==6.and.m2==1) then
   wigd=Cos(t*0.5_id)**7*7.322909257938404906286509527879665713489932169733_id*(-2._id + Cos(t)*14._id)*Sin(t*0.5_id)**5&
&
endif

if(l==7.and.m1==6.and.m2==2) then
   wigd=Cos(t*0.5_id)**8*5.9791303715506990105649203780978083069955431425686_id*(-4._id + Cos(t)*14._id)*Sin(t*0.5_id)**&
&4
endif

if(l==7.and.m1==6.and.m2==3) then
   wigd=Cos(t*0.5_id)**9*4.2278836313219407245431294806941488436657699679026_id*(-6._id + Cos(t)*14._id)*Sin(t*0.5_id)**&
&3
endif

if(l==7.and.m1==6.and.m2==4) then
   wigd=Cos(t*0.5_id)**10*2.5495097567963924150141120545113909947818854730498_id*(-8._id + Cos(t)*14._id)*Sin(t*0.5_id)*&
&*2
endif

if(l==7.and.m1==6.and.m2==5) then
   wigd=Cos(t*0.5_id)**11*1.2747548783981962075070560272556954973909427365249_id*(-10._id + Cos(t)*14._id)*Sin(t*0.5_id)&
&
endif

if(l==7.and.m1==6.and.m2==6) then
   wigd=Cos(t*0.5_id)**12*0.5_id*(-12._id + Cos(t)*14._id)
endif

if(l==7.and.m1==6.and.m2==7) then
   wigd=Cos(t*0.5_id)**13*-3.7416573867739413855837487323165493017560198077787_id*Sin(t*0.5_id)
endif

if(l==7.and.m1==7.and.m2==-7) then
   wigd=Sin(t*0.5_id)**14
endif

if(l==7.and.m1==7.and.m2==-6) then
   wigd=Cos(t*0.5_id)*3.7416573867739413855837487323165493017560198077787_id*Sin(t*0.5_id)**13
endif

if(l==7.and.m1==7.and.m2==-5) then
   wigd=Cos(t*0.5_id)**2*9.5393920141694564915262158602322654025462342525055_id*Sin(t*0.5_id)**12
endif

if(l==7.and.m1==7.and.m2==-4) then
   wigd=Cos(t*0.5_id)**3*19.078784028338912983052431720464530805092468505011_id*Sin(t*0.5_id)**11
endif

if(l==7.and.m1==7.and.m2==-3) then
   wigd=Cos(t*0.5_id)**4*31.63858403911274914310629158480098308708005351898_id*Sin(t*0.5_id)**10
endif

if(l==7.and.m1==7.and.m2==-2) then
   wigd=Cos(t*0.5_id)**5*44.743714642394187341373897459828840314950565642978_id*Sin(t*0.5_id)**9
endif

if(l==7.and.m1==7.and.m2==-1) then
   wigd=Cos(t*0.5_id)**6*54.799635035281028776516933082406334274967746248937_id*Sin(t*0.5_id)**8
endif

if(l==7.and.m1==7.and.m2==0) then
   wigd=Cos(t*0.5_id)**7*58.583274063507239250292076223037325707919457357862_id*Sin(t*0.5_id)**7
endif

if(l==7.and.m1==7.and.m2==1) then
   wigd=Cos(t*0.5_id)**8*54.799635035281028776516933082406334274967746248937_id*Sin(t*0.5_id)**6
endif

if(l==7.and.m1==7.and.m2==2) then
   wigd=Cos(t*0.5_id)**9*44.743714642394187341373897459828840314950565642978_id*Sin(t*0.5_id)**5
endif

if(l==7.and.m1==7.and.m2==3) then
   wigd=Cos(t*0.5_id)**10*31.63858403911274914310629158480098308708005351898_id*Sin(t*0.5_id)**4
endif

if(l==7.and.m1==7.and.m2==4) then
   wigd=Cos(t*0.5_id)**11*19.078784028338912983052431720464530805092468505011_id*Sin(t*0.5_id)**3
endif

if(l==7.and.m1==7.and.m2==5) then
   wigd=Cos(t*0.5_id)**12*9.5393920141694564915262158602322654025462342525055_id*Sin(t*0.5_id)**2
endif

if(l==7.and.m1==7.and.m2==6) then
   wigd=Cos(t*0.5_id)**13*3.7416573867739413855837487323165493017560198077787_id*Sin(t*0.5_id)
endif

if(l==7.and.m1==7.and.m2==7) then
   wigd=Cos(t*0.5_id)**14
endif

if(l==8.and.m1==-8.and.m2==-8) then
   wigd=Cos(t*0.5_id)**16
endif

if(l==8.and.m1==-8.and.m2==-7) then
   wigd=Cos(t*0.5_id)**15*-4._id*Sin(t*0.5_id)
endif

if(l==8.and.m1==-8.and.m2==-6) then
   wigd=Cos(t*0.5_id)**14*10.95445115010332226913939565601604267905489389996_id*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==-8.and.m2==-5) then
   wigd=Cos(t*0.5_id)**13*-23.664319132398464170269313166246468193662004923177_id*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==-8.and.m2==-4) then
   wigd=Cos(t*0.5_id)**12*42.661458015403083501772783042687091232665683613366_id*Sin(t*0.5_id)**4
endif

if(l==8.and.m1==-8.and.m2==-3) then
   wigd=Cos(t*0.5_id)**11*-66.09084656743322424726970525046488165135351934318_id*Sin(t*0.5_id)**5
endif

if(l==8.and.m1==-8.and.m2==-2) then
   wigd=Cos(t*0.5_id)**10*89.48742928478837468274779491965768062990113128596_id*Sin(t*0.5_id)**6
endif

if(l==8.and.m1==-8.and.m2==-1) then
   wigd=Cos(t*0.5_id)**9*-106.95793565696750114142397403237585587699286318683_id*Sin(t*0.5_id)**7
endif

if(l==8.and.m1==-8.and.m2==0) then
   wigd=Cos(t*0.5_id)**8*113.44602240713422209833732731669398460285660783012_id*Sin(t*0.5_id)**8
endif

if(l==8.and.m1==-8.and.m2==1) then
   wigd=Cos(t*0.5_id)**7*-106.95793565696750114142397403237585587699286318683_id*Sin(t*0.5_id)**9
endif

if(l==8.and.m1==-8.and.m2==2) then
   wigd=Cos(t*0.5_id)**6*89.48742928478837468274779491965768062990113128596_id*Sin(t*0.5_id)**10
endif

if(l==8.and.m1==-8.and.m2==3) then
   wigd=Cos(t*0.5_id)**5*-66.09084656743322424726970525046488165135351934318_id*Sin(t*0.5_id)**11
endif

if(l==8.and.m1==-8.and.m2==4) then
   wigd=Cos(t*0.5_id)**4*42.661458015403083501772783042687091232665683613366_id*Sin(t*0.5_id)**12
endif

if(l==8.and.m1==-8.and.m2==5) then
   wigd=Cos(t*0.5_id)**3*-23.664319132398464170269313166246468193662004923177_id*Sin(t*0.5_id)**13
endif

if(l==8.and.m1==-8.and.m2==6) then
   wigd=Cos(t*0.5_id)**2*10.95445115010332226913939565601604267905489389996_id*Sin(t*0.5_id)**14
endif

if(l==8.and.m1==-8.and.m2==7) then
   wigd=Cos(t*0.5_id)*-4._id*Sin(t*0.5_id)**15
endif

if(l==8.and.m1==-8.and.m2==8) then
   wigd=Sin(t*0.5_id)**16
endif

if(l==8.and.m1==-7.and.m2==-8) then
   wigd=Cos(t*0.5_id)**15*4._id*Sin(t*0.5_id)
endif

if(l==8.and.m1==-7.and.m2==-7) then
   wigd=Cos(t*0.5_id)**14*0.5_id*(-14._id + Cos(t)*16._id)
endif

if(l==8.and.m1==-7.and.m2==-6) then
   wigd=Cos(t*0.5_id)**13*-1.369306393762915283642424457002005334881861737495_id*(-12._id + Cos(t)*16._id)*Sin(t*0.5_id)&
&
endif

if(l==8.and.m1==-7.and.m2==-5) then
   wigd=Cos(t*0.5_id)**12*2.9580398915498080212836641457808085242077506153972_id*(-10._id + Cos(t)*16._id)*Sin(t*0.5_id)&
&**2
endif

if(l==8.and.m1==-7.and.m2==-4) then
   wigd=Cos(t*0.5_id)**11*-5.3326822519253854377215978803358864040832104516708_id*(-8._id + Cos(t)*16._id)*Sin(t*0.5_id)&
&**3
endif

if(l==8.and.m1==-7.and.m2==-3) then
   wigd=Cos(t*0.5_id)**10*8.261355820929153030908713156308110206419189917898_id*(-6._id + Cos(t)*16._id)*Sin(t*0.5_id)**&
&4
endif

if(l==8.and.m1==-7.and.m2==-2) then
   wigd=Cos(t*0.5_id)**9*-11.185928660598546835343474364957210078737641410745_id*(-4._id + Cos(t)*16._id)*Sin(t*0.5_id)*&
&*5
endif

if(l==8.and.m1==-7.and.m2==-1) then
   wigd=Cos(t*0.5_id)**8*13.369741957120937642677996754046981984624107898354_id*(-2._id + Cos(t)*16._id)*Sin(t*0.5_id)**&
&6
endif

if(l==8.and.m1==-7.and.m2==0) then
   wigd=Cos(t)*Cos(t*0.5_id)**7*-226.89204481426844419667465463338796920571321566023_id*Sin(t*0.5_id)**7
endif

if(l==8.and.m1==-7.and.m2==1) then
   wigd=Cos(t*0.5_id)**6*13.369741957120937642677996754046981984624107898354_id*(2._id + Cos(t)*16._id)*Sin(t*0.5_id)**8&
&
endif

if(l==8.and.m1==-7.and.m2==2) then
   wigd=Cos(t*0.5_id)**5*-11.185928660598546835343474364957210078737641410745_id*(4._id + Cos(t)*16._id)*Sin(t*0.5_id)**&
&9
endif

if(l==8.and.m1==-7.and.m2==3) then
   wigd=Cos(t*0.5_id)**4*8.261355820929153030908713156308110206419189917898_id*(6._id + Cos(t)*16._id)*Sin(t*0.5_id)**10&
&
endif

if(l==8.and.m1==-7.and.m2==4) then
   wigd=Cos(t*0.5_id)**3*-5.3326822519253854377215978803358864040832104516708_id*(8._id + Cos(t)*16._id)*Sin(t*0.5_id)**&
&11
endif

if(l==8.and.m1==-7.and.m2==5) then
   wigd=Cos(t*0.5_id)**2*2.9580398915498080212836641457808085242077506153972_id*(10._id + Cos(t)*16._id)*Sin(t*0.5_id)**&
&12
endif

if(l==8.and.m1==-7.and.m2==6) then
   wigd=Cos(t*0.5_id)*-1.369306393762915283642424457002005334881861737495_id*(12._id + Cos(t)*16._id)*Sin(t*0.5_id)**13
endif

if(l==8.and.m1==-7.and.m2==7) then
   wigd=0.5_id*(14._id + Cos(t)*16._id)*Sin(t*0.5_id)**14
endif

if(l==8.and.m1==-7.and.m2==8) then
   wigd=Cos(t*0.5_id)*-4._id*Sin(t*0.5_id)**15
endif

if(l==8.and.m1==-6.and.m2==-8) then
   wigd=Cos(t*0.5_id)**14*10.95445115010332226913939565601604267905489389996_id*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==-6.and.m2==-7) then
   wigd=Cos(t*0.5_id)**13*1.369306393762915283642424457002005334881861737495_id*(-12._id + Cos(t)*16._id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==-6.and.m2==-6) then
   wigd=Cos(t*0.5_id)**12*(1._id + (Cos(t) + -1._id)*15._id + (Cos(t) + -1._id)**2*30._id)
endif

if(l==8.and.m1==-6.and.m2==-5) then
   wigd=Cos(t*0.5_id)**11*-2.1602468994692867436553224786959988859017347690194_id*(3._id + (Cos(t) + -1._id)*22.5_id + (&
&Cos(t) + -1._id)**2*30._id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==-6.and.m2==-4) then
   wigd=Cos(t*0.5_id)**10*3.8944404818493075368873949642026809176529333861249_id*(6._id + (Cos(t) + -1._id)*30._id + (Co&
&s(t) + -1._id)**2*30._id)*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==-6.and.m2==-3) then
   wigd=Cos(t*0.5_id)**9*-6.0332412515993424345033528849187730050984461505421_id*(10._id + (Cos(t) + -1._id)**2*30._id +&
& (Cos(t) + -1._id)*37.5_id)*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==-6.and.m2==-2) then
   wigd=Cos(t*0.5_id)**8*8.1690472720711644400267750682082849910494805874321_id*(15._id + (Cos(t) + -1._id)**2*30._id + &
&(Cos(t) + -1._id)*45._id)*Sin(t*0.5_id)**4
endif

if(l==8.and.m1==-6.and.m2==-1) then
   wigd=Cos(t*0.5_id)**7*-9.763879010584539875048679370506220951319909559644_id*(21._id + (Cos(t) + -1._id)**2*30._id + &
&(Cos(t) + -1._id)*52.5_id)*Sin(t*0.5_id)**5
endif

if(l==8.and.m1==-6.and.m2==0) then
   wigd=Cos(t*0.5_id)**6*10.356157588603989566078588172024067830946506245733_id*(28._id + (Cos(t) + -1._id)**2*30._id + &
&(Cos(t) + -1._id)*60._id)*Sin(t*0.5_id)**6
endif

if(l==8.and.m1==-6.and.m2==1) then
   wigd=Cos(t*0.5_id)**5*-9.763879010584539875048679370506220951319909559644_id*((Cos(t) + -1._id)**2*30._id + 36._id + &
&(Cos(t) + -1._id)*67.5_id)*Sin(t*0.5_id)**7
endif

if(l==8.and.m1==-6.and.m2==2) then
   wigd=Cos(t*0.5_id)**4*8.1690472720711644400267750682082849910494805874321_id*((Cos(t) + -1._id)**2*30._id + 45._id + &
&(Cos(t) + -1._id)*75._id)*Sin(t*0.5_id)**8
endif

if(l==8.and.m1==-6.and.m2==3) then
   wigd=Cos(t*0.5_id)**3*-6.0332412515993424345033528849187730050984461505421_id*((Cos(t) + -1._id)**2*30._id + 55._id +&
& (Cos(t) + -1._id)*82.5_id)*Sin(t*0.5_id)**9
endif

if(l==8.and.m1==-6.and.m2==4) then
   wigd=Cos(t*0.5_id)**2*3.8944404818493075368873949642026809176529333861249_id*((Cos(t) + -1._id)**2*30._id + 66._id + &
&(Cos(t) + -1._id)*90._id)*Sin(t*0.5_id)**10
endif

if(l==8.and.m1==-6.and.m2==5) then
   wigd=Cos(t*0.5_id)*-2.1602468994692867436553224786959988859017347690194_id*((Cos(t) + -1._id)**2*30._id + 78._id + (C&
&os(t) + -1._id)*97.5_id)*Sin(t*0.5_id)**11
endif

if(l==8.and.m1==-6.and.m2==6) then
   wigd=((Cos(t) + -1._id)**2*30._id + 91._id + (Cos(t) + -1._id)*105._id)*Sin(t*0.5_id)**12
endif

if(l==8.and.m1==-6.and.m2==7) then
   wigd=Cos(t*0.5_id)*-1.369306393762915283642424457002005334881861737495_id*(12._id + Cos(t)*16._id)*Sin(t*0.5_id)**13
endif

if(l==8.and.m1==-6.and.m2==8) then
   wigd=Cos(t*0.5_id)**2*10.95445115010332226913939565601604267905489389996_id*Sin(t*0.5_id)**14
endif

if(l==8.and.m1==-5.and.m2==-8) then
   wigd=Cos(t*0.5_id)**13*23.664319132398464170269313166246468193662004923177_id*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==-5.and.m2==-7) then
   wigd=Cos(t*0.5_id)**12*2.9580398915498080212836641457808085242077506153972_id*(-10._id + Cos(t)*16._id)*Sin(t*0.5_id)&
&**2
endif

if(l==8.and.m1==-5.and.m2==-6) then
   wigd=Cos(t*0.5_id)**11*2.1602468994692867436553224786959988859017347690194_id*(3._id + (Cos(t) + -1._id)*22.5_id + (C&
&os(t) + -1._id)**2*30._id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==-5.and.m2==-5) then
   wigd=Cos(t*0.5_id)**10*(1._id + (Cos(t) + -1._id)*21._id + (Cos(t) + -1._id)**3*70._id + (Cos(t) + -1._id)**2*78.75_i&
&d)
endif

if(l==8.and.m1==-5.and.m2==-4) then
   wigd=Cos(t*0.5_id)**9*-1.8027756377319946465596106337352479731256482869226_id*(4._id + (Cos(t) + -1._id)*42._id + (Co&
&s(t) + -1._id)**3*70._id + (Cos(t) + -1._id)**2*105._id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==-5.and.m2==-3) then
   wigd=Cos(t*0.5_id)**8*2.792848008753788233976784908217275204128035944785_id*(10._id + (Cos(t) + -1._id)*70._id + (Cos&
&(t) + -1._id)**3*70._id + (Cos(t) + -1._id)**2*131.25_id)*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==-5.and.m2==-2) then
   wigd=Cos(t*0.5_id)**7*-3.7815340802378074032779109105564661534285535943372_id*(20._id + (Cos(t) + -1._id)**3*70._id +&
& (Cos(t) + -1._id)*105._id + (Cos(t) + -1._id)**2*157.5_id)*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==-5.and.m2==-1) then
   wigd=Cos(t*0.5_id)**6*4.51979771987324987758661308354299758386857907414_id*(35._id + (Cos(t) + -1._id)**3*70._id + (C&
&os(t) + -1._id)*147._id + (Cos(t) + -1._id)**2*183.75_id)*Sin(t*0.5_id)**4
endif

if(l==8.and.m1==-5.and.m2==0) then
   wigd=Cos(t*0.5_id)**5*-4.7939694259708057865757747278388043194589891760334_id*(56._id + (Cos(t) + -1._id)**3*70._id +&
& (Cos(t) + -1._id)*196._id + (Cos(t) + -1._id)**2*210._id)*Sin(t*0.5_id)**5
endif

if(l==8.and.m1==-5.and.m2==1) then
   wigd=Cos(t*0.5_id)**4*4.51979771987324987758661308354299758386857907414_id*((Cos(t) + -1._id)**3*70._id + 84._id + (C&
&os(t) + -1._id)**2*236.25_id + (Cos(t) + -1._id)*252._id)*Sin(t*0.5_id)**6
endif

if(l==8.and.m1==-5.and.m2==2) then
   wigd=Cos(t*0.5_id)**3*-3.7815340802378074032779109105564661534285535943372_id*((Cos(t) + -1._id)**3*70._id + 120._id &
&+ (Cos(t) + -1._id)**2*262.5_id + (Cos(t) + -1._id)*315._id)*Sin(t*0.5_id)**7
endif

if(l==8.and.m1==-5.and.m2==3) then
   wigd=Cos(t*0.5_id)**2*2.792848008753788233976784908217275204128035944785_id*((Cos(t) + -1._id)**3*70._id + 165._id + &
&(Cos(t) + -1._id)**2*288.75_id + (Cos(t) + -1._id)*385._id)*Sin(t*0.5_id)**8
endif

if(l==8.and.m1==-5.and.m2==4) then
   wigd=Cos(t*0.5_id)*-1.8027756377319946465596106337352479731256482869226_id*((Cos(t) + -1._id)**3*70._id + 220._id + (&
&Cos(t) + -1._id)**2*315._id + (Cos(t) + -1._id)*462._id)*Sin(t*0.5_id)**9
endif

if(l==8.and.m1==-5.and.m2==5) then
   wigd=((Cos(t) + -1._id)**3*70._id + 286._id + (Cos(t) + -1._id)**2*341.25_id + (Cos(t) + -1._id)*546._id)*Sin(t*0.5_i&
&d)**10
endif

if(l==8.and.m1==-5.and.m2==6) then
   wigd=Cos(t*0.5_id)*-2.1602468994692867436553224786959988859017347690194_id*((Cos(t) + -1._id)**2*30._id + 78._id + (C&
&os(t) + -1._id)*97.5_id)*Sin(t*0.5_id)**11
endif

if(l==8.and.m1==-5.and.m2==7) then
   wigd=Cos(t*0.5_id)**2*2.9580398915498080212836641457808085242077506153972_id*(10._id + Cos(t)*16._id)*Sin(t*0.5_id)**&
&12
endif

if(l==8.and.m1==-5.and.m2==8) then
   wigd=Cos(t*0.5_id)**3*-23.664319132398464170269313166246468193662004923177_id*Sin(t*0.5_id)**13
endif

if(l==8.and.m1==-4.and.m2==-8) then
   wigd=Cos(t*0.5_id)**12*42.661458015403083501772783042687091232665683613366_id*Sin(t*0.5_id)**4
endif

if(l==8.and.m1==-4.and.m2==-7) then
   wigd=Cos(t*0.5_id)**11*5.3326822519253854377215978803358864040832104516708_id*(-8._id + Cos(t)*16._id)*Sin(t*0.5_id)*&
&*3
endif

if(l==8.and.m1==-4.and.m2==-6) then
   wigd=Cos(t*0.5_id)**10*3.8944404818493075368873949642026809176529333861249_id*(6._id + (Cos(t) + -1._id)*30._id + (Co&
&s(t) + -1._id)**2*30._id)*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==-4.and.m2==-5) then
   wigd=Cos(t*0.5_id)**9*1.8027756377319946465596106337352479731256482869226_id*(4._id + (Cos(t) + -1._id)*42._id + (Cos&
&(t) + -1._id)**3*70._id + (Cos(t) + -1._id)**2*105._id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==-4.and.m2==-4) then
   wigd=Cos(t*0.5_id)**8*(1._id + (Cos(t) + -1._id)*26._id + (Cos(t) + -1._id)**4*113.75_id + (Cos(t) + -1._id)**2*136.5&
&_id + (Cos(t) + -1._id)**3*227.5_id)
endif

if(l==8.and.m1==-4.and.m2==-3) then
   wigd=Cos(t*0.5_id)**7*-1.5491933384829667540717061599129598443331686821166_id*(5._id + (Cos(t) + -1._id)*65._id + (Co&
&s(t) + -1._id)**4*113.75_id + (Cos(t) + -1._id)**2*227.5_id + (Cos(t) + -1._id)**3*284.375_id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==-4.and.m2==-2) then
   wigd=Cos(t*0.5_id)**6*2.097617696340303093982907027359875196950543715363_id*(15._id + (Cos(t) + -1._id)**4*113.75_id &
&+ (Cos(t) + -1._id)*130._id + (Cos(t) + -1._id)**2*341.25_id + (Cos(t) + -1._id)**3*341.25_id)*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==-4.and.m2==-1) then
   wigd=Cos(t*0.5_id)**5*-2.5071326821120348744018252306903741602502505221584_id*(35._id + (Cos(t) + -1._id)**4*113.75_i&
&d + (Cos(t) + -1._id)*227.5_id + (Cos(t) + -1._id)**3*398.125_id + (Cos(t) + -1._id)**2*477.75_id)*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==-4.and.m2==0) then
   wigd=Cos(t*0.5_id)**4*2.6592157812837549848856947047391733642548422694418_id*(70._id + (Cos(t) + -1._id)**4*113.75_id&
& + (Cos(t) + -1._id)*364._id + (Cos(t) + -1._id)**3*455._id + (Cos(t) + -1._id)**2*637._id)*Sin(t*0.5_id)**4
endif

if(l==8.and.m1==-4.and.m2==1) then
   wigd=Cos(t*0.5_id)**3*-2.5071326821120348744018252306903741602502505221584_id*((Cos(t) + -1._id)**4*113.75_id + 126._&
&id + (Cos(t) + -1._id)**3*511.875_id + (Cos(t) + -1._id)*546._id + (Cos(t) + -1._id)**2*819._id)*Sin(t*0.5_id)**5
endif

if(l==8.and.m1==-4.and.m2==2) then
   wigd=Cos(t*0.5_id)**2*2.097617696340303093982907027359875196950543715363_id*((Cos(t) + -1._id)**4*113.75_id + 210._id&
& + (Cos(t) + -1._id)**3*568.75_id + (Cos(t) + -1._id)*780._id + (Cos(t) + -1._id)**2*1023.75_id)*Sin(t*0.5_id)**6
endif

if(l==8.and.m1==-4.and.m2==3) then
   wigd=Cos(t*0.5_id)*-1.5491933384829667540717061599129598443331686821166_id*((Cos(t) + -1._id)**4*113.75_id + 330._id &
&+ (Cos(t) + -1._id)**3*625.625_id + (Cos(t) + -1._id)*1072.5_id + (Cos(t) + -1._id)**2*1251.25_id)*Sin(t*0.5_id)**7
endif

if(l==8.and.m1==-4.and.m2==4) then
   wigd=((Cos(t) + -1._id)**4*113.75_id + 495._id + (Cos(t) + -1._id)**3*682.5_id + (Cos(t) + -1._id)*1430._id + (Cos(t)&
& + -1._id)**2*1501.5_id)*Sin(t*0.5_id)**8
endif

if(l==8.and.m1==-4.and.m2==5) then
   wigd=Cos(t*0.5_id)*-1.8027756377319946465596106337352479731256482869226_id*((Cos(t) + -1._id)**3*70._id + 220._id + (&
&Cos(t) + -1._id)**2*315._id + (Cos(t) + -1._id)*462._id)*Sin(t*0.5_id)**9
endif

if(l==8.and.m1==-4.and.m2==6) then
   wigd=Cos(t*0.5_id)**2*3.8944404818493075368873949642026809176529333861249_id*((Cos(t) + -1._id)**2*30._id + 66._id + &
&(Cos(t) + -1._id)*90._id)*Sin(t*0.5_id)**10
endif

if(l==8.and.m1==-4.and.m2==7) then
   wigd=Cos(t*0.5_id)**3*-5.3326822519253854377215978803358864040832104516708_id*(8._id + Cos(t)*16._id)*Sin(t*0.5_id)**&
&11
endif

if(l==8.and.m1==-4.and.m2==8) then
   wigd=Cos(t*0.5_id)**4*42.661458015403083501772783042687091232665683613366_id*Sin(t*0.5_id)**12
endif

if(l==8.and.m1==-3.and.m2==-8) then
   wigd=Cos(t*0.5_id)**11*66.09084656743322424726970525046488165135351934318_id*Sin(t*0.5_id)**5
endif

if(l==8.and.m1==-3.and.m2==-7) then
   wigd=Cos(t*0.5_id)**10*8.261355820929153030908713156308110206419189917898_id*(-6._id + Cos(t)*16._id)*Sin(t*0.5_id)**&
&4
endif

if(l==8.and.m1==-3.and.m2==-6) then
   wigd=Cos(t*0.5_id)**9*6.0332412515993424345033528849187730050984461505421_id*(10._id + (Cos(t) + -1._id)**2*30._id + &
&(Cos(t) + -1._id)*37.5_id)*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==-3.and.m2==-5) then
   wigd=Cos(t*0.5_id)**8*2.792848008753788233976784908217275204128035944785_id*(10._id + (Cos(t) + -1._id)*70._id + (Cos&
&(t) + -1._id)**3*70._id + (Cos(t) + -1._id)**2*131.25_id)*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==-3.and.m2==-4) then
   wigd=Cos(t*0.5_id)**7*1.5491933384829667540717061599129598443331686821166_id*(5._id + (Cos(t) + -1._id)*65._id + (Cos&
&(t) + -1._id)**4*113.75_id + (Cos(t) + -1._id)**2*227.5_id + (Cos(t) + -1._id)**3*284.375_id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==-3.and.m2==-3) then
   wigd=Cos(t*0.5_id)**6*(1._id + (Cos(t) + -1._id)*30._id + (Cos(t) + -1._id)**5*136.5_id + (Cos(t) + -1._id)**2*195._i&
&d + (Cos(t) + -1._id)**4*426.5625_id + (Cos(t) + -1._id)**3*455._id)
endif

if(l==8.and.m1==-3.and.m2==-2) then
   wigd=Cos(t*0.5_id)**5*-1.3540064007726600600766472613776733914173673477876_id*(6._id + (Cos(t) + -1._id)*90._id + (Co&
&s(t) + -1._id)**5*136.5_id + (Cos(t) + -1._id)**2*390._id + (Cos(t) + -1._id)**4*511.875_id + (Cos(t) + -1._id)**3*682.5_id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==-3.and.m2==-1) then
   wigd=Cos(t*0.5_id)**4*1.61834718742537413773071441128209229367911210829573_id*(21._id + (Cos(t) + -1._id)**5*136.5_id&
& + (Cos(t) + -1._id)*210._id + (Cos(t) + -1._id)**4*597.1875_id + (Cos(t) + -1._id)**2*682.5_id + (Cos(t) + -1._id)**3*955.5_id)*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==-3.and.m2==0) then
   wigd=Cos(t*0.5_id)**3*-1.716516405813987968530033407550486175377456109921_id*(56._id + (Cos(t) + -1._id)**5*136.5_id &
&+ (Cos(t) + -1._id)*420._id + (Cos(t) + -1._id)**4*682.5_id + (Cos(t) + -1._id)**2*1092._id + (Cos(t) + -1._id)**3*1274._id)*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==-3.and.m2==1) then
   wigd=Cos(t*0.5_id)**2*1.61834718742537413773071441128209229367911210829573_id*(126._id + (Cos(t) + -1._id)**5*136.5_i&
&d + (Cos(t) + -1._id)*756._id + (Cos(t) + -1._id)**4*767.8125_id + (Cos(t) + -1._id)**2*1638._id + (Cos(t) + -1._id)**3*1638._id)*Sin(t*0.5_id)**4
endif

if(l==8.and.m1==-3.and.m2==2) then
   wigd=Cos(t*0.5_id)*-1.3540064007726600600766472613776733914173673477876_id*((Cos(t) + -1._id)**5*136.5_id + 252._id +&
& (Cos(t) + -1._id)**4*853.125_id + (Cos(t) + -1._id)*1260._id + (Cos(t) + -1._id)**3*2047.5_id + (Cos(t) + -1._id)**2*2340._id)*Sin(t*0.5_id)**5
endif

if(l==8.and.m1==-3.and.m2==3) then
   wigd=((Cos(t) + -1._id)**5*136.5_id + 462._id + (Cos(t) + -1._id)**4*938.4375_id + (Cos(t) + -1._id)*1980._id + (Cos(&
&t) + -1._id)**3*2502.5_id + (Cos(t) + -1._id)**2*3217.5_id)*Sin(t*0.5_id)**6
endif

if(l==8.and.m1==-3.and.m2==4) then
   wigd=Cos(t*0.5_id)*-1.5491933384829667540717061599129598443331686821166_id*((Cos(t) + -1._id)**4*113.75_id + 330._id &
&+ (Cos(t) + -1._id)**3*625.625_id + (Cos(t) + -1._id)*1072.5_id + (Cos(t) + -1._id)**2*1251.25_id)*Sin(t*0.5_id)**7
endif

if(l==8.and.m1==-3.and.m2==5) then
   wigd=Cos(t*0.5_id)**2*2.792848008753788233976784908217275204128035944785_id*((Cos(t) + -1._id)**3*70._id + 165._id + &
&(Cos(t) + -1._id)**2*288.75_id + (Cos(t) + -1._id)*385._id)*Sin(t*0.5_id)**8
endif

if(l==8.and.m1==-3.and.m2==6) then
   wigd=Cos(t*0.5_id)**3*-6.0332412515993424345033528849187730050984461505421_id*((Cos(t) + -1._id)**2*30._id + 55._id +&
& (Cos(t) + -1._id)*82.5_id)*Sin(t*0.5_id)**9
endif

if(l==8.and.m1==-3.and.m2==7) then
   wigd=Cos(t*0.5_id)**4*8.261355820929153030908713156308110206419189917898_id*(6._id + Cos(t)*16._id)*Sin(t*0.5_id)**10&
&
endif

if(l==8.and.m1==-3.and.m2==8) then
   wigd=Cos(t*0.5_id)**5*-66.09084656743322424726970525046488165135351934318_id*Sin(t*0.5_id)**11
endif

if(l==8.and.m1==-2.and.m2==-8) then
   wigd=Cos(t*0.5_id)**10*89.48742928478837468274779491965768062990113128596_id*Sin(t*0.5_id)**6
endif

if(l==8.and.m1==-2.and.m2==-7) then
   wigd=Cos(t*0.5_id)**9*11.185928660598546835343474364957210078737641410745_id*(-4._id + Cos(t)*16._id)*Sin(t*0.5_id)**&
&5
endif

if(l==8.and.m1==-2.and.m2==-6) then
   wigd=Cos(t*0.5_id)**8*8.1690472720711644400267750682082849910494805874321_id*(15._id + (Cos(t) + -1._id)**2*30._id + &
&(Cos(t) + -1._id)*45._id)*Sin(t*0.5_id)**4
endif

if(l==8.and.m1==-2.and.m2==-5) then
   wigd=Cos(t*0.5_id)**7*3.7815340802378074032779109105564661534285535943372_id*(20._id + (Cos(t) + -1._id)**3*70._id + &
&(Cos(t) + -1._id)*105._id + (Cos(t) + -1._id)**2*157.5_id)*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==-2.and.m2==-4) then
   wigd=Cos(t*0.5_id)**6*2.097617696340303093982907027359875196950543715363_id*(15._id + (Cos(t) + -1._id)**4*113.75_id &
&+ (Cos(t) + -1._id)*130._id + (Cos(t) + -1._id)**2*341.25_id + (Cos(t) + -1._id)**3*341.25_id)*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==-2.and.m2==-3) then
   wigd=Cos(t*0.5_id)**5*1.35400640077266006007664726137767339141736734778755_id*(6._id + (Cos(t) + -1._id)*90._id + (Co&
&s(t) + -1._id)**5*136.5_id + (Cos(t) + -1._id)**2*390._id + (Cos(t) + -1._id)**4*511.875_id + (Cos(t) + -1._id)**3*682.5_id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==-2.and.m2==-2) then
   wigd=Cos(t*0.5_id)**4*(1._id + (Cos(t) + -1._id)*33._id + (Cos(t) + -1._id)**6*125.125_id + (Cos(t) + -1._id)**2*247.&
&5_id + (Cos(t) + -1._id)**5*563.0625_id + (Cos(t) + -1._id)**3*715._id + (Cos(t) + -1._id)**4*938.4375_id)
endif

if(l==8.and.m1==-2.and.m2==-1) then
   wigd=Cos(t*0.5_id)**3*-1.1952286093343936399688171796931249848468790989981_id*(7._id + (Cos(t) + -1._id)*115.5_id + (&
&Cos(t) + -1._id)**6*125.125_id + (Cos(t) + -1._id)**2*577.5_id + (Cos(t) + -1._id)**5*656.90625_id + (Cos(t) + -1._id)**3*1251.25_id + (Cos(t) + -1._id)**4*1313.8125_id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==-2.and.m2==0) then
   wigd=Cos(t*0.5_id)**2*1.2677313820927748662644274910489179389461788351702_id*(28._id + (Cos(t) + -1._id)**6*125.125_i&
&d + (Cos(t) + -1._id)*308._id + (Cos(t) + -1._id)**5*750.75_id + (Cos(t) + -1._id)**2*1155._id + (Cos(t) + -1._id)**4*1751.75_id + (Cos(t) + -1._id)**3*2002._id)*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==-2.and.m2==1) then
   wigd=Cos(t*0.5_id)*-1.1952286093343936399688171796931249848468790989981_id*(84._id + (Cos(t) + -1._id)**6*125.125_id &
&+ (Cos(t) + -1._id)*693._id + (Cos(t) + -1._id)**5*844.59375_id + (Cos(t) + -1._id)**2*2079._id + (Cos(t) + -1._id)**4*2252.25_id + (Cos(t) + -1._id)**3*3003._id)*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==-2.and.m2==2) then
   wigd=((Cos(t) + -1._id)**6*125.125_id + 210._id + (Cos(t) + -1._id)**5*938.4375_id + (Cos(t) + -1._id)*1386._id + (Co&
&s(t) + -1._id)**4*2815.3125_id + (Cos(t) + -1._id)**2*3465._id + (Cos(t) + -1._id)**3*4290._id)*Sin(t*0.5_id)**4
endif

if(l==8.and.m1==-2.and.m2==3) then
   wigd=Cos(t*0.5_id)*-1.3540064007726600600766472613776733914173673477876_id*((Cos(t) + -1._id)**5*136.5_id + 252._id +&
& (Cos(t) + -1._id)**4*853.125_id + (Cos(t) + -1._id)*1260._id + (Cos(t) + -1._id)**3*2047.5_id + (Cos(t) + -1._id)**2*2340._id)*Sin(t*0.5_id)**5
endif

if(l==8.and.m1==-2.and.m2==4) then
   wigd=Cos(t*0.5_id)**2*2.097617696340303093982907027359875196950543715363_id*((Cos(t) + -1._id)**4*113.75_id + 210._id&
& + (Cos(t) + -1._id)**3*568.75_id + (Cos(t) + -1._id)*780._id + (Cos(t) + -1._id)**2*1023.75_id)*Sin(t*0.5_id)**6
endif

if(l==8.and.m1==-2.and.m2==5) then
   wigd=Cos(t*0.5_id)**3*-3.7815340802378074032779109105564661534285535943372_id*((Cos(t) + -1._id)**3*70._id + 120._id &
&+ (Cos(t) + -1._id)**2*262.5_id + (Cos(t) + -1._id)*315._id)*Sin(t*0.5_id)**7
endif

if(l==8.and.m1==-2.and.m2==6) then
   wigd=Cos(t*0.5_id)**4*8.1690472720711644400267750682082849910494805874321_id*((Cos(t) + -1._id)**2*30._id + 45._id + &
&(Cos(t) + -1._id)*75._id)*Sin(t*0.5_id)**8
endif

if(l==8.and.m1==-2.and.m2==7) then
   wigd=Cos(t*0.5_id)**5*-11.185928660598546835343474364957210078737641410745_id*(4._id + Cos(t)*16._id)*Sin(t*0.5_id)**&
&9
endif

if(l==8.and.m1==-2.and.m2==8) then
   wigd=Cos(t*0.5_id)**6*89.48742928478837468274779491965768062990113128596_id*Sin(t*0.5_id)**10
endif

if(l==8.and.m1==-1.and.m2==-8) then
   wigd=Cos(t*0.5_id)**9*106.95793565696750114142397403237585587699286318683_id*Sin(t*0.5_id)**7
endif

if(l==8.and.m1==-1.and.m2==-7) then
   wigd=Cos(t*0.5_id)**8*13.369741957120937642677996754046981984624107898354_id*(-2._id + Cos(t)*16._id)*Sin(t*0.5_id)**&
&6
endif

if(l==8.and.m1==-1.and.m2==-6) then
   wigd=Cos(t*0.5_id)**7*9.7638790105845398750486793705062209513199095596437_id*(21._id + (Cos(t) + -1._id)**2*30._id + &
&(Cos(t) + -1._id)*52.5_id)*Sin(t*0.5_id)**5
endif

if(l==8.and.m1==-1.and.m2==-5) then
   wigd=Cos(t*0.5_id)**6*4.51979771987324987758661308354299758386857907414_id*(35._id + (Cos(t) + -1._id)**3*70._id + (C&
&os(t) + -1._id)*147._id + (Cos(t) + -1._id)**2*183.75_id)*Sin(t*0.5_id)**4
endif

if(l==8.and.m1==-1.and.m2==-4) then
   wigd=Cos(t*0.5_id)**5*2.5071326821120348744018252306903741602502505221584_id*(35._id + (Cos(t) + -1._id)**4*113.75_id&
& + (Cos(t) + -1._id)*227.5_id + (Cos(t) + -1._id)**3*398.125_id + (Cos(t) + -1._id)**2*477.75_id)*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==-1.and.m2==-3) then
   wigd=Cos(t*0.5_id)**4*1.61834718742537413773071441128209229367911210829573_id*(21._id + (Cos(t) + -1._id)**5*136.5_id&
& + (Cos(t) + -1._id)*210._id + (Cos(t) + -1._id)**4*597.1875_id + (Cos(t) + -1._id)**2*682.5_id + (Cos(t) + -1._id)**3*955.5_id)*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==-1.and.m2==-2) then
   wigd=Cos(t*0.5_id)**3*1.1952286093343936399688171796931249848468790989981_id*(7._id + (Cos(t) + -1._id)*115.5_id + (C&
&os(t) + -1._id)**6*125.125_id + (Cos(t) + -1._id)**2*577.5_id + (Cos(t) + -1._id)**5*656.90625_id + (Cos(t) + -1._id)**3*1251.25_id + (Cos(t) + -1._id)**4*1313.8125_id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==-1.and.m2==-1) then
   wigd=Cos(t*0.5_id)**2*(1._id + (Cos(t) + -1._id)*35._id + (Cos(t) + -1._id)**7*89.375_id + (Cos(t) + -1._id)**2*288.7&
&5_id + (Cos(t) + -1._id)**6*547.421875_id + (Cos(t) + -1._id)**3*962.5_id + (Cos(t) + -1._id)**5*1313.8125_id + (Cos(t) + -1._id)**4*1564.0625_id)
endif

if(l==8.and.m1==-1.and.m2==0) then
   wigd=Cos(t*0.5_id)*-1.0606601717798212866012665431572735589272539065327_id*(8._id + (Cos(t) + -1._id)**7*89.375_id + &
&(Cos(t) + -1._id)*140._id + (Cos(t) + -1._id)**6*625.625_id + (Cos(t) + -1._id)**2*770._id + (Cos(t) + -1._id)**5*1751.75_id + (Cos(t) + -1._id)**3*1925._id + (Cos(t) + -1._id)**4*2502.5_id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==-1.and.m2==1) then
   wigd=(36._id + (Cos(t) + -1._id)**7*89.375_id + (Cos(t) + -1._id)*420._id + (Cos(t) + -1._id)**6*703.828125_id + (Cos&
&(t) + -1._id)**2*1732.5_id + (Cos(t) + -1._id)**5*2252.25_id + (Cos(t) + -1._id)**3*3465._id + (Cos(t) + -1._id)**4*3753.75_id)*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==-1.and.m2==2) then
   wigd=Cos(t*0.5_id)*-1.1952286093343936399688171796931249848468790989981_id*(84._id + (Cos(t) + -1._id)**6*125.125_id &
&+ (Cos(t) + -1._id)*693._id + (Cos(t) + -1._id)**5*844.59375_id + (Cos(t) + -1._id)**2*2079._id + (Cos(t) + -1._id)**4*2252.25_id + (Cos(t) + -1._id)**3*3003._id)*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==-1.and.m2==3) then
   wigd=Cos(t*0.5_id)**2*1.61834718742537413773071441128209229367911210829573_id*(126._id + (Cos(t) + -1._id)**5*136.5_i&
&d + (Cos(t) + -1._id)*756._id + (Cos(t) + -1._id)**4*767.8125_id + (Cos(t) + -1._id)**2*1638._id + (Cos(t) + -1._id)**3*1638._id)*Sin(t*0.5_id)**4
endif

if(l==8.and.m1==-1.and.m2==4) then
   wigd=Cos(t*0.5_id)**3*-2.5071326821120348744018252306903741602502505221584_id*((Cos(t) + -1._id)**4*113.75_id + 126._&
&id + (Cos(t) + -1._id)**3*511.875_id + (Cos(t) + -1._id)*546._id + (Cos(t) + -1._id)**2*819._id)*Sin(t*0.5_id)**5
endif

if(l==8.and.m1==-1.and.m2==5) then
   wigd=Cos(t*0.5_id)**4*4.51979771987324987758661308354299758386857907414_id*((Cos(t) + -1._id)**3*70._id + 84._id + (C&
&os(t) + -1._id)**2*236.25_id + (Cos(t) + -1._id)*252._id)*Sin(t*0.5_id)**6
endif

if(l==8.and.m1==-1.and.m2==6) then
   wigd=Cos(t*0.5_id)**5*-9.763879010584539875048679370506220951319909559644_id*((Cos(t) + -1._id)**2*30._id + 36._id + &
&(Cos(t) + -1._id)*67.5_id)*Sin(t*0.5_id)**7
endif

if(l==8.and.m1==-1.and.m2==7) then
   wigd=Cos(t*0.5_id)**6*13.369741957120937642677996754046981984624107898354_id*(2._id + Cos(t)*16._id)*Sin(t*0.5_id)**8&
&
endif

if(l==8.and.m1==-1.and.m2==8) then
   wigd=Cos(t*0.5_id)**7*-106.95793565696750114142397403237585587699286318683_id*Sin(t*0.5_id)**9
endif

if(l==8.and.m1==0.and.m2==-8) then
   wigd=Cos(t*0.5_id)**8*113.44602240713422209833732731669398460285660783012_id*Sin(t*0.5_id)**8
endif

if(l==8.and.m1==0.and.m2==-7) then
   wigd=Cos(t)*Cos(t*0.5_id)**7*226.89204481426844419667465463338796920571321566023_id*Sin(t*0.5_id)**7
endif

if(l==8.and.m1==0.and.m2==-6) then
   wigd=Cos(t*0.5_id)**6*10.356157588603989566078588172024067830946506245733_id*(28._id + (Cos(t) + -1._id)**2*30._id + &
&(Cos(t) + -1._id)*60._id)*Sin(t*0.5_id)**6
endif

if(l==8.and.m1==0.and.m2==-5) then
   wigd=Cos(t*0.5_id)**5*4.7939694259708057865757747278388043194589891760334_id*(56._id + (Cos(t) + -1._id)**3*70._id + &
&(Cos(t) + -1._id)*196._id + (Cos(t) + -1._id)**2*210._id)*Sin(t*0.5_id)**5
endif

if(l==8.and.m1==0.and.m2==-4) then
   wigd=Cos(t*0.5_id)**4*2.6592157812837549848856947047391733642548422694418_id*(70._id + (Cos(t) + -1._id)**4*113.75_id&
& + (Cos(t) + -1._id)*364._id + (Cos(t) + -1._id)**3*455._id + (Cos(t) + -1._id)**2*637._id)*Sin(t*0.5_id)**4
endif

if(l==8.and.m1==0.and.m2==-3) then
   wigd=Cos(t*0.5_id)**3*1.716516405813987968530033407550486175377456109921_id*(56._id + (Cos(t) + -1._id)**5*136.5_id +&
& (Cos(t) + -1._id)*420._id + (Cos(t) + -1._id)**4*682.5_id + (Cos(t) + -1._id)**2*1092._id + (Cos(t) + -1._id)**3*1274._id)*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==0.and.m2==-2) then
   wigd=Cos(t*0.5_id)**2*1.2677313820927748662644274910489179389461788351702_id*(28._id + (Cos(t) + -1._id)**6*125.125_i&
&d + (Cos(t) + -1._id)*308._id + (Cos(t) + -1._id)**5*750.75_id + (Cos(t) + -1._id)**2*1155._id + (Cos(t) + -1._id)**4*1751.75_id + (Cos(t) + -1._id)**3*2002._id)*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==0.and.m2==-1) then
   wigd=Cos(t*0.5_id)*1.0606601717798212866012665431572735589272539065327_id*(8._id + (Cos(t) + -1._id)**7*89.375_id + (&
&Cos(t) + -1._id)*140._id + (Cos(t) + -1._id)**6*625.625_id + (Cos(t) + -1._id)**2*770._id + (Cos(t) + -1._id)**5*1751.75_id + (Cos(t) + -1._id)**3*1925._id + (Cos(t) + -1._id)**4*2502.5_id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==0.and.m2==0) then
   wigd=0.0078125_id*(Cos(t)**6*-12012._id + Cos(t)**2*-1260._id + 35._id + Cos(t)**8*6435._id + Cos(t)**4*6930._id)
endif

if(l==8.and.m1==0.and.m2==1) then
   wigd=Cos(t*0.5_id)*-1.0606601717798212866012665431572735589272539065327_id*(8._id + (Cos(t) + -1._id)**7*89.375_id + &
&(Cos(t) + -1._id)*140._id + (Cos(t) + -1._id)**6*625.625_id + (Cos(t) + -1._id)**2*770._id + (Cos(t) + -1._id)**5*1751.75_id + (Cos(t) + -1._id)**3*1925._id + (Cos(t) + -1._id)**4*2502.5_id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==0.and.m2==2) then
   wigd=Cos(t*0.5_id)**2*1.2677313820927748662644274910489179389461788351702_id*(28._id + (Cos(t) + -1._id)**6*125.125_i&
&d + (Cos(t) + -1._id)*308._id + (Cos(t) + -1._id)**5*750.75_id + (Cos(t) + -1._id)**2*1155._id + (Cos(t) + -1._id)**4*1751.75_id + (Cos(t) + -1._id)**3*2002._id)*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==0.and.m2==3) then
   wigd=Cos(t*0.5_id)**3*-1.716516405813987968530033407550486175377456109921_id*(56._id + (Cos(t) + -1._id)**5*136.5_id &
&+ (Cos(t) + -1._id)*420._id + (Cos(t) + -1._id)**4*682.5_id + (Cos(t) + -1._id)**2*1092._id + (Cos(t) + -1._id)**3*1274._id)*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==0.and.m2==4) then
   wigd=Cos(t*0.5_id)**4*2.6592157812837549848856947047391733642548422694418_id*(70._id + (Cos(t) + -1._id)**4*113.75_id&
& + (Cos(t) + -1._id)*364._id + (Cos(t) + -1._id)**3*455._id + (Cos(t) + -1._id)**2*637._id)*Sin(t*0.5_id)**4
endif

if(l==8.and.m1==0.and.m2==5) then
   wigd=Cos(t*0.5_id)**5*-4.7939694259708057865757747278388043194589891760334_id*(56._id + (Cos(t) + -1._id)**3*70._id +&
& (Cos(t) + -1._id)*196._id + (Cos(t) + -1._id)**2*210._id)*Sin(t*0.5_id)**5
endif

if(l==8.and.m1==0.and.m2==6) then
   wigd=Cos(t*0.5_id)**6*10.356157588603989566078588172024067830946506245733_id*(28._id + (Cos(t) + -1._id)**2*30._id + &
&(Cos(t) + -1._id)*60._id)*Sin(t*0.5_id)**6
endif

if(l==8.and.m1==0.and.m2==7) then
   wigd=Cos(t)*Cos(t*0.5_id)**7*-226.89204481426844419667465463338796920571321566023_id*Sin(t*0.5_id)**7
endif

if(l==8.and.m1==0.and.m2==8) then
   wigd=Cos(t*0.5_id)**8*113.44602240713422209833732731669398460285660783012_id*Sin(t*0.5_id)**8
endif

if(l==8.and.m1==1.and.m2==-8) then
   wigd=Cos(t*0.5_id)**7*106.95793565696750114142397403237585587699286318683_id*Sin(t*0.5_id)**9
endif

if(l==8.and.m1==1.and.m2==-7) then
   wigd=Cos(t*0.5_id)**6*13.369741957120937642677996754046981984624107898354_id*(2._id + Cos(t)*16._id)*Sin(t*0.5_id)**8&
&
endif

if(l==8.and.m1==1.and.m2==-6) then
   wigd=Cos(t*0.5_id)**5*9.7638790105845398750486793705062209513199095596437_id*((Cos(t) + -1._id)**2*30._id + 36._id + &
&(Cos(t) + -1._id)*67.5_id)*Sin(t*0.5_id)**7
endif

if(l==8.and.m1==1.and.m2==-5) then
   wigd=Cos(t*0.5_id)**4*4.51979771987324987758661308354299758386857907414_id*((Cos(t) + -1._id)**3*70._id + 84._id + (C&
&os(t) + -1._id)**2*236.25_id + (Cos(t) + -1._id)*252._id)*Sin(t*0.5_id)**6
endif

if(l==8.and.m1==1.and.m2==-4) then
   wigd=Cos(t*0.5_id)**3*2.5071326821120348744018252306903741602502505221584_id*((Cos(t) + -1._id)**4*113.75_id + 126._i&
&d + (Cos(t) + -1._id)**3*511.875_id + (Cos(t) + -1._id)*546._id + (Cos(t) + -1._id)**2*819._id)*Sin(t*0.5_id)**5
endif

if(l==8.and.m1==1.and.m2==-3) then
   wigd=Cos(t*0.5_id)**2*1.61834718742537413773071441128209229367911210829573_id*(126._id + (Cos(t) + -1._id)**5*136.5_i&
&d + (Cos(t) + -1._id)*756._id + (Cos(t) + -1._id)**4*767.8125_id + (Cos(t) + -1._id)**2*1638._id + (Cos(t) + -1._id)**3*1638._id)*Sin(t*0.5_id)**4
endif

if(l==8.and.m1==1.and.m2==-2) then
   wigd=Cos(t*0.5_id)*1.1952286093343936399688171796931249848468790989981_id*(84._id + (Cos(t) + -1._id)**6*125.125_id +&
& (Cos(t) + -1._id)*693._id + (Cos(t) + -1._id)**5*844.59375_id + (Cos(t) + -1._id)**2*2079._id + (Cos(t) + -1._id)**4*2252.25_id + (Cos(t) + -1._id)**3*3003._id)*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==1.and.m2==-1) then
   wigd=(36._id + (Cos(t) + -1._id)**7*89.375_id + (Cos(t) + -1._id)*420._id + (Cos(t) + -1._id)**6*703.828125_id + (Cos&
&(t) + -1._id)**2*1732.5_id + (Cos(t) + -1._id)**5*2252.25_id + (Cos(t) + -1._id)**3*3465._id + (Cos(t) + -1._id)**4*3753.75_id)*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==1.and.m2==0) then
   wigd=Cos(t*0.5_id)*1.0606601717798212866012665431572735589272539065327_id*(8._id + (Cos(t) + -1._id)**7*89.375_id + (&
&Cos(t) + -1._id)*140._id + (Cos(t) + -1._id)**6*625.625_id + (Cos(t) + -1._id)**2*770._id + (Cos(t) + -1._id)**5*1751.75_id + (Cos(t) + -1._id)**3*1925._id + (Cos(t) + -1._id)**4*2502.5_id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==1.and.m2==1) then
   wigd=Cos(t*0.5_id)**2*(1._id + (Cos(t) + -1._id)*35._id + (Cos(t) + -1._id)**7*89.375_id + (Cos(t) + -1._id)**2*288.7&
&5_id + (Cos(t) + -1._id)**6*547.421875_id + (Cos(t) + -1._id)**3*962.5_id + (Cos(t) + -1._id)**5*1313.8125_id + (Cos(t) + -1._id)**4*1564.0625_id)
endif

if(l==8.and.m1==1.and.m2==2) then
   wigd=Cos(t*0.5_id)**3*-1.1952286093343936399688171796931249848468790989981_id*(7._id + (Cos(t) + -1._id)*115.5_id + (&
&Cos(t) + -1._id)**6*125.125_id + (Cos(t) + -1._id)**2*577.5_id + (Cos(t) + -1._id)**5*656.90625_id + (Cos(t) + -1._id)**3*1251.25_id + (Cos(t) + -1._id)**4*1313.8125_id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==1.and.m2==3) then
   wigd=Cos(t*0.5_id)**4*1.61834718742537413773071441128209229367911210829573_id*(21._id + (Cos(t) + -1._id)**5*136.5_id&
& + (Cos(t) + -1._id)*210._id + (Cos(t) + -1._id)**4*597.1875_id + (Cos(t) + -1._id)**2*682.5_id + (Cos(t) + -1._id)**3*955.5_id)*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==1.and.m2==4) then
   wigd=Cos(t*0.5_id)**5*-2.5071326821120348744018252306903741602502505221584_id*(35._id + (Cos(t) + -1._id)**4*113.75_i&
&d + (Cos(t) + -1._id)*227.5_id + (Cos(t) + -1._id)**3*398.125_id + (Cos(t) + -1._id)**2*477.75_id)*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==1.and.m2==5) then
   wigd=Cos(t*0.5_id)**6*4.51979771987324987758661308354299758386857907414_id*(35._id + (Cos(t) + -1._id)**3*70._id + (C&
&os(t) + -1._id)*147._id + (Cos(t) + -1._id)**2*183.75_id)*Sin(t*0.5_id)**4
endif

if(l==8.and.m1==1.and.m2==6) then
   wigd=Cos(t*0.5_id)**7*-9.763879010584539875048679370506220951319909559644_id*(21._id + (Cos(t) + -1._id)**2*30._id + &
&(Cos(t) + -1._id)*52.5_id)*Sin(t*0.5_id)**5
endif

if(l==8.and.m1==1.and.m2==7) then
   wigd=Cos(t*0.5_id)**8*13.369741957120937642677996754046981984624107898354_id*(-2._id + Cos(t)*16._id)*Sin(t*0.5_id)**&
&6
endif

if(l==8.and.m1==1.and.m2==8) then
   wigd=Cos(t*0.5_id)**9*-106.95793565696750114142397403237585587699286318683_id*Sin(t*0.5_id)**7
endif

if(l==8.and.m1==2.and.m2==-8) then
   wigd=Cos(t*0.5_id)**6*89.48742928478837468274779491965768062990113128596_id*Sin(t*0.5_id)**10
endif

if(l==8.and.m1==2.and.m2==-7) then
   wigd=Cos(t*0.5_id)**5*11.185928660598546835343474364957210078737641410745_id*(4._id + Cos(t)*16._id)*Sin(t*0.5_id)**9&
&
endif

if(l==8.and.m1==2.and.m2==-6) then
   wigd=Cos(t*0.5_id)**4*8.1690472720711644400267750682082849910494805874321_id*((Cos(t) + -1._id)**2*30._id + 45._id + &
&(Cos(t) + -1._id)*75._id)*Sin(t*0.5_id)**8
endif

if(l==8.and.m1==2.and.m2==-5) then
   wigd=Cos(t*0.5_id)**3*3.7815340802378074032779109105564661534285535943372_id*((Cos(t) + -1._id)**3*70._id + 120._id +&
& (Cos(t) + -1._id)**2*262.5_id + (Cos(t) + -1._id)*315._id)*Sin(t*0.5_id)**7
endif

if(l==8.and.m1==2.and.m2==-4) then
   wigd=Cos(t*0.5_id)**2*2.097617696340303093982907027359875196950543715363_id*((Cos(t) + -1._id)**4*113.75_id + 210._id&
& + (Cos(t) + -1._id)**3*568.75_id + (Cos(t) + -1._id)*780._id + (Cos(t) + -1._id)**2*1023.75_id)*Sin(t*0.5_id)**6
endif

if(l==8.and.m1==2.and.m2==-3) then
   wigd=Cos(t*0.5_id)*1.35400640077266006007664726137767339141736734778755_id*((Cos(t) + -1._id)**5*136.5_id + 252._id +&
& (Cos(t) + -1._id)**4*853.125_id + (Cos(t) + -1._id)*1260._id + (Cos(t) + -1._id)**3*2047.5_id + (Cos(t) + -1._id)**2*2340._id)*Sin(t*0.5_id)**5
endif

if(l==8.and.m1==2.and.m2==-2) then
   wigd=((Cos(t) + -1._id)**6*125.125_id + 210._id + (Cos(t) + -1._id)**5*938.4375_id + (Cos(t) + -1._id)*1386._id + (Co&
&s(t) + -1._id)**4*2815.3125_id + (Cos(t) + -1._id)**2*3465._id + (Cos(t) + -1._id)**3*4290._id)*Sin(t*0.5_id)**4
endif

if(l==8.and.m1==2.and.m2==-1) then
   wigd=Cos(t*0.5_id)*1.1952286093343936399688171796931249848468790989981_id*(84._id + (Cos(t) + -1._id)**6*125.125_id +&
& (Cos(t) + -1._id)*693._id + (Cos(t) + -1._id)**5*844.59375_id + (Cos(t) + -1._id)**2*2079._id + (Cos(t) + -1._id)**4*2252.25_id + (Cos(t) + -1._id)**3*3003._id)*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==2.and.m2==0) then
   wigd=Cos(t*0.5_id)**2*1.2677313820927748662644274910489179389461788351702_id*(28._id + (Cos(t) + -1._id)**6*125.125_i&
&d + (Cos(t) + -1._id)*308._id + (Cos(t) + -1._id)**5*750.75_id + (Cos(t) + -1._id)**2*1155._id + (Cos(t) + -1._id)**4*1751.75_id + (Cos(t) + -1._id)**3*2002._id)*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==2.and.m2==1) then
   wigd=Cos(t*0.5_id)**3*1.1952286093343936399688171796931249848468790989981_id*(7._id + (Cos(t) + -1._id)*115.5_id + (C&
&os(t) + -1._id)**6*125.125_id + (Cos(t) + -1._id)**2*577.5_id + (Cos(t) + -1._id)**5*656.90625_id + (Cos(t) + -1._id)**3*1251.25_id + (Cos(t) + -1._id)**4*1313.8125_id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==2.and.m2==2) then
   wigd=Cos(t*0.5_id)**4*(1._id + (Cos(t) + -1._id)*33._id + (Cos(t) + -1._id)**6*125.125_id + (Cos(t) + -1._id)**2*247.&
&5_id + (Cos(t) + -1._id)**5*563.0625_id + (Cos(t) + -1._id)**3*715._id + (Cos(t) + -1._id)**4*938.4375_id)
endif

if(l==8.and.m1==2.and.m2==3) then
   wigd=Cos(t*0.5_id)**5*-1.3540064007726600600766472613776733914173673477876_id*(6._id + (Cos(t) + -1._id)*90._id + (Co&
&s(t) + -1._id)**5*136.5_id + (Cos(t) + -1._id)**2*390._id + (Cos(t) + -1._id)**4*511.875_id + (Cos(t) + -1._id)**3*682.5_id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==2.and.m2==4) then
   wigd=Cos(t*0.5_id)**6*2.097617696340303093982907027359875196950543715363_id*(15._id + (Cos(t) + -1._id)**4*113.75_id &
&+ (Cos(t) + -1._id)*130._id + (Cos(t) + -1._id)**2*341.25_id + (Cos(t) + -1._id)**3*341.25_id)*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==2.and.m2==5) then
   wigd=Cos(t*0.5_id)**7*-3.7815340802378074032779109105564661534285535943372_id*(20._id + (Cos(t) + -1._id)**3*70._id +&
& (Cos(t) + -1._id)*105._id + (Cos(t) + -1._id)**2*157.5_id)*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==2.and.m2==6) then
   wigd=Cos(t*0.5_id)**8*8.1690472720711644400267750682082849910494805874321_id*(15._id + (Cos(t) + -1._id)**2*30._id + &
&(Cos(t) + -1._id)*45._id)*Sin(t*0.5_id)**4
endif

if(l==8.and.m1==2.and.m2==7) then
   wigd=Cos(t*0.5_id)**9*-11.185928660598546835343474364957210078737641410745_id*(-4._id + Cos(t)*16._id)*Sin(t*0.5_id)*&
&*5
endif

if(l==8.and.m1==2.and.m2==8) then
   wigd=Cos(t*0.5_id)**10*89.48742928478837468274779491965768062990113128596_id*Sin(t*0.5_id)**6
endif

if(l==8.and.m1==3.and.m2==-8) then
   wigd=Cos(t*0.5_id)**5*66.09084656743322424726970525046488165135351934318_id*Sin(t*0.5_id)**11
endif

if(l==8.and.m1==3.and.m2==-7) then
   wigd=Cos(t*0.5_id)**4*8.261355820929153030908713156308110206419189917898_id*(6._id + Cos(t)*16._id)*Sin(t*0.5_id)**10&
&
endif

if(l==8.and.m1==3.and.m2==-6) then
   wigd=Cos(t*0.5_id)**3*6.0332412515993424345033528849187730050984461505421_id*((Cos(t) + -1._id)**2*30._id + 55._id + &
&(Cos(t) + -1._id)*82.5_id)*Sin(t*0.5_id)**9
endif

if(l==8.and.m1==3.and.m2==-5) then
   wigd=Cos(t*0.5_id)**2*2.792848008753788233976784908217275204128035944785_id*((Cos(t) + -1._id)**3*70._id + 165._id + &
&(Cos(t) + -1._id)**2*288.75_id + (Cos(t) + -1._id)*385._id)*Sin(t*0.5_id)**8
endif

if(l==8.and.m1==3.and.m2==-4) then
   wigd=Cos(t*0.5_id)*1.5491933384829667540717061599129598443331686821166_id*((Cos(t) + -1._id)**4*113.75_id + 330._id +&
& (Cos(t) + -1._id)**3*625.625_id + (Cos(t) + -1._id)*1072.5_id + (Cos(t) + -1._id)**2*1251.25_id)*Sin(t*0.5_id)**7
endif

if(l==8.and.m1==3.and.m2==-3) then
   wigd=((Cos(t) + -1._id)**5*136.5_id + 462._id + (Cos(t) + -1._id)**4*938.4375_id + (Cos(t) + -1._id)*1980._id + (Cos(&
&t) + -1._id)**3*2502.5_id + (Cos(t) + -1._id)**2*3217.5_id)*Sin(t*0.5_id)**6
endif

if(l==8.and.m1==3.and.m2==-2) then
   wigd=Cos(t*0.5_id)*1.35400640077266006007664726137767339141736734778755_id*((Cos(t) + -1._id)**5*136.5_id + 252._id +&
& (Cos(t) + -1._id)**4*853.125_id + (Cos(t) + -1._id)*1260._id + (Cos(t) + -1._id)**3*2047.5_id + (Cos(t) + -1._id)**2*2340._id)*Sin(t*0.5_id)**5
endif

if(l==8.and.m1==3.and.m2==-1) then
   wigd=Cos(t*0.5_id)**2*1.61834718742537413773071441128209229367911210829573_id*(126._id + (Cos(t) + -1._id)**5*136.5_i&
&d + (Cos(t) + -1._id)*756._id + (Cos(t) + -1._id)**4*767.8125_id + (Cos(t) + -1._id)**2*1638._id + (Cos(t) + -1._id)**3*1638._id)*Sin(t*0.5_id)**4
endif

if(l==8.and.m1==3.and.m2==0) then
   wigd=Cos(t*0.5_id)**3*1.716516405813987968530033407550486175377456109921_id*(56._id + (Cos(t) + -1._id)**5*136.5_id +&
& (Cos(t) + -1._id)*420._id + (Cos(t) + -1._id)**4*682.5_id + (Cos(t) + -1._id)**2*1092._id + (Cos(t) + -1._id)**3*1274._id)*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==3.and.m2==1) then
   wigd=Cos(t*0.5_id)**4*1.61834718742537413773071441128209229367911210829573_id*(21._id + (Cos(t) + -1._id)**5*136.5_id&
& + (Cos(t) + -1._id)*210._id + (Cos(t) + -1._id)**4*597.1875_id + (Cos(t) + -1._id)**2*682.5_id + (Cos(t) + -1._id)**3*955.5_id)*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==3.and.m2==2) then
   wigd=Cos(t*0.5_id)**5*1.35400640077266006007664726137767339141736734778755_id*(6._id + (Cos(t) + -1._id)*90._id + (Co&
&s(t) + -1._id)**5*136.5_id + (Cos(t) + -1._id)**2*390._id + (Cos(t) + -1._id)**4*511.875_id + (Cos(t) + -1._id)**3*682.5_id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==3.and.m2==3) then
   wigd=Cos(t*0.5_id)**6*(1._id + (Cos(t) + -1._id)*30._id + (Cos(t) + -1._id)**5*136.5_id + (Cos(t) + -1._id)**2*195._i&
&d + (Cos(t) + -1._id)**4*426.5625_id + (Cos(t) + -1._id)**3*455._id)
endif

if(l==8.and.m1==3.and.m2==4) then
   wigd=Cos(t*0.5_id)**7*-1.5491933384829667540717061599129598443331686821166_id*(5._id + (Cos(t) + -1._id)*65._id + (Co&
&s(t) + -1._id)**4*113.75_id + (Cos(t) + -1._id)**2*227.5_id + (Cos(t) + -1._id)**3*284.375_id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==3.and.m2==5) then
   wigd=Cos(t*0.5_id)**8*2.792848008753788233976784908217275204128035944785_id*(10._id + (Cos(t) + -1._id)*70._id + (Cos&
&(t) + -1._id)**3*70._id + (Cos(t) + -1._id)**2*131.25_id)*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==3.and.m2==6) then
   wigd=Cos(t*0.5_id)**9*-6.0332412515993424345033528849187730050984461505421_id*(10._id + (Cos(t) + -1._id)**2*30._id +&
& (Cos(t) + -1._id)*37.5_id)*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==3.and.m2==7) then
   wigd=Cos(t*0.5_id)**10*8.261355820929153030908713156308110206419189917898_id*(-6._id + Cos(t)*16._id)*Sin(t*0.5_id)**&
&4
endif

if(l==8.and.m1==3.and.m2==8) then
   wigd=Cos(t*0.5_id)**11*-66.09084656743322424726970525046488165135351934318_id*Sin(t*0.5_id)**5
endif

if(l==8.and.m1==4.and.m2==-8) then
   wigd=Cos(t*0.5_id)**4*42.661458015403083501772783042687091232665683613366_id*Sin(t*0.5_id)**12
endif

if(l==8.and.m1==4.and.m2==-7) then
   wigd=Cos(t*0.5_id)**3*5.3326822519253854377215978803358864040832104516708_id*(8._id + Cos(t)*16._id)*Sin(t*0.5_id)**1&
&1
endif

if(l==8.and.m1==4.and.m2==-6) then
   wigd=Cos(t*0.5_id)**2*3.8944404818493075368873949642026809176529333861249_id*((Cos(t) + -1._id)**2*30._id + 66._id + &
&(Cos(t) + -1._id)*90._id)*Sin(t*0.5_id)**10
endif

if(l==8.and.m1==4.and.m2==-5) then
   wigd=Cos(t*0.5_id)*1.8027756377319946465596106337352479731256482869226_id*((Cos(t) + -1._id)**3*70._id + 220._id + (C&
&os(t) + -1._id)**2*315._id + (Cos(t) + -1._id)*462._id)*Sin(t*0.5_id)**9
endif

if(l==8.and.m1==4.and.m2==-4) then
   wigd=((Cos(t) + -1._id)**4*113.75_id + 495._id + (Cos(t) + -1._id)**3*682.5_id + (Cos(t) + -1._id)*1430._id + (Cos(t)&
& + -1._id)**2*1501.5_id)*Sin(t*0.5_id)**8
endif

if(l==8.and.m1==4.and.m2==-3) then
   wigd=Cos(t*0.5_id)*1.5491933384829667540717061599129598443331686821166_id*((Cos(t) + -1._id)**4*113.75_id + 330._id +&
& (Cos(t) + -1._id)**3*625.625_id + (Cos(t) + -1._id)*1072.5_id + (Cos(t) + -1._id)**2*1251.25_id)*Sin(t*0.5_id)**7
endif

if(l==8.and.m1==4.and.m2==-2) then
   wigd=Cos(t*0.5_id)**2*2.097617696340303093982907027359875196950543715363_id*((Cos(t) + -1._id)**4*113.75_id + 210._id&
& + (Cos(t) + -1._id)**3*568.75_id + (Cos(t) + -1._id)*780._id + (Cos(t) + -1._id)**2*1023.75_id)*Sin(t*0.5_id)**6
endif

if(l==8.and.m1==4.and.m2==-1) then
   wigd=Cos(t*0.5_id)**3*2.5071326821120348744018252306903741602502505221584_id*((Cos(t) + -1._id)**4*113.75_id + 126._i&
&d + (Cos(t) + -1._id)**3*511.875_id + (Cos(t) + -1._id)*546._id + (Cos(t) + -1._id)**2*819._id)*Sin(t*0.5_id)**5
endif

if(l==8.and.m1==4.and.m2==0) then
   wigd=Cos(t*0.5_id)**4*2.6592157812837549848856947047391733642548422694418_id*(70._id + (Cos(t) + -1._id)**4*113.75_id&
& + (Cos(t) + -1._id)*364._id + (Cos(t) + -1._id)**3*455._id + (Cos(t) + -1._id)**2*637._id)*Sin(t*0.5_id)**4
endif

if(l==8.and.m1==4.and.m2==1) then
   wigd=Cos(t*0.5_id)**5*2.5071326821120348744018252306903741602502505221584_id*(35._id + (Cos(t) + -1._id)**4*113.75_id&
& + (Cos(t) + -1._id)*227.5_id + (Cos(t) + -1._id)**3*398.125_id + (Cos(t) + -1._id)**2*477.75_id)*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==4.and.m2==2) then
   wigd=Cos(t*0.5_id)**6*2.097617696340303093982907027359875196950543715363_id*(15._id + (Cos(t) + -1._id)**4*113.75_id &
&+ (Cos(t) + -1._id)*130._id + (Cos(t) + -1._id)**2*341.25_id + (Cos(t) + -1._id)**3*341.25_id)*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==4.and.m2==3) then
   wigd=Cos(t*0.5_id)**7*1.5491933384829667540717061599129598443331686821166_id*(5._id + (Cos(t) + -1._id)*65._id + (Cos&
&(t) + -1._id)**4*113.75_id + (Cos(t) + -1._id)**2*227.5_id + (Cos(t) + -1._id)**3*284.375_id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==4.and.m2==4) then
   wigd=Cos(t*0.5_id)**8*(1._id + (Cos(t) + -1._id)*26._id + (Cos(t) + -1._id)**4*113.75_id + (Cos(t) + -1._id)**2*136.5&
&_id + (Cos(t) + -1._id)**3*227.5_id)
endif

if(l==8.and.m1==4.and.m2==5) then
   wigd=Cos(t*0.5_id)**9*-1.8027756377319946465596106337352479731256482869226_id*(4._id + (Cos(t) + -1._id)*42._id + (Co&
&s(t) + -1._id)**3*70._id + (Cos(t) + -1._id)**2*105._id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==4.and.m2==6) then
   wigd=Cos(t*0.5_id)**10*3.8944404818493075368873949642026809176529333861249_id*(6._id + (Cos(t) + -1._id)*30._id + (Co&
&s(t) + -1._id)**2*30._id)*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==4.and.m2==7) then
   wigd=Cos(t*0.5_id)**11*-5.3326822519253854377215978803358864040832104516708_id*(-8._id + Cos(t)*16._id)*Sin(t*0.5_id)&
&**3
endif

if(l==8.and.m1==4.and.m2==8) then
   wigd=Cos(t*0.5_id)**12*42.661458015403083501772783042687091232665683613366_id*Sin(t*0.5_id)**4
endif

if(l==8.and.m1==5.and.m2==-8) then
   wigd=Cos(t*0.5_id)**3*23.664319132398464170269313166246468193662004923177_id*Sin(t*0.5_id)**13
endif

if(l==8.and.m1==5.and.m2==-7) then
   wigd=Cos(t*0.5_id)**2*2.9580398915498080212836641457808085242077506153972_id*(10._id + Cos(t)*16._id)*Sin(t*0.5_id)**&
&12
endif

if(l==8.and.m1==5.and.m2==-6) then
   wigd=Cos(t*0.5_id)*2.1602468994692867436553224786959988859017347690194_id*((Cos(t) + -1._id)**2*30._id + 78._id + (Co&
&s(t) + -1._id)*97.5_id)*Sin(t*0.5_id)**11
endif

if(l==8.and.m1==5.and.m2==-5) then
   wigd=((Cos(t) + -1._id)**3*70._id + 286._id + (Cos(t) + -1._id)**2*341.25_id + (Cos(t) + -1._id)*546._id)*Sin(t*0.5_i&
&d)**10
endif

if(l==8.and.m1==5.and.m2==-4) then
   wigd=Cos(t*0.5_id)*1.8027756377319946465596106337352479731256482869226_id*((Cos(t) + -1._id)**3*70._id + 220._id + (C&
&os(t) + -1._id)**2*315._id + (Cos(t) + -1._id)*462._id)*Sin(t*0.5_id)**9
endif

if(l==8.and.m1==5.and.m2==-3) then
   wigd=Cos(t*0.5_id)**2*2.792848008753788233976784908217275204128035944785_id*((Cos(t) + -1._id)**3*70._id + 165._id + &
&(Cos(t) + -1._id)**2*288.75_id + (Cos(t) + -1._id)*385._id)*Sin(t*0.5_id)**8
endif

if(l==8.and.m1==5.and.m2==-2) then
   wigd=Cos(t*0.5_id)**3*3.7815340802378074032779109105564661534285535943372_id*((Cos(t) + -1._id)**3*70._id + 120._id +&
& (Cos(t) + -1._id)**2*262.5_id + (Cos(t) + -1._id)*315._id)*Sin(t*0.5_id)**7
endif

if(l==8.and.m1==5.and.m2==-1) then
   wigd=Cos(t*0.5_id)**4*4.51979771987324987758661308354299758386857907414_id*((Cos(t) + -1._id)**3*70._id + 84._id + (C&
&os(t) + -1._id)**2*236.25_id + (Cos(t) + -1._id)*252._id)*Sin(t*0.5_id)**6
endif

if(l==8.and.m1==5.and.m2==0) then
   wigd=Cos(t*0.5_id)**5*4.7939694259708057865757747278388043194589891760334_id*(56._id + (Cos(t) + -1._id)**3*70._id + &
&(Cos(t) + -1._id)*196._id + (Cos(t) + -1._id)**2*210._id)*Sin(t*0.5_id)**5
endif

if(l==8.and.m1==5.and.m2==1) then
   wigd=Cos(t*0.5_id)**6*4.51979771987324987758661308354299758386857907414_id*(35._id + (Cos(t) + -1._id)**3*70._id + (C&
&os(t) + -1._id)*147._id + (Cos(t) + -1._id)**2*183.75_id)*Sin(t*0.5_id)**4
endif

if(l==8.and.m1==5.and.m2==2) then
   wigd=Cos(t*0.5_id)**7*3.7815340802378074032779109105564661534285535943372_id*(20._id + (Cos(t) + -1._id)**3*70._id + &
&(Cos(t) + -1._id)*105._id + (Cos(t) + -1._id)**2*157.5_id)*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==5.and.m2==3) then
   wigd=Cos(t*0.5_id)**8*2.792848008753788233976784908217275204128035944785_id*(10._id + (Cos(t) + -1._id)*70._id + (Cos&
&(t) + -1._id)**3*70._id + (Cos(t) + -1._id)**2*131.25_id)*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==5.and.m2==4) then
   wigd=Cos(t*0.5_id)**9*1.8027756377319946465596106337352479731256482869226_id*(4._id + (Cos(t) + -1._id)*42._id + (Cos&
&(t) + -1._id)**3*70._id + (Cos(t) + -1._id)**2*105._id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==5.and.m2==5) then
   wigd=Cos(t*0.5_id)**10*(1._id + (Cos(t) + -1._id)*21._id + (Cos(t) + -1._id)**3*70._id + (Cos(t) + -1._id)**2*78.75_i&
&d)
endif

if(l==8.and.m1==5.and.m2==6) then
   wigd=Cos(t*0.5_id)**11*-2.1602468994692867436553224786959988859017347690194_id*(3._id + (Cos(t) + -1._id)*22.5_id + (&
&Cos(t) + -1._id)**2*30._id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==5.and.m2==7) then
   wigd=Cos(t*0.5_id)**12*2.9580398915498080212836641457808085242077506153972_id*(-10._id + Cos(t)*16._id)*Sin(t*0.5_id)&
&**2
endif

if(l==8.and.m1==5.and.m2==8) then
   wigd=Cos(t*0.5_id)**13*-23.664319132398464170269313166246468193662004923177_id*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==6.and.m2==-8) then
   wigd=Cos(t*0.5_id)**2*10.95445115010332226913939565601604267905489389996_id*Sin(t*0.5_id)**14
endif

if(l==8.and.m1==6.and.m2==-7) then
   wigd=Cos(t*0.5_id)*1.369306393762915283642424457002005334881861737495_id*(12._id + Cos(t)*16._id)*Sin(t*0.5_id)**13
endif

if(l==8.and.m1==6.and.m2==-6) then
   wigd=((Cos(t) + -1._id)**2*30._id + 91._id + (Cos(t) + -1._id)*105._id)*Sin(t*0.5_id)**12
endif

if(l==8.and.m1==6.and.m2==-5) then
   wigd=Cos(t*0.5_id)*2.1602468994692867436553224786959988859017347690194_id*((Cos(t) + -1._id)**2*30._id + 78._id + (Co&
&s(t) + -1._id)*97.5_id)*Sin(t*0.5_id)**11
endif

if(l==8.and.m1==6.and.m2==-4) then
   wigd=Cos(t*0.5_id)**2*3.8944404818493075368873949642026809176529333861249_id*((Cos(t) + -1._id)**2*30._id + 66._id + &
&(Cos(t) + -1._id)*90._id)*Sin(t*0.5_id)**10
endif

if(l==8.and.m1==6.and.m2==-3) then
   wigd=Cos(t*0.5_id)**3*6.0332412515993424345033528849187730050984461505421_id*((Cos(t) + -1._id)**2*30._id + 55._id + &
&(Cos(t) + -1._id)*82.5_id)*Sin(t*0.5_id)**9
endif

if(l==8.and.m1==6.and.m2==-2) then
   wigd=Cos(t*0.5_id)**4*8.1690472720711644400267750682082849910494805874321_id*((Cos(t) + -1._id)**2*30._id + 45._id + &
&(Cos(t) + -1._id)*75._id)*Sin(t*0.5_id)**8
endif

if(l==8.and.m1==6.and.m2==-1) then
   wigd=Cos(t*0.5_id)**5*9.7638790105845398750486793705062209513199095596437_id*((Cos(t) + -1._id)**2*30._id + 36._id + &
&(Cos(t) + -1._id)*67.5_id)*Sin(t*0.5_id)**7
endif

if(l==8.and.m1==6.and.m2==0) then
   wigd=Cos(t*0.5_id)**6*10.356157588603989566078588172024067830946506245733_id*(28._id + (Cos(t) + -1._id)**2*30._id + &
&(Cos(t) + -1._id)*60._id)*Sin(t*0.5_id)**6
endif

if(l==8.and.m1==6.and.m2==1) then
   wigd=Cos(t*0.5_id)**7*9.7638790105845398750486793705062209513199095596437_id*(21._id + (Cos(t) + -1._id)**2*30._id + &
&(Cos(t) + -1._id)*52.5_id)*Sin(t*0.5_id)**5
endif

if(l==8.and.m1==6.and.m2==2) then
   wigd=Cos(t*0.5_id)**8*8.1690472720711644400267750682082849910494805874321_id*(15._id + (Cos(t) + -1._id)**2*30._id + &
&(Cos(t) + -1._id)*45._id)*Sin(t*0.5_id)**4
endif

if(l==8.and.m1==6.and.m2==3) then
   wigd=Cos(t*0.5_id)**9*6.0332412515993424345033528849187730050984461505421_id*(10._id + (Cos(t) + -1._id)**2*30._id + &
&(Cos(t) + -1._id)*37.5_id)*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==6.and.m2==4) then
   wigd=Cos(t*0.5_id)**10*3.8944404818493075368873949642026809176529333861249_id*(6._id + (Cos(t) + -1._id)*30._id + (Co&
&s(t) + -1._id)**2*30._id)*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==6.and.m2==5) then
   wigd=Cos(t*0.5_id)**11*2.1602468994692867436553224786959988859017347690194_id*(3._id + (Cos(t) + -1._id)*22.5_id + (C&
&os(t) + -1._id)**2*30._id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==6.and.m2==6) then
   wigd=Cos(t*0.5_id)**12*(1._id + (Cos(t) + -1._id)*15._id + (Cos(t) + -1._id)**2*30._id)
endif

if(l==8.and.m1==6.and.m2==7) then
   wigd=Cos(t*0.5_id)**13*-1.369306393762915283642424457002005334881861737495_id*(-12._id + Cos(t)*16._id)*Sin(t*0.5_id)&
&
endif

if(l==8.and.m1==6.and.m2==8) then
   wigd=Cos(t*0.5_id)**14*10.95445115010332226913939565601604267905489389996_id*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==7.and.m2==-8) then
   wigd=Cos(t*0.5_id)*4._id*Sin(t*0.5_id)**15
endif

if(l==8.and.m1==7.and.m2==-7) then
   wigd=0.5_id*(14._id + Cos(t)*16._id)*Sin(t*0.5_id)**14
endif

if(l==8.and.m1==7.and.m2==-6) then
   wigd=Cos(t*0.5_id)*1.369306393762915283642424457002005334881861737495_id*(12._id + Cos(t)*16._id)*Sin(t*0.5_id)**13
endif

if(l==8.and.m1==7.and.m2==-5) then
   wigd=Cos(t*0.5_id)**2*2.9580398915498080212836641457808085242077506153972_id*(10._id + Cos(t)*16._id)*Sin(t*0.5_id)**&
&12
endif

if(l==8.and.m1==7.and.m2==-4) then
   wigd=Cos(t*0.5_id)**3*5.3326822519253854377215978803358864040832104516708_id*(8._id + Cos(t)*16._id)*Sin(t*0.5_id)**1&
&1
endif

if(l==8.and.m1==7.and.m2==-3) then
   wigd=Cos(t*0.5_id)**4*8.261355820929153030908713156308110206419189917898_id*(6._id + Cos(t)*16._id)*Sin(t*0.5_id)**10&
&
endif

if(l==8.and.m1==7.and.m2==-2) then
   wigd=Cos(t*0.5_id)**5*11.185928660598546835343474364957210078737641410745_id*(4._id + Cos(t)*16._id)*Sin(t*0.5_id)**9&
&
endif

if(l==8.and.m1==7.and.m2==-1) then
   wigd=Cos(t*0.5_id)**6*13.369741957120937642677996754046981984624107898354_id*(2._id + Cos(t)*16._id)*Sin(t*0.5_id)**8&
&
endif

if(l==8.and.m1==7.and.m2==0) then
   wigd=Cos(t)*Cos(t*0.5_id)**7*226.89204481426844419667465463338796920571321566023_id*Sin(t*0.5_id)**7
endif

if(l==8.and.m1==7.and.m2==1) then
   wigd=Cos(t*0.5_id)**8*13.369741957120937642677996754046981984624107898354_id*(-2._id + Cos(t)*16._id)*Sin(t*0.5_id)**&
&6
endif

if(l==8.and.m1==7.and.m2==2) then
   wigd=Cos(t*0.5_id)**9*11.185928660598546835343474364957210078737641410745_id*(-4._id + Cos(t)*16._id)*Sin(t*0.5_id)**&
&5
endif

if(l==8.and.m1==7.and.m2==3) then
   wigd=Cos(t*0.5_id)**10*8.261355820929153030908713156308110206419189917898_id*(-6._id + Cos(t)*16._id)*Sin(t*0.5_id)**&
&4
endif

if(l==8.and.m1==7.and.m2==4) then
   wigd=Cos(t*0.5_id)**11*5.3326822519253854377215978803358864040832104516708_id*(-8._id + Cos(t)*16._id)*Sin(t*0.5_id)*&
&*3
endif

if(l==8.and.m1==7.and.m2==5) then
   wigd=Cos(t*0.5_id)**12*2.9580398915498080212836641457808085242077506153972_id*(-10._id + Cos(t)*16._id)*Sin(t*0.5_id)&
&**2
endif

if(l==8.and.m1==7.and.m2==6) then
   wigd=Cos(t*0.5_id)**13*1.369306393762915283642424457002005334881861737495_id*(-12._id + Cos(t)*16._id)*Sin(t*0.5_id)
endif

if(l==8.and.m1==7.and.m2==7) then
   wigd=Cos(t*0.5_id)**14*0.5_id*(-14._id + Cos(t)*16._id)
endif

if(l==8.and.m1==7.and.m2==8) then
   wigd=Cos(t*0.5_id)**15*-4._id*Sin(t*0.5_id)
endif

if(l==8.and.m1==8.and.m2==-8) then
   wigd=Sin(t*0.5_id)**16
endif

if(l==8.and.m1==8.and.m2==-7) then
   wigd=Cos(t*0.5_id)*4._id*Sin(t*0.5_id)**15
endif

if(l==8.and.m1==8.and.m2==-6) then
   wigd=Cos(t*0.5_id)**2*10.95445115010332226913939565601604267905489389996_id*Sin(t*0.5_id)**14
endif

if(l==8.and.m1==8.and.m2==-5) then
   wigd=Cos(t*0.5_id)**3*23.664319132398464170269313166246468193662004923177_id*Sin(t*0.5_id)**13
endif

if(l==8.and.m1==8.and.m2==-4) then
   wigd=Cos(t*0.5_id)**4*42.661458015403083501772783042687091232665683613366_id*Sin(t*0.5_id)**12
endif

if(l==8.and.m1==8.and.m2==-3) then
   wigd=Cos(t*0.5_id)**5*66.09084656743322424726970525046488165135351934318_id*Sin(t*0.5_id)**11
endif

if(l==8.and.m1==8.and.m2==-2) then
   wigd=Cos(t*0.5_id)**6*89.48742928478837468274779491965768062990113128596_id*Sin(t*0.5_id)**10
endif

if(l==8.and.m1==8.and.m2==-1) then
   wigd=Cos(t*0.5_id)**7*106.95793565696750114142397403237585587699286318683_id*Sin(t*0.5_id)**9
endif

if(l==8.and.m1==8.and.m2==0) then
   wigd=Cos(t*0.5_id)**8*113.44602240713422209833732731669398460285660783012_id*Sin(t*0.5_id)**8
endif

if(l==8.and.m1==8.and.m2==1) then
   wigd=Cos(t*0.5_id)**9*106.95793565696750114142397403237585587699286318683_id*Sin(t*0.5_id)**7
endif

if(l==8.and.m1==8.and.m2==2) then
   wigd=Cos(t*0.5_id)**10*89.48742928478837468274779491965768062990113128596_id*Sin(t*0.5_id)**6
endif

if(l==8.and.m1==8.and.m2==3) then
   wigd=Cos(t*0.5_id)**11*66.09084656743322424726970525046488165135351934318_id*Sin(t*0.5_id)**5
endif

if(l==8.and.m1==8.and.m2==4) then
   wigd=Cos(t*0.5_id)**12*42.661458015403083501772783042687091232665683613366_id*Sin(t*0.5_id)**4
endif

if(l==8.and.m1==8.and.m2==5) then
   wigd=Cos(t*0.5_id)**13*23.664319132398464170269313166246468193662004923177_id*Sin(t*0.5_id)**3
endif

if(l==8.and.m1==8.and.m2==6) then
   wigd=Cos(t*0.5_id)**14*10.95445115010332226913939565601604267905489389996_id*Sin(t*0.5_id)**2
endif

if(l==8.and.m1==8.and.m2==7) then
   wigd=Cos(t*0.5_id)**15*4._id*Sin(t*0.5_id)
endif

if(l==8.and.m1==8.and.m2==8) then
   wigd=Cos(t*0.5_id)**16
endif

if(l==9.and.m1==-9.and.m2==-9) then
   wigd=Cos(t*0.5_id)**18
endif

if(l==9.and.m1==-9.and.m2==-8) then
   wigd=Cos(t*0.5_id)**17*-4.2426406871192851464050661726290942357090156261308_id*Sin(t*0.5_id)
endif

if(l==9.and.m1==-9.and.m2==-7) then
   wigd=Cos(t*0.5_id)**16*12.369316876852981649464229567922231075441597676121_id*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==-9.and.m2==-6) then
   wigd=Cos(t*0.5_id)**15*-28.565713714171399991997599245469061115064684639611_id*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==-9.and.m2==-5) then
   wigd=Cos(t*0.5_id)**14*55.317266743757323860013645690576758943480830292335_id*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==-9.and.m2==-4) then
   wigd=Cos(t*0.5_id)**13*-92.56349172324907493447268322253272345128247145089_id*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==-9.and.m2==-3) then
   wigd=Cos(t*0.5_id)**12*136.24977064200878866740367808992632648136146844452_id*Sin(t*0.5_id)**6
endif

if(l==9.and.m1==-9.and.m2==-2) then
   wigd=Cos(t*0.5_id)**11*-178.39282496782206627134409856717675981738782616465_id*Sin(t*0.5_id)**7
endif

if(l==9.and.m1==-9.and.m2==-1) then
   wigd=Cos(t*0.5_id)**10*209.18412941712380298700321806359345076313674734586_id*Sin(t*0.5_id)**8
endif

if(l==9.and.m1==-9.and.m2==0) then
   wigd=Cos(t*0.5_id)**9*-220.49943310584723587202231639809357727806572375952_id*Sin(t*0.5_id)**9
endif

if(l==9.and.m1==-9.and.m2==1) then
   wigd=Cos(t*0.5_id)**8*209.18412941712380298700321806359345076313674734586_id*Sin(t*0.5_id)**10
endif

if(l==9.and.m1==-9.and.m2==2) then
   wigd=Cos(t*0.5_id)**7*-178.39282496782206627134409856717675981738782616465_id*Sin(t*0.5_id)**11
endif

if(l==9.and.m1==-9.and.m2==3) then
   wigd=Cos(t*0.5_id)**6*136.24977064200878866740367808992632648136146844452_id*Sin(t*0.5_id)**12
endif

if(l==9.and.m1==-9.and.m2==4) then
   wigd=Cos(t*0.5_id)**5*-92.56349172324907493447268322253272345128247145089_id*Sin(t*0.5_id)**13
endif

if(l==9.and.m1==-9.and.m2==5) then
   wigd=Cos(t*0.5_id)**4*55.317266743757323860013645690576758943480830292335_id*Sin(t*0.5_id)**14
endif

if(l==9.and.m1==-9.and.m2==6) then
   wigd=Cos(t*0.5_id)**3*-28.565713714171399991997599245469061115064684639611_id*Sin(t*0.5_id)**15
endif

if(l==9.and.m1==-9.and.m2==7) then
   wigd=Cos(t*0.5_id)**2*12.369316876852981649464229567922231075441597676121_id*Sin(t*0.5_id)**16
endif

if(l==9.and.m1==-9.and.m2==8) then
   wigd=Cos(t*0.5_id)*-4.2426406871192851464050661726290942357090156261308_id*Sin(t*0.5_id)**17
endif

if(l==9.and.m1==-9.and.m2==9) then
   wigd=Sin(t*0.5_id)**18
endif

if(l==9.and.m1==-8.and.m2==-9) then
   wigd=Cos(t*0.5_id)**17*4.2426406871192851464050661726290942357090156261308_id*Sin(t*0.5_id)
endif

if(l==9.and.m1==-8.and.m2==-8) then
   wigd=Cos(t*0.5_id)**16*0.5_id*(-16._id + Cos(t)*18._id)
endif

if(l==9.and.m1==-8.and.m2==-7) then
   wigd=Cos(t*0.5_id)**15*-1.4577379737113251177185382193863957691303495837215_id*(-14._id + Cos(t)*18._id)*Sin(t*0.5_id&
&)
endif

if(l==9.and.m1==-8.and.m2==-6) then
   wigd=Cos(t*0.5_id)**14*3.3665016461206926511211286390232002368679663214932_id*(-12._id + Cos(t)*18._id)*Sin(t*0.5_id)&
&**2
endif

if(l==9.and.m1==-8.and.m2==-5) then
   wigd=Cos(t*0.5_id)**13*-6.5192024052026487145829715574291844165280937789101_id*(-10._id + Cos(t)*18._id)*Sin(t*0.5_id&
&)**3
endif

if(l==9.and.m1==-8.and.m2==-4) then
   wigd=Cos(t*0.5_id)**12*10.9087121146357144115021544873729018728051337026946_id*(-8._id + Cos(t)*18._id)*Sin(t*0.5_id)&
&**4
endif

if(l==9.and.m1==-8.and.m2==-3) then
   wigd=Cos(t*0.5_id)**11*-16.057189459346032556961112792299908291321915539858_id*(-6._id + Cos(t)*18._id)*Sin(t*0.5_id)&
&**5
endif

if(l==9.and.m1==-8.and.m2==-2) then
   wigd=Cos(t*0.5_id)**10*(-4._id + Cos(t)*18._id)*21.023796041628638288419857057495807938240784503821_id*Sin(t*0.5_id)*&
&*6
endif

if(l==9.and.m1==-8.and.m2==-1) then
   wigd=Cos(t*0.5_id)**9*-24.652586071242100015888052988157175076326378239492_id*(-2._id + Cos(t)*18._id)*Sin(t*0.5_id)*&
&*7
endif

if(l==9.and.m1==-8.and.m2==0) then
   wigd=Cos(t)*Cos(t*0.5_id)**8*467.74993319080228383869644347115012625778024243715_id*Sin(t*0.5_id)**8
endif

if(l==9.and.m1==-8.and.m2==1) then
   wigd=Cos(t*0.5_id)**7*-24.652586071242100015888052988157175076326378239492_id*(2._id + Cos(t)*18._id)*Sin(t*0.5_id)**&
&9
endif

if(l==9.and.m1==-8.and.m2==2) then
   wigd=Cos(t*0.5_id)**6*(4._id + Cos(t)*18._id)*21.023796041628638288419857057495807938240784503821_id*Sin(t*0.5_id)**1&
&0
endif

if(l==9.and.m1==-8.and.m2==3) then
   wigd=Cos(t*0.5_id)**5*-16.057189459346032556961112792299908291321915539858_id*(6._id + Cos(t)*18._id)*Sin(t*0.5_id)**&
&11
endif

if(l==9.and.m1==-8.and.m2==4) then
   wigd=Cos(t*0.5_id)**4*10.9087121146357144115021544873729018728051337026946_id*(8._id + Cos(t)*18._id)*Sin(t*0.5_id)**&
&12
endif

if(l==9.and.m1==-8.and.m2==5) then
   wigd=Cos(t*0.5_id)**3*-6.5192024052026487145829715574291844165280937789101_id*(10._id + Cos(t)*18._id)*Sin(t*0.5_id)*&
&*13
endif

if(l==9.and.m1==-8.and.m2==6) then
   wigd=Cos(t*0.5_id)**2*3.3665016461206926511211286390232002368679663214932_id*(12._id + Cos(t)*18._id)*Sin(t*0.5_id)**&
&14
endif

if(l==9.and.m1==-8.and.m2==7) then
   wigd=Cos(t*0.5_id)*-1.4577379737113251177185382193863957691303495837215_id*(14._id + Cos(t)*18._id)*Sin(t*0.5_id)**15&
&
endif

if(l==9.and.m1==-8.and.m2==8) then
   wigd=0.5_id*(16._id + Cos(t)*18._id)*Sin(t*0.5_id)**16
endif

if(l==9.and.m1==-8.and.m2==9) then
   wigd=Cos(t*0.5_id)*-4.2426406871192851464050661726290942357090156261308_id*Sin(t*0.5_id)**17
endif

if(l==9.and.m1==-7.and.m2==-9) then
   wigd=Cos(t*0.5_id)**16*12.369316876852981649464229567922231075441597676121_id*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==-7.and.m2==-8) then
   wigd=Cos(t*0.5_id)**15*1.4577379737113251177185382193863957691303495837215_id*(-14._id + Cos(t)*18._id)*Sin(t*0.5_id)&
&
endif

if(l==9.and.m1==-7.and.m2==-7) then
   wigd=Cos(t*0.5_id)**14*(1._id + (Cos(t) + -1._id)*17._id + (Cos(t) + -1._id)**2*38.25_id)
endif

if(l==9.and.m1==-7.and.m2==-6) then
   wigd=Cos(t*0.5_id)**13*-2.3094010767585030580365951220078298225904070050805_id*(3._id + (Cos(t) + -1._id)*25.5_id + (&
&Cos(t) + -1._id)**2*38.25_id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==-7.and.m2==-5) then
   wigd=Cos(t*0.5_id)**12*4.4721359549995793928183473374625524708812367192231_id*(6._id + (Cos(t) + -1._id)*34._id + (Co&
&s(t) + -1._id)**2*38.25_id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==-7.and.m2==-4) then
   wigd=Cos(t*0.5_id)**11*-7.483314773547882771167497464633098603512039615557_id*(10._id + (Cos(t) + -1._id)**2*38.25_id&
& + (Cos(t) + -1._id)*42.5_id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==-7.and.m2==-3) then
   wigd=Cos(t*0.5_id)**10*11.015141094572204041211617541744146941892253223863_id*(15._id + (Cos(t) + -1._id)**2*38.25_id&
& + (Cos(t) + -1._id)*51._id)*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==-7.and.m2==-2) then
   wigd=Cos(t*0.5_id)**9*-14.422205101855957172476885069881983785005186295381_id*(21._id + (Cos(t) + -1._id)**2*38.25_id&
& + (Cos(t) + -1._id)*59.5_id)*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==-7.and.m2==-1) then
   wigd=Cos(t*0.5_id)**8*16.9115345252877628981725179227765953746630798716105_id*(28._id + (Cos(t) + -1._id)**2*38.25_id&
& + (Cos(t) + -1._id)*68._id)*Sin(t*0.5_id)**6
endif

if(l==9.and.m1==-7.and.m2==0) then
   wigd=Cos(t*0.5_id)**7*-17.826322609494583523570662338729309312832143864472_id*(36._id + (Cos(t) + -1._id)**2*38.25_id&
& + (Cos(t) + -1._id)*76.5_id)*Sin(t*0.5_id)**7
endif

if(l==9.and.m1==-7.and.m2==1) then
   wigd=Cos(t*0.5_id)**6*16.9115345252877628981725179227765953746630798716105_id*((Cos(t) + -1._id)**2*38.25_id + 45._id&
& + (Cos(t) + -1._id)*85._id)*Sin(t*0.5_id)**8
endif

if(l==9.and.m1==-7.and.m2==2) then
   wigd=Cos(t*0.5_id)**5*-14.422205101855957172476885069881983785005186295381_id*((Cos(t) + -1._id)**2*38.25_id + 55._id&
& + (Cos(t) + -1._id)*93.5_id)*Sin(t*0.5_id)**9
endif

if(l==9.and.m1==-7.and.m2==3) then
   wigd=Cos(t*0.5_id)**4*11.015141094572204041211617541744146941892253223863_id*((Cos(t) + -1._id)**2*38.25_id + 66._id &
&+ (Cos(t) + -1._id)*102._id)*Sin(t*0.5_id)**10
endif

if(l==9.and.m1==-7.and.m2==4) then
   wigd=Cos(t*0.5_id)**3*-7.483314773547882771167497464633098603512039615557_id*((Cos(t) + -1._id)**2*38.25_id + 78._id &
&+ (Cos(t) + -1._id)*110.5_id)*Sin(t*0.5_id)**11
endif

if(l==9.and.m1==-7.and.m2==5) then
   wigd=Cos(t*0.5_id)**2*4.4721359549995793928183473374625524708812367192231_id*((Cos(t) + -1._id)**2*38.25_id + 91._id &
&+ (Cos(t) + -1._id)*119._id)*Sin(t*0.5_id)**12
endif

if(l==9.and.m1==-7.and.m2==6) then
   wigd=Cos(t*0.5_id)*-2.3094010767585030580365951220078298225904070050805_id*((Cos(t) + -1._id)**2*38.25_id + 105._id +&
& (Cos(t) + -1._id)*127.5_id)*Sin(t*0.5_id)**13
endif

if(l==9.and.m1==-7.and.m2==7) then
   wigd=((Cos(t) + -1._id)**2*38.25_id + 120._id + (Cos(t) + -1._id)*136._id)*Sin(t*0.5_id)**14
endif

if(l==9.and.m1==-7.and.m2==8) then
   wigd=Cos(t*0.5_id)*-1.4577379737113251177185382193863957691303495837215_id*(14._id + Cos(t)*18._id)*Sin(t*0.5_id)**15&
&
endif

if(l==9.and.m1==-7.and.m2==9) then
   wigd=Cos(t*0.5_id)**2*12.369316876852981649464229567922231075441597676121_id*Sin(t*0.5_id)**16
endif

if(l==9.and.m1==-6.and.m2==-9) then
   wigd=Cos(t*0.5_id)**15*28.565713714171399991997599245469061115064684639611_id*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==-6.and.m2==-8) then
   wigd=Cos(t*0.5_id)**14*3.3665016461206926511211286390232002368679663214932_id*(-12._id + Cos(t)*18._id)*Sin(t*0.5_id)&
&**2
endif

if(l==9.and.m1==-6.and.m2==-7) then
   wigd=Cos(t*0.5_id)**13*2.3094010767585030580365951220078298225904070050805_id*(3._id + (Cos(t) + -1._id)*25.5_id + (C&
&os(t) + -1._id)**2*38.25_id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==-6.and.m2==-6) then
   wigd=Cos(t*0.5_id)**12*(1._id + (Cos(t) + -1._id)*24._id + (Cos(t) + -1._id)**2*102._id + (Cos(t) + -1._id)**3*102._i&
&d)
endif

if(l==9.and.m1==-6.and.m2==-5) then
   wigd=Cos(t*0.5_id)**11*-1.9364916731037084425896326998911998054164608526458_id*(4._id + (Cos(t) + -1._id)*48._id + (C&
&os(t) + -1._id)**3*102._id + (Cos(t) + -1._id)**2*136._id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==-6.and.m2==-4) then
   wigd=Cos(t*0.5_id)**10*3.2403703492039301154829837180439983288526021535292_id*(10._id + (Cos(t) + -1._id)*80._id + (C&
&os(t) + -1._id)**3*102._id + (Cos(t) + -1._id)**2*170._id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==-6.and.m2==-3) then
   wigd=Cos(t*0.5_id)**9*-4.7696960070847282457631079301161327012731171262527_id*(20._id + (Cos(t) + -1._id)**3*102._id &
&+ (Cos(t) + -1._id)*120._id + (Cos(t) + -1._id)**2*204._id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==-6.and.m2==-2) then
   wigd=Cos(t*0.5_id)**8*6.2449979983983982058468931209397944610729599779917_id*(35._id + (Cos(t) + -1._id)**3*102._id +&
& (Cos(t) + -1._id)*168._id + (Cos(t) + -1._id)**2*238._id)*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==-6.and.m2==-1) then
   wigd=Cos(t*0.5_id)**7*-7.322909257938404906286509527879665713489932169733_id*(56._id + (Cos(t) + -1._id)**3*102._id +&
& (Cos(t) + -1._id)*224._id + (Cos(t) + -1._id)**2*272._id)*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==-6.and.m2==0) then
   wigd=Cos(t*0.5_id)**6*7.719024117939607353441468160313305396623552881689_id*(84._id + (Cos(t) + -1._id)**3*102._id + &
&(Cos(t) + -1._id)*288._id + (Cos(t) + -1._id)**2*306._id)*Sin(t*0.5_id)**6
endif

if(l==9.and.m1==-6.and.m2==1) then
   wigd=Cos(t*0.5_id)**5*-7.322909257938404906286509527879665713489932169733_id*((Cos(t) + -1._id)**3*102._id + 120._id &
&+ (Cos(t) + -1._id)**2*340._id + (Cos(t) + -1._id)*360._id)*Sin(t*0.5_id)**7
endif

if(l==9.and.m1==-6.and.m2==2) then
   wigd=Cos(t*0.5_id)**4*6.2449979983983982058468931209397944610729599779917_id*((Cos(t) + -1._id)**3*102._id + 165._id &
&+ (Cos(t) + -1._id)**2*374._id + (Cos(t) + -1._id)*440._id)*Sin(t*0.5_id)**8
endif

if(l==9.and.m1==-6.and.m2==3) then
   wigd=Cos(t*0.5_id)**3*-4.7696960070847282457631079301161327012731171262527_id*((Cos(t) + -1._id)**3*102._id + 220._id&
& + (Cos(t) + -1._id)**2*408._id + (Cos(t) + -1._id)*528._id)*Sin(t*0.5_id)**9
endif

if(l==9.and.m1==-6.and.m2==4) then
   wigd=Cos(t*0.5_id)**2*3.2403703492039301154829837180439983288526021535292_id*((Cos(t) + -1._id)**3*102._id + 286._id &
&+ (Cos(t) + -1._id)**2*442._id + (Cos(t) + -1._id)*624._id)*Sin(t*0.5_id)**10
endif

if(l==9.and.m1==-6.and.m2==5) then
   wigd=Cos(t*0.5_id)*-1.9364916731037084425896326998911998054164608526458_id*((Cos(t) + -1._id)**3*102._id + 364._id + &
&(Cos(t) + -1._id)**2*476._id + (Cos(t) + -1._id)*728._id)*Sin(t*0.5_id)**11
endif

if(l==9.and.m1==-6.and.m2==6) then
   wigd=((Cos(t) + -1._id)**3*102._id + 455._id + (Cos(t) + -1._id)**2*510._id + (Cos(t) + -1._id)*840._id)*Sin(t*0.5_id&
&)**12
endif

if(l==9.and.m1==-6.and.m2==7) then
   wigd=Cos(t*0.5_id)*-2.3094010767585030580365951220078298225904070050805_id*((Cos(t) + -1._id)**2*38.25_id + 105._id +&
& (Cos(t) + -1._id)*127.5_id)*Sin(t*0.5_id)**13
endif

if(l==9.and.m1==-6.and.m2==8) then
   wigd=Cos(t*0.5_id)**2*3.3665016461206926511211286390232002368679663214932_id*(12._id + Cos(t)*18._id)*Sin(t*0.5_id)**&
&14
endif

if(l==9.and.m1==-6.and.m2==9) then
   wigd=Cos(t*0.5_id)**3*-28.565713714171399991997599245469061115064684639611_id*Sin(t*0.5_id)**15
endif

if(l==9.and.m1==-5.and.m2==-9) then
   wigd=Cos(t*0.5_id)**14*55.317266743757323860013645690576758943480830292335_id*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==-5.and.m2==-8) then
   wigd=Cos(t*0.5_id)**13*6.5192024052026487145829715574291844165280937789101_id*(-10._id + Cos(t)*18._id)*Sin(t*0.5_id)&
&**3
endif

if(l==9.and.m1==-5.and.m2==-7) then
   wigd=Cos(t*0.5_id)**12*4.4721359549995793928183473374625524708812367192231_id*(6._id + (Cos(t) + -1._id)*34._id + (Co&
&s(t) + -1._id)**2*38.25_id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==-5.and.m2==-6) then
   wigd=Cos(t*0.5_id)**11*1.9364916731037084425896326998911998054164608526458_id*(4._id + (Cos(t) + -1._id)*48._id + (Co&
&s(t) + -1._id)**3*102._id + (Cos(t) + -1._id)**2*136._id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==-5.and.m2==-5) then
   wigd=Cos(t*0.5_id)**10*(1._id + (Cos(t) + -1._id)*30._id + (Cos(t) + -1._id)**2*180._id + (Cos(t) + -1._id)**4*191.25&
&_id + (Cos(t) + -1._id)**3*340._id)
endif

if(l==9.and.m1==-5.and.m2==-4) then
   wigd=Cos(t*0.5_id)**9*-1.6733200530681510959563440515703749787856307385973_id*(5._id + (Cos(t) + -1._id)*75._id + (Co&
&s(t) + -1._id)**4*191.25_id + (Cos(t) + -1._id)**2*300._id + (Cos(t) + -1._id)**3*425._id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==-5.and.m2==-3) then
   wigd=Cos(t*0.5_id)**8*2.4630604269214887994423741066531428037677316148781_id*(15._id + (Cos(t) + -1._id)*150._id + (C&
&os(t) + -1._id)**4*191.25_id + (Cos(t) + -1._id)**2*450._id + (Cos(t) + -1._id)**3*510._id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==-5.and.m2==-2) then
   wigd=Cos(t*0.5_id)**7*-3.2249030993194198609466452921215084524537585222434_id*(35._id + (Cos(t) + -1._id)**4*191.25_i&
&d + (Cos(t) + -1._id)*262.5_id + (Cos(t) + -1._id)**3*595._id + (Cos(t) + -1._id)**2*630._id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==-5.and.m2==-1) then
   wigd=Cos(t*0.5_id)**6*3.7815340802378074032779109105564661534285535943372_id*(70._id + (Cos(t) + -1._id)**4*191.25_id&
& + (Cos(t) + -1._id)*420._id + (Cos(t) + -1._id)**3*680._id + (Cos(t) + -1._id)**2*840._id)*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==-5.and.m2==0) then
   wigd=Cos(t*0.5_id)**5*-3.9860869143671326737099469187318722046636954283791_id*(126._id + (Cos(t) + -1._id)**4*191.25_&
&id + (Cos(t) + -1._id)*630._id + (Cos(t) + -1._id)**3*765._id + (Cos(t) + -1._id)**2*1080._id)*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==-5.and.m2==1) then
   wigd=Cos(t*0.5_id)**4*3.7815340802378074032779109105564661534285535943372_id*((Cos(t) + -1._id)**4*191.25_id + 210._i&
&d + (Cos(t) + -1._id)**3*850._id + (Cos(t) + -1._id)*900._id + (Cos(t) + -1._id)**2*1350._id)*Sin(t*0.5_id)**6
endif

if(l==9.and.m1==-5.and.m2==2) then
   wigd=Cos(t*0.5_id)**3*-3.2249030993194198609466452921215084524537585222434_id*((Cos(t) + -1._id)**4*191.25_id + 330._&
&id + (Cos(t) + -1._id)**3*935._id + (Cos(t) + -1._id)*1237.5_id + (Cos(t) + -1._id)**2*1650._id)*Sin(t*0.5_id)**7
endif

if(l==9.and.m1==-5.and.m2==3) then
   wigd=Cos(t*0.5_id)**2*2.4630604269214887994423741066531428037677316148781_id*((Cos(t) + -1._id)**4*191.25_id + 495._i&
&d + (Cos(t) + -1._id)**3*1020._id + (Cos(t) + -1._id)*1650._id + (Cos(t) + -1._id)**2*1980._id)*Sin(t*0.5_id)**8
endif

if(l==9.and.m1==-5.and.m2==4) then
   wigd=Cos(t*0.5_id)*-1.6733200530681510959563440515703749787856307385973_id*((Cos(t) + -1._id)**4*191.25_id + 715._id &
&+ (Cos(t) + -1._id)**3*1105._id + (Cos(t) + -1._id)*2145._id + (Cos(t) + -1._id)**2*2340._id)*Sin(t*0.5_id)**9
endif

if(l==9.and.m1==-5.and.m2==5) then
   wigd=((Cos(t) + -1._id)**4*191.25_id + 1001._id + (Cos(t) + -1._id)**3*1190._id + (Cos(t) + -1._id)*2730._id + (Cos(t&
&) + -1._id)**2*2730._id)*Sin(t*0.5_id)**10
endif

if(l==9.and.m1==-5.and.m2==6) then
   wigd=Cos(t*0.5_id)*-1.9364916731037084425896326998911998054164608526458_id*((Cos(t) + -1._id)**3*102._id + 364._id + &
&(Cos(t) + -1._id)**2*476._id + (Cos(t) + -1._id)*728._id)*Sin(t*0.5_id)**11
endif

if(l==9.and.m1==-5.and.m2==7) then
   wigd=Cos(t*0.5_id)**2*4.4721359549995793928183473374625524708812367192231_id*((Cos(t) + -1._id)**2*38.25_id + 91._id &
&+ (Cos(t) + -1._id)*119._id)*Sin(t*0.5_id)**12
endif

if(l==9.and.m1==-5.and.m2==8) then
   wigd=Cos(t*0.5_id)**3*-6.5192024052026487145829715574291844165280937789101_id*(10._id + Cos(t)*18._id)*Sin(t*0.5_id)*&
&*13
endif

if(l==9.and.m1==-5.and.m2==9) then
   wigd=Cos(t*0.5_id)**4*55.317266743757323860013645690576758943480830292335_id*Sin(t*0.5_id)**14
endif

if(l==9.and.m1==-4.and.m2==-9) then
   wigd=Cos(t*0.5_id)**13*92.56349172324907493447268322253272345128247145089_id*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==-4.and.m2==-8) then
   wigd=Cos(t*0.5_id)**12*10.9087121146357144115021544873729018728051337026946_id*(-8._id + Cos(t)*18._id)*Sin(t*0.5_id)&
&**4
endif

if(l==9.and.m1==-4.and.m2==-7) then
   wigd=Cos(t*0.5_id)**11*7.483314773547882771167497464633098603512039615557_id*(10._id + (Cos(t) + -1._id)**2*38.25_id &
&+ (Cos(t) + -1._id)*42.5_id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==-4.and.m2==-6) then
   wigd=Cos(t*0.5_id)**10*3.2403703492039301154829837180439983288526021535292_id*(10._id + (Cos(t) + -1._id)*80._id + (C&
&os(t) + -1._id)**3*102._id + (Cos(t) + -1._id)**2*170._id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==-4.and.m2==-5) then
   wigd=Cos(t*0.5_id)**9*1.67332005306815109595634405157037497878563073859734_id*(5._id + (Cos(t) + -1._id)*75._id + (Co&
&s(t) + -1._id)**4*191.25_id + (Cos(t) + -1._id)**2*300._id + (Cos(t) + -1._id)**3*425._id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==-4.and.m2==-4) then
   wigd=Cos(t*0.5_id)**8*(1._id + (Cos(t) + -1._id)*35._id + (Cos(t) + -1._id)**2*262.5_id + (Cos(t) + -1._id)**5*267.75&
&_id + (Cos(t) + -1._id)**3*700._id + (Cos(t) + -1._id)**4*743.75_id)
endif

if(l==9.and.m1==-4.and.m2==-3) then
   wigd=Cos(t*0.5_id)**7*-1.4719601443879744757940071211598756606957732468219_id*(6._id + (Cos(t) + -1._id)*105._id + (C&
&os(t) + -1._id)**5*267.75_id + (Cos(t) + -1._id)**2*525._id + (Cos(t) + -1._id)**4*892.5_id + (Cos(t) + -1._id)**3*1050._id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==-4.and.m2==-2) then
   wigd=Cos(t*0.5_id)**6*1.92724822331886306650718651592793499606505258574536_id*(21._id + (Cos(t) + -1._id)*245._id + (&
&Cos(t) + -1._id)**5*267.75_id + (Cos(t) + -1._id)**2*918.75_id + (Cos(t) + -1._id)**4*1041.25_id + (Cos(t) + -1._id)**3*1470._id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==-4.and.m2==-1) then
   wigd=Cos(t*0.5_id)**5*-2.25989885993662493879330654177149879193428953707_id*(56._id + (Cos(t) + -1._id)**5*267.75_id &
&+ (Cos(t) + -1._id)*490._id + (Cos(t) + -1._id)**4*1190._id + (Cos(t) + -1._id)**2*1470._id + (Cos(t) + -1._id)**3*1960._id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==-4.and.m2==0) then
   wigd=Cos(t*0.5_id)**4*2.3821425596725261067220435421431184595895000743328_id*(126._id + (Cos(t) + -1._id)**5*267.75_i&
&d + (Cos(t) + -1._id)*882._id + (Cos(t) + -1._id)**4*1338.75_id + (Cos(t) + -1._id)**2*2205._id + (Cos(t) + -1._id)**3*2520._id)*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==-4.and.m2==1) then
   wigd=Cos(t*0.5_id)**3*-2.25989885993662493879330654177149879193428953707_id*(252._id + (Cos(t) + -1._id)**5*267.75_id&
& + (Cos(t) + -1._id)*1470._id + (Cos(t) + -1._id)**4*1487.5_id + (Cos(t) + -1._id)**2*3150._id + (Cos(t) + -1._id)**3*3150._id)*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==-4.and.m2==2) then
   wigd=Cos(t*0.5_id)**2*1.92724822331886306650718651592793499606505258574536_id*((Cos(t) + -1._id)**5*267.75_id + 462._&
&id + (Cos(t) + -1._id)**4*1636.25_id + (Cos(t) + -1._id)*2310._id + (Cos(t) + -1._id)**3*3850._id + (Cos(t) + -1._id)**2*4331.25_id)*Sin(t*0.5_id)**6
endif

if(l==9.and.m1==-4.and.m2==3) then
   wigd=Cos(t*0.5_id)*-1.4719601443879744757940071211598756606957732468219_id*((Cos(t) + -1._id)**5*267.75_id + 792._id &
&+ (Cos(t) + -1._id)**4*1785._id + (Cos(t) + -1._id)*3465._id + (Cos(t) + -1._id)**3*4620._id + (Cos(t) + -1._id)**2*5775._id)*Sin(t*0.5_id)**7
endif

if(l==9.and.m1==-4.and.m2==4) then
   wigd=((Cos(t) + -1._id)**5*267.75_id + 1287._id + (Cos(t) + -1._id)**4*1933.75_id + (Cos(t) + -1._id)*5005._id + (Cos&
&(t) + -1._id)**3*5460._id + (Cos(t) + -1._id)**2*7507.5_id)*Sin(t*0.5_id)**8
endif

if(l==9.and.m1==-4.and.m2==5) then
   wigd=Cos(t*0.5_id)*-1.6733200530681510959563440515703749787856307385973_id*((Cos(t) + -1._id)**4*191.25_id + 715._id &
&+ (Cos(t) + -1._id)**3*1105._id + (Cos(t) + -1._id)*2145._id + (Cos(t) + -1._id)**2*2340._id)*Sin(t*0.5_id)**9
endif

if(l==9.and.m1==-4.and.m2==6) then
   wigd=Cos(t*0.5_id)**2*3.2403703492039301154829837180439983288526021535292_id*((Cos(t) + -1._id)**3*102._id + 286._id &
&+ (Cos(t) + -1._id)**2*442._id + (Cos(t) + -1._id)*624._id)*Sin(t*0.5_id)**10
endif

if(l==9.and.m1==-4.and.m2==7) then
   wigd=Cos(t*0.5_id)**3*-7.483314773547882771167497464633098603512039615557_id*((Cos(t) + -1._id)**2*38.25_id + 78._id &
&+ (Cos(t) + -1._id)*110.5_id)*Sin(t*0.5_id)**11
endif

if(l==9.and.m1==-4.and.m2==8) then
   wigd=Cos(t*0.5_id)**4*10.9087121146357144115021544873729018728051337026946_id*(8._id + Cos(t)*18._id)*Sin(t*0.5_id)**&
&12
endif

if(l==9.and.m1==-4.and.m2==9) then
   wigd=Cos(t*0.5_id)**5*-92.56349172324907493447268322253272345128247145089_id*Sin(t*0.5_id)**13
endif

if(l==9.and.m1==-3.and.m2==-9) then
   wigd=Cos(t*0.5_id)**12*136.24977064200878866740367808992632648136146844452_id*Sin(t*0.5_id)**6
endif

if(l==9.and.m1==-3.and.m2==-8) then
   wigd=Cos(t*0.5_id)**11*16.0571894593460325569611127922999082913219155398578_id*(-6._id + Cos(t)*18._id)*Sin(t*0.5_id)&
&**5
endif

if(l==9.and.m1==-3.and.m2==-7) then
   wigd=Cos(t*0.5_id)**10*11.015141094572204041211617541744146941892253223863_id*(15._id + (Cos(t) + -1._id)**2*38.25_id&
& + (Cos(t) + -1._id)*51._id)*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==-3.and.m2==-6) then
   wigd=Cos(t*0.5_id)**9*4.7696960070847282457631079301161327012731171262527_id*(20._id + (Cos(t) + -1._id)**3*102._id +&
& (Cos(t) + -1._id)*120._id + (Cos(t) + -1._id)**2*204._id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==-3.and.m2==-5) then
   wigd=Cos(t*0.5_id)**8*2.4630604269214887994423741066531428037677316148781_id*(15._id + (Cos(t) + -1._id)*150._id + (C&
&os(t) + -1._id)**4*191.25_id + (Cos(t) + -1._id)**2*450._id + (Cos(t) + -1._id)**3*510._id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==-3.and.m2==-4) then
   wigd=Cos(t*0.5_id)**7*1.4719601443879744757940071211598756606957732468219_id*(6._id + (Cos(t) + -1._id)*105._id + (Co&
&s(t) + -1._id)**5*267.75_id + (Cos(t) + -1._id)**2*525._id + (Cos(t) + -1._id)**4*892.5_id + (Cos(t) + -1._id)**3*1050._id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==-3.and.m2==-3) then
   wigd=Cos(t*0.5_id)**6*(1._id + (Cos(t) + -1._id)*39._id + (Cos(t) + -1._id)**6*290.0625_id + (Cos(t) + -1._id)**2*341&
&.25_id + (Cos(t) + -1._id)**3*1137.5_id + (Cos(t) + -1._id)**5*1160.25_id + (Cos(t) + -1._id)**4*1706.25_id)
endif

if(l==9.and.m1==-3.and.m2==-2) then
   wigd=Cos(t*0.5_id)**5*-1.3093073414159542875965849124937167111384161647908_id*(7._id + (Cos(t) + -1._id)*136.5_id + (&
&Cos(t) + -1._id)**6*290.0625_id + (Cos(t) + -1._id)**2*796.25_id + (Cos(t) + -1._id)**5*1353.625_id + (Cos(t) + -1._id)**3*1990.625_id + (Cos(t) + -1._id)**4*2388.75_id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==-3.and.m2==-1) then
   wigd=Cos(t*0.5_id)**4*1.53529894715747693137457218396297949971170131913177_id*(28._id + (Cos(t) + -1._id)**6*290.0625&
&_id + (Cos(t) + -1._id)*364._id + (Cos(t) + -1._id)**5*1547._id + (Cos(t) + -1._id)**2*1592.5_id + (Cos(t) + -1._id)**3*3185._id + (Cos(t) + -1._id)**4*3185._id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==-3.and.m2==0) then
   wigd=Cos(t*0.5_id)**3*-1.6183471874253741377307144112820922936791121082957_id*(84._id + (Cos(t) + -1._id)**6*290.0625&
&_id + (Cos(t) + -1._id)*819._id + (Cos(t) + -1._id)**5*1740.375_id + (Cos(t) + -1._id)**2*2866.5_id + (Cos(t) + -1._id)**4*4095._id + (Cos(t) + -1._id)**3*4777.5_id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==-3.and.m2==1) then
   wigd=Cos(t*0.5_id)**2*1.53529894715747693137457218396297949971170131913177_id*(210._id + (Cos(t) + -1._id)**6*290.062&
&5_id + (Cos(t) + -1._id)*1638._id + (Cos(t) + -1._id)**5*1933.75_id + (Cos(t) + -1._id)**2*4777.5_id + (Cos(t) + -1._id)**4*5118.75_id + (Cos(t) + -1._id)**3*6825._id)*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==-3.and.m2==2) then
   wigd=Cos(t*0.5_id)*-1.3093073414159542875965849124937167111384161647908_id*((Cos(t) + -1._id)**6*290.0625_id + 462._i&
&d + (Cos(t) + -1._id)**5*2127.125_id + (Cos(t) + -1._id)*3003._id + (Cos(t) + -1._id)**4*6256.25_id + (Cos(t) + -1._id)**2*7507.5_id + (Cos(t) + -1._id)**3*9384.375_id)*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==-3.and.m2==3) then
   wigd=((Cos(t) + -1._id)**6*290.0625_id + 924._id + (Cos(t) + -1._id)**5*2320.5_id + (Cos(t) + -1._id)*5148._id + (Cos&
&(t) + -1._id)**4*7507.5_id + (Cos(t) + -1._id)**2*11261.25_id + (Cos(t) + -1._id)**3*12512.5_id)*Sin(t*0.5_id)**6
endif

if(l==9.and.m1==-3.and.m2==4) then
   wigd=Cos(t*0.5_id)*-1.4719601443879744757940071211598756606957732468219_id*((Cos(t) + -1._id)**5*267.75_id + 792._id &
&+ (Cos(t) + -1._id)**4*1785._id + (Cos(t) + -1._id)*3465._id + (Cos(t) + -1._id)**3*4620._id + (Cos(t) + -1._id)**2*5775._id)*Sin(t*0.5_id)**7
endif

if(l==9.and.m1==-3.and.m2==5) then
   wigd=Cos(t*0.5_id)**2*2.4630604269214887994423741066531428037677316148781_id*((Cos(t) + -1._id)**4*191.25_id + 495._i&
&d + (Cos(t) + -1._id)**3*1020._id + (Cos(t) + -1._id)*1650._id + (Cos(t) + -1._id)**2*1980._id)*Sin(t*0.5_id)**8
endif

if(l==9.and.m1==-3.and.m2==6) then
   wigd=Cos(t*0.5_id)**3*-4.7696960070847282457631079301161327012731171262527_id*((Cos(t) + -1._id)**3*102._id + 220._id&
& + (Cos(t) + -1._id)**2*408._id + (Cos(t) + -1._id)*528._id)*Sin(t*0.5_id)**9
endif

if(l==9.and.m1==-3.and.m2==7) then
   wigd=Cos(t*0.5_id)**4*11.015141094572204041211617541744146941892253223863_id*((Cos(t) + -1._id)**2*38.25_id + 66._id &
&+ (Cos(t) + -1._id)*102._id)*Sin(t*0.5_id)**10
endif

if(l==9.and.m1==-3.and.m2==8) then
   wigd=Cos(t*0.5_id)**5*-16.057189459346032556961112792299908291321915539858_id*(6._id + Cos(t)*18._id)*Sin(t*0.5_id)**&
&11
endif

if(l==9.and.m1==-3.and.m2==9) then
   wigd=Cos(t*0.5_id)**6*136.24977064200878866740367808992632648136146844452_id*Sin(t*0.5_id)**12
endif

if(l==9.and.m1==-2.and.m2==-9) then
   wigd=Cos(t*0.5_id)**11*178.39282496782206627134409856717675981738782616465_id*Sin(t*0.5_id)**7
endif

if(l==9.and.m1==-2.and.m2==-8) then
   wigd=Cos(t*0.5_id)**10*(-4._id + Cos(t)*18._id)*21.023796041628638288419857057495807938240784503821_id*Sin(t*0.5_id)*&
&*6
endif

if(l==9.and.m1==-2.and.m2==-7) then
   wigd=Cos(t*0.5_id)**9*14.422205101855957172476885069881983785005186295381_id*(21._id + (Cos(t) + -1._id)**2*38.25_id &
&+ (Cos(t) + -1._id)*59.5_id)*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==-2.and.m2==-6) then
   wigd=Cos(t*0.5_id)**8*6.2449979983983982058468931209397944610729599779917_id*(35._id + (Cos(t) + -1._id)**3*102._id +&
& (Cos(t) + -1._id)*168._id + (Cos(t) + -1._id)**2*238._id)*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==-2.and.m2==-5) then
   wigd=Cos(t*0.5_id)**7*3.2249030993194198609466452921215084524537585222434_id*(35._id + (Cos(t) + -1._id)**4*191.25_id&
& + (Cos(t) + -1._id)*262.5_id + (Cos(t) + -1._id)**3*595._id + (Cos(t) + -1._id)**2*630._id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==-2.and.m2==-4) then
   wigd=Cos(t*0.5_id)**6*1.92724822331886306650718651592793499606505258574536_id*(21._id + (Cos(t) + -1._id)*245._id + (&
&Cos(t) + -1._id)**5*267.75_id + (Cos(t) + -1._id)**2*918.75_id + (Cos(t) + -1._id)**4*1041.25_id + (Cos(t) + -1._id)**3*1470._id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==-2.and.m2==-3) then
   wigd=Cos(t*0.5_id)**5*1.3093073414159542875965849124937167111384161647908_id*(7._id + (Cos(t) + -1._id)*136.5_id + (C&
&os(t) + -1._id)**6*290.0625_id + (Cos(t) + -1._id)**2*796.25_id + (Cos(t) + -1._id)**5*1353.625_id + (Cos(t) + -1._id)**3*1990.625_id + (Cos(t) + -1._id)**4*2388.75_id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==-2.and.m2==-2) then
   wigd=Cos(t*0.5_id)**4*(1._id + (Cos(t) + -1._id)*42._id + (Cos(t) + -1._id)**7*248.625_id + (Cos(t) + -1._id)**2*409.&
&5_id + (Cos(t) + -1._id)**6*1353.625_id + (Cos(t) + -1._id)**3*1592.5_id + (Cos(t) + -1._id)**5*2866.5_id + (Cos(t) + -1._id)**4*2985.9375_id)
endif

if(l==9.and.m1==-2.and.m2==-1) then
   wigd=Cos(t*0.5_id)**3*-1.1726039399558573886414075283861165701470570883529_id*(8._id + (Cos(t) + -1._id)*168._id + (C&
&os(t) + -1._id)**7*248.625_id + (Cos(t) + -1._id)**2*1092._id + (Cos(t) + -1._id)**6*1547._id + (Cos(t) + -1._id)**3*3185._id + (Cos(t) + -1._id)**5*3822._id + (Cos(t) + -1._id)**4*4777.5_id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==-2.and.m2==0) then
   wigd=Cos(t*0.5_id)**2*1.2360330811826104914518995734667855101633173865163_id*(36._id + (Cos(t) + -1._id)**7*248.625_i&
&d + (Cos(t) + -1._id)*504._id + (Cos(t) + -1._id)**6*1740.375_id + (Cos(t) + -1._id)**2*2457._id + (Cos(t) + -1._id)**5*4914._id + (Cos(t) + -1._id)**3*5733._id + (Cos(t) + -1._id)**4*7166.25_id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==-2.and.m2==1) then
   wigd=Cos(t*0.5_id)*-1.1726039399558573886414075283861165701470570883529_id*(120._id + (Cos(t) + -1._id)**7*248.625_id&
& + (Cos(t) + -1._id)*1260._id + (Cos(t) + -1._id)**6*1933.75_id + (Cos(t) + -1._id)**2*4914._id + (Cos(t) + -1._id)**5*6142.5_id + (Cos(t) + -1._id)**3*9555._id + (Cos(t) + -1._id)**4*10237.5_id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==-2.and.m2==2) then
   wigd=((Cos(t) + -1._id)**7*248.625_id + 330._id + (Cos(t) + -1._id)**6*2127.125_id + (Cos(t) + -1._id)*2772._id + (Co&
&s(t) + -1._id)**5*7507.5_id + (Cos(t) + -1._id)**2*9009._id + (Cos(t) + -1._id)**4*14076.5625_id + (Cos(t) + -1._id)**3*15015._id)*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==-2.and.m2==3) then
   wigd=Cos(t*0.5_id)*-1.3093073414159542875965849124937167111384161647908_id*((Cos(t) + -1._id)**6*290.0625_id + 462._i&
&d + (Cos(t) + -1._id)**5*2127.125_id + (Cos(t) + -1._id)*3003._id + (Cos(t) + -1._id)**4*6256.25_id + (Cos(t) + -1._id)**2*7507.5_id + (Cos(t) + -1._id)**3*9384.375_id)*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==-2.and.m2==4) then
   wigd=Cos(t*0.5_id)**2*1.92724822331886306650718651592793499606505258574536_id*((Cos(t) + -1._id)**5*267.75_id + 462._&
&id + (Cos(t) + -1._id)**4*1636.25_id + (Cos(t) + -1._id)*2310._id + (Cos(t) + -1._id)**3*3850._id + (Cos(t) + -1._id)**2*4331.25_id)*Sin(t*0.5_id)**6
endif

if(l==9.and.m1==-2.and.m2==5) then
   wigd=Cos(t*0.5_id)**3*-3.2249030993194198609466452921215084524537585222434_id*((Cos(t) + -1._id)**4*191.25_id + 330._&
&id + (Cos(t) + -1._id)**3*935._id + (Cos(t) + -1._id)*1237.5_id + (Cos(t) + -1._id)**2*1650._id)*Sin(t*0.5_id)**7
endif

if(l==9.and.m1==-2.and.m2==6) then
   wigd=Cos(t*0.5_id)**4*6.2449979983983982058468931209397944610729599779917_id*((Cos(t) + -1._id)**3*102._id + 165._id &
&+ (Cos(t) + -1._id)**2*374._id + (Cos(t) + -1._id)*440._id)*Sin(t*0.5_id)**8
endif

if(l==9.and.m1==-2.and.m2==7) then
   wigd=Cos(t*0.5_id)**5*-14.422205101855957172476885069881983785005186295381_id*((Cos(t) + -1._id)**2*38.25_id + 55._id&
& + (Cos(t) + -1._id)*93.5_id)*Sin(t*0.5_id)**9
endif

if(l==9.and.m1==-2.and.m2==8) then
   wigd=Cos(t*0.5_id)**6*(4._id + Cos(t)*18._id)*21.023796041628638288419857057495807938240784503821_id*Sin(t*0.5_id)**1&
&0
endif

if(l==9.and.m1==-2.and.m2==9) then
   wigd=Cos(t*0.5_id)**7*-178.39282496782206627134409856717675981738782616465_id*Sin(t*0.5_id)**11
endif

if(l==9.and.m1==-1.and.m2==-9) then
   wigd=Cos(t*0.5_id)**10*209.18412941712380298700321806359345076313674734586_id*Sin(t*0.5_id)**8
endif

if(l==9.and.m1==-1.and.m2==-8) then
   wigd=Cos(t*0.5_id)**9*(-2._id + Cos(t)*18._id)*24.652586071242100015888052988157175076326378239492_id*Sin(t*0.5_id)**&
&7
endif

if(l==9.and.m1==-1.and.m2==-7) then
   wigd=Cos(t*0.5_id)**8*16.9115345252877628981725179227765953746630798716105_id*(28._id + (Cos(t) + -1._id)**2*38.25_id&
& + (Cos(t) + -1._id)*68._id)*Sin(t*0.5_id)**6
endif

if(l==9.and.m1==-1.and.m2==-6) then
   wigd=Cos(t*0.5_id)**7*7.322909257938404906286509527879665713489932169733_id*(56._id + (Cos(t) + -1._id)**3*102._id + &
&(Cos(t) + -1._id)*224._id + (Cos(t) + -1._id)**2*272._id)*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==-1.and.m2==-5) then
   wigd=Cos(t*0.5_id)**6*3.7815340802378074032779109105564661534285535943372_id*(70._id + (Cos(t) + -1._id)**4*191.25_id&
& + (Cos(t) + -1._id)*420._id + (Cos(t) + -1._id)**3*680._id + (Cos(t) + -1._id)**2*840._id)*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==-1.and.m2==-4) then
   wigd=Cos(t*0.5_id)**5*2.25989885993662493879330654177149879193428953707_id*(56._id + (Cos(t) + -1._id)**5*267.75_id +&
& (Cos(t) + -1._id)*490._id + (Cos(t) + -1._id)**4*1190._id + (Cos(t) + -1._id)**2*1470._id + (Cos(t) + -1._id)**3*1960._id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==-1.and.m2==-3) then
   wigd=Cos(t*0.5_id)**4*1.53529894715747693137457218396297949971170131913177_id*(28._id + (Cos(t) + -1._id)**6*290.0625&
&_id + (Cos(t) + -1._id)*364._id + (Cos(t) + -1._id)**5*1547._id + (Cos(t) + -1._id)**2*1592.5_id + (Cos(t) + -1._id)**3*3185._id + (Cos(t) + -1._id)**4*3185._id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==-1.and.m2==-2) then
   wigd=Cos(t*0.5_id)**3*1.1726039399558573886414075283861165701470570883529_id*(8._id + (Cos(t) + -1._id)*168._id + (Co&
&s(t) + -1._id)**7*248.625_id + (Cos(t) + -1._id)**2*1092._id + (Cos(t) + -1._id)**6*1547._id + (Cos(t) + -1._id)**3*3185._id + (Cos(t) + -1._id)**5*3822._id + (Cos(t) + -1._id)**4*4777.5_id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==-1.and.m2==-1) then
   wigd=Cos(t*0.5_id)**2*(1._id + (Cos(t) + -1._id)*44._id + (Cos(t) + -1._id)**8*170.9296875_id + (Cos(t) + -1._id)**2*&
&462._id + (Cos(t) + -1._id)**7*1215.5_id + (Cos(t) + -1._id)**3*2002._id + (Cos(t) + -1._id)**6*3503.5_id + (Cos(t) + -1._id)**4*4379.375_id + (Cos(t) + -1._id)**5*5255.25_id)
endif

if(l==9.and.m1==-1.and.m2==0) then
   wigd=Cos(t*0.5_id)*-1.0540925533894597773329645148109061779065183797751_id*(9._id + (Cos(t) + -1._id)**8*170.9296875_&
&id + (Cos(t) + -1._id)*198._id + (Cos(t) + -1._id)**7*1367.4375_id + (Cos(t) + -1._id)**2*1386._id + (Cos(t) + -1._id)**3*4504.5_id + (Cos(t) + -1._id)**6*4504.5_id + (Cos(t) + -1._id)**4*7882.875_id + (Cos(t) + -1._id)**5*7882.875_id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==-1.and.m2==1) then
   wigd=(45._id + (Cos(t) + -1._id)**8*170.9296875_id + (Cos(t) + -1._id)*660._id + (Cos(t) + -1._id)**7*1519.375_id + (&
&Cos(t) + -1._id)**2*3465._id + (Cos(t) + -1._id)**6*5630.625_id + (Cos(t) + -1._id)**3*9009._id + (Cos(t) + -1._id)**5*11261.25_id + (Cos(t) + -1._id)**4*13138.125_id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==-1.and.m2==2) then
   wigd=Cos(t*0.5_id)*-1.1726039399558573886414075283861165701470570883529_id*(120._id + (Cos(t) + -1._id)**7*248.625_id&
& + (Cos(t) + -1._id)*1260._id + (Cos(t) + -1._id)**6*1933.75_id + (Cos(t) + -1._id)**2*4914._id + (Cos(t) + -1._id)**5*6142.5_id + (Cos(t) + -1._id)**3*9555._id + (Cos(t) + -1._id)**4*10237.5_id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==-1.and.m2==3) then
   wigd=Cos(t*0.5_id)**2*1.53529894715747693137457218396297949971170131913177_id*(210._id + (Cos(t) + -1._id)**6*290.062&
&5_id + (Cos(t) + -1._id)*1638._id + (Cos(t) + -1._id)**5*1933.75_id + (Cos(t) + -1._id)**2*4777.5_id + (Cos(t) + -1._id)**4*5118.75_id + (Cos(t) + -1._id)**3*6825._id)*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==-1.and.m2==4) then
   wigd=Cos(t*0.5_id)**3*-2.25989885993662493879330654177149879193428953707_id*(252._id + (Cos(t) + -1._id)**5*267.75_id&
& + (Cos(t) + -1._id)*1470._id + (Cos(t) + -1._id)**4*1487.5_id + (Cos(t) + -1._id)**2*3150._id + (Cos(t) + -1._id)**3*3150._id)*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==-1.and.m2==5) then
   wigd=Cos(t*0.5_id)**4*3.7815340802378074032779109105564661534285535943372_id*((Cos(t) + -1._id)**4*191.25_id + 210._i&
&d + (Cos(t) + -1._id)**3*850._id + (Cos(t) + -1._id)*900._id + (Cos(t) + -1._id)**2*1350._id)*Sin(t*0.5_id)**6
endif

if(l==9.and.m1==-1.and.m2==6) then
   wigd=Cos(t*0.5_id)**5*-7.322909257938404906286509527879665713489932169733_id*((Cos(t) + -1._id)**3*102._id + 120._id &
&+ (Cos(t) + -1._id)**2*340._id + (Cos(t) + -1._id)*360._id)*Sin(t*0.5_id)**7
endif

if(l==9.and.m1==-1.and.m2==7) then
   wigd=Cos(t*0.5_id)**6*16.9115345252877628981725179227765953746630798716105_id*((Cos(t) + -1._id)**2*38.25_id + 45._id&
& + (Cos(t) + -1._id)*85._id)*Sin(t*0.5_id)**8
endif

if(l==9.and.m1==-1.and.m2==8) then
   wigd=Cos(t*0.5_id)**7*-24.652586071242100015888052988157175076326378239492_id*(2._id + Cos(t)*18._id)*Sin(t*0.5_id)**&
&9
endif

if(l==9.and.m1==-1.and.m2==9) then
   wigd=Cos(t*0.5_id)**8*209.18412941712380298700321806359345076313674734586_id*Sin(t*0.5_id)**10
endif

if(l==9.and.m1==0.and.m2==-9) then
   wigd=Cos(t*0.5_id)**9*220.49943310584723587202231639809357727806572375952_id*Sin(t*0.5_id)**9
endif

if(l==9.and.m1==0.and.m2==-8) then
   wigd=Cos(t)*Cos(t*0.5_id)**8*467.74993319080228383869644347115012625778024243715_id*Sin(t*0.5_id)**8
endif

if(l==9.and.m1==0.and.m2==-7) then
   wigd=Cos(t*0.5_id)**7*17.826322609494583523570662338729309312832143864472_id*(36._id + (Cos(t) + -1._id)**2*38.25_id &
&+ (Cos(t) + -1._id)*76.5_id)*Sin(t*0.5_id)**7
endif

if(l==9.and.m1==0.and.m2==-6) then
   wigd=Cos(t*0.5_id)**6*7.719024117939607353441468160313305396623552881689_id*(84._id + (Cos(t) + -1._id)**3*102._id + &
&(Cos(t) + -1._id)*288._id + (Cos(t) + -1._id)**2*306._id)*Sin(t*0.5_id)**6
endif

if(l==9.and.m1==0.and.m2==-5) then
   wigd=Cos(t*0.5_id)**5*3.9860869143671326737099469187318722046636954283791_id*(126._id + (Cos(t) + -1._id)**4*191.25_i&
&d + (Cos(t) + -1._id)*630._id + (Cos(t) + -1._id)**3*765._id + (Cos(t) + -1._id)**2*1080._id)*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==0.and.m2==-4) then
   wigd=Cos(t*0.5_id)**4*2.3821425596725261067220435421431184595895000743328_id*(126._id + (Cos(t) + -1._id)**5*267.75_i&
&d + (Cos(t) + -1._id)*882._id + (Cos(t) + -1._id)**4*1338.75_id + (Cos(t) + -1._id)**2*2205._id + (Cos(t) + -1._id)**3*2520._id)*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==0.and.m2==-3) then
   wigd=Cos(t*0.5_id)**3*1.61834718742537413773071441128209229367911210829573_id*(84._id + (Cos(t) + -1._id)**6*290.0625&
&_id + (Cos(t) + -1._id)*819._id + (Cos(t) + -1._id)**5*1740.375_id + (Cos(t) + -1._id)**2*2866.5_id + (Cos(t) + -1._id)**4*4095._id + (Cos(t) + -1._id)**3*4777.5_id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==0.and.m2==-2) then
   wigd=Cos(t*0.5_id)**2*1.2360330811826104914518995734667855101633173865163_id*(36._id + (Cos(t) + -1._id)**7*248.625_i&
&d + (Cos(t) + -1._id)*504._id + (Cos(t) + -1._id)**6*1740.375_id + (Cos(t) + -1._id)**2*2457._id + (Cos(t) + -1._id)**5*4914._id + (Cos(t) + -1._id)**3*5733._id + (Cos(t) + -1._id)**4*7166.25_id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==0.and.m2==-1) then
   wigd=Cos(t*0.5_id)*1.0540925533894597773329645148109061779065183797751_id*(9._id + (Cos(t) + -1._id)**8*170.9296875_i&
&d + (Cos(t) + -1._id)*198._id + (Cos(t) + -1._id)**7*1367.4375_id + (Cos(t) + -1._id)**2*1386._id + (Cos(t) + -1._id)**3*4504.5_id + (Cos(t) + -1._id)**6*4504.5_id + (Cos(t) + -1._id)**4*7882.875_id + (Cos(t) + -1._id)**5*7882.875_id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==0.and.m2==0) then
   wigd=0.0078125_id*(Cos(t)**7*-25740._id + Cos(t)**3*-4620._id + Cos(t)*315._id + Cos(t)**9*12155._id + Cos(t)**5*1801&
&8._id)
endif

if(l==9.and.m1==0.and.m2==1) then
   wigd=Cos(t*0.5_id)*-1.0540925533894597773329645148109061779065183797751_id*(9._id + (Cos(t) + -1._id)**8*170.9296875_&
&id + (Cos(t) + -1._id)*198._id + (Cos(t) + -1._id)**7*1367.4375_id + (Cos(t) + -1._id)**2*1386._id + (Cos(t) + -1._id)**3*4504.5_id + (Cos(t) + -1._id)**6*4504.5_id + (Cos(t) + -1._id)**4*7882.875_id + (Cos(t) + -1._id)**5*7882.875_id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==0.and.m2==2) then
   wigd=Cos(t*0.5_id)**2*1.2360330811826104914518995734667855101633173865163_id*(36._id + (Cos(t) + -1._id)**7*248.625_i&
&d + (Cos(t) + -1._id)*504._id + (Cos(t) + -1._id)**6*1740.375_id + (Cos(t) + -1._id)**2*2457._id + (Cos(t) + -1._id)**5*4914._id + (Cos(t) + -1._id)**3*5733._id + (Cos(t) + -1._id)**4*7166.25_id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==0.and.m2==3) then
   wigd=Cos(t*0.5_id)**3*-1.6183471874253741377307144112820922936791121082957_id*(84._id + (Cos(t) + -1._id)**6*290.0625&
&_id + (Cos(t) + -1._id)*819._id + (Cos(t) + -1._id)**5*1740.375_id + (Cos(t) + -1._id)**2*2866.5_id + (Cos(t) + -1._id)**4*4095._id + (Cos(t) + -1._id)**3*4777.5_id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==0.and.m2==4) then
   wigd=Cos(t*0.5_id)**4*2.3821425596725261067220435421431184595895000743328_id*(126._id + (Cos(t) + -1._id)**5*267.75_i&
&d + (Cos(t) + -1._id)*882._id + (Cos(t) + -1._id)**4*1338.75_id + (Cos(t) + -1._id)**2*2205._id + (Cos(t) + -1._id)**3*2520._id)*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==0.and.m2==5) then
   wigd=Cos(t*0.5_id)**5*-3.9860869143671326737099469187318722046636954283791_id*(126._id + (Cos(t) + -1._id)**4*191.25_&
&id + (Cos(t) + -1._id)*630._id + (Cos(t) + -1._id)**3*765._id + (Cos(t) + -1._id)**2*1080._id)*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==0.and.m2==6) then
   wigd=Cos(t*0.5_id)**6*7.719024117939607353441468160313305396623552881689_id*(84._id + (Cos(t) + -1._id)**3*102._id + &
&(Cos(t) + -1._id)*288._id + (Cos(t) + -1._id)**2*306._id)*Sin(t*0.5_id)**6
endif

if(l==9.and.m1==0.and.m2==7) then
   wigd=Cos(t*0.5_id)**7*-17.826322609494583523570662338729309312832143864472_id*(36._id + (Cos(t) + -1._id)**2*38.25_id&
& + (Cos(t) + -1._id)*76.5_id)*Sin(t*0.5_id)**7
endif

if(l==9.and.m1==0.and.m2==8) then
   wigd=Cos(t)*Cos(t*0.5_id)**8*467.74993319080228383869644347115012625778024243715_id*Sin(t*0.5_id)**8
endif

if(l==9.and.m1==0.and.m2==9) then
   wigd=Cos(t*0.5_id)**9*-220.49943310584723587202231639809357727806572375952_id*Sin(t*0.5_id)**9
endif

if(l==9.and.m1==1.and.m2==-9) then
   wigd=Cos(t*0.5_id)**8*209.18412941712380298700321806359345076313674734586_id*Sin(t*0.5_id)**10
endif

if(l==9.and.m1==1.and.m2==-8) then
   wigd=Cos(t*0.5_id)**7*(2._id + Cos(t)*18._id)*24.652586071242100015888052988157175076326378239492_id*Sin(t*0.5_id)**9&
&
endif

if(l==9.and.m1==1.and.m2==-7) then
   wigd=Cos(t*0.5_id)**6*16.9115345252877628981725179227765953746630798716105_id*((Cos(t) + -1._id)**2*38.25_id + 45._id&
& + (Cos(t) + -1._id)*85._id)*Sin(t*0.5_id)**8
endif

if(l==9.and.m1==1.and.m2==-6) then
   wigd=Cos(t*0.5_id)**5*7.322909257938404906286509527879665713489932169733_id*((Cos(t) + -1._id)**3*102._id + 120._id +&
& (Cos(t) + -1._id)**2*340._id + (Cos(t) + -1._id)*360._id)*Sin(t*0.5_id)**7
endif

if(l==9.and.m1==1.and.m2==-5) then
   wigd=Cos(t*0.5_id)**4*3.7815340802378074032779109105564661534285535943372_id*((Cos(t) + -1._id)**4*191.25_id + 210._i&
&d + (Cos(t) + -1._id)**3*850._id + (Cos(t) + -1._id)*900._id + (Cos(t) + -1._id)**2*1350._id)*Sin(t*0.5_id)**6
endif

if(l==9.and.m1==1.and.m2==-4) then
   wigd=Cos(t*0.5_id)**3*2.25989885993662493879330654177149879193428953707_id*(252._id + (Cos(t) + -1._id)**5*267.75_id &
&+ (Cos(t) + -1._id)*1470._id + (Cos(t) + -1._id)**4*1487.5_id + (Cos(t) + -1._id)**2*3150._id + (Cos(t) + -1._id)**3*3150._id)*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==1.and.m2==-3) then
   wigd=Cos(t*0.5_id)**2*1.53529894715747693137457218396297949971170131913177_id*(210._id + (Cos(t) + -1._id)**6*290.062&
&5_id + (Cos(t) + -1._id)*1638._id + (Cos(t) + -1._id)**5*1933.75_id + (Cos(t) + -1._id)**2*4777.5_id + (Cos(t) + -1._id)**4*5118.75_id + (Cos(t) + -1._id)**3*6825._id)*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==1.and.m2==-2) then
   wigd=Cos(t*0.5_id)*1.1726039399558573886414075283861165701470570883529_id*(120._id + (Cos(t) + -1._id)**7*248.625_id &
&+ (Cos(t) + -1._id)*1260._id + (Cos(t) + -1._id)**6*1933.75_id + (Cos(t) + -1._id)**2*4914._id + (Cos(t) + -1._id)**5*6142.5_id + (Cos(t) + -1._id)**3*9555._id + (Cos(t) + -1._id)**4*10237.5_id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==1.and.m2==-1) then
   wigd=(45._id + (Cos(t) + -1._id)**8*170.9296875_id + (Cos(t) + -1._id)*660._id + (Cos(t) + -1._id)**7*1519.375_id + (&
&Cos(t) + -1._id)**2*3465._id + (Cos(t) + -1._id)**6*5630.625_id + (Cos(t) + -1._id)**3*9009._id + (Cos(t) + -1._id)**5*11261.25_id + (Cos(t) + -1._id)**4*13138.125_id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==1.and.m2==0) then
   wigd=Cos(t*0.5_id)*1.0540925533894597773329645148109061779065183797751_id*(9._id + (Cos(t) + -1._id)**8*170.9296875_i&
&d + (Cos(t) + -1._id)*198._id + (Cos(t) + -1._id)**7*1367.4375_id + (Cos(t) + -1._id)**2*1386._id + (Cos(t) + -1._id)**3*4504.5_id + (Cos(t) + -1._id)**6*4504.5_id + (Cos(t) + -1._id)**4*7882.875_id + (Cos(t) + -1._id)**5*7882.875_id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==1.and.m2==1) then
   wigd=Cos(t*0.5_id)**2*(1._id + (Cos(t) + -1._id)*44._id + (Cos(t) + -1._id)**8*170.9296875_id + (Cos(t) + -1._id)**2*&
&462._id + (Cos(t) + -1._id)**7*1215.5_id + (Cos(t) + -1._id)**3*2002._id + (Cos(t) + -1._id)**6*3503.5_id + (Cos(t) + -1._id)**4*4379.375_id + (Cos(t) + -1._id)**5*5255.25_id)
endif

if(l==9.and.m1==1.and.m2==2) then
   wigd=Cos(t*0.5_id)**3*-1.1726039399558573886414075283861165701470570883529_id*(8._id + (Cos(t) + -1._id)*168._id + (C&
&os(t) + -1._id)**7*248.625_id + (Cos(t) + -1._id)**2*1092._id + (Cos(t) + -1._id)**6*1547._id + (Cos(t) + -1._id)**3*3185._id + (Cos(t) + -1._id)**5*3822._id + (Cos(t) + -1._id)**4*4777.5_id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==1.and.m2==3) then
   wigd=Cos(t*0.5_id)**4*1.53529894715747693137457218396297949971170131913177_id*(28._id + (Cos(t) + -1._id)**6*290.0625&
&_id + (Cos(t) + -1._id)*364._id + (Cos(t) + -1._id)**5*1547._id + (Cos(t) + -1._id)**2*1592.5_id + (Cos(t) + -1._id)**3*3185._id + (Cos(t) + -1._id)**4*3185._id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==1.and.m2==4) then
   wigd=Cos(t*0.5_id)**5*-2.25989885993662493879330654177149879193428953707_id*(56._id + (Cos(t) + -1._id)**5*267.75_id &
&+ (Cos(t) + -1._id)*490._id + (Cos(t) + -1._id)**4*1190._id + (Cos(t) + -1._id)**2*1470._id + (Cos(t) + -1._id)**3*1960._id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==1.and.m2==5) then
   wigd=Cos(t*0.5_id)**6*3.7815340802378074032779109105564661534285535943372_id*(70._id + (Cos(t) + -1._id)**4*191.25_id&
& + (Cos(t) + -1._id)*420._id + (Cos(t) + -1._id)**3*680._id + (Cos(t) + -1._id)**2*840._id)*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==1.and.m2==6) then
   wigd=Cos(t*0.5_id)**7*-7.322909257938404906286509527879665713489932169733_id*(56._id + (Cos(t) + -1._id)**3*102._id +&
& (Cos(t) + -1._id)*224._id + (Cos(t) + -1._id)**2*272._id)*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==1.and.m2==7) then
   wigd=Cos(t*0.5_id)**8*16.9115345252877628981725179227765953746630798716105_id*(28._id + (Cos(t) + -1._id)**2*38.25_id&
& + (Cos(t) + -1._id)*68._id)*Sin(t*0.5_id)**6
endif

if(l==9.and.m1==1.and.m2==8) then
   wigd=Cos(t*0.5_id)**9*-24.652586071242100015888052988157175076326378239492_id*(-2._id + Cos(t)*18._id)*Sin(t*0.5_id)*&
&*7
endif

if(l==9.and.m1==1.and.m2==9) then
   wigd=Cos(t*0.5_id)**10*209.18412941712380298700321806359345076313674734586_id*Sin(t*0.5_id)**8
endif

if(l==9.and.m1==2.and.m2==-9) then
   wigd=Cos(t*0.5_id)**7*178.39282496782206627134409856717675981738782616465_id*Sin(t*0.5_id)**11
endif

if(l==9.and.m1==2.and.m2==-8) then
   wigd=Cos(t*0.5_id)**6*(4._id + Cos(t)*18._id)*21.023796041628638288419857057495807938240784503821_id*Sin(t*0.5_id)**1&
&0
endif

if(l==9.and.m1==2.and.m2==-7) then
   wigd=Cos(t*0.5_id)**5*14.422205101855957172476885069881983785005186295381_id*((Cos(t) + -1._id)**2*38.25_id + 55._id &
&+ (Cos(t) + -1._id)*93.5_id)*Sin(t*0.5_id)**9
endif

if(l==9.and.m1==2.and.m2==-6) then
   wigd=Cos(t*0.5_id)**4*6.2449979983983982058468931209397944610729599779917_id*((Cos(t) + -1._id)**3*102._id + 165._id &
&+ (Cos(t) + -1._id)**2*374._id + (Cos(t) + -1._id)*440._id)*Sin(t*0.5_id)**8
endif

if(l==9.and.m1==2.and.m2==-5) then
   wigd=Cos(t*0.5_id)**3*3.2249030993194198609466452921215084524537585222434_id*((Cos(t) + -1._id)**4*191.25_id + 330._i&
&d + (Cos(t) + -1._id)**3*935._id + (Cos(t) + -1._id)*1237.5_id + (Cos(t) + -1._id)**2*1650._id)*Sin(t*0.5_id)**7
endif

if(l==9.and.m1==2.and.m2==-4) then
   wigd=Cos(t*0.5_id)**2*1.92724822331886306650718651592793499606505258574536_id*((Cos(t) + -1._id)**5*267.75_id + 462._&
&id + (Cos(t) + -1._id)**4*1636.25_id + (Cos(t) + -1._id)*2310._id + (Cos(t) + -1._id)**3*3850._id + (Cos(t) + -1._id)**2*4331.25_id)*Sin(t*0.5_id)**6
endif

if(l==9.and.m1==2.and.m2==-3) then
   wigd=Cos(t*0.5_id)*1.3093073414159542875965849124937167111384161647908_id*((Cos(t) + -1._id)**6*290.0625_id + 462._id&
& + (Cos(t) + -1._id)**5*2127.125_id + (Cos(t) + -1._id)*3003._id + (Cos(t) + -1._id)**4*6256.25_id + (Cos(t) + -1._id)**2*7507.5_id + (Cos(t) + -1._id)**3*9384.375_id)*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==2.and.m2==-2) then
   wigd=((Cos(t) + -1._id)**7*248.625_id + 330._id + (Cos(t) + -1._id)**6*2127.125_id + (Cos(t) + -1._id)*2772._id + (Co&
&s(t) + -1._id)**5*7507.5_id + (Cos(t) + -1._id)**2*9009._id + (Cos(t) + -1._id)**4*14076.5625_id + (Cos(t) + -1._id)**3*15015._id)*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==2.and.m2==-1) then
   wigd=Cos(t*0.5_id)*1.1726039399558573886414075283861165701470570883529_id*(120._id + (Cos(t) + -1._id)**7*248.625_id &
&+ (Cos(t) + -1._id)*1260._id + (Cos(t) + -1._id)**6*1933.75_id + (Cos(t) + -1._id)**2*4914._id + (Cos(t) + -1._id)**5*6142.5_id + (Cos(t) + -1._id)**3*9555._id + (Cos(t) + -1._id)**4*10237.5_id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==2.and.m2==0) then
   wigd=Cos(t*0.5_id)**2*1.2360330811826104914518995734667855101633173865163_id*(36._id + (Cos(t) + -1._id)**7*248.625_i&
&d + (Cos(t) + -1._id)*504._id + (Cos(t) + -1._id)**6*1740.375_id + (Cos(t) + -1._id)**2*2457._id + (Cos(t) + -1._id)**5*4914._id + (Cos(t) + -1._id)**3*5733._id + (Cos(t) + -1._id)**4*7166.25_id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==2.and.m2==1) then
   wigd=Cos(t*0.5_id)**3*1.1726039399558573886414075283861165701470570883529_id*(8._id + (Cos(t) + -1._id)*168._id + (Co&
&s(t) + -1._id)**7*248.625_id + (Cos(t) + -1._id)**2*1092._id + (Cos(t) + -1._id)**6*1547._id + (Cos(t) + -1._id)**3*3185._id + (Cos(t) + -1._id)**5*3822._id + (Cos(t) + -1._id)**4*4777.5_id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==2.and.m2==2) then
   wigd=Cos(t*0.5_id)**4*(1._id + (Cos(t) + -1._id)*42._id + (Cos(t) + -1._id)**7*248.625_id + (Cos(t) + -1._id)**2*409.&
&5_id + (Cos(t) + -1._id)**6*1353.625_id + (Cos(t) + -1._id)**3*1592.5_id + (Cos(t) + -1._id)**5*2866.5_id + (Cos(t) + -1._id)**4*2985.9375_id)
endif

if(l==9.and.m1==2.and.m2==3) then
   wigd=Cos(t*0.5_id)**5*-1.3093073414159542875965849124937167111384161647908_id*(7._id + (Cos(t) + -1._id)*136.5_id + (&
&Cos(t) + -1._id)**6*290.0625_id + (Cos(t) + -1._id)**2*796.25_id + (Cos(t) + -1._id)**5*1353.625_id + (Cos(t) + -1._id)**3*1990.625_id + (Cos(t) + -1._id)**4*2388.75_id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==2.and.m2==4) then
   wigd=Cos(t*0.5_id)**6*1.92724822331886306650718651592793499606505258574536_id*(21._id + (Cos(t) + -1._id)*245._id + (&
&Cos(t) + -1._id)**5*267.75_id + (Cos(t) + -1._id)**2*918.75_id + (Cos(t) + -1._id)**4*1041.25_id + (Cos(t) + -1._id)**3*1470._id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==2.and.m2==5) then
   wigd=Cos(t*0.5_id)**7*-3.2249030993194198609466452921215084524537585222434_id*(35._id + (Cos(t) + -1._id)**4*191.25_i&
&d + (Cos(t) + -1._id)*262.5_id + (Cos(t) + -1._id)**3*595._id + (Cos(t) + -1._id)**2*630._id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==2.and.m2==6) then
   wigd=Cos(t*0.5_id)**8*6.2449979983983982058468931209397944610729599779917_id*(35._id + (Cos(t) + -1._id)**3*102._id +&
& (Cos(t) + -1._id)*168._id + (Cos(t) + -1._id)**2*238._id)*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==2.and.m2==7) then
   wigd=Cos(t*0.5_id)**9*-14.422205101855957172476885069881983785005186295381_id*(21._id + (Cos(t) + -1._id)**2*38.25_id&
& + (Cos(t) + -1._id)*59.5_id)*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==2.and.m2==8) then
   wigd=Cos(t*0.5_id)**10*(-4._id + Cos(t)*18._id)*21.023796041628638288419857057495807938240784503821_id*Sin(t*0.5_id)*&
&*6
endif

if(l==9.and.m1==2.and.m2==9) then
   wigd=Cos(t*0.5_id)**11*-178.39282496782206627134409856717675981738782616465_id*Sin(t*0.5_id)**7
endif

if(l==9.and.m1==3.and.m2==-9) then
   wigd=Cos(t*0.5_id)**6*136.24977064200878866740367808992632648136146844452_id*Sin(t*0.5_id)**12
endif

if(l==9.and.m1==3.and.m2==-8) then
   wigd=Cos(t*0.5_id)**5*16.0571894593460325569611127922999082913219155398578_id*(6._id + Cos(t)*18._id)*Sin(t*0.5_id)**&
&11
endif

if(l==9.and.m1==3.and.m2==-7) then
   wigd=Cos(t*0.5_id)**4*11.015141094572204041211617541744146941892253223863_id*((Cos(t) + -1._id)**2*38.25_id + 66._id &
&+ (Cos(t) + -1._id)*102._id)*Sin(t*0.5_id)**10
endif

if(l==9.and.m1==3.and.m2==-6) then
   wigd=Cos(t*0.5_id)**3*4.7696960070847282457631079301161327012731171262527_id*((Cos(t) + -1._id)**3*102._id + 220._id &
&+ (Cos(t) + -1._id)**2*408._id + (Cos(t) + -1._id)*528._id)*Sin(t*0.5_id)**9
endif

if(l==9.and.m1==3.and.m2==-5) then
   wigd=Cos(t*0.5_id)**2*2.4630604269214887994423741066531428037677316148781_id*((Cos(t) + -1._id)**4*191.25_id + 495._i&
&d + (Cos(t) + -1._id)**3*1020._id + (Cos(t) + -1._id)*1650._id + (Cos(t) + -1._id)**2*1980._id)*Sin(t*0.5_id)**8
endif

if(l==9.and.m1==3.and.m2==-4) then
   wigd=Cos(t*0.5_id)*1.4719601443879744757940071211598756606957732468219_id*((Cos(t) + -1._id)**5*267.75_id + 792._id +&
& (Cos(t) + -1._id)**4*1785._id + (Cos(t) + -1._id)*3465._id + (Cos(t) + -1._id)**3*4620._id + (Cos(t) + -1._id)**2*5775._id)*Sin(t*0.5_id)**7
endif

if(l==9.and.m1==3.and.m2==-3) then
   wigd=((Cos(t) + -1._id)**6*290.0625_id + 924._id + (Cos(t) + -1._id)**5*2320.5_id + (Cos(t) + -1._id)*5148._id + (Cos&
&(t) + -1._id)**4*7507.5_id + (Cos(t) + -1._id)**2*11261.25_id + (Cos(t) + -1._id)**3*12512.5_id)*Sin(t*0.5_id)**6
endif

if(l==9.and.m1==3.and.m2==-2) then
   wigd=Cos(t*0.5_id)*1.3093073414159542875965849124937167111384161647908_id*((Cos(t) + -1._id)**6*290.0625_id + 462._id&
& + (Cos(t) + -1._id)**5*2127.125_id + (Cos(t) + -1._id)*3003._id + (Cos(t) + -1._id)**4*6256.25_id + (Cos(t) + -1._id)**2*7507.5_id + (Cos(t) + -1._id)**3*9384.375_id)*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==3.and.m2==-1) then
   wigd=Cos(t*0.5_id)**2*1.53529894715747693137457218396297949971170131913177_id*(210._id + (Cos(t) + -1._id)**6*290.062&
&5_id + (Cos(t) + -1._id)*1638._id + (Cos(t) + -1._id)**5*1933.75_id + (Cos(t) + -1._id)**2*4777.5_id + (Cos(t) + -1._id)**4*5118.75_id + (Cos(t) + -1._id)**3*6825._id)*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==3.and.m2==0) then
   wigd=Cos(t*0.5_id)**3*1.61834718742537413773071441128209229367911210829573_id*(84._id + (Cos(t) + -1._id)**6*290.0625&
&_id + (Cos(t) + -1._id)*819._id + (Cos(t) + -1._id)**5*1740.375_id + (Cos(t) + -1._id)**2*2866.5_id + (Cos(t) + -1._id)**4*4095._id + (Cos(t) + -1._id)**3*4777.5_id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==3.and.m2==1) then
   wigd=Cos(t*0.5_id)**4*1.53529894715747693137457218396297949971170131913177_id*(28._id + (Cos(t) + -1._id)**6*290.0625&
&_id + (Cos(t) + -1._id)*364._id + (Cos(t) + -1._id)**5*1547._id + (Cos(t) + -1._id)**2*1592.5_id + (Cos(t) + -1._id)**3*3185._id + (Cos(t) + -1._id)**4*3185._id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==3.and.m2==2) then
   wigd=Cos(t*0.5_id)**5*1.3093073414159542875965849124937167111384161647908_id*(7._id + (Cos(t) + -1._id)*136.5_id + (C&
&os(t) + -1._id)**6*290.0625_id + (Cos(t) + -1._id)**2*796.25_id + (Cos(t) + -1._id)**5*1353.625_id + (Cos(t) + -1._id)**3*1990.625_id + (Cos(t) + -1._id)**4*2388.75_id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==3.and.m2==3) then
   wigd=Cos(t*0.5_id)**6*(1._id + (Cos(t) + -1._id)*39._id + (Cos(t) + -1._id)**6*290.0625_id + (Cos(t) + -1._id)**2*341&
&.25_id + (Cos(t) + -1._id)**3*1137.5_id + (Cos(t) + -1._id)**5*1160.25_id + (Cos(t) + -1._id)**4*1706.25_id)
endif

if(l==9.and.m1==3.and.m2==4) then
   wigd=Cos(t*0.5_id)**7*-1.4719601443879744757940071211598756606957732468219_id*(6._id + (Cos(t) + -1._id)*105._id + (C&
&os(t) + -1._id)**5*267.75_id + (Cos(t) + -1._id)**2*525._id + (Cos(t) + -1._id)**4*892.5_id + (Cos(t) + -1._id)**3*1050._id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==3.and.m2==5) then
   wigd=Cos(t*0.5_id)**8*2.4630604269214887994423741066531428037677316148781_id*(15._id + (Cos(t) + -1._id)*150._id + (C&
&os(t) + -1._id)**4*191.25_id + (Cos(t) + -1._id)**2*450._id + (Cos(t) + -1._id)**3*510._id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==3.and.m2==6) then
   wigd=Cos(t*0.5_id)**9*-4.7696960070847282457631079301161327012731171262527_id*(20._id + (Cos(t) + -1._id)**3*102._id &
&+ (Cos(t) + -1._id)*120._id + (Cos(t) + -1._id)**2*204._id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==3.and.m2==7) then
   wigd=Cos(t*0.5_id)**10*11.015141094572204041211617541744146941892253223863_id*(15._id + (Cos(t) + -1._id)**2*38.25_id&
& + (Cos(t) + -1._id)*51._id)*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==3.and.m2==8) then
   wigd=Cos(t*0.5_id)**11*-16.057189459346032556961112792299908291321915539858_id*(-6._id + Cos(t)*18._id)*Sin(t*0.5_id)&
&**5
endif

if(l==9.and.m1==3.and.m2==9) then
   wigd=Cos(t*0.5_id)**12*136.24977064200878866740367808992632648136146844452_id*Sin(t*0.5_id)**6
endif

if(l==9.and.m1==4.and.m2==-9) then
   wigd=Cos(t*0.5_id)**5*92.56349172324907493447268322253272345128247145089_id*Sin(t*0.5_id)**13
endif

if(l==9.and.m1==4.and.m2==-8) then
   wigd=Cos(t*0.5_id)**4*10.9087121146357144115021544873729018728051337026946_id*(8._id + Cos(t)*18._id)*Sin(t*0.5_id)**&
&12
endif

if(l==9.and.m1==4.and.m2==-7) then
   wigd=Cos(t*0.5_id)**3*7.483314773547882771167497464633098603512039615557_id*((Cos(t) + -1._id)**2*38.25_id + 78._id +&
& (Cos(t) + -1._id)*110.5_id)*Sin(t*0.5_id)**11
endif

if(l==9.and.m1==4.and.m2==-6) then
   wigd=Cos(t*0.5_id)**2*3.2403703492039301154829837180439983288526021535292_id*((Cos(t) + -1._id)**3*102._id + 286._id &
&+ (Cos(t) + -1._id)**2*442._id + (Cos(t) + -1._id)*624._id)*Sin(t*0.5_id)**10
endif

if(l==9.and.m1==4.and.m2==-5) then
   wigd=Cos(t*0.5_id)*1.67332005306815109595634405157037497878563073859734_id*((Cos(t) + -1._id)**4*191.25_id + 715._id &
&+ (Cos(t) + -1._id)**3*1105._id + (Cos(t) + -1._id)*2145._id + (Cos(t) + -1._id)**2*2340._id)*Sin(t*0.5_id)**9
endif

if(l==9.and.m1==4.and.m2==-4) then
   wigd=((Cos(t) + -1._id)**5*267.75_id + 1287._id + (Cos(t) + -1._id)**4*1933.75_id + (Cos(t) + -1._id)*5005._id + (Cos&
&(t) + -1._id)**3*5460._id + (Cos(t) + -1._id)**2*7507.5_id)*Sin(t*0.5_id)**8
endif

if(l==9.and.m1==4.and.m2==-3) then
   wigd=Cos(t*0.5_id)*1.4719601443879744757940071211598756606957732468219_id*((Cos(t) + -1._id)**5*267.75_id + 792._id +&
& (Cos(t) + -1._id)**4*1785._id + (Cos(t) + -1._id)*3465._id + (Cos(t) + -1._id)**3*4620._id + (Cos(t) + -1._id)**2*5775._id)*Sin(t*0.5_id)**7
endif

if(l==9.and.m1==4.and.m2==-2) then
   wigd=Cos(t*0.5_id)**2*1.92724822331886306650718651592793499606505258574536_id*((Cos(t) + -1._id)**5*267.75_id + 462._&
&id + (Cos(t) + -1._id)**4*1636.25_id + (Cos(t) + -1._id)*2310._id + (Cos(t) + -1._id)**3*3850._id + (Cos(t) + -1._id)**2*4331.25_id)*Sin(t*0.5_id)**6
endif

if(l==9.and.m1==4.and.m2==-1) then
   wigd=Cos(t*0.5_id)**3*2.25989885993662493879330654177149879193428953707_id*(252._id + (Cos(t) + -1._id)**5*267.75_id &
&+ (Cos(t) + -1._id)*1470._id + (Cos(t) + -1._id)**4*1487.5_id + (Cos(t) + -1._id)**2*3150._id + (Cos(t) + -1._id)**3*3150._id)*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==4.and.m2==0) then
   wigd=Cos(t*0.5_id)**4*2.3821425596725261067220435421431184595895000743328_id*(126._id + (Cos(t) + -1._id)**5*267.75_i&
&d + (Cos(t) + -1._id)*882._id + (Cos(t) + -1._id)**4*1338.75_id + (Cos(t) + -1._id)**2*2205._id + (Cos(t) + -1._id)**3*2520._id)*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==4.and.m2==1) then
   wigd=Cos(t*0.5_id)**5*2.25989885993662493879330654177149879193428953707_id*(56._id + (Cos(t) + -1._id)**5*267.75_id +&
& (Cos(t) + -1._id)*490._id + (Cos(t) + -1._id)**4*1190._id + (Cos(t) + -1._id)**2*1470._id + (Cos(t) + -1._id)**3*1960._id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==4.and.m2==2) then
   wigd=Cos(t*0.5_id)**6*1.92724822331886306650718651592793499606505258574536_id*(21._id + (Cos(t) + -1._id)*245._id + (&
&Cos(t) + -1._id)**5*267.75_id + (Cos(t) + -1._id)**2*918.75_id + (Cos(t) + -1._id)**4*1041.25_id + (Cos(t) + -1._id)**3*1470._id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==4.and.m2==3) then
   wigd=Cos(t*0.5_id)**7*1.4719601443879744757940071211598756606957732468219_id*(6._id + (Cos(t) + -1._id)*105._id + (Co&
&s(t) + -1._id)**5*267.75_id + (Cos(t) + -1._id)**2*525._id + (Cos(t) + -1._id)**4*892.5_id + (Cos(t) + -1._id)**3*1050._id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==4.and.m2==4) then
   wigd=Cos(t*0.5_id)**8*(1._id + (Cos(t) + -1._id)*35._id + (Cos(t) + -1._id)**2*262.5_id + (Cos(t) + -1._id)**5*267.75&
&_id + (Cos(t) + -1._id)**3*700._id + (Cos(t) + -1._id)**4*743.75_id)
endif

if(l==9.and.m1==4.and.m2==5) then
   wigd=Cos(t*0.5_id)**9*-1.6733200530681510959563440515703749787856307385973_id*(5._id + (Cos(t) + -1._id)*75._id + (Co&
&s(t) + -1._id)**4*191.25_id + (Cos(t) + -1._id)**2*300._id + (Cos(t) + -1._id)**3*425._id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==4.and.m2==6) then
   wigd=Cos(t*0.5_id)**10*3.2403703492039301154829837180439983288526021535292_id*(10._id + (Cos(t) + -1._id)*80._id + (C&
&os(t) + -1._id)**3*102._id + (Cos(t) + -1._id)**2*170._id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==4.and.m2==7) then
   wigd=Cos(t*0.5_id)**11*-7.483314773547882771167497464633098603512039615557_id*(10._id + (Cos(t) + -1._id)**2*38.25_id&
& + (Cos(t) + -1._id)*42.5_id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==4.and.m2==8) then
   wigd=Cos(t*0.5_id)**12*10.9087121146357144115021544873729018728051337026946_id*(-8._id + Cos(t)*18._id)*Sin(t*0.5_id)&
&**4
endif

if(l==9.and.m1==4.and.m2==9) then
   wigd=Cos(t*0.5_id)**13*-92.56349172324907493447268322253272345128247145089_id*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==5.and.m2==-9) then
   wigd=Cos(t*0.5_id)**4*55.317266743757323860013645690576758943480830292335_id*Sin(t*0.5_id)**14
endif

if(l==9.and.m1==5.and.m2==-8) then
   wigd=Cos(t*0.5_id)**3*6.5192024052026487145829715574291844165280937789101_id*(10._id + Cos(t)*18._id)*Sin(t*0.5_id)**&
&13
endif

if(l==9.and.m1==5.and.m2==-7) then
   wigd=Cos(t*0.5_id)**2*4.4721359549995793928183473374625524708812367192231_id*((Cos(t) + -1._id)**2*38.25_id + 91._id &
&+ (Cos(t) + -1._id)*119._id)*Sin(t*0.5_id)**12
endif

if(l==9.and.m1==5.and.m2==-6) then
   wigd=Cos(t*0.5_id)*1.9364916731037084425896326998911998054164608526458_id*((Cos(t) + -1._id)**3*102._id + 364._id + (&
&Cos(t) + -1._id)**2*476._id + (Cos(t) + -1._id)*728._id)*Sin(t*0.5_id)**11
endif

if(l==9.and.m1==5.and.m2==-5) then
   wigd=((Cos(t) + -1._id)**4*191.25_id + 1001._id + (Cos(t) + -1._id)**3*1190._id + (Cos(t) + -1._id)*2730._id + (Cos(t&
&) + -1._id)**2*2730._id)*Sin(t*0.5_id)**10
endif

if(l==9.and.m1==5.and.m2==-4) then
   wigd=Cos(t*0.5_id)*1.67332005306815109595634405157037497878563073859734_id*((Cos(t) + -1._id)**4*191.25_id + 715._id &
&+ (Cos(t) + -1._id)**3*1105._id + (Cos(t) + -1._id)*2145._id + (Cos(t) + -1._id)**2*2340._id)*Sin(t*0.5_id)**9
endif

if(l==9.and.m1==5.and.m2==-3) then
   wigd=Cos(t*0.5_id)**2*2.4630604269214887994423741066531428037677316148781_id*((Cos(t) + -1._id)**4*191.25_id + 495._i&
&d + (Cos(t) + -1._id)**3*1020._id + (Cos(t) + -1._id)*1650._id + (Cos(t) + -1._id)**2*1980._id)*Sin(t*0.5_id)**8
endif

if(l==9.and.m1==5.and.m2==-2) then
   wigd=Cos(t*0.5_id)**3*3.2249030993194198609466452921215084524537585222434_id*((Cos(t) + -1._id)**4*191.25_id + 330._i&
&d + (Cos(t) + -1._id)**3*935._id + (Cos(t) + -1._id)*1237.5_id + (Cos(t) + -1._id)**2*1650._id)*Sin(t*0.5_id)**7
endif

if(l==9.and.m1==5.and.m2==-1) then
   wigd=Cos(t*0.5_id)**4*3.7815340802378074032779109105564661534285535943372_id*((Cos(t) + -1._id)**4*191.25_id + 210._i&
&d + (Cos(t) + -1._id)**3*850._id + (Cos(t) + -1._id)*900._id + (Cos(t) + -1._id)**2*1350._id)*Sin(t*0.5_id)**6
endif

if(l==9.and.m1==5.and.m2==0) then
   wigd=Cos(t*0.5_id)**5*3.9860869143671326737099469187318722046636954283791_id*(126._id + (Cos(t) + -1._id)**4*191.25_i&
&d + (Cos(t) + -1._id)*630._id + (Cos(t) + -1._id)**3*765._id + (Cos(t) + -1._id)**2*1080._id)*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==5.and.m2==1) then
   wigd=Cos(t*0.5_id)**6*3.7815340802378074032779109105564661534285535943372_id*(70._id + (Cos(t) + -1._id)**4*191.25_id&
& + (Cos(t) + -1._id)*420._id + (Cos(t) + -1._id)**3*680._id + (Cos(t) + -1._id)**2*840._id)*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==5.and.m2==2) then
   wigd=Cos(t*0.5_id)**7*3.2249030993194198609466452921215084524537585222434_id*(35._id + (Cos(t) + -1._id)**4*191.25_id&
& + (Cos(t) + -1._id)*262.5_id + (Cos(t) + -1._id)**3*595._id + (Cos(t) + -1._id)**2*630._id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==5.and.m2==3) then
   wigd=Cos(t*0.5_id)**8*2.4630604269214887994423741066531428037677316148781_id*(15._id + (Cos(t) + -1._id)*150._id + (C&
&os(t) + -1._id)**4*191.25_id + (Cos(t) + -1._id)**2*450._id + (Cos(t) + -1._id)**3*510._id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==5.and.m2==4) then
   wigd=Cos(t*0.5_id)**9*1.67332005306815109595634405157037497878563073859734_id*(5._id + (Cos(t) + -1._id)*75._id + (Co&
&s(t) + -1._id)**4*191.25_id + (Cos(t) + -1._id)**2*300._id + (Cos(t) + -1._id)**3*425._id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==5.and.m2==5) then
   wigd=Cos(t*0.5_id)**10*(1._id + (Cos(t) + -1._id)*30._id + (Cos(t) + -1._id)**2*180._id + (Cos(t) + -1._id)**4*191.25&
&_id + (Cos(t) + -1._id)**3*340._id)
endif

if(l==9.and.m1==5.and.m2==6) then
   wigd=Cos(t*0.5_id)**11*-1.9364916731037084425896326998911998054164608526458_id*(4._id + (Cos(t) + -1._id)*48._id + (C&
&os(t) + -1._id)**3*102._id + (Cos(t) + -1._id)**2*136._id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==5.and.m2==7) then
   wigd=Cos(t*0.5_id)**12*4.4721359549995793928183473374625524708812367192231_id*(6._id + (Cos(t) + -1._id)*34._id + (Co&
&s(t) + -1._id)**2*38.25_id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==5.and.m2==8) then
   wigd=Cos(t*0.5_id)**13*-6.5192024052026487145829715574291844165280937789101_id*(-10._id + Cos(t)*18._id)*Sin(t*0.5_id&
&)**3
endif

if(l==9.and.m1==5.and.m2==9) then
   wigd=Cos(t*0.5_id)**14*55.317266743757323860013645690576758943480830292335_id*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==6.and.m2==-9) then
   wigd=Cos(t*0.5_id)**3*28.565713714171399991997599245469061115064684639611_id*Sin(t*0.5_id)**15
endif

if(l==9.and.m1==6.and.m2==-8) then
   wigd=Cos(t*0.5_id)**2*3.3665016461206926511211286390232002368679663214932_id*(12._id + Cos(t)*18._id)*Sin(t*0.5_id)**&
&14
endif

if(l==9.and.m1==6.and.m2==-7) then
   wigd=Cos(t*0.5_id)*2.3094010767585030580365951220078298225904070050805_id*((Cos(t) + -1._id)**2*38.25_id + 105._id + &
&(Cos(t) + -1._id)*127.5_id)*Sin(t*0.5_id)**13
endif

if(l==9.and.m1==6.and.m2==-6) then
   wigd=((Cos(t) + -1._id)**3*102._id + 455._id + (Cos(t) + -1._id)**2*510._id + (Cos(t) + -1._id)*840._id)*Sin(t*0.5_id&
&)**12
endif

if(l==9.and.m1==6.and.m2==-5) then
   wigd=Cos(t*0.5_id)*1.9364916731037084425896326998911998054164608526458_id*((Cos(t) + -1._id)**3*102._id + 364._id + (&
&Cos(t) + -1._id)**2*476._id + (Cos(t) + -1._id)*728._id)*Sin(t*0.5_id)**11
endif

if(l==9.and.m1==6.and.m2==-4) then
   wigd=Cos(t*0.5_id)**2*3.2403703492039301154829837180439983288526021535292_id*((Cos(t) + -1._id)**3*102._id + 286._id &
&+ (Cos(t) + -1._id)**2*442._id + (Cos(t) + -1._id)*624._id)*Sin(t*0.5_id)**10
endif

if(l==9.and.m1==6.and.m2==-3) then
   wigd=Cos(t*0.5_id)**3*4.7696960070847282457631079301161327012731171262527_id*((Cos(t) + -1._id)**3*102._id + 220._id &
&+ (Cos(t) + -1._id)**2*408._id + (Cos(t) + -1._id)*528._id)*Sin(t*0.5_id)**9
endif

if(l==9.and.m1==6.and.m2==-2) then
   wigd=Cos(t*0.5_id)**4*6.2449979983983982058468931209397944610729599779917_id*((Cos(t) + -1._id)**3*102._id + 165._id &
&+ (Cos(t) + -1._id)**2*374._id + (Cos(t) + -1._id)*440._id)*Sin(t*0.5_id)**8
endif

if(l==9.and.m1==6.and.m2==-1) then
   wigd=Cos(t*0.5_id)**5*7.322909257938404906286509527879665713489932169733_id*((Cos(t) + -1._id)**3*102._id + 120._id +&
& (Cos(t) + -1._id)**2*340._id + (Cos(t) + -1._id)*360._id)*Sin(t*0.5_id)**7
endif

if(l==9.and.m1==6.and.m2==0) then
   wigd=Cos(t*0.5_id)**6*7.719024117939607353441468160313305396623552881689_id*(84._id + (Cos(t) + -1._id)**3*102._id + &
&(Cos(t) + -1._id)*288._id + (Cos(t) + -1._id)**2*306._id)*Sin(t*0.5_id)**6
endif

if(l==9.and.m1==6.and.m2==1) then
   wigd=Cos(t*0.5_id)**7*7.322909257938404906286509527879665713489932169733_id*(56._id + (Cos(t) + -1._id)**3*102._id + &
&(Cos(t) + -1._id)*224._id + (Cos(t) + -1._id)**2*272._id)*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==6.and.m2==2) then
   wigd=Cos(t*0.5_id)**8*6.2449979983983982058468931209397944610729599779917_id*(35._id + (Cos(t) + -1._id)**3*102._id +&
& (Cos(t) + -1._id)*168._id + (Cos(t) + -1._id)**2*238._id)*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==6.and.m2==3) then
   wigd=Cos(t*0.5_id)**9*4.7696960070847282457631079301161327012731171262527_id*(20._id + (Cos(t) + -1._id)**3*102._id +&
& (Cos(t) + -1._id)*120._id + (Cos(t) + -1._id)**2*204._id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==6.and.m2==4) then
   wigd=Cos(t*0.5_id)**10*3.2403703492039301154829837180439983288526021535292_id*(10._id + (Cos(t) + -1._id)*80._id + (C&
&os(t) + -1._id)**3*102._id + (Cos(t) + -1._id)**2*170._id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==6.and.m2==5) then
   wigd=Cos(t*0.5_id)**11*1.9364916731037084425896326998911998054164608526458_id*(4._id + (Cos(t) + -1._id)*48._id + (Co&
&s(t) + -1._id)**3*102._id + (Cos(t) + -1._id)**2*136._id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==6.and.m2==6) then
   wigd=Cos(t*0.5_id)**12*(1._id + (Cos(t) + -1._id)*24._id + (Cos(t) + -1._id)**2*102._id + (Cos(t) + -1._id)**3*102._i&
&d)
endif

if(l==9.and.m1==6.and.m2==7) then
   wigd=Cos(t*0.5_id)**13*-2.3094010767585030580365951220078298225904070050805_id*(3._id + (Cos(t) + -1._id)*25.5_id + (&
&Cos(t) + -1._id)**2*38.25_id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==6.and.m2==8) then
   wigd=Cos(t*0.5_id)**14*3.3665016461206926511211286390232002368679663214932_id*(-12._id + Cos(t)*18._id)*Sin(t*0.5_id)&
&**2
endif

if(l==9.and.m1==6.and.m2==9) then
   wigd=Cos(t*0.5_id)**15*-28.565713714171399991997599245469061115064684639611_id*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==7.and.m2==-9) then
   wigd=Cos(t*0.5_id)**2*12.369316876852981649464229567922231075441597676121_id*Sin(t*0.5_id)**16
endif

if(l==9.and.m1==7.and.m2==-8) then
   wigd=Cos(t*0.5_id)*1.4577379737113251177185382193863957691303495837215_id*(14._id + Cos(t)*18._id)*Sin(t*0.5_id)**15
endif

if(l==9.and.m1==7.and.m2==-7) then
   wigd=((Cos(t) + -1._id)**2*38.25_id + 120._id + (Cos(t) + -1._id)*136._id)*Sin(t*0.5_id)**14
endif

if(l==9.and.m1==7.and.m2==-6) then
   wigd=Cos(t*0.5_id)*2.3094010767585030580365951220078298225904070050805_id*((Cos(t) + -1._id)**2*38.25_id + 105._id + &
&(Cos(t) + -1._id)*127.5_id)*Sin(t*0.5_id)**13
endif

if(l==9.and.m1==7.and.m2==-5) then
   wigd=Cos(t*0.5_id)**2*4.4721359549995793928183473374625524708812367192231_id*((Cos(t) + -1._id)**2*38.25_id + 91._id &
&+ (Cos(t) + -1._id)*119._id)*Sin(t*0.5_id)**12
endif

if(l==9.and.m1==7.and.m2==-4) then
   wigd=Cos(t*0.5_id)**3*7.483314773547882771167497464633098603512039615557_id*((Cos(t) + -1._id)**2*38.25_id + 78._id +&
& (Cos(t) + -1._id)*110.5_id)*Sin(t*0.5_id)**11
endif

if(l==9.and.m1==7.and.m2==-3) then
   wigd=Cos(t*0.5_id)**4*11.015141094572204041211617541744146941892253223863_id*((Cos(t) + -1._id)**2*38.25_id + 66._id &
&+ (Cos(t) + -1._id)*102._id)*Sin(t*0.5_id)**10
endif

if(l==9.and.m1==7.and.m2==-2) then
   wigd=Cos(t*0.5_id)**5*14.422205101855957172476885069881983785005186295381_id*((Cos(t) + -1._id)**2*38.25_id + 55._id &
&+ (Cos(t) + -1._id)*93.5_id)*Sin(t*0.5_id)**9
endif

if(l==9.and.m1==7.and.m2==-1) then
   wigd=Cos(t*0.5_id)**6*16.9115345252877628981725179227765953746630798716105_id*((Cos(t) + -1._id)**2*38.25_id + 45._id&
& + (Cos(t) + -1._id)*85._id)*Sin(t*0.5_id)**8
endif

if(l==9.and.m1==7.and.m2==0) then
   wigd=Cos(t*0.5_id)**7*17.826322609494583523570662338729309312832143864472_id*(36._id + (Cos(t) + -1._id)**2*38.25_id &
&+ (Cos(t) + -1._id)*76.5_id)*Sin(t*0.5_id)**7
endif

if(l==9.and.m1==7.and.m2==1) then
   wigd=Cos(t*0.5_id)**8*16.9115345252877628981725179227765953746630798716105_id*(28._id + (Cos(t) + -1._id)**2*38.25_id&
& + (Cos(t) + -1._id)*68._id)*Sin(t*0.5_id)**6
endif

if(l==9.and.m1==7.and.m2==2) then
   wigd=Cos(t*0.5_id)**9*14.422205101855957172476885069881983785005186295381_id*(21._id + (Cos(t) + -1._id)**2*38.25_id &
&+ (Cos(t) + -1._id)*59.5_id)*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==7.and.m2==3) then
   wigd=Cos(t*0.5_id)**10*11.015141094572204041211617541744146941892253223863_id*(15._id + (Cos(t) + -1._id)**2*38.25_id&
& + (Cos(t) + -1._id)*51._id)*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==7.and.m2==4) then
   wigd=Cos(t*0.5_id)**11*7.483314773547882771167497464633098603512039615557_id*(10._id + (Cos(t) + -1._id)**2*38.25_id &
&+ (Cos(t) + -1._id)*42.5_id)*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==7.and.m2==5) then
   wigd=Cos(t*0.5_id)**12*4.4721359549995793928183473374625524708812367192231_id*(6._id + (Cos(t) + -1._id)*34._id + (Co&
&s(t) + -1._id)**2*38.25_id)*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==7.and.m2==6) then
   wigd=Cos(t*0.5_id)**13*2.3094010767585030580365951220078298225904070050805_id*(3._id + (Cos(t) + -1._id)*25.5_id + (C&
&os(t) + -1._id)**2*38.25_id)*Sin(t*0.5_id)
endif

if(l==9.and.m1==7.and.m2==7) then
   wigd=Cos(t*0.5_id)**14*(1._id + (Cos(t) + -1._id)*17._id + (Cos(t) + -1._id)**2*38.25_id)
endif

if(l==9.and.m1==7.and.m2==8) then
   wigd=Cos(t*0.5_id)**15*-1.4577379737113251177185382193863957691303495837215_id*(-14._id + Cos(t)*18._id)*Sin(t*0.5_id&
&)
endif

if(l==9.and.m1==7.and.m2==9) then
   wigd=Cos(t*0.5_id)**16*12.369316876852981649464229567922231075441597676121_id*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==8.and.m2==-9) then
   wigd=Cos(t*0.5_id)*4.2426406871192851464050661726290942357090156261308_id*Sin(t*0.5_id)**17
endif

if(l==9.and.m1==8.and.m2==-8) then
   wigd=0.5_id*(16._id + Cos(t)*18._id)*Sin(t*0.5_id)**16
endif

if(l==9.and.m1==8.and.m2==-7) then
   wigd=Cos(t*0.5_id)*1.4577379737113251177185382193863957691303495837215_id*(14._id + Cos(t)*18._id)*Sin(t*0.5_id)**15
endif

if(l==9.and.m1==8.and.m2==-6) then
   wigd=Cos(t*0.5_id)**2*3.3665016461206926511211286390232002368679663214932_id*(12._id + Cos(t)*18._id)*Sin(t*0.5_id)**&
&14
endif

if(l==9.and.m1==8.and.m2==-5) then
   wigd=Cos(t*0.5_id)**3*6.5192024052026487145829715574291844165280937789101_id*(10._id + Cos(t)*18._id)*Sin(t*0.5_id)**&
&13
endif

if(l==9.and.m1==8.and.m2==-4) then
   wigd=Cos(t*0.5_id)**4*10.9087121146357144115021544873729018728051337026946_id*(8._id + Cos(t)*18._id)*Sin(t*0.5_id)**&
&12
endif

if(l==9.and.m1==8.and.m2==-3) then
   wigd=Cos(t*0.5_id)**5*16.0571894593460325569611127922999082913219155398578_id*(6._id + Cos(t)*18._id)*Sin(t*0.5_id)**&
&11
endif

if(l==9.and.m1==8.and.m2==-2) then
   wigd=Cos(t*0.5_id)**6*(4._id + Cos(t)*18._id)*21.023796041628638288419857057495807938240784503821_id*Sin(t*0.5_id)**1&
&0
endif

if(l==9.and.m1==8.and.m2==-1) then
   wigd=Cos(t*0.5_id)**7*(2._id + Cos(t)*18._id)*24.652586071242100015888052988157175076326378239492_id*Sin(t*0.5_id)**9&
&
endif

if(l==9.and.m1==8.and.m2==0) then
   wigd=Cos(t)*Cos(t*0.5_id)**8*467.74993319080228383869644347115012625778024243715_id*Sin(t*0.5_id)**8
endif

if(l==9.and.m1==8.and.m2==1) then
   wigd=Cos(t*0.5_id)**9*(-2._id + Cos(t)*18._id)*24.652586071242100015888052988157175076326378239492_id*Sin(t*0.5_id)**&
&7
endif

if(l==9.and.m1==8.and.m2==2) then
   wigd=Cos(t*0.5_id)**10*(-4._id + Cos(t)*18._id)*21.023796041628638288419857057495807938240784503821_id*Sin(t*0.5_id)*&
&*6
endif

if(l==9.and.m1==8.and.m2==3) then
   wigd=Cos(t*0.5_id)**11*16.0571894593460325569611127922999082913219155398578_id*(-6._id + Cos(t)*18._id)*Sin(t*0.5_id)&
&**5
endif

if(l==9.and.m1==8.and.m2==4) then
   wigd=Cos(t*0.5_id)**12*10.9087121146357144115021544873729018728051337026946_id*(-8._id + Cos(t)*18._id)*Sin(t*0.5_id)&
&**4
endif

if(l==9.and.m1==8.and.m2==5) then
   wigd=Cos(t*0.5_id)**13*6.5192024052026487145829715574291844165280937789101_id*(-10._id + Cos(t)*18._id)*Sin(t*0.5_id)&
&**3
endif

if(l==9.and.m1==8.and.m2==6) then
   wigd=Cos(t*0.5_id)**14*3.3665016461206926511211286390232002368679663214932_id*(-12._id + Cos(t)*18._id)*Sin(t*0.5_id)&
&**2
endif

if(l==9.and.m1==8.and.m2==7) then
   wigd=Cos(t*0.5_id)**15*1.4577379737113251177185382193863957691303495837215_id*(-14._id + Cos(t)*18._id)*Sin(t*0.5_id)&
&
endif

if(l==9.and.m1==8.and.m2==8) then
   wigd=Cos(t*0.5_id)**16*0.5_id*(-16._id + Cos(t)*18._id)
endif

if(l==9.and.m1==8.and.m2==9) then
   wigd=Cos(t*0.5_id)**17*-4.2426406871192851464050661726290942357090156261308_id*Sin(t*0.5_id)
endif

if(l==9.and.m1==9.and.m2==-9) then
   wigd=Sin(t*0.5_id)**18
endif

if(l==9.and.m1==9.and.m2==-8) then
   wigd=Cos(t*0.5_id)*4.2426406871192851464050661726290942357090156261308_id*Sin(t*0.5_id)**17
endif

if(l==9.and.m1==9.and.m2==-7) then
   wigd=Cos(t*0.5_id)**2*12.369316876852981649464229567922231075441597676121_id*Sin(t*0.5_id)**16
endif

if(l==9.and.m1==9.and.m2==-6) then
   wigd=Cos(t*0.5_id)**3*28.565713714171399991997599245469061115064684639611_id*Sin(t*0.5_id)**15
endif

if(l==9.and.m1==9.and.m2==-5) then
   wigd=Cos(t*0.5_id)**4*55.317266743757323860013645690576758943480830292335_id*Sin(t*0.5_id)**14
endif

if(l==9.and.m1==9.and.m2==-4) then
   wigd=Cos(t*0.5_id)**5*92.56349172324907493447268322253272345128247145089_id*Sin(t*0.5_id)**13
endif

if(l==9.and.m1==9.and.m2==-3) then
   wigd=Cos(t*0.5_id)**6*136.24977064200878866740367808992632648136146844452_id*Sin(t*0.5_id)**12
endif

if(l==9.and.m1==9.and.m2==-2) then
   wigd=Cos(t*0.5_id)**7*178.39282496782206627134409856717675981738782616465_id*Sin(t*0.5_id)**11
endif

if(l==9.and.m1==9.and.m2==-1) then
   wigd=Cos(t*0.5_id)**8*209.18412941712380298700321806359345076313674734586_id*Sin(t*0.5_id)**10
endif

if(l==9.and.m1==9.and.m2==0) then
   wigd=Cos(t*0.5_id)**9*220.49943310584723587202231639809357727806572375952_id*Sin(t*0.5_id)**9
endif

if(l==9.and.m1==9.and.m2==1) then
   wigd=Cos(t*0.5_id)**10*209.18412941712380298700321806359345076313674734586_id*Sin(t*0.5_id)**8
endif

if(l==9.and.m1==9.and.m2==2) then
   wigd=Cos(t*0.5_id)**11*178.39282496782206627134409856717675981738782616465_id*Sin(t*0.5_id)**7
endif

if(l==9.and.m1==9.and.m2==3) then
   wigd=Cos(t*0.5_id)**12*136.24977064200878866740367808992632648136146844452_id*Sin(t*0.5_id)**6
endif

if(l==9.and.m1==9.and.m2==4) then
   wigd=Cos(t*0.5_id)**13*92.56349172324907493447268322253272345128247145089_id*Sin(t*0.5_id)**5
endif

if(l==9.and.m1==9.and.m2==5) then
   wigd=Cos(t*0.5_id)**14*55.317266743757323860013645690576758943480830292335_id*Sin(t*0.5_id)**4
endif

if(l==9.and.m1==9.and.m2==6) then
   wigd=Cos(t*0.5_id)**15*28.565713714171399991997599245469061115064684639611_id*Sin(t*0.5_id)**3
endif

if(l==9.and.m1==9.and.m2==7) then
   wigd=Cos(t*0.5_id)**16*12.369316876852981649464229567922231075441597676121_id*Sin(t*0.5_id)**2
endif

if(l==9.and.m1==9.and.m2==8) then
   wigd=Cos(t*0.5_id)**17*4.2426406871192851464050661726290942357090156261308_id*Sin(t*0.5_id)
endif

if(l==9.and.m1==9.and.m2==9) then
   wigd=Cos(t*0.5_id)**18
endif

if(l==10.and.m1==-10.and.m2==-10) then
   wigd=Cos(t*0.5_id)**20
endif

if(l==10.and.m1==-10.and.m2==-9) then
   wigd=Cos(t*0.5_id)**19*-4.4721359549995793928183473374625524708812367192231_id*Sin(t*0.5_id)
endif

if(l==10.and.m1==-10.and.m2==-8) then
   wigd=Cos(t*0.5_id)**18*13.784048752090221767955912552934175427198163558399_id*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==-10.and.m2==-7) then
   wigd=Cos(t*0.5_id)**17*-33.76388603226826436623377881904421997769695431525_id*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==-10.and.m2==-6) then
   wigd=Cos(t*0.5_id)**16*69.606034221179416364101895080211399709692904715228_id*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==-10.and.m2==-5) then
   wigd=Cos(t*0.5_id)**15*-124.51505933018704545197426688644896998829430818083_id*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==-10.and.m2==-4) then
   wigd=Cos(t*0.5_id)**14*196.87559523719540998400115941877306606266144591633_id*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==-10.and.m2==-3) then
   wigd=Cos(t*0.5_id)**13*-278.42413688471766545640758032084559883877161886091_id*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==-10.and.m2==-2) then
   wigd=Cos(t*0.5_id)**12*354.92252675760100307924766676167204916378679569317_id*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==-10.and.m2==-1) then
   wigd=Cos(t*0.5_id)**11*-409.82923272992618480080474681969664551446288789265_id*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==-10.and.m2==0) then
   wigd=Cos(t*0.5_id)**10*429.83252552593085495728450959347453264579252131871_id*Sin(t*0.5_id)**10
endif

if(l==10.and.m1==-10.and.m2==1) then
   wigd=Cos(t*0.5_id)**9*-409.82923272992618480080474681969664551446288789265_id*Sin(t*0.5_id)**11
endif

if(l==10.and.m1==-10.and.m2==2) then
   wigd=Cos(t*0.5_id)**8*354.92252675760100307924766676167204916378679569317_id*Sin(t*0.5_id)**12
endif

if(l==10.and.m1==-10.and.m2==3) then
   wigd=Cos(t*0.5_id)**7*-278.42413688471766545640758032084559883877161886091_id*Sin(t*0.5_id)**13
endif

if(l==10.and.m1==-10.and.m2==4) then
   wigd=Cos(t*0.5_id)**6*196.87559523719540998400115941877306606266144591633_id*Sin(t*0.5_id)**14
endif

if(l==10.and.m1==-10.and.m2==5) then
   wigd=Cos(t*0.5_id)**5*-124.51505933018704545197426688644896998829430818083_id*Sin(t*0.5_id)**15
endif

if(l==10.and.m1==-10.and.m2==6) then
   wigd=Cos(t*0.5_id)**4*69.606034221179416364101895080211399709692904715228_id*Sin(t*0.5_id)**16
endif

if(l==10.and.m1==-10.and.m2==7) then
   wigd=Cos(t*0.5_id)**3*-33.76388603226826436623377881904421997769695431525_id*Sin(t*0.5_id)**17
endif

if(l==10.and.m1==-10.and.m2==8) then
   wigd=Cos(t*0.5_id)**2*13.784048752090221767955912552934175427198163558399_id*Sin(t*0.5_id)**18
endif

if(l==10.and.m1==-10.and.m2==9) then
   wigd=Cos(t*0.5_id)*-4.4721359549995793928183473374625524708812367192231_id*Sin(t*0.5_id)**19
endif

if(l==10.and.m1==-10.and.m2==10) then
   wigd=Sin(t*0.5_id)**20
endif

if(l==10.and.m1==-9.and.m2==-10) then
   wigd=Cos(t*0.5_id)**19*4.4721359549995793928183473374625524708812367192231_id*Sin(t*0.5_id)
endif

if(l==10.and.m1==-9.and.m2==-9) then
   wigd=Cos(t*0.5_id)**18*0.5_id*(-18._id + Cos(t)*20._id)
endif

if(l==10.and.m1==-9.and.m2==-8) then
   wigd=Cos(t*0.5_id)**17*-1.5411035007422441125625480953635610563089060058611_id*(-16._id + Cos(t)*20._id)*Sin(t*0.5_id&
&)
endif

if(l==10.and.m1==-9.and.m2==-7) then
   wigd=Cos(t*0.5_id)**16*3.7749172176353748486183424034730585291110973523117_id*(-14._id + Cos(t)*20._id)*Sin(t*0.5_id)&
&**2
endif

if(l==10.and.m1==-9.and.m2==-6) then
   wigd=Cos(t*0.5_id)**15*-7.782191208136690340748391680403060624268394261302_id*(-12._id + Cos(t)*20._id)*Sin(t*0.5_id)&
&**3
endif

if(l==10.and.m1==-9.and.m2==-5) then
   wigd=Cos(t*0.5_id)**14*13.9212068442358832728203790160422799419385809430456_id*(-10._id + Cos(t)*20._id)*Sin(t*0.5_id&
&)**4
endif

if(l==10.and.m1==-9.and.m2==-4) then
   wigd=Cos(t*0.5_id)**13*-22.011360703055138476529216317127335502674181544772_id*(-8._id + Cos(t)*20._id)*Sin(t*0.5_id)&
&**5
endif

if(l==10.and.m1==-9.and.m2==-3) then
   wigd=Cos(t*0.5_id)**12*(-6._id + Cos(t)*20._id)*31.128764832546761362993566721612242497073577045207_id*Sin(t*0.5_id)*&
&*6
endif

if(l==10.and.m1==-9.and.m2==-2) then
   wigd=Cos(t*0.5_id)**11*-39.681544828799193311277116215239828697071216149532_id*(-4._id + Cos(t)*20._id)*Sin(t*0.5_id)&
&**7
endif

if(l==10.and.m1==-9.and.m2==-1) then
   wigd=Cos(t*0.5_id)**10*(-2._id + Cos(t)*20._id)*45.820301177534832960627900345328910787301645567449_id*Sin(t*0.5_id)*&
&*8
endif

if(l==10.and.m1==-9.and.m2==0) then
   wigd=Cos(t)*Cos(t*0.5_id)**9*-961.13474601639493532560896830571407682331375986_id*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==-9.and.m2==1) then
   wigd=Cos(t*0.5_id)**8*(2._id + Cos(t)*20._id)*45.820301177534832960627900345328910787301645567449_id*Sin(t*0.5_id)**1&
&0
endif

if(l==10.and.m1==-9.and.m2==2) then
   wigd=Cos(t*0.5_id)**7*-39.681544828799193311277116215239828697071216149532_id*(4._id + Cos(t)*20._id)*Sin(t*0.5_id)**&
&11
endif

if(l==10.and.m1==-9.and.m2==3) then
   wigd=Cos(t*0.5_id)**6*(6._id + Cos(t)*20._id)*31.128764832546761362993566721612242497073577045207_id*Sin(t*0.5_id)**1&
&2
endif

if(l==10.and.m1==-9.and.m2==4) then
   wigd=Cos(t*0.5_id)**5*-22.011360703055138476529216317127335502674181544772_id*(8._id + Cos(t)*20._id)*Sin(t*0.5_id)**&
&13
endif

if(l==10.and.m1==-9.and.m2==5) then
   wigd=Cos(t*0.5_id)**4*13.9212068442358832728203790160422799419385809430456_id*(10._id + Cos(t)*20._id)*Sin(t*0.5_id)*&
&*14
endif

if(l==10.and.m1==-9.and.m2==6) then
   wigd=Cos(t*0.5_id)**3*-7.782191208136690340748391680403060624268394261302_id*(12._id + Cos(t)*20._id)*Sin(t*0.5_id)**&
&15
endif

if(l==10.and.m1==-9.and.m2==7) then
   wigd=Cos(t*0.5_id)**2*3.7749172176353748486183424034730585291110973523117_id*(14._id + Cos(t)*20._id)*Sin(t*0.5_id)**&
&16
endif

if(l==10.and.m1==-9.and.m2==8) then
   wigd=Cos(t*0.5_id)*-1.5411035007422441125625480953635610563089060058611_id*(16._id + Cos(t)*20._id)*Sin(t*0.5_id)**17&
&
endif

if(l==10.and.m1==-9.and.m2==9) then
   wigd=0.5_id*(18._id + Cos(t)*20._id)*Sin(t*0.5_id)**18
endif

if(l==10.and.m1==-9.and.m2==10) then
   wigd=Cos(t*0.5_id)*-4.4721359549995793928183473374625524708812367192231_id*Sin(t*0.5_id)**19
endif

if(l==10.and.m1==-8.and.m2==-10) then
   wigd=Cos(t*0.5_id)**18*13.784048752090221767955912552934175427198163558399_id*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==-8.and.m2==-9) then
   wigd=Cos(t*0.5_id)**17*1.5411035007422441125625480953635610563089060058611_id*(-16._id + Cos(t)*20._id)*Sin(t*0.5_id)&
&
endif

if(l==10.and.m1==-8.and.m2==-8) then
   wigd=Cos(t*0.5_id)**16*(1._id + (Cos(t) + -1._id)*19._id + (Cos(t) + -1._id)**2*47.5_id)
endif

if(l==10.and.m1==-8.and.m2==-7) then
   wigd=Cos(t*0.5_id)**15*-2.4494897427831780981972840747058913919659474806567_id*(3._id + (Cos(t) + -1._id)*28.5_id + (&
&Cos(t) + -1._id)**2*47.5_id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==-8.and.m2==-6) then
   wigd=Cos(t*0.5_id)**14*5.0497524691810389766816929585348003553019494822398_id*(6._id + (Cos(t) + -1._id)*38._id + (Co&
&s(t) + -1._id)**2*47.5_id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==-8.and.m2==-5) then
   wigd=Cos(t*0.5_id)**13*-9.033271832508971939888199943876461581439881615896_id*(10._id + (Cos(t) + -1._id)*47.5_id + (&
&Cos(t) + -1._id)**2*47.5_id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==-8.and.m2==-4) then
   wigd=Cos(t*0.5_id)**12*14.282856857085699995998799622734530557532342319805_id*(15._id + (Cos(t) + -1._id)**2*47.5_id &
&+ (Cos(t) + -1._id)*57._id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==-8.and.m2==-3) then
   wigd=Cos(t*0.5_id)**11*-20.199009876724155906726771834139201421207797928959_id*(21._id + (Cos(t) + -1._id)**2*47.5_id&
& + (Cos(t) + -1._id)*66.5_id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==-8.and.m2==-2) then
   wigd=Cos(t*0.5_id)**10*25.748786379167465530841591578687506418644843651043_id*(28._id + (Cos(t) + -1._id)**2*47.5_id &
&+ (Cos(t) + -1._id)*76._id)*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==-8.and.m2==-1) then
   wigd=Cos(t*0.5_id)**9*-29.732137494637011045224016427862793302897971027442_id*(36._id + (Cos(t) + -1._id)**2*47.5_id &
&+ (Cos(t) + -1._id)*85.5_id)*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==-8.and.m2==0) then
   wigd=Cos(t*0.5_id)**8*31.183328879386818922579762898076675083852016162476_id*(45._id + (Cos(t) + -1._id)**2*47.5_id +&
& (Cos(t) + -1._id)*95._id)*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==-8.and.m2==1) then
   wigd=Cos(t*0.5_id)**7*-29.732137494637011045224016427862793302897971027442_id*((Cos(t) + -1._id)**2*47.5_id + 55._id &
&+ (Cos(t) + -1._id)*104.5_id)*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==-8.and.m2==2) then
   wigd=Cos(t*0.5_id)**6*25.748786379167465530841591578687506418644843651043_id*((Cos(t) + -1._id)**2*47.5_id + 66._id +&
& (Cos(t) + -1._id)*114._id)*Sin(t*0.5_id)**10
endif

if(l==10.and.m1==-8.and.m2==3) then
   wigd=Cos(t*0.5_id)**5*-20.199009876724155906726771834139201421207797928959_id*((Cos(t) + -1._id)**2*47.5_id + 78._id &
&+ (Cos(t) + -1._id)*123.5_id)*Sin(t*0.5_id)**11
endif

if(l==10.and.m1==-8.and.m2==4) then
   wigd=Cos(t*0.5_id)**4*14.282856857085699995998799622734530557532342319805_id*((Cos(t) + -1._id)**2*47.5_id + 91._id +&
& (Cos(t) + -1._id)*133._id)*Sin(t*0.5_id)**12
endif

if(l==10.and.m1==-8.and.m2==5) then
   wigd=Cos(t*0.5_id)**3*-9.033271832508971939888199943876461581439881615896_id*((Cos(t) + -1._id)**2*47.5_id + 105._id &
&+ (Cos(t) + -1._id)*142.5_id)*Sin(t*0.5_id)**13
endif

if(l==10.and.m1==-8.and.m2==6) then
   wigd=Cos(t*0.5_id)**2*5.0497524691810389766816929585348003553019494822398_id*((Cos(t) + -1._id)**2*47.5_id + 120._id &
&+ (Cos(t) + -1._id)*152._id)*Sin(t*0.5_id)**14
endif

if(l==10.and.m1==-8.and.m2==7) then
   wigd=Cos(t*0.5_id)*-2.4494897427831780981972840747058913919659474806567_id*((Cos(t) + -1._id)**2*47.5_id + 136._id + &
&(Cos(t) + -1._id)*161.5_id)*Sin(t*0.5_id)**15
endif

if(l==10.and.m1==-8.and.m2==8) then
   wigd=((Cos(t) + -1._id)**2*47.5_id + 153._id + (Cos(t) + -1._id)*171._id)*Sin(t*0.5_id)**16
endif

if(l==10.and.m1==-8.and.m2==9) then
   wigd=Cos(t*0.5_id)*-1.5411035007422441125625480953635610563089060058611_id*(16._id + Cos(t)*20._id)*Sin(t*0.5_id)**17&
&
endif

if(l==10.and.m1==-8.and.m2==10) then
   wigd=Cos(t*0.5_id)**2*13.784048752090221767955912552934175427198163558399_id*Sin(t*0.5_id)**18
endif

if(l==10.and.m1==-7.and.m2==-10) then
   wigd=Cos(t*0.5_id)**17*33.76388603226826436623377881904421997769695431525_id*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==-7.and.m2==-9) then
   wigd=Cos(t*0.5_id)**16*3.7749172176353748486183424034730585291110973523117_id*(-14._id + Cos(t)*20._id)*Sin(t*0.5_id)&
&**2
endif

if(l==10.and.m1==-7.and.m2==-8) then
   wigd=Cos(t*0.5_id)**15*2.4494897427831780981972840747058913919659474806567_id*(3._id + (Cos(t) + -1._id)*28.5_id + (C&
&os(t) + -1._id)**2*47.5_id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==-7.and.m2==-7) then
   wigd=Cos(t*0.5_id)**14*(1._id + (Cos(t) + -1._id)*27._id + (Cos(t) + -1._id)**2*128.25_id + (Cos(t) + -1._id)**3*142.&
&5_id)
endif

if(l==10.and.m1==-7.and.m2==-6) then
   wigd=Cos(t*0.5_id)**13*-2.0615528128088302749107049279870385125735996126868_id*(4._id + (Cos(t) + -1._id)*54._id + (C&
&os(t) + -1._id)**3*142.5_id + (Cos(t) + -1._id)**2*171._id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==-7.and.m2==-5) then
   wigd=Cos(t*0.5_id)**12*3.687817782917154924000909712705117262898722019489_id*(10._id + (Cos(t) + -1._id)*90._id + (Co&
&s(t) + -1._id)**3*142.5_id + (Cos(t) + -1._id)**2*213.75_id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==-7.and.m2==-4) then
   wigd=Cos(t*0.5_id)**11*-5.830951894845300470874152877545583076521398334886_id*(20._id + (Cos(t) + -1._id)*135._id + (&
&Cos(t) + -1._id)**3*142.5_id + (Cos(t) + -1._id)**2*256.5_id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==-7.and.m2==-3) then
   wigd=Cos(t*0.5_id)**10*8.246211251235321099642819711948154050294398450747_id*(35._id + (Cos(t) + -1._id)**3*142.5_id &
&+ (Cos(t) + -1._id)*189._id + (Cos(t) + -1._id)**2*299.25_id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==-7.and.m2==-2) then
   wigd=Cos(t*0.5_id)**9*-10.51189802081431914420992852874790396912039225191_id*(56._id + (Cos(t) + -1._id)**3*142.5_id &
&+ (Cos(t) + -1._id)*252._id + (Cos(t) + -1._id)**2*342._id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==-7.and.m2==-1) then
   wigd=Cos(t*0.5_id)**8*12.138094304022082911201150512925122493534377641012_id*(84._id + (Cos(t) + -1._id)**3*142.5_id &
&+ (Cos(t) + -1._id)*324._id + (Cos(t) + -1._id)**2*384.75_id)*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==-7.and.m2==0) then
   wigd=Cos(t*0.5_id)**7*-12.730540705982078068014872638317794132486379784447_id*(120._id + (Cos(t) + -1._id)**3*142.5_i&
&d + (Cos(t) + -1._id)*405._id + (Cos(t) + -1._id)**2*427.5_id)*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==-7.and.m2==1) then
   wigd=Cos(t*0.5_id)**6*12.138094304022082911201150512925122493534377641012_id*((Cos(t) + -1._id)**3*142.5_id + 165._id&
& + (Cos(t) + -1._id)**2*470.25_id + (Cos(t) + -1._id)*495._id)*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==-7.and.m2==2) then
   wigd=Cos(t*0.5_id)**5*-10.51189802081431914420992852874790396912039225191_id*((Cos(t) + -1._id)**3*142.5_id + 220._id&
& + (Cos(t) + -1._id)**2*513._id + (Cos(t) + -1._id)*594._id)*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==-7.and.m2==3) then
   wigd=Cos(t*0.5_id)**4*8.246211251235321099642819711948154050294398450747_id*((Cos(t) + -1._id)**3*142.5_id + 286._id &
&+ (Cos(t) + -1._id)**2*555.75_id + (Cos(t) + -1._id)*702._id)*Sin(t*0.5_id)**10
endif

if(l==10.and.m1==-7.and.m2==4) then
   wigd=Cos(t*0.5_id)**3*-5.830951894845300470874152877545583076521398334886_id*((Cos(t) + -1._id)**3*142.5_id + 364._id&
& + (Cos(t) + -1._id)**2*598.5_id + (Cos(t) + -1._id)*819._id)*Sin(t*0.5_id)**11
endif

if(l==10.and.m1==-7.and.m2==5) then
   wigd=Cos(t*0.5_id)**2*3.687817782917154924000909712705117262898722019489_id*((Cos(t) + -1._id)**3*142.5_id + 455._id &
&+ (Cos(t) + -1._id)**2*641.25_id + (Cos(t) + -1._id)*945._id)*Sin(t*0.5_id)**12
endif

if(l==10.and.m1==-7.and.m2==6) then
   wigd=Cos(t*0.5_id)*-2.0615528128088302749107049279870385125735996126868_id*((Cos(t) + -1._id)**3*142.5_id + 560._id +&
& (Cos(t) + -1._id)**2*684._id + (Cos(t) + -1._id)*1080._id)*Sin(t*0.5_id)**13
endif

if(l==10.and.m1==-7.and.m2==7) then
   wigd=((Cos(t) + -1._id)**3*142.5_id + 680._id + (Cos(t) + -1._id)**2*726.75_id + (Cos(t) + -1._id)*1224._id)*Sin(t*0.&
&5_id)**14
endif

if(l==10.and.m1==-7.and.m2==8) then
   wigd=Cos(t*0.5_id)*-2.4494897427831780981972840747058913919659474806567_id*((Cos(t) + -1._id)**2*47.5_id + 136._id + &
&(Cos(t) + -1._id)*161.5_id)*Sin(t*0.5_id)**15
endif

if(l==10.and.m1==-7.and.m2==9) then
   wigd=Cos(t*0.5_id)**2*3.7749172176353748486183424034730585291110973523117_id*(14._id + Cos(t)*20._id)*Sin(t*0.5_id)**&
&16
endif

if(l==10.and.m1==-7.and.m2==10) then
   wigd=Cos(t*0.5_id)**3*-33.76388603226826436623377881904421997769695431525_id*Sin(t*0.5_id)**17
endif

if(l==10.and.m1==-6.and.m2==-10) then
   wigd=Cos(t*0.5_id)**16*69.606034221179416364101895080211399709692904715228_id*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==-6.and.m2==-9) then
   wigd=Cos(t*0.5_id)**15*7.782191208136690340748391680403060624268394261302_id*(-12._id + Cos(t)*20._id)*Sin(t*0.5_id)*&
&*3
endif

if(l==10.and.m1==-6.and.m2==-8) then
   wigd=Cos(t*0.5_id)**14*5.0497524691810389766816929585348003553019494822398_id*(6._id + (Cos(t) + -1._id)*38._id + (Co&
&s(t) + -1._id)**2*47.5_id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==-6.and.m2==-7) then
   wigd=Cos(t*0.5_id)**13*2.0615528128088302749107049279870385125735996126868_id*(4._id + (Cos(t) + -1._id)*54._id + (Co&
&s(t) + -1._id)**3*142.5_id + (Cos(t) + -1._id)**2*171._id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==-6.and.m2==-6) then
   wigd=Cos(t*0.5_id)**12*(1._id + (Cos(t) + -1._id)*34._id + (Cos(t) + -1._id)**2*229.5_id + (Cos(t) + -1._id)**4*302.8&
&125_id + (Cos(t) + -1._id)**3*484.5_id)
endif

if(l==10.and.m1==-6.and.m2==-5) then
   wigd=Cos(t*0.5_id)**11*-1.7888543819998317571273389349850209883524946876892_id*(5._id + (Cos(t) + -1._id)*85._id + (C&
&os(t) + -1._id)**4*302.8125_id + (Cos(t) + -1._id)**2*382.5_id + (Cos(t) + -1._id)**3*605.625_id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==-6.and.m2==-4) then
   wigd=Cos(t*0.5_id)**10*2.8284271247461900976033774484193961571393437507539_id*(15._id + (Cos(t) + -1._id)*170._id + (&
&Cos(t) + -1._id)**4*302.8125_id + (Cos(t) + -1._id)**2*573.75_id + (Cos(t) + -1._id)**3*726.75_id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==-6.and.m2==-3) then
   wigd=Cos(t*0.5_id)**9*-4._id*(35._id + (Cos(t) + -1._id)*297.5_id + (Cos(t) + -1._id)**4*302.8125_id + (Cos(t) + -1._&
&id)**2*803.25_id + (Cos(t) + -1._id)**3*847.875_id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==-6.and.m2==-2) then
   wigd=Cos(t*0.5_id)**8*5.0990195135927848300282241090227819895637709460996_id*(70._id + (Cos(t) + -1._id)**4*302.8125_&
&id + (Cos(t) + -1._id)*476._id + (Cos(t) + -1._id)**3*969._id + (Cos(t) + -1._id)**2*1071._id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==-6.and.m2==-1) then
   wigd=Cos(t*0.5_id)**7*-5.8878405775518979031760284846395026427830929872876_id*(126._id + (Cos(t) + -1._id)**4*302.812&
&5_id + (Cos(t) + -1._id)*714._id + (Cos(t) + -1._id)**3*1090.125_id + (Cos(t) + -1._id)**2*1377._id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==-6.and.m2==0) then
   wigd=Cos(t*0.5_id)**6*6.1752192943516858827531745282506443172988423053509_id*(210._id + (Cos(t) + -1._id)**4*302.8125&
&_id + (Cos(t) + -1._id)*1020._id + (Cos(t) + -1._id)**3*1211.25_id + (Cos(t) + -1._id)**2*1721.25_id)*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==-6.and.m2==1) then
   wigd=Cos(t*0.5_id)**5*-5.8878405775518979031760284846395026427830929872876_id*((Cos(t) + -1._id)**4*302.8125_id + 330&
&._id + (Cos(t) + -1._id)**3*1332.375_id + (Cos(t) + -1._id)*1402.5_id + (Cos(t) + -1._id)**2*2103.75_id)*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==-6.and.m2==2) then
   wigd=Cos(t*0.5_id)**4*5.0990195135927848300282241090227819895637709460996_id*((Cos(t) + -1._id)**4*302.8125_id + 495.&
&_id + (Cos(t) + -1._id)**3*1453.5_id + (Cos(t) + -1._id)*1870._id + (Cos(t) + -1._id)**2*2524.5_id)*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==-6.and.m2==3) then
   wigd=Cos(t*0.5_id)**3*-4._id*((Cos(t) + -1._id)**4*302.8125_id + 715._id + (Cos(t) + -1._id)**3*1574.625_id + (Cos(t)&
& + -1._id)*2431._id + (Cos(t) + -1._id)**2*2983.5_id)*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==-6.and.m2==4) then
   wigd=Cos(t*0.5_id)**2*2.8284271247461900976033774484193961571393437507539_id*((Cos(t) + -1._id)**4*302.8125_id + 1001&
&._id + (Cos(t) + -1._id)**3*1695.75_id + (Cos(t) + -1._id)*3094._id + (Cos(t) + -1._id)**2*3480.75_id)*Sin(t*0.5_id)**10
endif

if(l==10.and.m1==-6.and.m2==5) then
   wigd=Cos(t*0.5_id)*-1.7888543819998317571273389349850209883524946876892_id*((Cos(t) + -1._id)**4*302.8125_id + 1365._&
&id + (Cos(t) + -1._id)**3*1816.875_id + (Cos(t) + -1._id)*3867.5_id + (Cos(t) + -1._id)**2*4016.25_id)*Sin(t*0.5_id)**11
endif

if(l==10.and.m1==-6.and.m2==6) then
   wigd=((Cos(t) + -1._id)**4*302.8125_id + 1820._id + (Cos(t) + -1._id)**3*1938._id + (Cos(t) + -1._id)**2*4590._id + (&
&Cos(t) + -1._id)*4760._id)*Sin(t*0.5_id)**12
endif

if(l==10.and.m1==-6.and.m2==7) then
   wigd=Cos(t*0.5_id)*-2.0615528128088302749107049279870385125735996126868_id*((Cos(t) + -1._id)**3*142.5_id + 560._id +&
& (Cos(t) + -1._id)**2*684._id + (Cos(t) + -1._id)*1080._id)*Sin(t*0.5_id)**13
endif

if(l==10.and.m1==-6.and.m2==8) then
   wigd=Cos(t*0.5_id)**2*5.0497524691810389766816929585348003553019494822398_id*((Cos(t) + -1._id)**2*47.5_id + 120._id &
&+ (Cos(t) + -1._id)*152._id)*Sin(t*0.5_id)**14
endif

if(l==10.and.m1==-6.and.m2==9) then
   wigd=Cos(t*0.5_id)**3*-7.782191208136690340748391680403060624268394261302_id*(12._id + Cos(t)*20._id)*Sin(t*0.5_id)**&
&15
endif

if(l==10.and.m1==-6.and.m2==10) then
   wigd=Cos(t*0.5_id)**4*69.606034221179416364101895080211399709692904715228_id*Sin(t*0.5_id)**16
endif

if(l==10.and.m1==-5.and.m2==-10) then
   wigd=Cos(t*0.5_id)**15*124.51505933018704545197426688644896998829430818083_id*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==-5.and.m2==-9) then
   wigd=Cos(t*0.5_id)**14*13.9212068442358832728203790160422799419385809430456_id*(-10._id + Cos(t)*20._id)*Sin(t*0.5_id&
&)**4
endif

if(l==10.and.m1==-5.and.m2==-8) then
   wigd=Cos(t*0.5_id)**13*9.033271832508971939888199943876461581439881615896_id*(10._id + (Cos(t) + -1._id)*47.5_id + (C&
&os(t) + -1._id)**2*47.5_id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==-5.and.m2==-7) then
   wigd=Cos(t*0.5_id)**12*3.687817782917154924000909712705117262898722019489_id*(10._id + (Cos(t) + -1._id)*90._id + (Co&
&s(t) + -1._id)**3*142.5_id + (Cos(t) + -1._id)**2*213.75_id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==-5.and.m2==-6) then
   wigd=Cos(t*0.5_id)**11*1.7888543819998317571273389349850209883524946876892_id*(5._id + (Cos(t) + -1._id)*85._id + (Co&
&s(t) + -1._id)**4*302.8125_id + (Cos(t) + -1._id)**2*382.5_id + (Cos(t) + -1._id)**3*605.625_id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==-5.and.m2==-5) then
   wigd=Cos(t*0.5_id)**10*(1._id + (Cos(t) + -1._id)*40._id + (Cos(t) + -1._id)**2*340._id + (Cos(t) + -1._id)**5*484.5_&
&id + (Cos(t) + -1._id)**3*1020._id + (Cos(t) + -1._id)**4*1211.25_id)
endif

if(l==10.and.m1==-5.and.m2==-4) then
   wigd=Cos(t*0.5_id)**9*-1.5811388300841896659994467722163592668597775696626_id*(6._id + (Cos(t) + -1._id)*120._id + (C&
&os(t) + -1._id)**5*484.5_id + (Cos(t) + -1._id)**2*680._id + (Cos(t) + -1._id)**4*1453.5_id + (Cos(t) + -1._id)**3*1530._id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==-5.and.m2==-3) then
   wigd=Cos(t*0.5_id)**8*2.2360679774997896964091736687312762354406183596115_id*(21._id + (Cos(t) + -1._id)*280._id + (C&
&os(t) + -1._id)**5*484.5_id + (Cos(t) + -1._id)**2*1190._id + (Cos(t) + -1._id)**4*1695.75_id + (Cos(t) + -1._id)**3*2142._id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==-5.and.m2==-2) then
   wigd=Cos(t*0.5_id)**7*-2.850438562747844947840122563916886197690013277291_id*(56._id + (Cos(t) + -1._id)**5*484.5_id &
&+ (Cos(t) + -1._id)*560._id + (Cos(t) + -1._id)**2*1904._id + (Cos(t) + -1._id)**4*1938._id + (Cos(t) + -1._id)**3*2856._id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==-5.and.m2==-1) then
   wigd=Cos(t*0.5_id)**6*3.2914029430219165029064101739538498443475436446803_id*(126._id + (Cos(t) + -1._id)**5*484.5_id&
& + (Cos(t) + -1._id)*1008._id + (Cos(t) + -1._id)**4*2180.25_id + (Cos(t) + -1._id)**2*2856._id + (Cos(t) + -1._id)**3*3672._id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==-5.and.m2==0) then
   wigd=Cos(t*0.5_id)**5*-3.4520525295346631886928627240080226103155020819111_id*(252._id + (Cos(t) + -1._id)**5*484.5_i&
&d + (Cos(t) + -1._id)*1680._id + (Cos(t) + -1._id)**4*2422.5_id + (Cos(t) + -1._id)**2*4080._id + (Cos(t) + -1._id)**3*4590._id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==-5.and.m2==1) then
   wigd=Cos(t*0.5_id)**4*3.2914029430219165029064101739538498443475436446803_id*(462._id + (Cos(t) + -1._id)**5*484.5_id&
& + (Cos(t) + -1._id)*2640._id + (Cos(t) + -1._id)**4*2664.75_id + (Cos(t) + -1._id)**2*5610._id + (Cos(t) + -1._id)**3*5610._id)*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==-5.and.m2==2) then
   wigd=Cos(t*0.5_id)**3*-2.850438562747844947840122563916886197690013277291_id*((Cos(t) + -1._id)**5*484.5_id + 792._id&
& + (Cos(t) + -1._id)**4*2907._id + (Cos(t) + -1._id)*3960._id + (Cos(t) + -1._id)**3*6732._id + (Cos(t) + -1._id)**2*7480._id)*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==-5.and.m2==3) then
   wigd=Cos(t*0.5_id)**2*2.2360679774997896964091736687312762354406183596115_id*((Cos(t) + -1._id)**5*484.5_id + 1287._i&
&d + (Cos(t) + -1._id)**4*3149.25_id + (Cos(t) + -1._id)*5720._id + (Cos(t) + -1._id)**3*7956._id + (Cos(t) + -1._id)**2*9724._id)*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==-5.and.m2==4) then
   wigd=Cos(t*0.5_id)*-1.5811388300841896659994467722163592668597775696626_id*((Cos(t) + -1._id)**5*484.5_id + 2002._id &
&+ (Cos(t) + -1._id)**4*3391.5_id + (Cos(t) + -1._id)*8008._id + (Cos(t) + -1._id)**3*9282._id + (Cos(t) + -1._id)**2*12376._id)*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==-5.and.m2==5) then
   wigd=((Cos(t) + -1._id)**5*484.5_id + 3003._id + (Cos(t) + -1._id)**4*3633.75_id + (Cos(t) + -1._id)**3*10710._id + (&
&Cos(t) + -1._id)*10920._id + (Cos(t) + -1._id)**2*15470._id)*Sin(t*0.5_id)**10
endif

if(l==10.and.m1==-5.and.m2==6) then
   wigd=Cos(t*0.5_id)*-1.7888543819998317571273389349850209883524946876892_id*((Cos(t) + -1._id)**4*302.8125_id + 1365._&
&id + (Cos(t) + -1._id)**3*1816.875_id + (Cos(t) + -1._id)*3867.5_id + (Cos(t) + -1._id)**2*4016.25_id)*Sin(t*0.5_id)**11
endif

if(l==10.and.m1==-5.and.m2==7) then
   wigd=Cos(t*0.5_id)**2*3.687817782917154924000909712705117262898722019489_id*((Cos(t) + -1._id)**3*142.5_id + 455._id &
&+ (Cos(t) + -1._id)**2*641.25_id + (Cos(t) + -1._id)*945._id)*Sin(t*0.5_id)**12
endif

if(l==10.and.m1==-5.and.m2==8) then
   wigd=Cos(t*0.5_id)**3*-9.033271832508971939888199943876461581439881615896_id*((Cos(t) + -1._id)**2*47.5_id + 105._id &
&+ (Cos(t) + -1._id)*142.5_id)*Sin(t*0.5_id)**13
endif

if(l==10.and.m1==-5.and.m2==9) then
   wigd=Cos(t*0.5_id)**4*13.9212068442358832728203790160422799419385809430456_id*(10._id + Cos(t)*20._id)*Sin(t*0.5_id)*&
&*14
endif

if(l==10.and.m1==-5.and.m2==10) then
   wigd=Cos(t*0.5_id)**5*-124.51505933018704545197426688644896998829430818083_id*Sin(t*0.5_id)**15
endif

if(l==10.and.m1==-4.and.m2==-10) then
   wigd=Cos(t*0.5_id)**14*196.87559523719540998400115941877306606266144591633_id*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==-4.and.m2==-9) then
   wigd=Cos(t*0.5_id)**13*(-8._id + Cos(t)*20._id)*22.011360703055138476529216317127335502674181544772_id*Sin(t*0.5_id)*&
&*5
endif

if(l==10.and.m1==-4.and.m2==-8) then
   wigd=Cos(t*0.5_id)**12*14.282856857085699995998799622734530557532342319805_id*(15._id + (Cos(t) + -1._id)**2*47.5_id &
&+ (Cos(t) + -1._id)*57._id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==-4.and.m2==-7) then
   wigd=Cos(t*0.5_id)**11*5.830951894845300470874152877545583076521398334886_id*(20._id + (Cos(t) + -1._id)*135._id + (C&
&os(t) + -1._id)**3*142.5_id + (Cos(t) + -1._id)**2*256.5_id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==-4.and.m2==-6) then
   wigd=Cos(t*0.5_id)**10*2.8284271247461900976033774484193961571393437507539_id*(15._id + (Cos(t) + -1._id)*170._id + (&
&Cos(t) + -1._id)**4*302.8125_id + (Cos(t) + -1._id)**2*573.75_id + (Cos(t) + -1._id)**3*726.75_id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==-4.and.m2==-5) then
   wigd=Cos(t*0.5_id)**9*1.58113883008418966599944677221635926685977756966261_id*(6._id + (Cos(t) + -1._id)*120._id + (C&
&os(t) + -1._id)**5*484.5_id + (Cos(t) + -1._id)**2*680._id + (Cos(t) + -1._id)**4*1453.5_id + (Cos(t) + -1._id)**3*1530._id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==-4.and.m2==-4) then
   wigd=Cos(t*0.5_id)**8*(1._id + (Cos(t) + -1._id)*45._id + (Cos(t) + -1._id)**2*450._id + (Cos(t) + -1._id)**6*605.625&
&_id + (Cos(t) + -1._id)**3*1700._id + (Cos(t) + -1._id)**5*2180.25_id + (Cos(t) + -1._id)**4*2868.75_id)
endif

if(l==10.and.m1==-4.and.m2==-3) then
   wigd=Cos(t*0.5_id)**7*-1.4142135623730950488016887242096980785696718753769_id*(7._id + (Cos(t) + -1._id)*157.5_id + (&
&Cos(t) + -1._id)**6*605.625_id + (Cos(t) + -1._id)**2*1050._id + (Cos(t) + -1._id)**5*2543.625_id + (Cos(t) + -1._id)**3*2975._id + (Cos(t) + -1._id)**4*4016.25_id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==-4.and.m2==-2) then
   wigd=Cos(t*0.5_id)**6*1.8027756377319946465596106337352479731256482869226_id*(28._id + (Cos(t) + -1._id)*420._id + (C&
&os(t) + -1._id)**6*605.625_id + (Cos(t) + -1._id)**2*2100._id + (Cos(t) + -1._id)**5*2907._id + (Cos(t) + -1._id)**3*4760._id + (Cos(t) + -1._id)**4*5355._id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==-4.and.m2==-1) then
   wigd=Cos(t*0.5_id)**5*-2.0816659994661327352822977069799314870243199926639_id*(84._id + (Cos(t) + -1._id)**6*605.625_&
&id + (Cos(t) + -1._id)*945._id + (Cos(t) + -1._id)**5*3270.375_id + (Cos(t) + -1._id)**2*3780._id + (Cos(t) + -1._id)**4*6885._id + (Cos(t) + -1._id)**3*7140._id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==-4.and.m2==0) then
   wigd=Cos(t*0.5_id)**4*2.1832697191750419792351883418789116087243562822933_id*(210._id + (Cos(t) + -1._id)**6*605.625_&
&id + (Cos(t) + -1._id)*1890._id + (Cos(t) + -1._id)**5*3633.75_id + (Cos(t) + -1._id)**2*6300._id + (Cos(t) + -1._id)**4*8606.25_id + (Cos(t) + -1._id)**3*10200._id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==-4.and.m2==1) then
   wigd=Cos(t*0.5_id)**3*-2.0816659994661327352822977069799314870243199926639_id*(462._id + (Cos(t) + -1._id)**6*605.625&
&_id + (Cos(t) + -1._id)*3465._id + (Cos(t) + -1._id)**5*3997.125_id + (Cos(t) + -1._id)**2*9900._id + (Cos(t) + -1._id)**4*10518.75_id + (Cos(t) + -1._id)**3*14025._id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==-4.and.m2==2) then
   wigd=Cos(t*0.5_id)**2*1.8027756377319946465596106337352479731256482869226_id*((Cos(t) + -1._id)**6*605.625_id + 924._&
&id + (Cos(t) + -1._id)**5*4360.5_id + (Cos(t) + -1._id)*5940._id + (Cos(t) + -1._id)**4*12622.5_id + (Cos(t) + -1._id)**2*14850._id + (Cos(t) + -1._id)**3*18700._id)*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==-4.and.m2==3) then
   wigd=Cos(t*0.5_id)*-1.4142135623730950488016887242096980785696718753769_id*((Cos(t) + -1._id)**6*605.625_id + 1716._i&
&d + (Cos(t) + -1._id)**5*4723.875_id + (Cos(t) + -1._id)*9652.5_id + (Cos(t) + -1._id)**4*14917.5_id + (Cos(t) + -1._id)**2*21450._id + (Cos(t) + -1._id)**3*24310._id)*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==-4.and.m2==4) then
   wigd=((Cos(t) + -1._id)**6*605.625_id + 3003._id + (Cos(t) + -1._id)**5*5087.25_id + (Cos(t) + -1._id)*15015._id + (C&
&os(t) + -1._id)**4*17403.75_id + (Cos(t) + -1._id)**2*30030._id + (Cos(t) + -1._id)**3*30940._id)*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==-4.and.m2==5) then
   wigd=Cos(t*0.5_id)*-1.5811388300841896659994467722163592668597775696626_id*((Cos(t) + -1._id)**5*484.5_id + 2002._id &
&+ (Cos(t) + -1._id)**4*3391.5_id + (Cos(t) + -1._id)*8008._id + (Cos(t) + -1._id)**3*9282._id + (Cos(t) + -1._id)**2*12376._id)*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==-4.and.m2==6) then
   wigd=Cos(t*0.5_id)**2*2.8284271247461900976033774484193961571393437507539_id*((Cos(t) + -1._id)**4*302.8125_id + 1001&
&._id + (Cos(t) + -1._id)**3*1695.75_id + (Cos(t) + -1._id)*3094._id + (Cos(t) + -1._id)**2*3480.75_id)*Sin(t*0.5_id)**10
endif

if(l==10.and.m1==-4.and.m2==7) then
   wigd=Cos(t*0.5_id)**3*-5.830951894845300470874152877545583076521398334886_id*((Cos(t) + -1._id)**3*142.5_id + 364._id&
& + (Cos(t) + -1._id)**2*598.5_id + (Cos(t) + -1._id)*819._id)*Sin(t*0.5_id)**11
endif

if(l==10.and.m1==-4.and.m2==8) then
   wigd=Cos(t*0.5_id)**4*14.282856857085699995998799622734530557532342319805_id*((Cos(t) + -1._id)**2*47.5_id + 91._id +&
& (Cos(t) + -1._id)*133._id)*Sin(t*0.5_id)**12
endif

if(l==10.and.m1==-4.and.m2==9) then
   wigd=Cos(t*0.5_id)**5*-22.011360703055138476529216317127335502674181544772_id*(8._id + Cos(t)*20._id)*Sin(t*0.5_id)**&
&13
endif

if(l==10.and.m1==-4.and.m2==10) then
   wigd=Cos(t*0.5_id)**6*196.87559523719540998400115941877306606266144591633_id*Sin(t*0.5_id)**14
endif

if(l==10.and.m1==-3.and.m2==-10) then
   wigd=Cos(t*0.5_id)**13*278.42413688471766545640758032084559883877161886091_id*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==-3.and.m2==-9) then
   wigd=Cos(t*0.5_id)**12*(-6._id + Cos(t)*20._id)*31.128764832546761362993566721612242497073577045207_id*Sin(t*0.5_id)*&
&*6
endif

if(l==10.and.m1==-3.and.m2==-8) then
   wigd=Cos(t*0.5_id)**11*20.199009876724155906726771834139201421207797928959_id*(21._id + (Cos(t) + -1._id)**2*47.5_id &
&+ (Cos(t) + -1._id)*66.5_id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==-3.and.m2==-7) then
   wigd=Cos(t*0.5_id)**10*8.246211251235321099642819711948154050294398450747_id*(35._id + (Cos(t) + -1._id)**3*142.5_id &
&+ (Cos(t) + -1._id)*189._id + (Cos(t) + -1._id)**2*299.25_id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==-3.and.m2==-6) then
   wigd=Cos(t*0.5_id)**9*4._id*(35._id + (Cos(t) + -1._id)*297.5_id + (Cos(t) + -1._id)**4*302.8125_id + (Cos(t) + -1._i&
&d)**2*803.25_id + (Cos(t) + -1._id)**3*847.875_id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==-3.and.m2==-5) then
   wigd=Cos(t*0.5_id)**8*2.2360679774997896964091736687312762354406183596115_id*(21._id + (Cos(t) + -1._id)*280._id + (C&
&os(t) + -1._id)**5*484.5_id + (Cos(t) + -1._id)**2*1190._id + (Cos(t) + -1._id)**4*1695.75_id + (Cos(t) + -1._id)**3*2142._id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==-3.and.m2==-4) then
   wigd=Cos(t*0.5_id)**7*1.41421356237309504880168872420969807856967187537695_id*(7._id + (Cos(t) + -1._id)*157.5_id + (&
&Cos(t) + -1._id)**6*605.625_id + (Cos(t) + -1._id)**2*1050._id + (Cos(t) + -1._id)**5*2543.625_id + (Cos(t) + -1._id)**3*2975._id + (Cos(t) + -1._id)**4*4016.25_id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==-3.and.m2==-3) then
   wigd=Cos(t*0.5_id)**6*(1._id + (Cos(t) + -1._id)*49._id + (Cos(t) + -1._id)**2*551.25_id + (Cos(t) + -1._id)**7*605.6&
&25_id + (Cos(t) + -1._id)**3*2450._id + (Cos(t) + -1._id)**6*2967.5625_id + (Cos(t) + -1._id)**4*5206.25_id + (Cos(t) + -1._id)**5*5622.75_id)
endif

if(l==10.and.m1==-3.and.m2==-2) then
   wigd=Cos(t*0.5_id)**5*-1.2747548783981962075070560272556954973909427365249_id*(8._id + (Cos(t) + -1._id)*196._id + (C&
&os(t) + -1._id)**7*605.625_id + (Cos(t) + -1._id)**2*1470._id + (Cos(t) + -1._id)**6*3391.5_id + (Cos(t) + -1._id)**3*4900._id + (Cos(t) + -1._id)**5*7497._id + (Cos(t) + -1._id)**4*8330._id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==-3.and.m2==-1) then
   wigd=Cos(t*0.5_id)**4*1.4719601443879744757940071211598756606957732468219_id*(36._id + (Cos(t) + -1._id)*588._id + (C&
&os(t) + -1._id)**7*605.625_id + (Cos(t) + -1._id)**2*3307.5_id + (Cos(t) + -1._id)**6*3815.4375_id + (Cos(t) + -1._id)**3*8820._id + (Cos(t) + -1._id)**5*9639._id + (Cos(t) + -1._id)**4*12495._id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==-3.and.m2==0) then
   wigd=Cos(t*0.5_id)**3*-1.5438048235879214706882936320626610793247105763377_id*(120._id + (Cos(t) + -1._id)**7*605.625&
&_id + (Cos(t) + -1._id)*1470._id + (Cos(t) + -1._id)**6*4239.375_id + (Cos(t) + -1._id)**2*6615._id + (Cos(t) + -1._id)**5*12048.75_id + (Cos(t) + -1._id)**3*14700._id + (Cos(t) + -1._id)**4*17850._id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==-3.and.m2==1) then
   wigd=Cos(t*0.5_id)**2*1.4719601443879744757940071211598756606957732468219_id*(330._id + (Cos(t) + -1._id)**7*605.625_&
&id + (Cos(t) + -1._id)*3234._id + (Cos(t) + -1._id)**6*4663.3125_id + (Cos(t) + -1._id)**2*12127.5_id + (Cos(t) + -1._id)**5*14726.25_id + (Cos(t) + -1._id)**3*23100._id + (Cos(t) + -1._id)**4*24543.75_id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==-3.and.m2==2) then
   wigd=Cos(t*0.5_id)*-1.2747548783981962075070560272556954973909427365249_id*((Cos(t) + -1._id)**7*605.625_id + 792._id&
& + (Cos(t) + -1._id)**6*5087.25_id + (Cos(t) + -1._id)*6468._id + (Cos(t) + -1._id)**5*17671.5_id + (Cos(t) + -1._id)**2*20790._id + (Cos(t) + -1._id)**4*32725._id + (Cos(t) + -1._id)**3*34650._id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==-3.and.m2==3) then
   wigd=((Cos(t) + -1._id)**7*605.625_id + 1716._id + (Cos(t) + -1._id)**6*5511.1875_id + (Cos(t) + -1._id)*12012._id + &
&(Cos(t) + -1._id)**5*20884.5_id + (Cos(t) + -1._id)**2*33783.75_id + (Cos(t) + -1._id)**4*42542.5_id + (Cos(t) + -1._id)**3*50050._id)*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==-3.and.m2==4) then
   wigd=Cos(t*0.5_id)*-1.4142135623730950488016887242096980785696718753769_id*((Cos(t) + -1._id)**6*605.625_id + 1716._i&
&d + (Cos(t) + -1._id)**5*4723.875_id + (Cos(t) + -1._id)*9652.5_id + (Cos(t) + -1._id)**4*14917.5_id + (Cos(t) + -1._id)**2*21450._id + (Cos(t) + -1._id)**3*24310._id)*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==-3.and.m2==5) then
   wigd=Cos(t*0.5_id)**2*2.2360679774997896964091736687312762354406183596115_id*((Cos(t) + -1._id)**5*484.5_id + 1287._i&
&d + (Cos(t) + -1._id)**4*3149.25_id + (Cos(t) + -1._id)*5720._id + (Cos(t) + -1._id)**3*7956._id + (Cos(t) + -1._id)**2*9724._id)*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==-3.and.m2==6) then
   wigd=Cos(t*0.5_id)**3*-4._id*((Cos(t) + -1._id)**4*302.8125_id + 715._id + (Cos(t) + -1._id)**3*1574.625_id + (Cos(t)&
& + -1._id)*2431._id + (Cos(t) + -1._id)**2*2983.5_id)*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==-3.and.m2==7) then
   wigd=Cos(t*0.5_id)**4*8.246211251235321099642819711948154050294398450747_id*((Cos(t) + -1._id)**3*142.5_id + 286._id &
&+ (Cos(t) + -1._id)**2*555.75_id + (Cos(t) + -1._id)*702._id)*Sin(t*0.5_id)**10
endif

if(l==10.and.m1==-3.and.m2==8) then
   wigd=Cos(t*0.5_id)**5*-20.199009876724155906726771834139201421207797928959_id*((Cos(t) + -1._id)**2*47.5_id + 78._id &
&+ (Cos(t) + -1._id)*123.5_id)*Sin(t*0.5_id)**11
endif

if(l==10.and.m1==-3.and.m2==9) then
   wigd=Cos(t*0.5_id)**6*(6._id + Cos(t)*20._id)*31.128764832546761362993566721612242497073577045207_id*Sin(t*0.5_id)**1&
&2
endif

if(l==10.and.m1==-3.and.m2==10) then
   wigd=Cos(t*0.5_id)**7*-278.42413688471766545640758032084559883877161886091_id*Sin(t*0.5_id)**13
endif

if(l==10.and.m1==-2.and.m2==-10) then
   wigd=Cos(t*0.5_id)**12*354.92252675760100307924766676167204916378679569317_id*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==-2.and.m2==-9) then
   wigd=Cos(t*0.5_id)**11*(-4._id + Cos(t)*20._id)*39.681544828799193311277116215239828697071216149532_id*Sin(t*0.5_id)*&
&*7
endif

if(l==10.and.m1==-2.and.m2==-8) then
   wigd=Cos(t*0.5_id)**10*25.748786379167465530841591578687506418644843651043_id*(28._id + (Cos(t) + -1._id)**2*47.5_id &
&+ (Cos(t) + -1._id)*76._id)*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==-2.and.m2==-7) then
   wigd=Cos(t*0.5_id)**9*10.5118980208143191442099285287479039691203922519104_id*(56._id + (Cos(t) + -1._id)**3*142.5_id&
& + (Cos(t) + -1._id)*252._id + (Cos(t) + -1._id)**2*342._id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==-2.and.m2==-6) then
   wigd=Cos(t*0.5_id)**8*5.0990195135927848300282241090227819895637709460996_id*(70._id + (Cos(t) + -1._id)**4*302.8125_&
&id + (Cos(t) + -1._id)*476._id + (Cos(t) + -1._id)**3*969._id + (Cos(t) + -1._id)**2*1071._id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==-2.and.m2==-5) then
   wigd=Cos(t*0.5_id)**7*2.850438562747844947840122563916886197690013277291_id*(56._id + (Cos(t) + -1._id)**5*484.5_id +&
& (Cos(t) + -1._id)*560._id + (Cos(t) + -1._id)**2*1904._id + (Cos(t) + -1._id)**4*1938._id + (Cos(t) + -1._id)**3*2856._id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==-2.and.m2==-4) then
   wigd=Cos(t*0.5_id)**6*1.8027756377319946465596106337352479731256482869226_id*(28._id + (Cos(t) + -1._id)*420._id + (C&
&os(t) + -1._id)**6*605.625_id + (Cos(t) + -1._id)**2*2100._id + (Cos(t) + -1._id)**5*2907._id + (Cos(t) + -1._id)**3*4760._id + (Cos(t) + -1._id)**4*5355._id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==-2.and.m2==-3) then
   wigd=Cos(t*0.5_id)**5*1.2747548783981962075070560272556954973909427365249_id*(8._id + (Cos(t) + -1._id)*196._id + (Co&
&s(t) + -1._id)**7*605.625_id + (Cos(t) + -1._id)**2*1470._id + (Cos(t) + -1._id)**6*3391.5_id + (Cos(t) + -1._id)**3*4900._id + (Cos(t) + -1._id)**5*7497._id + (Cos(t) + -1._id)**4*8330._id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==-2.and.m2==-2) then
   wigd=Cos(t*0.5_id)**4*(1._id + (Cos(t) + -1._id)*52._id + (Cos(t) + -1._id)**8*492.0703125_id + (Cos(t) + -1._id)**2*&
&637._id + (Cos(t) + -1._id)**7*3149.25_id + (Cos(t) + -1._id)**3*3185._id + (Cos(t) + -1._id)**4*7962.5_id + (Cos(t) + -1._id)**6*8121.75_id + (Cos(t) + -1._id)**5*10829._id)
endif

if(l==10.and.m1==-2.and.m2==-1) then
   wigd=Cos(t*0.5_id)**3*-1.1547005383792515290182975610039149112952035025403_id*(9._id + (Cos(t) + -1._id)*234._id + (C&
&os(t) + -1._id)**8*492.0703125_id + (Cos(t) + -1._id)**2*1911._id + (Cos(t) + -1._id)**7*3542.90625_id + (Cos(t) + -1._id)**3*7166.25_id + (Cos(t) + -1._id)**6*10442.25_id + (Cos(t) + -1._id)**4*14332.5_id + (Cos(t) + -1._id)**5*16243.5_id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==-2.and.m2==0) then
   wigd=Cos(t*0.5_id)**2*1.21106014163899666616901312388727747967789433339067_id*(45._id + (Cos(t) + -1._id)**8*492.0703&
&125_id + (Cos(t) + -1._id)*780._id + (Cos(t) + -1._id)**7*3936.5625_id + (Cos(t) + -1._id)**2*4777.5_id + (Cos(t) + -1._id)**6*13052.8125_id + (Cos(t) + -1._id)**3*14332.5_id + (Cos(t) + -1._id)**5*23205._id + (Cos(t) + -1._id)**4*23887.5_id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==-2.and.m2==1) then
   wigd=Cos(t*0.5_id)*-1.1547005383792515290182975610039149112952035025403_id*(165._id + (Cos(t) + -1._id)**8*492.070312&
&5_id + (Cos(t) + -1._id)*2145._id + (Cos(t) + -1._id)**7*4330.21875_id + (Cos(t) + -1._id)**2*10510.5_id + &
&(Cos(t) + -1._id)**6*15953.4375_id + (Cos(t) + -1._id)**3*26276.25_id + (Cos(t) + -1._id)**5*31906.875_id + (Cos(t) + -1._id)**4*37537.5_id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==-2.and.m2==2) then
   wigd=((Cos(t) + -1._id)**8*492.0703125_id + 495._id + (Cos(t) + -1._id)**7*4723.875_id + (Cos(t) + -1._id)*5148._id +&
& (Cos(t) + -1._id)**6*19144.125_id + (Cos(t) + -1._id)**2*21021._id + (Cos(t) + -1._id)**5*42542.5_id + (Cos(t) + -1._id)**3*45045._id + (Cos(t) + -1._id)**4*56306.25_id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==-2.and.m2==3) then
   wigd=Cos(t*0.5_id)*-1.2747548783981962075070560272556954973909427365249_id*((Cos(t) + -1._id)**7*605.625_id + 792._id&
& + (Cos(t) + -1._id)**6*5087.25_id + (Cos(t) + -1._id)*6468._id + (Cos(t) + -1._id)**5*17671.5_id + (Cos(t) + -1._id)**2*20790._id + (Cos(t) + -1._id)**4*32725._id + (Cos(t) + -1._id)**3*34650._id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==-2.and.m2==4) then
   wigd=Cos(t*0.5_id)**2*1.8027756377319946465596106337352479731256482869226_id*((Cos(t) + -1._id)**6*605.625_id + 924._&
&id + (Cos(t) + -1._id)**5*4360.5_id + (Cos(t) + -1._id)*5940._id + (Cos(t) + -1._id)**4*12622.5_id + (Cos(t) + -1._id)**2*14850._id + (Cos(t) + -1._id)**3*18700._id)*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==-2.and.m2==5) then
   wigd=Cos(t*0.5_id)**3*-2.850438562747844947840122563916886197690013277291_id*((Cos(t) + -1._id)**5*484.5_id + 792._id&
& + (Cos(t) + -1._id)**4*2907._id + (Cos(t) + -1._id)*3960._id + (Cos(t) + -1._id)**3*6732._id + (Cos(t) + -1._id)**2*7480._id)*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==-2.and.m2==6) then
   wigd=Cos(t*0.5_id)**4*5.0990195135927848300282241090227819895637709460996_id*((Cos(t) + -1._id)**4*302.8125_id + 495.&
&_id + (Cos(t) + -1._id)**3*1453.5_id + (Cos(t) + -1._id)*1870._id + (Cos(t) + -1._id)**2*2524.5_id)*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==-2.and.m2==7) then
   wigd=Cos(t*0.5_id)**5*-10.51189802081431914420992852874790396912039225191_id*((Cos(t) + -1._id)**3*142.5_id + 220._id&
& + (Cos(t) + -1._id)**2*513._id + (Cos(t) + -1._id)*594._id)*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==-2.and.m2==8) then
   wigd=Cos(t*0.5_id)**6*25.748786379167465530841591578687506418644843651043_id*((Cos(t) + -1._id)**2*47.5_id + 66._id +&
& (Cos(t) + -1._id)*114._id)*Sin(t*0.5_id)**10
endif

if(l==10.and.m1==-2.and.m2==9) then
   wigd=Cos(t*0.5_id)**7*-39.681544828799193311277116215239828697071216149532_id*(4._id + Cos(t)*20._id)*Sin(t*0.5_id)**&
&11
endif

if(l==10.and.m1==-2.and.m2==10) then
   wigd=Cos(t*0.5_id)**8*354.92252675760100307924766676167204916378679569317_id*Sin(t*0.5_id)**12
endif

if(l==10.and.m1==-1.and.m2==-10) then
   wigd=Cos(t*0.5_id)**11*409.82923272992618480080474681969664551446288789265_id*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==-1.and.m2==-9) then
   wigd=Cos(t*0.5_id)**10*(-2._id + Cos(t)*20._id)*45.820301177534832960627900345328910787301645567449_id*Sin(t*0.5_id)*&
&*8
endif

if(l==10.and.m1==-1.and.m2==-8) then
   wigd=Cos(t*0.5_id)**9*29.732137494637011045224016427862793302897971027442_id*(36._id + (Cos(t) + -1._id)**2*47.5_id +&
& (Cos(t) + -1._id)*85.5_id)*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==-1.and.m2==-7) then
   wigd=Cos(t*0.5_id)**8*12.138094304022082911201150512925122493534377641012_id*(84._id + (Cos(t) + -1._id)**3*142.5_id &
&+ (Cos(t) + -1._id)*324._id + (Cos(t) + -1._id)**2*384.75_id)*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==-1.and.m2==-6) then
   wigd=Cos(t*0.5_id)**7*5.8878405775518979031760284846395026427830929872876_id*(126._id + (Cos(t) + -1._id)**4*302.8125&
&_id + (Cos(t) + -1._id)*714._id + (Cos(t) + -1._id)**3*1090.125_id + (Cos(t) + -1._id)**2*1377._id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==-1.and.m2==-5) then
   wigd=Cos(t*0.5_id)**6*3.2914029430219165029064101739538498443475436446803_id*(126._id + (Cos(t) + -1._id)**5*484.5_id&
& + (Cos(t) + -1._id)*1008._id + (Cos(t) + -1._id)**4*2180.25_id + (Cos(t) + -1._id)**2*2856._id + (Cos(t) + -1._id)**3*3672._id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==-1.and.m2==-4) then
   wigd=Cos(t*0.5_id)**5*2.0816659994661327352822977069799314870243199926639_id*(84._id + (Cos(t) + -1._id)**6*605.625_i&
&d + (Cos(t) + -1._id)*945._id + (Cos(t) + -1._id)**5*3270.375_id + (Cos(t) + -1._id)**2*3780._id + (Cos(t) + -1._id)**4*6885._id + (Cos(t) + -1._id)**3*7140._id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==-1.and.m2==-3) then
   wigd=Cos(t*0.5_id)**4*1.4719601443879744757940071211598756606957732468219_id*(36._id + (Cos(t) + -1._id)*588._id + (C&
&os(t) + -1._id)**7*605.625_id + (Cos(t) + -1._id)**2*3307.5_id + (Cos(t) + -1._id)**6*3815.4375_id + (Cos(t) + -1._id)**3*8820._id + (Cos(t) + -1._id)**5*9639._id + (Cos(t) + -1._id)**4*12495._id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==-1.and.m2==-2) then
   wigd=Cos(t*0.5_id)**3*1.1547005383792515290182975610039149112952035025403_id*(9._id + (Cos(t) + -1._id)*234._id + (Co&
&s(t) + -1._id)**8*492.0703125_id + (Cos(t) + -1._id)**2*1911._id + (Cos(t) + -1._id)**7*3542.90625_id + (Cos(t) + -1._id)**3*7166.25_id + (Cos(t) + -1._id)**6*10442.25_id + (Cos(t) + -1._id)**4*14332.5_id + (Cos(t) + -1._id)**5*16243.5_id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==-1.and.m2==-1) then
   wigd=Cos(t*0.5_id)**2*(1._id + (Cos(t) + -1._id)*54._id + (Cos(t) + -1._id)**9*328.046875_id + (Cos(t) + -1._id)**2*7&
&02._id + (Cos(t) + -1._id)**8*2657.1796875_id + (Cos(t) + -1._id)**3*3822._id + (Cos(t) + -1._id)**7*8950.5_id + (Cos(t) + -1._id)**4*10749.375_id + (Cos(t) + -1._id)**6*16243.5_id + (Cos(t) + -1._id)**5*17199._id)
endif

if(l==10.and.m1==-1.and.m2==0) then
   wigd=Cos(t*0.5_id)*-1.0488088481701515469914535136799375984752718576815_id*(10._id + (Cos(t) + -1._id)*270._id + (Cos&
&(t) + -1._id)**9*328.046875_id + (Cos(t) + -1._id)**2*2340._id + (Cos(t) + -1._id)**8*2952.421875_id + &
&(Cos(t) + -1._id)**3*9555._id + (Cos(t) + -1._id)**7*11188.125_id + (Cos(t) + -1._id)**4*21498.75_id + (Cos(t) + -1._id)**6*23205._id + (Cos(t) + -1._id)**5*28665._id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==-1.and.m2==1) then
   wigd=(55._id + (Cos(t) + -1._id)**9*328.046875_id + (Cos(t) + -1._id)*990._id + (Cos(t) + -1._id)**8*3247.6640625_id &
&+ (Cos(t) + -1._id)**2*6435._id + (Cos(t) + -1._id)**7*13674.375_id + (Cos(t) + -1._id)**3*21021._id + (Cos(t) + -1._id)**6*31906.875_id + (Cos(t) + -1._id)**4*39414.375_id + (Cos(t) + -1._id)**5*45045._id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==-1.and.m2==2) then
   wigd=Cos(t*0.5_id)*-1.1547005383792515290182975610039149112952035025403_id*(165._id + (Cos(t) + -1._id)**8*492.070312&
&5_id + (Cos(t) + -1._id)*2145._id + (Cos(t) + -1._id)**7*4330.21875_id + (Cos(t) + -1._id)**2*10510.5_id + &
&(Cos(t) + -1._id)**6*15953.4375_id + (Cos(t) + -1._id)**3*26276.25_id + (Cos(t) + -1._id)**5*31906.875_id + (Cos(t) + -1._id)**4*37537.5_id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==-1.and.m2==3) then
   wigd=Cos(t*0.5_id)**2*1.4719601443879744757940071211598756606957732468219_id*(330._id + (Cos(t) + -1._id)**7*605.625_&
&id + (Cos(t) + -1._id)*3234._id + (Cos(t) + -1._id)**6*4663.3125_id + (Cos(t) + -1._id)**2*12127.5_id + (Cos(t) + -1._id)**5*14726.25_id + (Cos(t) + -1._id)**3*23100._id + (Cos(t) + -1._id)**4*24543.75_id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==-1.and.m2==4) then
   wigd=Cos(t*0.5_id)**3*-2.0816659994661327352822977069799314870243199926639_id*(462._id + (Cos(t) + -1._id)**6*605.625&
&_id + (Cos(t) + -1._id)*3465._id + (Cos(t) + -1._id)**5*3997.125_id + (Cos(t) + -1._id)**2*9900._id + (Cos(t) + -1._id)**4*10518.75_id + (Cos(t) + -1._id)**3*14025._id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==-1.and.m2==5) then
   wigd=Cos(t*0.5_id)**4*3.2914029430219165029064101739538498443475436446803_id*(462._id + (Cos(t) + -1._id)**5*484.5_id&
& + (Cos(t) + -1._id)*2640._id + (Cos(t) + -1._id)**4*2664.75_id + (Cos(t) + -1._id)**2*5610._id + (Cos(t) + -1._id)**3*5610._id)*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==-1.and.m2==6) then
   wigd=Cos(t*0.5_id)**5*-5.8878405775518979031760284846395026427830929872876_id*((Cos(t) + -1._id)**4*302.8125_id + 330&
&._id + (Cos(t) + -1._id)**3*1332.375_id + (Cos(t) + -1._id)*1402.5_id + (Cos(t) + -1._id)**2*2103.75_id)*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==-1.and.m2==7) then
   wigd=Cos(t*0.5_id)**6*12.138094304022082911201150512925122493534377641012_id*((Cos(t) + -1._id)**3*142.5_id + 165._id&
& + (Cos(t) + -1._id)**2*470.25_id + (Cos(t) + -1._id)*495._id)*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==-1.and.m2==8) then
   wigd=Cos(t*0.5_id)**7*-29.732137494637011045224016427862793302897971027442_id*((Cos(t) + -1._id)**2*47.5_id + 55._id &
&+ (Cos(t) + -1._id)*104.5_id)*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==-1.and.m2==9) then
   wigd=Cos(t*0.5_id)**8*(2._id + Cos(t)*20._id)*45.820301177534832960627900345328910787301645567449_id*Sin(t*0.5_id)**1&
&0
endif

if(l==10.and.m1==-1.and.m2==10) then
   wigd=Cos(t*0.5_id)**9*-409.82923272992618480080474681969664551446288789265_id*Sin(t*0.5_id)**11
endif

if(l==10.and.m1==0.and.m2==-10) then
   wigd=Cos(t*0.5_id)**10*429.83252552593085495728450959347453264579252131871_id*Sin(t*0.5_id)**10
endif

if(l==10.and.m1==0.and.m2==-9) then
   wigd=Cos(t)*Cos(t*0.5_id)**9*961.13474601639493532560896830571407682331375986_id*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==0.and.m2==-8) then
   wigd=Cos(t*0.5_id)**8*31.183328879386818922579762898076675083852016162476_id*(45._id + (Cos(t) + -1._id)**2*47.5_id +&
& (Cos(t) + -1._id)*95._id)*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==0.and.m2==-7) then
   wigd=Cos(t*0.5_id)**7*12.7305407059820780680148726383177941324863797844472_id*(120._id + (Cos(t) + -1._id)**3*142.5_i&
&d + (Cos(t) + -1._id)*405._id + (Cos(t) + -1._id)**2*427.5_id)*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==0.and.m2==-6) then
   wigd=Cos(t*0.5_id)**6*6.1752192943516858827531745282506443172988423053509_id*(210._id + (Cos(t) + -1._id)**4*302.8125&
&_id + (Cos(t) + -1._id)*1020._id + (Cos(t) + -1._id)**3*1211.25_id + (Cos(t) + -1._id)**2*1721.25_id)*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==0.and.m2==-5) then
   wigd=Cos(t*0.5_id)**5*3.4520525295346631886928627240080226103155020819111_id*(252._id + (Cos(t) + -1._id)**5*484.5_id&
& + (Cos(t) + -1._id)*1680._id + (Cos(t) + -1._id)**4*2422.5_id + (Cos(t) + -1._id)**2*4080._id + (Cos(t) + -1._id)**3*4590._id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==0.and.m2==-4) then
   wigd=Cos(t*0.5_id)**4*2.1832697191750419792351883418789116087243562822933_id*(210._id + (Cos(t) + -1._id)**6*605.625_&
&id + (Cos(t) + -1._id)*1890._id + (Cos(t) + -1._id)**5*3633.75_id + (Cos(t) + -1._id)**2*6300._id + (Cos(t) + -1._id)**4*8606.25_id + (Cos(t) + -1._id)**3*10200._id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==0.and.m2==-3) then
   wigd=Cos(t*0.5_id)**3*1.5438048235879214706882936320626610793247105763377_id*(120._id + (Cos(t) + -1._id)**7*605.625_&
&id + (Cos(t) + -1._id)*1470._id + (Cos(t) + -1._id)**6*4239.375_id + (Cos(t) + -1._id)**2*6615._id + (Cos(t) + -1._id)**5*12048.75_id + (Cos(t) + -1._id)**3*14700._id + (Cos(t) + -1._id)**4*17850._id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==0.and.m2==-2) then
   wigd=Cos(t*0.5_id)**2*1.21106014163899666616901312388727747967789433339067_id*(45._id + (Cos(t) + -1._id)**8*492.0703&
&125_id + (Cos(t) + -1._id)*780._id + (Cos(t) + -1._id)**7*3936.5625_id + (Cos(t) + -1._id)**2*4777.5_id + (Cos(t) + -1._id)**6*13052.8125_id + (Cos(t) + -1._id)**3*14332.5_id + (Cos(t) + -1._id)**5*23205._id + (Cos(t) + -1._id)**4*23887.5_id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==0.and.m2==-1) then
   wigd=Cos(t*0.5_id)*1.0488088481701515469914535136799375984752718576815_id*(10._id + (Cos(t) + -1._id)*270._id + (Cos(&
&t) + -1._id)**9*328.046875_id + (Cos(t) + -1._id)**2*2340._id + (Cos(t) + -1._id)**8*2952.421875_id + &
&(Cos(t) + -1._id)**3*9555._id + (Cos(t) + -1._id)**7*11188.125_id + (Cos(t) + -1._id)**4*21498.75_id + (Cos(t) + -1._id)**6*23205._id + (Cos(t) + -1._id)**5*28665._id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==0.and.m2==0) then
   wigd=0.00390625_id*(Cos(t)**8*-109395._id + Cos(t)**4*-30030._id + -63._id + Cos(t)**2*3465._id + Cos(t)**10*46189._i&
&d + Cos(t)**6*90090._id)
endif

if(l==10.and.m1==0.and.m2==1) then
   wigd=Cos(t*0.5_id)*-1.0488088481701515469914535136799375984752718576815_id*(10._id + (Cos(t) + -1._id)*270._id + (Cos&
&(t) + -1._id)**9*328.046875_id + (Cos(t) + -1._id)**2*2340._id + (Cos(t) + -1._id)**8*2952.421875_id + &
&(Cos(t) + -1._id)**3*9555._id + (Cos(t) + -1._id)**7*11188.125_id + (Cos(t) + -1._id)**4*21498.75_id + (Cos(t) + -1._id)**6*23205._id + (Cos(t) + -1._id)**5*28665._id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==0.and.m2==2) then
   wigd=Cos(t*0.5_id)**2*1.21106014163899666616901312388727747967789433339067_id*(45._id + (Cos(t) + -1._id)**8*492.0703&
&125_id + (Cos(t) + -1._id)*780._id + (Cos(t) + -1._id)**7*3936.5625_id + (Cos(t) + -1._id)**2*4777.5_id + (Cos(t) + -1._id)**6*13052.8125_id + (Cos(t) + -1._id)**3*14332.5_id + (Cos(t) + -1._id)**5*23205._id + (Cos(t) + -1._id)**4*23887.5_id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==0.and.m2==3) then
   wigd=Cos(t*0.5_id)**3*-1.5438048235879214706882936320626610793247105763377_id*(120._id + (Cos(t) + -1._id)**7*605.625&
&_id + (Cos(t) + -1._id)*1470._id + (Cos(t) + -1._id)**6*4239.375_id + (Cos(t) + -1._id)**2*6615._id + (Cos(t) + -1._id)**5*12048.75_id + (Cos(t) + -1._id)**3*14700._id + (Cos(t) + -1._id)**4*17850._id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==0.and.m2==4) then
   wigd=Cos(t*0.5_id)**4*2.1832697191750419792351883418789116087243562822933_id*(210._id + (Cos(t) + -1._id)**6*605.625_&
&id + (Cos(t) + -1._id)*1890._id + (Cos(t) + -1._id)**5*3633.75_id + (Cos(t) + -1._id)**2*6300._id + (Cos(t) + -1._id)**4*8606.25_id + (Cos(t) + -1._id)**3*10200._id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==0.and.m2==5) then
   wigd=Cos(t*0.5_id)**5*-3.4520525295346631886928627240080226103155020819111_id*(252._id + (Cos(t) + -1._id)**5*484.5_i&
&d + (Cos(t) + -1._id)*1680._id + (Cos(t) + -1._id)**4*2422.5_id + (Cos(t) + -1._id)**2*4080._id + (Cos(t) + -1._id)**3*4590._id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==0.and.m2==6) then
   wigd=Cos(t*0.5_id)**6*6.1752192943516858827531745282506443172988423053509_id*(210._id + (Cos(t) + -1._id)**4*302.8125&
&_id + (Cos(t) + -1._id)*1020._id + (Cos(t) + -1._id)**3*1211.25_id + (Cos(t) + -1._id)**2*1721.25_id)*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==0.and.m2==7) then
   wigd=Cos(t*0.5_id)**7*-12.730540705982078068014872638317794132486379784447_id*(120._id + (Cos(t) + -1._id)**3*142.5_i&
&d + (Cos(t) + -1._id)*405._id + (Cos(t) + -1._id)**2*427.5_id)*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==0.and.m2==8) then
   wigd=Cos(t*0.5_id)**8*31.183328879386818922579762898076675083852016162476_id*(45._id + (Cos(t) + -1._id)**2*47.5_id +&
& (Cos(t) + -1._id)*95._id)*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==0.and.m2==9) then
   wigd=Cos(t)*Cos(t*0.5_id)**9*-961.13474601639493532560896830571407682331375986_id*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==0.and.m2==10) then
   wigd=Cos(t*0.5_id)**10*429.83252552593085495728450959347453264579252131871_id*Sin(t*0.5_id)**10
endif

if(l==10.and.m1==1.and.m2==-10) then
   wigd=Cos(t*0.5_id)**9*409.82923272992618480080474681969664551446288789265_id*Sin(t*0.5_id)**11
endif

if(l==10.and.m1==1.and.m2==-9) then
   wigd=Cos(t*0.5_id)**8*(2._id + Cos(t)*20._id)*45.820301177534832960627900345328910787301645567449_id*Sin(t*0.5_id)**1&
&0
endif

if(l==10.and.m1==1.and.m2==-8) then
   wigd=Cos(t*0.5_id)**7*29.732137494637011045224016427862793302897971027442_id*((Cos(t) + -1._id)**2*47.5_id + 55._id +&
& (Cos(t) + -1._id)*104.5_id)*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==1.and.m2==-7) then
   wigd=Cos(t*0.5_id)**6*12.138094304022082911201150512925122493534377641012_id*((Cos(t) + -1._id)**3*142.5_id + 165._id&
& + (Cos(t) + -1._id)**2*470.25_id + (Cos(t) + -1._id)*495._id)*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==1.and.m2==-6) then
   wigd=Cos(t*0.5_id)**5*5.8878405775518979031760284846395026427830929872876_id*((Cos(t) + -1._id)**4*302.8125_id + 330.&
&_id + (Cos(t) + -1._id)**3*1332.375_id + (Cos(t) + -1._id)*1402.5_id + (Cos(t) + -1._id)**2*2103.75_id)*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==1.and.m2==-5) then
   wigd=Cos(t*0.5_id)**4*3.2914029430219165029064101739538498443475436446803_id*(462._id + (Cos(t) + -1._id)**5*484.5_id&
& + (Cos(t) + -1._id)*2640._id + (Cos(t) + -1._id)**4*2664.75_id + (Cos(t) + -1._id)**2*5610._id + (Cos(t) + -1._id)**3*5610._id)*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==1.and.m2==-4) then
   wigd=Cos(t*0.5_id)**3*2.0816659994661327352822977069799314870243199926639_id*(462._id + (Cos(t) + -1._id)**6*605.625_&
&id + (Cos(t) + -1._id)*3465._id + (Cos(t) + -1._id)**5*3997.125_id + (Cos(t) + -1._id)**2*9900._id + (Cos(t) + -1._id)**4*10518.75_id + (Cos(t) + -1._id)**3*14025._id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==1.and.m2==-3) then
   wigd=Cos(t*0.5_id)**2*1.4719601443879744757940071211598756606957732468219_id*(330._id + (Cos(t) + -1._id)**7*605.625_&
&id + (Cos(t) + -1._id)*3234._id + (Cos(t) + -1._id)**6*4663.3125_id + (Cos(t) + -1._id)**2*12127.5_id + (Cos(t) + -1._id)**5*14726.25_id + (Cos(t) + -1._id)**3*23100._id + (Cos(t) + -1._id)**4*24543.75_id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==1.and.m2==-2) then
   wigd=Cos(t*0.5_id)*1.1547005383792515290182975610039149112952035025403_id*(165._id + (Cos(t) + -1._id)**8*492.0703125&
&_id + (Cos(t) + -1._id)*2145._id + (Cos(t) + -1._id)**7*4330.21875_id + (Cos(t) + -1._id)**2*10510.5_id + (Cos(t) + -1._id)**6*15953.4375_id + (Cos(t) + -1._id)**3*26276.25_id + (Cos(t) + -1._id)**5*31906.875_id + (Cos(t) + -1._id)**4*37537.5_id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==1.and.m2==-1) then
   wigd=(55._id + (Cos(t) + -1._id)**9*328.046875_id + (Cos(t) + -1._id)*990._id + (Cos(t) + -1._id)**8*3247.6640625_id &
&+ (Cos(t) + -1._id)**2*6435._id + (Cos(t) + -1._id)**7*13674.375_id + (Cos(t) + -1._id)**3*21021._id + (Cos(t) + -1._id)**6*31906.875_id + (Cos(t) + -1._id)**4*39414.375_id + (Cos(t) + -1._id)**5*45045._id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==1.and.m2==0) then
   wigd=Cos(t*0.5_id)*1.0488088481701515469914535136799375984752718576815_id*(10._id + (Cos(t) + -1._id)*270._id + (Cos(&
&t) + -1._id)**9*328.046875_id + (Cos(t) + -1._id)**2*2340._id + (Cos(t) + -1._id)**8*2952.421875_id + (Cos(t) +&
& -1._id)**3*9555._id + (Cos(t) + -1._id)**7*11188.125_id + (Cos(t) + -1._id)**4*21498.75_id + (Cos(t) + -1._id)**6*23205._id + (Cos(t) + -1._id)**5*28665._id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==1.and.m2==1) then
   wigd=Cos(t*0.5_id)**2*(1._id + (Cos(t) + -1._id)*54._id + (Cos(t) + -1._id)**9*328.046875_id + (Cos(t) + -1._id)**2*7&
&02._id + (Cos(t) + -1._id)**8*2657.1796875_id + (Cos(t) + -1._id)**3*3822._id + (Cos(t) + -1._id)**7*8950.5_id + (Cos(t) + -1._id)**4*10749.375_id + (Cos(t) + -1._id)**6*16243.5_id + (Cos(t) + -1._id)**5*17199._id)
endif

if(l==10.and.m1==1.and.m2==2) then
   wigd=Cos(t*0.5_id)**3*-1.1547005383792515290182975610039149112952035025403_id*(9._id + (Cos(t) + -1._id)*234._id + (C&
&os(t) + -1._id)**8*492.0703125_id + (Cos(t) + -1._id)**2*1911._id + (Cos(t) + -1._id)**7*3542.90625_id + (Cos(t) + -1._id)**3*7166.25_id + (Cos(t) + -1._id)**6*10442.25_id + (Cos(t) + -1._id)**4*14332.5_id + (Cos(t) + -1._id)**5*16243.5_id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==1.and.m2==3) then
   wigd=Cos(t*0.5_id)**4*1.4719601443879744757940071211598756606957732468219_id*(36._id + (Cos(t) + -1._id)*588._id + (C&
&os(t) + -1._id)**7*605.625_id + (Cos(t) + -1._id)**2*3307.5_id + (Cos(t) + -1._id)**6*3815.4375_id + (Cos(t) + -1._id)**3*8820._id + (Cos(t) + -1._id)**5*9639._id + (Cos(t) + -1._id)**4*12495._id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==1.and.m2==4) then
   wigd=Cos(t*0.5_id)**5*-2.0816659994661327352822977069799314870243199926639_id*(84._id + (Cos(t) + -1._id)**6*605.625_&
&id + (Cos(t) + -1._id)*945._id + (Cos(t) + -1._id)**5*3270.375_id + (Cos(t) + -1._id)**2*3780._id + (Cos(t) + -1._id)**4*6885._id + (Cos(t) + -1._id)**3*7140._id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==1.and.m2==5) then
   wigd=Cos(t*0.5_id)**6*3.2914029430219165029064101739538498443475436446803_id*(126._id + (Cos(t) + -1._id)**5*484.5_id&
& + (Cos(t) + -1._id)*1008._id + (Cos(t) + -1._id)**4*2180.25_id + (Cos(t) + -1._id)**2*2856._id + (Cos(t) + -1._id)**3*3672._id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==1.and.m2==6) then
   wigd=Cos(t*0.5_id)**7*-5.8878405775518979031760284846395026427830929872876_id*(126._id + (Cos(t) + -1._id)**4*302.812&
&5_id + (Cos(t) + -1._id)*714._id + (Cos(t) + -1._id)**3*1090.125_id + (Cos(t) + -1._id)**2*1377._id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==1.and.m2==7) then
   wigd=Cos(t*0.5_id)**8*12.138094304022082911201150512925122493534377641012_id*(84._id + (Cos(t) + -1._id)**3*142.5_id &
&+ (Cos(t) + -1._id)*324._id + (Cos(t) + -1._id)**2*384.75_id)*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==1.and.m2==8) then
   wigd=Cos(t*0.5_id)**9*-29.732137494637011045224016427862793302897971027442_id*(36._id + (Cos(t) + -1._id)**2*47.5_id &
&+ (Cos(t) + -1._id)*85.5_id)*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==1.and.m2==9) then
   wigd=Cos(t*0.5_id)**10*(-2._id + Cos(t)*20._id)*45.820301177534832960627900345328910787301645567449_id*Sin(t*0.5_id)*&
&*8
endif

if(l==10.and.m1==1.and.m2==10) then
   wigd=Cos(t*0.5_id)**11*-409.82923272992618480080474681969664551446288789265_id*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==2.and.m2==-10) then
   wigd=Cos(t*0.5_id)**8*354.92252675760100307924766676167204916378679569317_id*Sin(t*0.5_id)**12
endif

if(l==10.and.m1==2.and.m2==-9) then
   wigd=Cos(t*0.5_id)**7*(4._id + Cos(t)*20._id)*39.681544828799193311277116215239828697071216149532_id*Sin(t*0.5_id)**1&
&1
endif

if(l==10.and.m1==2.and.m2==-8) then
   wigd=Cos(t*0.5_id)**6*25.748786379167465530841591578687506418644843651043_id*((Cos(t) + -1._id)**2*47.5_id + 66._id +&
& (Cos(t) + -1._id)*114._id)*Sin(t*0.5_id)**10
endif

if(l==10.and.m1==2.and.m2==-7) then
   wigd=Cos(t*0.5_id)**5*10.5118980208143191442099285287479039691203922519104_id*((Cos(t) + -1._id)**3*142.5_id + 220._i&
&d + (Cos(t) + -1._id)**2*513._id + (Cos(t) + -1._id)*594._id)*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==2.and.m2==-6) then
   wigd=Cos(t*0.5_id)**4*5.0990195135927848300282241090227819895637709460996_id*((Cos(t) + -1._id)**4*302.8125_id + 495.&
&_id + (Cos(t) + -1._id)**3*1453.5_id + (Cos(t) + -1._id)*1870._id + (Cos(t) + -1._id)**2*2524.5_id)*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==2.and.m2==-5) then
   wigd=Cos(t*0.5_id)**3*2.850438562747844947840122563916886197690013277291_id*((Cos(t) + -1._id)**5*484.5_id + 792._id &
&+ (Cos(t) + -1._id)**4*2907._id + (Cos(t) + -1._id)*3960._id + (Cos(t) + -1._id)**3*6732._id + (Cos(t) + -1._id)**2*7480._id)*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==2.and.m2==-4) then
   wigd=Cos(t*0.5_id)**2*1.8027756377319946465596106337352479731256482869226_id*((Cos(t) + -1._id)**6*605.625_id + 924._&
&id + (Cos(t) + -1._id)**5*4360.5_id + (Cos(t) + -1._id)*5940._id + (Cos(t) + -1._id)**4*12622.5_id + (Cos(t) + -1._id)**2*14850._id + (Cos(t) + -1._id)**3*18700._id)*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==2.and.m2==-3) then
   wigd=Cos(t*0.5_id)*1.2747548783981962075070560272556954973909427365249_id*((Cos(t) + -1._id)**7*605.625_id + 792._id &
&+ (Cos(t) + -1._id)**6*5087.25_id + (Cos(t) + -1._id)*6468._id + (Cos(t) + -1._id)**5*17671.5_id + (Cos(t) + -1._id)**2*20790._id + (Cos(t) + -1._id)**4*32725._id + (Cos(t) + -1._id)**3*34650._id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==2.and.m2==-2) then
   wigd=((Cos(t) + -1._id)**8*492.0703125_id + 495._id + (Cos(t) + -1._id)**7*4723.875_id + (Cos(t) + -1._id)*5148._id +&
& (Cos(t) + -1._id)**6*19144.125_id + (Cos(t) + -1._id)**2*21021._id + (Cos(t) + -1._id)**5*42542.5_id + (Cos(t) + -1._id)**3*45045._id + (Cos(t) + -1._id)**4*56306.25_id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==2.and.m2==-1) then
   wigd=Cos(t*0.5_id)*1.1547005383792515290182975610039149112952035025403_id*(165._id + (Cos(t) + -1._id)**8*492.0703125&
&_id + (Cos(t) + -1._id)*2145._id + (Cos(t) + -1._id)**7*4330.21875_id + (Cos(t) + -1._id)**2*10510.5_id + (Cos(t) + -1._id)**6*15953.4375_id + (Cos(t) + -1._id)**3*26276.25_id + (Cos(t) + -1._id)**5*31906.875_id + (Cos(t) + -1._id)**4*37537.5_id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==2.and.m2==0) then
   wigd=Cos(t*0.5_id)**2*1.21106014163899666616901312388727747967789433339067_id*(45._id + (Cos(t) + -1._id)**8*492.0703&
&125_id + (Cos(t) + -1._id)*780._id + (Cos(t) + -1._id)**7*3936.5625_id + (Cos(t) + -1._id)**2*4777.5_id + (Cos(t) + -1._id)**6*13052.8125_id + (Cos(t) + -1._id)**3*14332.5_id + (Cos(t) + -1._id)**5*23205._id + (Cos(t) + -1._id)**4*23887.5_id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==2.and.m2==1) then
   wigd=Cos(t*0.5_id)**3*1.1547005383792515290182975610039149112952035025403_id*(9._id + (Cos(t) + -1._id)*234._id + (Co&
&s(t) + -1._id)**8*492.0703125_id + (Cos(t) + -1._id)**2*1911._id + (Cos(t) + -1._id)**7*3542.90625_id + (Cos(t) + -1._id)**3*7166.25_id + (Cos(t) + -1._id)**6*10442.25_id + (Cos(t) + -1._id)**4*14332.5_id + (Cos(t) + -1._id)**5*16243.5_id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==2.and.m2==2) then
   wigd=Cos(t*0.5_id)**4*(1._id + (Cos(t) + -1._id)*52._id + (Cos(t) + -1._id)**8*492.0703125_id + (Cos(t) + -1._id)**2*&
&637._id + (Cos(t) + -1._id)**7*3149.25_id + (Cos(t) + -1._id)**3*3185._id + (Cos(t) + -1._id)**4*7962.5_id + (Cos(t) + -1._id)**6*8121.75_id + (Cos(t) + -1._id)**5*10829._id)
endif

if(l==10.and.m1==2.and.m2==3) then
   wigd=Cos(t*0.5_id)**5*-1.2747548783981962075070560272556954973909427365249_id*(8._id + (Cos(t) + -1._id)*196._id + (C&
&os(t) + -1._id)**7*605.625_id + (Cos(t) + -1._id)**2*1470._id + (Cos(t) + -1._id)**6*3391.5_id + (Cos(t) + -1._id)**3*4900._id + (Cos(t) + -1._id)**5*7497._id + (Cos(t) + -1._id)**4*8330._id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==2.and.m2==4) then
   wigd=Cos(t*0.5_id)**6*1.8027756377319946465596106337352479731256482869226_id*(28._id + (Cos(t) + -1._id)*420._id + (C&
&os(t) + -1._id)**6*605.625_id + (Cos(t) + -1._id)**2*2100._id + (Cos(t) + -1._id)**5*2907._id + (Cos(t) + -1._id)**3*4760._id + (Cos(t) + -1._id)**4*5355._id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==2.and.m2==5) then
   wigd=Cos(t*0.5_id)**7*-2.850438562747844947840122563916886197690013277291_id*(56._id + (Cos(t) + -1._id)**5*484.5_id &
&+ (Cos(t) + -1._id)*560._id + (Cos(t) + -1._id)**2*1904._id + (Cos(t) + -1._id)**4*1938._id + (Cos(t) + -1._id)**3*2856._id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==2.and.m2==6) then
   wigd=Cos(t*0.5_id)**8*5.0990195135927848300282241090227819895637709460996_id*(70._id + (Cos(t) + -1._id)**4*302.8125_&
&id + (Cos(t) + -1._id)*476._id + (Cos(t) + -1._id)**3*969._id + (Cos(t) + -1._id)**2*1071._id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==2.and.m2==7) then
   wigd=Cos(t*0.5_id)**9*-10.51189802081431914420992852874790396912039225191_id*(56._id + (Cos(t) + -1._id)**3*142.5_id &
&+ (Cos(t) + -1._id)*252._id + (Cos(t) + -1._id)**2*342._id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==2.and.m2==8) then
   wigd=Cos(t*0.5_id)**10*25.748786379167465530841591578687506418644843651043_id*(28._id + (Cos(t) + -1._id)**2*47.5_id &
&+ (Cos(t) + -1._id)*76._id)*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==2.and.m2==9) then
   wigd=Cos(t*0.5_id)**11*-39.681544828799193311277116215239828697071216149532_id*(-4._id + Cos(t)*20._id)*Sin(t*0.5_id)&
&**7
endif

if(l==10.and.m1==2.and.m2==10) then
   wigd=Cos(t*0.5_id)**12*354.92252675760100307924766676167204916378679569317_id*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==3.and.m2==-10) then
   wigd=Cos(t*0.5_id)**7*278.42413688471766545640758032084559883877161886091_id*Sin(t*0.5_id)**13
endif

if(l==10.and.m1==3.and.m2==-9) then
   wigd=Cos(t*0.5_id)**6*(6._id + Cos(t)*20._id)*31.128764832546761362993566721612242497073577045207_id*Sin(t*0.5_id)**1&
&2
endif

if(l==10.and.m1==3.and.m2==-8) then
   wigd=Cos(t*0.5_id)**5*20.199009876724155906726771834139201421207797928959_id*((Cos(t) + -1._id)**2*47.5_id + 78._id +&
& (Cos(t) + -1._id)*123.5_id)*Sin(t*0.5_id)**11
endif

if(l==10.and.m1==3.and.m2==-7) then
   wigd=Cos(t*0.5_id)**4*8.246211251235321099642819711948154050294398450747_id*((Cos(t) + -1._id)**3*142.5_id + 286._id &
&+ (Cos(t) + -1._id)**2*555.75_id + (Cos(t) + -1._id)*702._id)*Sin(t*0.5_id)**10
endif

if(l==10.and.m1==3.and.m2==-6) then
   wigd=Cos(t*0.5_id)**3*4._id*((Cos(t) + -1._id)**4*302.8125_id + 715._id + (Cos(t) + -1._id)**3*1574.625_id + (Cos(t) &
&+ -1._id)*2431._id + (Cos(t) + -1._id)**2*2983.5_id)*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==3.and.m2==-5) then
   wigd=Cos(t*0.5_id)**2*2.2360679774997896964091736687312762354406183596115_id*((Cos(t) + -1._id)**5*484.5_id + 1287._i&
&d + (Cos(t) + -1._id)**4*3149.25_id + (Cos(t) + -1._id)*5720._id + (Cos(t) + -1._id)**3*7956._id + (Cos(t) + -1._id)**2*9724._id)*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==3.and.m2==-4) then
   wigd=Cos(t*0.5_id)*1.41421356237309504880168872420969807856967187537695_id*((Cos(t) + -1._id)**6*605.625_id + 1716._i&
&d + (Cos(t) + -1._id)**5*4723.875_id + (Cos(t) + -1._id)*9652.5_id + (Cos(t) + -1._id)**4*14917.5_id + (Cos(t) + -1._id)**2*21450._id + (Cos(t) + -1._id)**3*24310._id)*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==3.and.m2==-3) then
   wigd=((Cos(t) + -1._id)**7*605.625_id + 1716._id + (Cos(t) + -1._id)**6*5511.1875_id + (Cos(t) + -1._id)*12012._id + &
&(Cos(t) + -1._id)**5*20884.5_id + (Cos(t) + -1._id)**2*33783.75_id + (Cos(t) + -1._id)**4*42542.5_id + (Cos(t) + -1._id)**3*50050._id)*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==3.and.m2==-2) then
   wigd=Cos(t*0.5_id)*1.2747548783981962075070560272556954973909427365249_id*((Cos(t) + -1._id)**7*605.625_id + 792._id &
&+ (Cos(t) + -1._id)**6*5087.25_id + (Cos(t) + -1._id)*6468._id + (Cos(t) + -1._id)**5*17671.5_id + (Cos(t) + -1._id)**2*20790._id + (Cos(t) + -1._id)**4*32725._id + (Cos(t) + -1._id)**3*34650._id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==3.and.m2==-1) then
   wigd=Cos(t*0.5_id)**2*1.4719601443879744757940071211598756606957732468219_id*(330._id + (Cos(t) + -1._id)**7*605.625_&
&id + (Cos(t) + -1._id)*3234._id + (Cos(t) + -1._id)**6*4663.3125_id + (Cos(t) + -1._id)**2*12127.5_id + (Cos(t) + -1._id)**5*14726.25_id + (Cos(t) + -1._id)**3*23100._id + (Cos(t) + -1._id)**4*24543.75_id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==3.and.m2==0) then
   wigd=Cos(t*0.5_id)**3*1.5438048235879214706882936320626610793247105763377_id*(120._id + (Cos(t) + -1._id)**7*605.625_&
&id + (Cos(t) + -1._id)*1470._id + (Cos(t) + -1._id)**6*4239.375_id + (Cos(t) + -1._id)**2*6615._id + (Cos(t) + -1._id)**5*12048.75_id + (Cos(t) + -1._id)**3*14700._id + (Cos(t) + -1._id)**4*17850._id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==3.and.m2==1) then
   wigd=Cos(t*0.5_id)**4*1.4719601443879744757940071211598756606957732468219_id*(36._id + (Cos(t) + -1._id)*588._id + (C&
&os(t) + -1._id)**7*605.625_id + (Cos(t) + -1._id)**2*3307.5_id + (Cos(t) + -1._id)**6*3815.4375_id + (Cos(t) + -1._id)**3*8820._id + (Cos(t) + -1._id)**5*9639._id + (Cos(t) + -1._id)**4*12495._id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==3.and.m2==2) then
   wigd=Cos(t*0.5_id)**5*1.2747548783981962075070560272556954973909427365249_id*(8._id + (Cos(t) + -1._id)*196._id + (Co&
&s(t) + -1._id)**7*605.625_id + (Cos(t) + -1._id)**2*1470._id + (Cos(t) + -1._id)**6*3391.5_id + (Cos(t) + -1._id)**3*4900._id + (Cos(t) + -1._id)**5*7497._id + (Cos(t) + -1._id)**4*8330._id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==3.and.m2==3) then
   wigd=Cos(t*0.5_id)**6*(1._id + (Cos(t) + -1._id)*49._id + (Cos(t) + -1._id)**2*551.25_id + (Cos(t) + -1._id)**7*605.6&
&25_id + (Cos(t) + -1._id)**3*2450._id + (Cos(t) + -1._id)**6*2967.5625_id + (Cos(t) + -1._id)**4*5206.25_id + (Cos(t) + -1._id)**5*5622.75_id)
endif

if(l==10.and.m1==3.and.m2==4) then
   wigd=Cos(t*0.5_id)**7*-1.4142135623730950488016887242096980785696718753769_id*(7._id + (Cos(t) + -1._id)*157.5_id + (&
&Cos(t) + -1._id)**6*605.625_id + (Cos(t) + -1._id)**2*1050._id + (Cos(t) + -1._id)**5*2543.625_id + (Cos(t) + -1._id)**3*2975._id + (Cos(t) + -1._id)**4*4016.25_id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==3.and.m2==5) then
   wigd=Cos(t*0.5_id)**8*2.2360679774997896964091736687312762354406183596115_id*(21._id + (Cos(t) + -1._id)*280._id + (C&
&os(t) + -1._id)**5*484.5_id + (Cos(t) + -1._id)**2*1190._id + (Cos(t) + -1._id)**4*1695.75_id + (Cos(t) + -1._id)**3*2142._id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==3.and.m2==6) then
   wigd=Cos(t*0.5_id)**9*-4._id*(35._id + (Cos(t) + -1._id)*297.5_id + (Cos(t) + -1._id)**4*302.8125_id + (Cos(t) + -1._&
&id)**2*803.25_id + (Cos(t) + -1._id)**3*847.875_id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==3.and.m2==7) then
   wigd=Cos(t*0.5_id)**10*8.246211251235321099642819711948154050294398450747_id*(35._id + (Cos(t) + -1._id)**3*142.5_id &
&+ (Cos(t) + -1._id)*189._id + (Cos(t) + -1._id)**2*299.25_id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==3.and.m2==8) then
   wigd=Cos(t*0.5_id)**11*-20.199009876724155906726771834139201421207797928959_id*(21._id + (Cos(t) + -1._id)**2*47.5_id&
& + (Cos(t) + -1._id)*66.5_id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==3.and.m2==9) then
   wigd=Cos(t*0.5_id)**12*(-6._id + Cos(t)*20._id)*31.128764832546761362993566721612242497073577045207_id*Sin(t*0.5_id)*&
&*6
endif

if(l==10.and.m1==3.and.m2==10) then
   wigd=Cos(t*0.5_id)**13*-278.42413688471766545640758032084559883877161886091_id*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==4.and.m2==-10) then
   wigd=Cos(t*0.5_id)**6*196.87559523719540998400115941877306606266144591633_id*Sin(t*0.5_id)**14
endif

if(l==10.and.m1==4.and.m2==-9) then
   wigd=Cos(t*0.5_id)**5*(8._id + Cos(t)*20._id)*22.011360703055138476529216317127335502674181544772_id*Sin(t*0.5_id)**1&
&3
endif

if(l==10.and.m1==4.and.m2==-8) then
   wigd=Cos(t*0.5_id)**4*14.282856857085699995998799622734530557532342319805_id*((Cos(t) + -1._id)**2*47.5_id + 91._id +&
& (Cos(t) + -1._id)*133._id)*Sin(t*0.5_id)**12
endif

if(l==10.and.m1==4.and.m2==-7) then
   wigd=Cos(t*0.5_id)**3*5.830951894845300470874152877545583076521398334886_id*((Cos(t) + -1._id)**3*142.5_id + 364._id &
&+ (Cos(t) + -1._id)**2*598.5_id + (Cos(t) + -1._id)*819._id)*Sin(t*0.5_id)**11
endif

if(l==10.and.m1==4.and.m2==-6) then
   wigd=Cos(t*0.5_id)**2*2.8284271247461900976033774484193961571393437507539_id*((Cos(t) + -1._id)**4*302.8125_id + 1001&
&._id + (Cos(t) + -1._id)**3*1695.75_id + (Cos(t) + -1._id)*3094._id + (Cos(t) + -1._id)**2*3480.75_id)*Sin(t*0.5_id)**10
endif

if(l==10.and.m1==4.and.m2==-5) then
   wigd=Cos(t*0.5_id)*1.58113883008418966599944677221635926685977756966261_id*((Cos(t) + -1._id)**5*484.5_id + 2002._id &
&+ (Cos(t) + -1._id)**4*3391.5_id + (Cos(t) + -1._id)*8008._id + (Cos(t) + -1._id)**3*9282._id + (Cos(t) + -1._id)**2*12376._id)*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==4.and.m2==-4) then
   wigd=((Cos(t) + -1._id)**6*605.625_id + 3003._id + (Cos(t) + -1._id)**5*5087.25_id + (Cos(t) + -1._id)*15015._id + (C&
&os(t) + -1._id)**4*17403.75_id + (Cos(t) + -1._id)**2*30030._id + (Cos(t) + -1._id)**3*30940._id)*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==4.and.m2==-3) then
   wigd=Cos(t*0.5_id)*1.41421356237309504880168872420969807856967187537695_id*((Cos(t) + -1._id)**6*605.625_id + 1716._i&
&d + (Cos(t) + -1._id)**5*4723.875_id + (Cos(t) + -1._id)*9652.5_id + (Cos(t) + -1._id)**4*14917.5_id + (Cos(t) + -1._id)**2*21450._id + (Cos(t) + -1._id)**3*24310._id)*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==4.and.m2==-2) then
   wigd=Cos(t*0.5_id)**2*1.8027756377319946465596106337352479731256482869226_id*((Cos(t) + -1._id)**6*605.625_id + 924._&
&id + (Cos(t) + -1._id)**5*4360.5_id + (Cos(t) + -1._id)*5940._id + (Cos(t) + -1._id)**4*12622.5_id + (Cos(t) + -1._id)**2*14850._id + (Cos(t) + -1._id)**3*18700._id)*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==4.and.m2==-1) then
   wigd=Cos(t*0.5_id)**3*2.0816659994661327352822977069799314870243199926639_id*(462._id + (Cos(t) + -1._id)**6*605.625_&
&id + (Cos(t) + -1._id)*3465._id + (Cos(t) + -1._id)**5*3997.125_id + (Cos(t) + -1._id)**2*9900._id + (Cos(t) + -1._id)**4*10518.75_id + (Cos(t) + -1._id)**3*14025._id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==4.and.m2==0) then
   wigd=Cos(t*0.5_id)**4*2.1832697191750419792351883418789116087243562822933_id*(210._id + (Cos(t) + -1._id)**6*605.625_&
&id + (Cos(t) + -1._id)*1890._id + (Cos(t) + -1._id)**5*3633.75_id + (Cos(t) + -1._id)**2*6300._id + (Cos(t) + -1._id)**4*8606.25_id + (Cos(t) + -1._id)**3*10200._id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==4.and.m2==1) then
   wigd=Cos(t*0.5_id)**5*2.0816659994661327352822977069799314870243199926639_id*(84._id + (Cos(t) + -1._id)**6*605.625_i&
&d + (Cos(t) + -1._id)*945._id + (Cos(t) + -1._id)**5*3270.375_id + (Cos(t) + -1._id)**2*3780._id + (Cos(t) + -1._id)**4*6885._id + (Cos(t) + -1._id)**3*7140._id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==4.and.m2==2) then
   wigd=Cos(t*0.5_id)**6*1.8027756377319946465596106337352479731256482869226_id*(28._id + (Cos(t) + -1._id)*420._id + (C&
&os(t) + -1._id)**6*605.625_id + (Cos(t) + -1._id)**2*2100._id + (Cos(t) + -1._id)**5*2907._id + (Cos(t) + -1._id)**3*4760._id + (Cos(t) + -1._id)**4*5355._id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==4.and.m2==3) then
   wigd=Cos(t*0.5_id)**7*1.41421356237309504880168872420969807856967187537695_id*(7._id + (Cos(t) + -1._id)*157.5_id + (&
&Cos(t) + -1._id)**6*605.625_id + (Cos(t) + -1._id)**2*1050._id + (Cos(t) + -1._id)**5*2543.625_id + (Cos(t) + -1._id)**3*2975._id + (Cos(t) + -1._id)**4*4016.25_id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==4.and.m2==4) then
   wigd=Cos(t*0.5_id)**8*(1._id + (Cos(t) + -1._id)*45._id + (Cos(t) + -1._id)**2*450._id + (Cos(t) + -1._id)**6*605.625&
&_id + (Cos(t) + -1._id)**3*1700._id + (Cos(t) + -1._id)**5*2180.25_id + (Cos(t) + -1._id)**4*2868.75_id)
endif

if(l==10.and.m1==4.and.m2==5) then
   wigd=Cos(t*0.5_id)**9*-1.5811388300841896659994467722163592668597775696626_id*(6._id + (Cos(t) + -1._id)*120._id + (C&
&os(t) + -1._id)**5*484.5_id + (Cos(t) + -1._id)**2*680._id + (Cos(t) + -1._id)**4*1453.5_id + (Cos(t) + -1._id)**3*1530._id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==4.and.m2==6) then
   wigd=Cos(t*0.5_id)**10*2.8284271247461900976033774484193961571393437507539_id*(15._id + (Cos(t) + -1._id)*170._id + (&
&Cos(t) + -1._id)**4*302.8125_id + (Cos(t) + -1._id)**2*573.75_id + (Cos(t) + -1._id)**3*726.75_id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==4.and.m2==7) then
   wigd=Cos(t*0.5_id)**11*-5.830951894845300470874152877545583076521398334886_id*(20._id + (Cos(t) + -1._id)*135._id + (&
&Cos(t) + -1._id)**3*142.5_id + (Cos(t) + -1._id)**2*256.5_id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==4.and.m2==8) then
   wigd=Cos(t*0.5_id)**12*14.282856857085699995998799622734530557532342319805_id*(15._id + (Cos(t) + -1._id)**2*47.5_id &
&+ (Cos(t) + -1._id)*57._id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==4.and.m2==9) then
   wigd=Cos(t*0.5_id)**13*-22.011360703055138476529216317127335502674181544772_id*(-8._id + Cos(t)*20._id)*Sin(t*0.5_id)&
&**5
endif

if(l==10.and.m1==4.and.m2==10) then
   wigd=Cos(t*0.5_id)**14*196.87559523719540998400115941877306606266144591633_id*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==5.and.m2==-10) then
   wigd=Cos(t*0.5_id)**5*124.51505933018704545197426688644896998829430818083_id*Sin(t*0.5_id)**15
endif

if(l==10.and.m1==5.and.m2==-9) then
   wigd=Cos(t*0.5_id)**4*13.9212068442358832728203790160422799419385809430456_id*(10._id + Cos(t)*20._id)*Sin(t*0.5_id)*&
&*14
endif

if(l==10.and.m1==5.and.m2==-8) then
   wigd=Cos(t*0.5_id)**3*9.033271832508971939888199943876461581439881615896_id*((Cos(t) + -1._id)**2*47.5_id + 105._id +&
& (Cos(t) + -1._id)*142.5_id)*Sin(t*0.5_id)**13
endif

if(l==10.and.m1==5.and.m2==-7) then
   wigd=Cos(t*0.5_id)**2*3.687817782917154924000909712705117262898722019489_id*((Cos(t) + -1._id)**3*142.5_id + 455._id &
&+ (Cos(t) + -1._id)**2*641.25_id + (Cos(t) + -1._id)*945._id)*Sin(t*0.5_id)**12
endif

if(l==10.and.m1==5.and.m2==-6) then
   wigd=Cos(t*0.5_id)*1.7888543819998317571273389349850209883524946876892_id*((Cos(t) + -1._id)**4*302.8125_id + 1365._i&
&d + (Cos(t) + -1._id)**3*1816.875_id + (Cos(t) + -1._id)*3867.5_id + (Cos(t) + -1._id)**2*4016.25_id)*Sin(t*0.5_id)**11
endif

if(l==10.and.m1==5.and.m2==-5) then
   wigd=((Cos(t) + -1._id)**5*484.5_id + 3003._id + (Cos(t) + -1._id)**4*3633.75_id + (Cos(t) + -1._id)**3*10710._id + (&
&Cos(t) + -1._id)*10920._id + (Cos(t) + -1._id)**2*15470._id)*Sin(t*0.5_id)**10
endif

if(l==10.and.m1==5.and.m2==-4) then
   wigd=Cos(t*0.5_id)*1.58113883008418966599944677221635926685977756966261_id*((Cos(t) + -1._id)**5*484.5_id + 2002._id &
&+ (Cos(t) + -1._id)**4*3391.5_id + (Cos(t) + -1._id)*8008._id + (Cos(t) + -1._id)**3*9282._id + (Cos(t) + -1._id)**2*12376._id)*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==5.and.m2==-3) then
   wigd=Cos(t*0.5_id)**2*2.2360679774997896964091736687312762354406183596115_id*((Cos(t) + -1._id)**5*484.5_id + 1287._i&
&d + (Cos(t) + -1._id)**4*3149.25_id + (Cos(t) + -1._id)*5720._id + (Cos(t) + -1._id)**3*7956._id + (Cos(t) + -1._id)**2*9724._id)*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==5.and.m2==-2) then
   wigd=Cos(t*0.5_id)**3*2.850438562747844947840122563916886197690013277291_id*((Cos(t) + -1._id)**5*484.5_id + 792._id &
&+ (Cos(t) + -1._id)**4*2907._id + (Cos(t) + -1._id)*3960._id + (Cos(t) + -1._id)**3*6732._id + (Cos(t) + -1._id)**2*7480._id)*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==5.and.m2==-1) then
   wigd=Cos(t*0.5_id)**4*3.2914029430219165029064101739538498443475436446803_id*(462._id + (Cos(t) + -1._id)**5*484.5_id&
& + (Cos(t) + -1._id)*2640._id + (Cos(t) + -1._id)**4*2664.75_id + (Cos(t) + -1._id)**2*5610._id + (Cos(t) + -1._id)**3*5610._id)*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==5.and.m2==0) then
   wigd=Cos(t*0.5_id)**5*3.4520525295346631886928627240080226103155020819111_id*(252._id + (Cos(t) + -1._id)**5*484.5_id&
& + (Cos(t) + -1._id)*1680._id + (Cos(t) + -1._id)**4*2422.5_id + (Cos(t) + -1._id)**2*4080._id + (Cos(t) + -1._id)**3*4590._id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==5.and.m2==1) then
   wigd=Cos(t*0.5_id)**6*3.2914029430219165029064101739538498443475436446803_id*(126._id + (Cos(t) + -1._id)**5*484.5_id&
& + (Cos(t) + -1._id)*1008._id + (Cos(t) + -1._id)**4*2180.25_id + (Cos(t) + -1._id)**2*2856._id + (Cos(t) + -1._id)**3*3672._id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==5.and.m2==2) then
   wigd=Cos(t*0.5_id)**7*2.850438562747844947840122563916886197690013277291_id*(56._id + (Cos(t) + -1._id)**5*484.5_id +&
& (Cos(t) + -1._id)*560._id + (Cos(t) + -1._id)**2*1904._id + (Cos(t) + -1._id)**4*1938._id + (Cos(t) + -1._id)**3*2856._id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==5.and.m2==3) then
   wigd=Cos(t*0.5_id)**8*2.2360679774997896964091736687312762354406183596115_id*(21._id + (Cos(t) + -1._id)*280._id + (C&
&os(t) + -1._id)**5*484.5_id + (Cos(t) + -1._id)**2*1190._id + (Cos(t) + -1._id)**4*1695.75_id + (Cos(t) + -1._id)**3*2142._id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==5.and.m2==4) then
   wigd=Cos(t*0.5_id)**9*1.58113883008418966599944677221635926685977756966261_id*(6._id + (Cos(t) + -1._id)*120._id + (C&
&os(t) + -1._id)**5*484.5_id + (Cos(t) + -1._id)**2*680._id + (Cos(t) + -1._id)**4*1453.5_id + (Cos(t) + -1._id)**3*1530._id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==5.and.m2==5) then
   wigd=Cos(t*0.5_id)**10*(1._id + (Cos(t) + -1._id)*40._id + (Cos(t) + -1._id)**2*340._id + (Cos(t) + -1._id)**5*484.5_&
&id + (Cos(t) + -1._id)**3*1020._id + (Cos(t) + -1._id)**4*1211.25_id)
endif

if(l==10.and.m1==5.and.m2==6) then
   wigd=Cos(t*0.5_id)**11*-1.7888543819998317571273389349850209883524946876892_id*(5._id + (Cos(t) + -1._id)*85._id + (C&
&os(t) + -1._id)**4*302.8125_id + (Cos(t) + -1._id)**2*382.5_id + (Cos(t) + -1._id)**3*605.625_id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==5.and.m2==7) then
   wigd=Cos(t*0.5_id)**12*3.687817782917154924000909712705117262898722019489_id*(10._id + (Cos(t) + -1._id)*90._id + (Co&
&s(t) + -1._id)**3*142.5_id + (Cos(t) + -1._id)**2*213.75_id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==5.and.m2==8) then
   wigd=Cos(t*0.5_id)**13*-9.033271832508971939888199943876461581439881615896_id*(10._id + (Cos(t) + -1._id)*47.5_id + (&
&Cos(t) + -1._id)**2*47.5_id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==5.and.m2==9) then
   wigd=Cos(t*0.5_id)**14*13.9212068442358832728203790160422799419385809430456_id*(-10._id + Cos(t)*20._id)*Sin(t*0.5_id&
&)**4
endif

if(l==10.and.m1==5.and.m2==10) then
   wigd=Cos(t*0.5_id)**15*-124.51505933018704545197426688644896998829430818083_id*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==6.and.m2==-10) then
   wigd=Cos(t*0.5_id)**4*69.606034221179416364101895080211399709692904715228_id*Sin(t*0.5_id)**16
endif

if(l==10.and.m1==6.and.m2==-9) then
   wigd=Cos(t*0.5_id)**3*7.782191208136690340748391680403060624268394261302_id*(12._id + Cos(t)*20._id)*Sin(t*0.5_id)**1&
&5
endif

if(l==10.and.m1==6.and.m2==-8) then
   wigd=Cos(t*0.5_id)**2*5.0497524691810389766816929585348003553019494822398_id*((Cos(t) + -1._id)**2*47.5_id + 120._id &
&+ (Cos(t) + -1._id)*152._id)*Sin(t*0.5_id)**14
endif

if(l==10.and.m1==6.and.m2==-7) then
   wigd=Cos(t*0.5_id)*2.0615528128088302749107049279870385125735996126868_id*((Cos(t) + -1._id)**3*142.5_id + 560._id + &
&(Cos(t) + -1._id)**2*684._id + (Cos(t) + -1._id)*1080._id)*Sin(t*0.5_id)**13
endif

if(l==10.and.m1==6.and.m2==-6) then
   wigd=((Cos(t) + -1._id)**4*302.8125_id + 1820._id + (Cos(t) + -1._id)**3*1938._id + (Cos(t) + -1._id)**2*4590._id + (&
&Cos(t) + -1._id)*4760._id)*Sin(t*0.5_id)**12
endif

if(l==10.and.m1==6.and.m2==-5) then
   wigd=Cos(t*0.5_id)*1.7888543819998317571273389349850209883524946876892_id*((Cos(t) + -1._id)**4*302.8125_id + 1365._i&
&d + (Cos(t) + -1._id)**3*1816.875_id + (Cos(t) + -1._id)*3867.5_id + (Cos(t) + -1._id)**2*4016.25_id)*Sin(t*0.5_id)**11
endif

if(l==10.and.m1==6.and.m2==-4) then
   wigd=Cos(t*0.5_id)**2*2.8284271247461900976033774484193961571393437507539_id*((Cos(t) + -1._id)**4*302.8125_id + 1001&
&._id + (Cos(t) + -1._id)**3*1695.75_id + (Cos(t) + -1._id)*3094._id + (Cos(t) + -1._id)**2*3480.75_id)*Sin(t*0.5_id)**10
endif

if(l==10.and.m1==6.and.m2==-3) then
   wigd=Cos(t*0.5_id)**3*4._id*((Cos(t) + -1._id)**4*302.8125_id + 715._id + (Cos(t) + -1._id)**3*1574.625_id + (Cos(t) &
&+ -1._id)*2431._id + (Cos(t) + -1._id)**2*2983.5_id)*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==6.and.m2==-2) then
   wigd=Cos(t*0.5_id)**4*5.0990195135927848300282241090227819895637709460996_id*((Cos(t) + -1._id)**4*302.8125_id + 495.&
&_id + (Cos(t) + -1._id)**3*1453.5_id + (Cos(t) + -1._id)*1870._id + (Cos(t) + -1._id)**2*2524.5_id)*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==6.and.m2==-1) then
   wigd=Cos(t*0.5_id)**5*5.8878405775518979031760284846395026427830929872876_id*((Cos(t) + -1._id)**4*302.8125_id + 330.&
&_id + (Cos(t) + -1._id)**3*1332.375_id + (Cos(t) + -1._id)*1402.5_id + (Cos(t) + -1._id)**2*2103.75_id)*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==6.and.m2==0) then
   wigd=Cos(t*0.5_id)**6*6.1752192943516858827531745282506443172988423053509_id*(210._id + (Cos(t) + -1._id)**4*302.8125&
&_id + (Cos(t) + -1._id)*1020._id + (Cos(t) + -1._id)**3*1211.25_id + (Cos(t) + -1._id)**2*1721.25_id)*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==6.and.m2==1) then
   wigd=Cos(t*0.5_id)**7*5.8878405775518979031760284846395026427830929872876_id*(126._id + (Cos(t) + -1._id)**4*302.8125&
&_id + (Cos(t) + -1._id)*714._id + (Cos(t) + -1._id)**3*1090.125_id + (Cos(t) + -1._id)**2*1377._id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==6.and.m2==2) then
   wigd=Cos(t*0.5_id)**8*5.0990195135927848300282241090227819895637709460996_id*(70._id + (Cos(t) + -1._id)**4*302.8125_&
&id + (Cos(t) + -1._id)*476._id + (Cos(t) + -1._id)**3*969._id + (Cos(t) + -1._id)**2*1071._id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==6.and.m2==3) then
   wigd=Cos(t*0.5_id)**9*4._id*(35._id + (Cos(t) + -1._id)*297.5_id + (Cos(t) + -1._id)**4*302.8125_id + (Cos(t) + -1._i&
&d)**2*803.25_id + (Cos(t) + -1._id)**3*847.875_id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==6.and.m2==4) then
   wigd=Cos(t*0.5_id)**10*2.8284271247461900976033774484193961571393437507539_id*(15._id + (Cos(t) + -1._id)*170._id + (&
&Cos(t) + -1._id)**4*302.8125_id + (Cos(t) + -1._id)**2*573.75_id + (Cos(t) + -1._id)**3*726.75_id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==6.and.m2==5) then
   wigd=Cos(t*0.5_id)**11*1.7888543819998317571273389349850209883524946876892_id*(5._id + (Cos(t) + -1._id)*85._id + (Co&
&s(t) + -1._id)**4*302.8125_id + (Cos(t) + -1._id)**2*382.5_id + (Cos(t) + -1._id)**3*605.625_id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==6.and.m2==6) then
   wigd=Cos(t*0.5_id)**12*(1._id + (Cos(t) + -1._id)*34._id + (Cos(t) + -1._id)**2*229.5_id + (Cos(t) + -1._id)**4*302.8&
&125_id + (Cos(t) + -1._id)**3*484.5_id)
endif

if(l==10.and.m1==6.and.m2==7) then
   wigd=Cos(t*0.5_id)**13*-2.0615528128088302749107049279870385125735996126868_id*(4._id + (Cos(t) + -1._id)*54._id + (C&
&os(t) + -1._id)**3*142.5_id + (Cos(t) + -1._id)**2*171._id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==6.and.m2==8) then
   wigd=Cos(t*0.5_id)**14*5.0497524691810389766816929585348003553019494822398_id*(6._id + (Cos(t) + -1._id)*38._id + (Co&
&s(t) + -1._id)**2*47.5_id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==6.and.m2==9) then
   wigd=Cos(t*0.5_id)**15*-7.782191208136690340748391680403060624268394261302_id*(-12._id + Cos(t)*20._id)*Sin(t*0.5_id)&
&**3
endif

if(l==10.and.m1==6.and.m2==10) then
   wigd=Cos(t*0.5_id)**16*69.606034221179416364101895080211399709692904715228_id*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==7.and.m2==-10) then
   wigd=Cos(t*0.5_id)**3*33.76388603226826436623377881904421997769695431525_id*Sin(t*0.5_id)**17
endif

if(l==10.and.m1==7.and.m2==-9) then
   wigd=Cos(t*0.5_id)**2*3.7749172176353748486183424034730585291110973523117_id*(14._id + Cos(t)*20._id)*Sin(t*0.5_id)**&
&16
endif

if(l==10.and.m1==7.and.m2==-8) then
   wigd=Cos(t*0.5_id)*2.4494897427831780981972840747058913919659474806567_id*((Cos(t) + -1._id)**2*47.5_id + 136._id + (&
&Cos(t) + -1._id)*161.5_id)*Sin(t*0.5_id)**15
endif

if(l==10.and.m1==7.and.m2==-7) then
   wigd=((Cos(t) + -1._id)**3*142.5_id + 680._id + (Cos(t) + -1._id)**2*726.75_id + (Cos(t) + -1._id)*1224._id)*Sin(t*0.&
&5_id)**14
endif

if(l==10.and.m1==7.and.m2==-6) then
   wigd=Cos(t*0.5_id)*2.0615528128088302749107049279870385125735996126868_id*((Cos(t) + -1._id)**3*142.5_id + 560._id + &
&(Cos(t) + -1._id)**2*684._id + (Cos(t) + -1._id)*1080._id)*Sin(t*0.5_id)**13
endif

if(l==10.and.m1==7.and.m2==-5) then
   wigd=Cos(t*0.5_id)**2*3.687817782917154924000909712705117262898722019489_id*((Cos(t) + -1._id)**3*142.5_id + 455._id &
&+ (Cos(t) + -1._id)**2*641.25_id + (Cos(t) + -1._id)*945._id)*Sin(t*0.5_id)**12
endif

if(l==10.and.m1==7.and.m2==-4) then
   wigd=Cos(t*0.5_id)**3*5.830951894845300470874152877545583076521398334886_id*((Cos(t) + -1._id)**3*142.5_id + 364._id &
&+ (Cos(t) + -1._id)**2*598.5_id + (Cos(t) + -1._id)*819._id)*Sin(t*0.5_id)**11
endif

if(l==10.and.m1==7.and.m2==-3) then
   wigd=Cos(t*0.5_id)**4*8.246211251235321099642819711948154050294398450747_id*((Cos(t) + -1._id)**3*142.5_id + 286._id &
&+ (Cos(t) + -1._id)**2*555.75_id + (Cos(t) + -1._id)*702._id)*Sin(t*0.5_id)**10
endif

if(l==10.and.m1==7.and.m2==-2) then
   wigd=Cos(t*0.5_id)**5*10.5118980208143191442099285287479039691203922519104_id*((Cos(t) + -1._id)**3*142.5_id + 220._i&
&d + (Cos(t) + -1._id)**2*513._id + (Cos(t) + -1._id)*594._id)*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==7.and.m2==-1) then
   wigd=Cos(t*0.5_id)**6*12.138094304022082911201150512925122493534377641012_id*((Cos(t) + -1._id)**3*142.5_id + 165._id&
& + (Cos(t) + -1._id)**2*470.25_id + (Cos(t) + -1._id)*495._id)*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==7.and.m2==0) then
   wigd=Cos(t*0.5_id)**7*12.7305407059820780680148726383177941324863797844472_id*(120._id + (Cos(t) + -1._id)**3*142.5_i&
&d + (Cos(t) + -1._id)*405._id + (Cos(t) + -1._id)**2*427.5_id)*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==7.and.m2==1) then
   wigd=Cos(t*0.5_id)**8*12.138094304022082911201150512925122493534377641012_id*(84._id + (Cos(t) + -1._id)**3*142.5_id &
&+ (Cos(t) + -1._id)*324._id + (Cos(t) + -1._id)**2*384.75_id)*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==7.and.m2==2) then
   wigd=Cos(t*0.5_id)**9*10.5118980208143191442099285287479039691203922519104_id*(56._id + (Cos(t) + -1._id)**3*142.5_id&
& + (Cos(t) + -1._id)*252._id + (Cos(t) + -1._id)**2*342._id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==7.and.m2==3) then
   wigd=Cos(t*0.5_id)**10*8.246211251235321099642819711948154050294398450747_id*(35._id + (Cos(t) + -1._id)**3*142.5_id &
&+ (Cos(t) + -1._id)*189._id + (Cos(t) + -1._id)**2*299.25_id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==7.and.m2==4) then
   wigd=Cos(t*0.5_id)**11*5.830951894845300470874152877545583076521398334886_id*(20._id + (Cos(t) + -1._id)*135._id + (C&
&os(t) + -1._id)**3*142.5_id + (Cos(t) + -1._id)**2*256.5_id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==7.and.m2==5) then
   wigd=Cos(t*0.5_id)**12*3.687817782917154924000909712705117262898722019489_id*(10._id + (Cos(t) + -1._id)*90._id + (Co&
&s(t) + -1._id)**3*142.5_id + (Cos(t) + -1._id)**2*213.75_id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==7.and.m2==6) then
   wigd=Cos(t*0.5_id)**13*2.0615528128088302749107049279870385125735996126868_id*(4._id + (Cos(t) + -1._id)*54._id + (Co&
&s(t) + -1._id)**3*142.5_id + (Cos(t) + -1._id)**2*171._id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==7.and.m2==7) then
   wigd=Cos(t*0.5_id)**14*(1._id + (Cos(t) + -1._id)*27._id + (Cos(t) + -1._id)**2*128.25_id + (Cos(t) + -1._id)**3*142.5_id)
endif

if(l==10.and.m1==7.and.m2==8) then
   wigd=Cos(t*0.5_id)**15*-2.4494897427831780981972840747058913919659474806567_id*(3._id + (Cos(t) + -1._id)*28.5_id + (Cos(t) + -1._id)**2*47.5_id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==7.and.m2==9) then
   wigd=Cos(t*0.5_id)**16*3.7749172176353748486183424034730585291110973523117_id*(-14._id + Cos(t)*20._id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==7.and.m2==10) then
   wigd=Cos(t*0.5_id)**17*-33.76388603226826436623377881904421997769695431525_id*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==8.and.m2==-10) then
   wigd=Cos(t*0.5_id)**2*13.784048752090221767955912552934175427198163558399_id*Sin(t*0.5_id)**18
endif

if(l==10.and.m1==8.and.m2==-9) then
   wigd=Cos(t*0.5_id)*1.5411035007422441125625480953635610563089060058611_id*(16._id + Cos(t)*20._id)*Sin(t*0.5_id)**17
endif

if(l==10.and.m1==8.and.m2==-8) then
   wigd=((Cos(t) + -1._id)**2*47.5_id + 153._id + (Cos(t) + -1._id)*171._id)*Sin(t*0.5_id)**16
endif

if(l==10.and.m1==8.and.m2==-7) then
   wigd=Cos(t*0.5_id)*2.4494897427831780981972840747058913919659474806567_id*((Cos(t) + -1._id)**2*47.5_id + 136._id + (Cos(t) + -1._id)*161.5_id)*Sin(t*0.5_id)**15
endif

if(l==10.and.m1==8.and.m2==-6) then
   wigd=Cos(t*0.5_id)**2*5.0497524691810389766816929585348003553019494822398_id*((Cos(t) + -1._id)**2*47.5_id + 120._id + (Cos(t) + -1._id)*152._id)*Sin(t*0.5_id)**14
endif

if(l==10.and.m1==8.and.m2==-5) then
   wigd=Cos(t*0.5_id)**3*9.033271832508971939888199943876461581439881615896_id*((Cos(t) + -1._id)**2*47.5_id + 105._id + (Cos(t) + -1._id)*142.5_id)*Sin(t*0.5_id)**13
endif

if(l==10.and.m1==8.and.m2==-4) then
   wigd=Cos(t*0.5_id)**4*14.282856857085699995998799622734530557532342319805_id*((Cos(t) + -1._id)**2*47.5_id + 91._id + (Cos(t) + -1._id)*133._id)*Sin(t*0.5_id)**12
endif

if(l==10.and.m1==8.and.m2==-3) then
   wigd=Cos(t*0.5_id)**5*20.199009876724155906726771834139201421207797928959_id*((Cos(t) + -1._id)**2*47.5_id + 78._id + (Cos(t) + -1._id)*123.5_id)*Sin(t*0.5_id)**11
endif

if(l==10.and.m1==8.and.m2==-2) then
   wigd=Cos(t*0.5_id)**6*25.748786379167465530841591578687506418644843651043_id*((Cos(t) + -1._id)**2*47.5_id + 66._id + (Cos(t) + -1._id)*114._id)*Sin(t*0.5_id)**10
endif

if(l==10.and.m1==8.and.m2==-1) then
   wigd=Cos(t*0.5_id)**7*29.732137494637011045224016427862793302897971027442_id*((Cos(t) + -1._id)**2*47.5_id + 55._id + (Cos(t) + -1._id)*104.5_id)*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==8.and.m2==0) then
   wigd=Cos(t*0.5_id)**8*31.183328879386818922579762898076675083852016162476_id*(45._id + (Cos(t) + -1._id)**2*47.5_id + (Cos(t) + -1._id)*95._id)*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==8.and.m2==1) then
   wigd=Cos(t*0.5_id)**9*29.732137494637011045224016427862793302897971027442_id*(36._id + (Cos(t) + -1._id)**2*47.5_id + (Cos(t) + -1._id)*85.5_id)*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==8.and.m2==2) then
   wigd=Cos(t*0.5_id)**10*25.748786379167465530841591578687506418644843651043_id*(28._id + (Cos(t) + -1._id)**2*47.5_id + (Cos(t) + -1._id)*76._id)*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==8.and.m2==3) then
   wigd=Cos(t*0.5_id)**11*20.199009876724155906726771834139201421207797928959_id*(21._id + (Cos(t) + -1._id)**2*47.5_id + (Cos(t) + -1._id)*66.5_id)*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==8.and.m2==4) then
   wigd=Cos(t*0.5_id)**12*14.282856857085699995998799622734530557532342319805_id*(15._id + (Cos(t) + -1._id)**2*47.5_id + (Cos(t) + -1._id)*57._id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==8.and.m2==5) then
   wigd=Cos(t*0.5_id)**13*9.033271832508971939888199943876461581439881615896_id*(10._id + (Cos(t) + -1._id)*47.5_id + (Cos(t) + -1._id)**2*47.5_id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==8.and.m2==6) then
   wigd=Cos(t*0.5_id)**14*5.0497524691810389766816929585348003553019494822398_id*(6._id + (Cos(t) + -1._id)*38._id + (Cos(t) + -1._id)**2*47.5_id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==8.and.m2==7) then
   wigd=Cos(t*0.5_id)**15*2.4494897427831780981972840747058913919659474806567_id*(3._id + (Cos(t) + -1._id)*28.5_id + (Cos(t) + -1._id)**2*47.5_id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==8.and.m2==8) then
   wigd=Cos(t*0.5_id)**16*(1._id + (Cos(t) + -1._id)*19._id + (Cos(t) + -1._id)**2*47.5_id)
endif

if(l==10.and.m1==8.and.m2==9) then
   wigd=Cos(t*0.5_id)**17*-1.5411035007422441125625480953635610563089060058611_id*(-16._id + Cos(t)*20._id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==8.and.m2==10) then
   wigd=Cos(t*0.5_id)**18*13.784048752090221767955912552934175427198163558399_id*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==9.and.m2==-10) then
   wigd=Cos(t*0.5_id)*4.4721359549995793928183473374625524708812367192231_id*Sin(t*0.5_id)**19
endif

if(l==10.and.m1==9.and.m2==-9) then
   wigd=0.5_id*(18._id + Cos(t)*20._id)*Sin(t*0.5_id)**18
endif

if(l==10.and.m1==9.and.m2==-8) then
   wigd=Cos(t*0.5_id)*1.5411035007422441125625480953635610563089060058611_id*(16._id + Cos(t)*20._id)*Sin(t*0.5_id)**17
endif

if(l==10.and.m1==9.and.m2==-7) then
   wigd=Cos(t*0.5_id)**2*3.7749172176353748486183424034730585291110973523117_id*(14._id + Cos(t)*20._id)*Sin(t*0.5_id)**16
endif

if(l==10.and.m1==9.and.m2==-6) then
   wigd=Cos(t*0.5_id)**3*7.782191208136690340748391680403060624268394261302_id*(12._id + Cos(t)*20._id)*Sin(t*0.5_id)**15
endif

if(l==10.and.m1==9.and.m2==-5) then
   wigd=Cos(t*0.5_id)**4*13.9212068442358832728203790160422799419385809430456_id*(10._id + Cos(t)*20._id)*Sin(t*0.5_id)**14
endif

if(l==10.and.m1==9.and.m2==-4) then
   wigd=Cos(t*0.5_id)**5*(8._id + Cos(t)*20._id)*22.011360703055138476529216317127335502674181544772_id*Sin(t*0.5_id)**13
endif

if(l==10.and.m1==9.and.m2==-3) then
   wigd=Cos(t*0.5_id)**6*(6._id + Cos(t)*20._id)*31.128764832546761362993566721612242497073577045207_id*Sin(t*0.5_id)**12
endif

if(l==10.and.m1==9.and.m2==-2) then
   wigd=Cos(t*0.5_id)**7*(4._id + Cos(t)*20._id)*39.681544828799193311277116215239828697071216149532_id*Sin(t*0.5_id)**11
endif

if(l==10.and.m1==9.and.m2==-1) then
   wigd=Cos(t*0.5_id)**8*(2._id + Cos(t)*20._id)*45.820301177534832960627900345328910787301645567449_id*Sin(t*0.5_id)**10
endif

if(l==10.and.m1==9.and.m2==0) then
   wigd=Cos(t)*Cos(t*0.5_id)**9*961.13474601639493532560896830571407682331375986_id*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==9.and.m2==1) then
   wigd=Cos(t*0.5_id)**10*(-2._id + Cos(t)*20._id)*45.820301177534832960627900345328910787301645567449_id*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==9.and.m2==2) then
   wigd=Cos(t*0.5_id)**11*(-4._id + Cos(t)*20._id)*39.681544828799193311277116215239828697071216149532_id*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==9.and.m2==3) then
   wigd=Cos(t*0.5_id)**12*(-6._id + Cos(t)*20._id)*31.128764832546761362993566721612242497073577045207_id*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==9.and.m2==4) then
   wigd=Cos(t*0.5_id)**13*(-8._id + Cos(t)*20._id)*22.011360703055138476529216317127335502674181544772_id*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==9.and.m2==5) then
   wigd=Cos(t*0.5_id)**14*13.9212068442358832728203790160422799419385809430456_id*(-10._id + Cos(t)*20._id)*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==9.and.m2==6) then
   wigd=Cos(t*0.5_id)**15*7.782191208136690340748391680403060624268394261302_id*(-12._id + Cos(t)*20._id)*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==9.and.m2==7) then
   wigd=Cos(t*0.5_id)**16*3.7749172176353748486183424034730585291110973523117_id*(-14._id + Cos(t)*20._id)*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==9.and.m2==8) then
   wigd=Cos(t*0.5_id)**17*1.5411035007422441125625480953635610563089060058611_id*(-16._id + Cos(t)*20._id)*Sin(t*0.5_id)
endif

if(l==10.and.m1==9.and.m2==9) then
   wigd=Cos(t*0.5_id)**18*0.5_id*(-18._id + Cos(t)*20._id)
endif

if(l==10.and.m1==9.and.m2==10) then
   wigd=Cos(t*0.5_id)**19*-4.4721359549995793928183473374625524708812367192231_id*Sin(t*0.5_id)
endif

if(l==10.and.m1==10.and.m2==-10) then
   wigd=Sin(t*0.5_id)**20
endif

if(l==10.and.m1==10.and.m2==-9) then
   wigd=Cos(t*0.5_id)*4.4721359549995793928183473374625524708812367192231_id*Sin(t*0.5_id)**19
endif

if(l==10.and.m1==10.and.m2==-8) then
   wigd=Cos(t*0.5_id)**2*13.784048752090221767955912552934175427198163558399_id*Sin(t*0.5_id)**18
endif

if(l==10.and.m1==10.and.m2==-7) then
   wigd=Cos(t*0.5_id)**3*33.76388603226826436623377881904421997769695431525_id*Sin(t*0.5_id)**17
endif

if(l==10.and.m1==10.and.m2==-6) then
   wigd=Cos(t*0.5_id)**4*69.606034221179416364101895080211399709692904715228_id*Sin(t*0.5_id)**16
endif

if(l==10.and.m1==10.and.m2==-5) then
   wigd=Cos(t*0.5_id)**5*124.51505933018704545197426688644896998829430818083_id*Sin(t*0.5_id)**15
endif

if(l==10.and.m1==10.and.m2==-4) then
   wigd=Cos(t*0.5_id)**6*196.87559523719540998400115941877306606266144591633_id*Sin(t*0.5_id)**14
endif

if(l==10.and.m1==10.and.m2==-3) then
   wigd=Cos(t*0.5_id)**7*278.42413688471766545640758032084559883877161886091_id*Sin(t*0.5_id)**13
endif

if(l==10.and.m1==10.and.m2==-2) then
   wigd=Cos(t*0.5_id)**8*354.92252675760100307924766676167204916378679569317_id*Sin(t*0.5_id)**12
endif

if(l==10.and.m1==10.and.m2==-1) then
   wigd=Cos(t*0.5_id)**9*409.82923272992618480080474681969664551446288789265_id*Sin(t*0.5_id)**11
endif

if(l==10.and.m1==10.and.m2==0) then
   wigd=Cos(t*0.5_id)**10*429.83252552593085495728450959347453264579252131871_id*Sin(t*0.5_id)**10
endif

if(l==10.and.m1==10.and.m2==1) then
   wigd=Cos(t*0.5_id)**11*409.82923272992618480080474681969664551446288789265_id*Sin(t*0.5_id)**9
endif

if(l==10.and.m1==10.and.m2==2) then
   wigd=Cos(t*0.5_id)**12*354.92252675760100307924766676167204916378679569317_id*Sin(t*0.5_id)**8
endif

if(l==10.and.m1==10.and.m2==3) then
   wigd=Cos(t*0.5_id)**13*278.42413688471766545640758032084559883877161886091_id*Sin(t*0.5_id)**7
endif

if(l==10.and.m1==10.and.m2==4) then
   wigd=Cos(t*0.5_id)**14*196.87559523719540998400115941877306606266144591633_id*Sin(t*0.5_id)**6
endif

if(l==10.and.m1==10.and.m2==5) then
   wigd=Cos(t*0.5_id)**15*124.51505933018704545197426688644896998829430818083_id*Sin(t*0.5_id)**5
endif

if(l==10.and.m1==10.and.m2==6) then
   wigd=Cos(t*0.5_id)**16*69.606034221179416364101895080211399709692904715228_id*Sin(t*0.5_id)**4
endif

if(l==10.and.m1==10.and.m2==7) then
   wigd=Cos(t*0.5_id)**17*33.76388603226826436623377881904421997769695431525_id*Sin(t*0.5_id)**3
endif

if(l==10.and.m1==10.and.m2==8) then
   wigd=Cos(t*0.5_id)**18*13.784048752090221767955912552934175427198163558399_id*Sin(t*0.5_id)**2
endif

if(l==10.and.m1==10.and.m2==9) then
   wigd=Cos(t*0.5_id)**19*4.4721359549995793928183473374625524708812367192231_id*Sin(t*0.5_id)
endif

if(l==10.and.m1==10.and.m2==10) then
   wigd=Cos(t*0.5_id)**20
endif

end function wigd

end module acc_wig
