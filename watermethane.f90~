module watermethane_mod
!Module providing a potential for water-methane dimer, written by C L Vaillant. Water-Methane interaction
!is taken from J Chem Phys 123, 134311 (2005).
!Energies are in hartree atomic units, distances in bohr lengths

!Atoms need to be in the order (H H Q D D T T O) (H H H H C M M M M), where
!Q, D, T and M are all partial charges required by the potential.

!From the paper, the rigid body coordinates that should be used for the coordsinirigid
!file are given as:
!H  0.0   1.45365   -1.12169
!H  0.0   -1.45365   -1.12169
!Q  0.0   0.0   -0.04490
!D  0.0   0.70785   0.34527
!D  0.0   -0.70785   0.34527
!T  0.60787   0.0   0.35218
!T  -0.60787   0.0   0.35218
!O  0.0	   0.0	 0.0
!H  0.0   0.0   2.07704
!H  1.95825   0.0   -0.69235
!H  -0.97913   1.69590   -0.69235
!H  -0.97913   -1.69590   -0.69235
!C  0.0   0.0   0.0
!M  0.0   0.0   1.03852
!M  0.97913   0.0   -0.34617
!M  -0.48956   0.84795   -0.34617
!m  -0.48956   -0.84795   -0.34617
implicit none

integer, parameter::             waterdof=7, methanedof=9
!Water charges, in a.u., in the order: H1 H2 Q D1 D2 T1 T2
double precision, parameter::    watercharge(waterdof)= (/0.494714d0, 0.494714d0, -1.830627d0, 0.420599d0, &
     0.420599d0, 0.0d0, 0.0d0/)
!Methane charges, in a.u., in the order: H1, H2, H3, H4, C, M1, M2, M3, M4
double precision, parameter::    methanecharge(methanedof)= (/0.279901d0,0.279901d0,0.279901d0,0.279901d0, 3.590472d0, &
     -1.177519d0,-1.177519d0,-1.177519d0,-1.177519d0/)
!Parameters in Angstroms and kcal/mol
double precision::               betaang(waterdof, methanedof), Aang0(waterdof, methanedof), Aang1(waterdof, methanedof), AangM(waterdof, methanedof)
double precision::               Cang6(waterdof, methanedof), Cang8(waterdof, methanedof), Cang10(waterdof, methanedof), deltaang6(waterdof, methanedof)
double precision::               deltaang8(waterdof, methanedof)
!Parameters in atomic units
double precision::               beta(waterdof, methanedof), A0(waterdof, methanedof), A1(waterdof, methanedof), AM(waterdof, methanedof)
double precision::               C6(waterdof, methanedof), C8(waterdof, methanedof), C10(waterdof, methanedof), delta6(waterdof, methanedof)
double precision::               delta8(waterdof, methanedof)
!Unit conversions
double precision, parameter::    autoang= 0.529177d0, autokcal= 627.503d0
!--------------------------------------------------------------------------------------------------------------------
!Input of data
!H1 and H2:
data betaang(1,:)/2.84808454d0,2.84808454d0,2.84808454d0,2.84808454d0,2.7971225d0,2.75581866d0,2.75581866d0,2.75581866d0,2.75581866d0/
data betaang(2,:)/2.84808454d0,2.84808454d0,2.84808454d0,2.84808454d0,2.7971225d0,2.75581866d0,2.75581866d0,2.75581866d0,2.75581866d0/
!Q:
data betaang(3,:)/2.86928398d0,2.86928398d0,2.86928398d0,2.86928398d0,2.3463075d0,2.31474866d0,2.31474866d0,2.31474866d0,2.31474866d0/
!D1 and D2:
data betaang(4,:)/5.71995231d0,5.71995231d0,5.71995231d0,5.71995231d0,3.16999754d0,2.35594058d0,2.35594058d0,2.35594058d0,2.35594058d0/
data betaang(5,:)/5.71995231d0,5.71995231d0,5.71995231d0,5.71995231d0,3.16999754d0,2.35594058d0,2.35594058d0,2.35594058d0,2.35594058d0/
!T1 and T2:
data betaang(6,:)/6.24776382d0,6.24776382d0,6.24776382d0,6.24776382d0,2.31915671d0,2.28762859d0,2.28762859d0,2.28762859d0,2.28762859d0/
data betaang(7,:)/6.24776382d0,6.24776382d0,6.24776382d0,6.24776382d0,2.31915671d0,2.28762859d0,2.28762859d0,2.28762859d0,2.28762859d0/

!H1 and H2:
data Aang0(1,:)/-752.765963d0,-752.765963d0,-752.765963d0,-752.765963d0,-40504.8858d0,5933.09667d0,5933.09667d0,5933.09667d0,5933.09667d0/
data Aang0(2,:)/-752.765963d0,-752.765963d0,-752.765963d0,-752.765963d0,-40504.8858d0,5933.09667d0,5933.09667d0,5933.09667d0,5933.09667d0/
!Q:
data Aang0(3,:)/4592.62807d0,4592.62807d0,4592.62807d0,4592.62807d0,43408.9282d0,-5121.6292d0,-5121.6292d0,-5121.6292d0,-5121.6292d0/
!D1 and D2:
data Aang0(4,:)/5367.76805d0,5367.76805d0,5367.76805d0,5367.76805d0,55943.4633d0,-2584.27027d0,-2584.27027d0,-2584.27027d0,-2584.27027d0/
data Aang0(5,:)/5367.76805d0,5367.76805d0,5367.76805d0,5367.76805d0,55943.4633d0,-2584.27027d0,-2584.27027d0,-2584.27027d0,-2584.27027d0/
!T1 and T2:
data Aang0(6,:)/1258.12101d0,1258.12101d0,1258.12101d0,1258.12101d0,-19777.5292d0,2979.79274d0,2979.79274d0,2979.79274d0,2979.79274d0/
data Aang0(7,:)/1258.12101d0,1258.12101d0,1258.12101d0,1258.12101d0,-19777.5292d0,2979.79274d0,2979.79274d0,2979.79274d0,2979.79274d0/

!H1 and H2:
data AangM(1,:)/908.685355d0,908.685355d0,908.685355d0,908.685355d0,26577.2947d0,-2622.41721d0,-2622.41721d0,-2622.41721d0,-2622.41721d0/
data AangM(2,:)/908.685355d0,908.685355d0,908.685355d0,908.685355d0,26577.2947d0,-2622.41721d0,-2622.41721d0,-2622.41721d0,-2622.41721d0/
!Q:
data AangM(3,:)/1252.30889d0,1252.30889d0,1252.30889d0,1252.30889d0,-23339.3574d0,8389.34399d0,8389.34399d0,8389.34399d0,8389.34399d0/
!D1 and D2:
data AangM(4,:)/-1430.20075d0,-1430.20075d0,-1430.20075d0,-1430.20075d0,-84638.8183d0,1242.92288d0,1242.92288d0,1242.92288d0,1242.92288d0/
data AangM(5,:)/-1430.20075d0,-1430.20075d0,-1430.20075d0,-1430.20075d0,-84638.8183d0,1242.92288d0,1242.92288d0,1242.92288d0,1242.92288d0/
!T1 and T2:
data AangM(6,:)/-162.796205d0,-162.796205d0,-162.796205d0,-162.796205d0,8863.16664d0,-668.609725d0,-668.609725d0,-668.609725d0,-668.609725d0/
data AangM(7,:)/-162.796205d0,-162.796205d0,-162.796205d0,-162.796205d0,8863.16664d0,-668.609725d0,-668.609725d0,-668.609725d0,-668.609725d0/

!H1 and H2:
data Aang1(1,:)/417.797177d0,417.797177d0,417.797177d0,417.797177d0,19352.3942d0,-3719.99887d0,-3719.99887d0,-3719.99887d0,-3719.99887d0/
data Aang1(2,:)/417.797177d0,417.797177d0,417.797177d0,417.797177d0,19352.3942d0,-3719.99887d0,-3719.99887d0,-3719.99887d0,-3719.99887d0/
!Q:
data Aang1(3,:)/-1789.69987d0,-1789.69987d0,-1789.69987d0,-1789.69987d0,-29232.8793d0,5058.10255d0,5058.10255d0,5058.10255d0,5058.10255d0/
!D1 and D2:
data Aang1(4,:)/-14542.3959d0,-14542.3959d0,-14542.3959d0,-14542.3959d0,-5391.49926d0,651.930524d0,651.930524d0,651.930524d0,651.930524d0/
data Aang1(5,:)/-14542.3959d0,-14542.3959d0,-14542.3959d0,-14542.3959d0,-5391.49926d0,651.930524d0,651.930524d0,651.930524d0,651.930524d0/
!T1 and T2:
data Aang1(6,:)/-9354.24387d0,-9354.24387d0,-9354.24387d0,-9354.24387d0,10463.6792d0,-2062.86173d0,-2062.86173d0,-2062.86173d0,-2062.86173d0/
data Aang1(7,:)/-9354.24387d0,-9354.24387d0,-9354.24387d0,-9354.24387d0,10463.6792d0,-2062.86173d0,-2062.86173d0,-2062.86173d0,-2062.86173d0/

!H1 and H2:
data Cang6(1,:)/-31.1396325d0,-31.1396325d0,-31.1396325d0,-31.1396325d0,-176.385261d0,26.1819133d0,26.1819133d0,26.1819133d0,26.1819133d0/
data Cang6(2,:)/-31.1396325d0,-31.1396325d0,-31.1396325d0,-31.1396325d0,-176.385261d0,26.1819133d0,26.1819133d0,26.1819133d0,26.1819133d0/
!Q:
data Cang6(3,:)/-1291.50705d0,-1291.50705d0,-1291.50705d0,-1291.50705d0,-23944.8325d0,7240.4699d0,7240.4699d0,7240.4699d0,7240.4699d0/
!D1 and D2:
data Cang6(4,:)/172.100547d0,172.100547d0,172.100547d0,172.100547d0,4688.06162d0,-1425.39796,-1425.39796,-1425.39796,-1425.39796/
data Cang6(5,:)/172.100547d0,172.100547d0,172.100547d0,172.100547d0,4688.06162d0,-1425.39796,-1425.39796,-1425.39796,-1425.39796/
!T1 and T2:
data Cang6(6,:)/methanedof*0.0d0/
data Cang6(7,:)/methanedof*0.0d0/

!H1 and H2:
data Cang8(1,:)/40.6973228d0,40.6973228d0,40.6973228d0,40.6973228d0,-470.183908d0,212.925753d0,212.925753d0,212.925753d0,212.925753d0/
data Cang8(2,:)/40.6973228d0,40.6973228d0,40.6973228d0,40.6973228d0,-470.183908d0,212.925753d0,212.925753d0,212.925753d0,212.925753d0/
!Q:
data Cang8(3,:)/7345.62345d0,7345.62345d0,7345.62345d0,7345.62345d0,132928.009d0,-46577.4854d0,-46577.4854d0,-46577.4854d0,-46577.4854d0/
!D1 and D2:
data Cang8(4,:)/-1195.7863d0,-1195.7863d0,-1195.7863d0,-1195.7863d0,-32880.3016d0,11355.4445d0,11355.4445d0,11355.4445d0,11355.4445d0/
data Cang8(5,:)/-1195.7863d0,-1195.7863d0,-1195.7863d0,-1195.7863d0,-32880.3016d0,11355.4445d0,11355.4445d0,11355.4445d0,11355.4445d0/
!T1 and T2:
data Cang8(6,:)/methanedof*0.0d0/
data Cang8(7,:)/methanedof*0.0d0/

!H1 and H2:
data Cang10(1,:)/-13.9555905d0,-13.9555905d0,-13.9555905d0,-13.9555905d0,334.00843d0,-620.561765d0,-620.561765d0,-620.561765d0,-620.561765d0/
data Cang10(2,:)/-13.9555905d0,-13.9555905d0,-13.9555905d0,-13.9555905d0,334.00843d0,-620.561765d0,-620.561765d0,-620.561765d0,-620.561765d0/
!Q:
data Cang10(3,:)/-12518.903d0,-12518.903d0,-12518.903d0,-12518.903d0,-119240.978d0,62124.298d0,62124.298d0,62124.298d0,62124.298d0/
!D1 and D2:
data Cang10(4,:)/1655.75062d0,1655.75062d0,1655.75062d0,1655.75062d0,-8388.62553,-13763.1195d0,-13763.1195d0,-13763.1195d0,-13763.1195d0/
data Cang10(5,:)/1655.75062d0,1655.75062d0,1655.75062d0,1655.75062d0,-8388.62553,-13763.1195d0,-13763.1195d0,-13763.1195d0,-13763.1195d0/
!T1 and T2:
data Cang10(6,:)/methanedof*0.0d0/
data Cang10(7,:)/methanedof*0.0d0/

!H1 and H2:
data deltaang6(1,:)/7.335799d0,7.335799d0,7.335799d0,7.335799d0,2.825277d0,1.943410d0,1.943410d0,1.943410d0,1.943410d0/
data deltaang6(2,:)/7.335799d0,7.335799d0,7.335799d0,7.335799d0,2.825277d0,1.943410d0,1.943410d0,1.943410d0,1.943410d0/
!Q:
data deltaang6(3,:)/4.341591d0,4.341591d0,4.341591d0,4.341591d0,4.288189d0,4.259787d0,4.259787d0,4.259787d0,4.259787d0/
!D1 and D2:
data deltaang6(4,:)/5.759895d0,5.759895d0,5.759895d0,5.759895d0,6.129260d0,3.737571d0,3.737571d0,3.737571d0,3.737571d0/
data deltaang6(5,:)/5.759895d0,5.759895d0,5.759895d0,5.759895d0,6.129260d0,3.737571d0,3.737571d0,3.737571d0,3.737571d0/
!T1 and T2:
data deltaang6(6,:)/methanedof*0.0d0/
data deltaang6(7,:)/methanedof*0.0d0/

!H1 and H2:
!H1, H2, H3, H4, C, M1, M2, M3, M4
data deltaang8(1,:)/1.2192d-2,1.2192d-2,1.2192d-2,1.2192d-2,31.106042,5.747d-3,5.747d-3,5.747d-3,5.747d-3/
data deltaang8(2,:)/1.2192d-2,1.2192d-2,1.2192d-2,1.2192d-2,31.106042,5.747d-3,5.747d-3,5.747d-3,5.747d-3/
!Q:
data deltaang8(3,:)/3.643903d0,3.643903d0,3.643903d0,3.643903d0,4.138380d0,4.368576d0,4.368576d0,4.368576d0,4.368576d0/
!D1 and D2:
data deltaang8(4,:)/4.415080d0,4.415080d0,4.415080d0,4.415080d0,3.962102d0,3.741763d0,3.741763d0,3.741763d0,3.741763d0/
data deltaang8(5,:)/4.415080d0,4.415080d0,4.415080d0,4.415080d0,3.962102d0,3.741763d0,3.741763d0,3.741763d0,3.741763d0/
!T1 and T2:
data deltaang8(6,:)/methanedof*0.0d0/
data deltaang8(7,:)/methanedof*0.0d0/

contains
!-------------------------------------------------------------------------------------------
subroutine clath_pot(x, grad, ereal, gradt, nwater)
implicit none

integer::              i,j,k, nwater
double precision::     x(nwater*24+3*methanedof), grad(nwater*24+3*methanedof), ereal, ewater, ewatermeth
double precision::     xwatermeth(51), tempgrad(51), xwater(nwater*9), gradwater(nwater*9)
logical::             gradt

ereal=0.0d0
grad(:)=0.0d0
do i=1, nwater
   xwatermeth(:)=0.0d0
   do j=1, waterdof+1
      do k=1, 3
         xwatermeth(3*(j-1) + k)= x(3*((waterdof+1)*(i-1) + j -1) +k)/autoang
      end do
   end do
   xwatermeth(25:51)= x(nwater*(waterdof+1)*3 + 1: nwater*(waterdof+1)*3 + methanedof*3)/autoang
   call wmrb(xwatermeth, tempgrad, ewatermeth, gradt)
   ereal= ereal+ ewatermeth*autokcal
   do j=1, waterdof+1
      do k=1, 3
         grad(3*((waterdof+1)*(i-1) + j -1) +k)= grad(3*((waterdof+1)*(i-1) + j -1) +k)&
              +tempgrad(3*(j-1) + k)*autokcal/autoang
      end do
   end do
   grad(nwater*(waterdof+1)*3 + 1: nwater*(waterdof+1)*3 + methanedof*3)= &
        grad(nwater*(waterdof+1)*3 + 1: nwater*(waterdof+1)*3 + methanedof*3)+&
        tempgrad(25:51)*autokcal/autoang
   do k=1,3
      xwater(3*(3*(i-1)+1 -1) + k)= x(3*((waterdof+1)*(i-1) + 8 -1) +k)
      xwater(3*(3*(i-1)+2 -1) + k)= x(3*((waterdof+1)*(i-1) + 1 -1) +k)
      xwater(3*(3*(i-1)+3 -1) + k)= x(3*((waterdof+1)*(i-1) + 2 -1) +k)
   end do
end do

if (gradt) then
  call mbpolenergygradient(nwater,ewater,xwater,gradwater)
else
  call mbpolenergy(nwater,ewater,xwater)
endif
do i=1,nwater
   do k=1,3
      grad(3*((waterdof+1)*(i-1) + 8 -1) +k)= grad(3*((waterdof+1)*(i-1) + 8 -1) +k) + &
           gradwater(3*(3*(i-1)+1 -1) + k)
      grad(3*((waterdof+1)*(i-1) + 1 -1) +k)= grad(3*((waterdof+1)*(i-1) + 1 -1) +k) + &
           gradwater(3*(3*(i-1)+2 -1) + k)
      grad(3*((waterdof+1)*(i-1) + 2 -1) +k)= grad(3*((waterdof+1)*(i-1) + 2 -1) +k) + &
           gradwater(3*(3*(i-1)+3 -1) + k)
   end do
end do


ereal= ereal+ewater
return
end subroutine clath_pot
!-------------------------------------------------------------------------------------------
subroutine wmrb(x, grad, ereal, gradt)
!function that calculates the potential from the water-methane interaction
!for rigid body
!atoms/sites need to be: H H Q D D T T O/ H H H H C M M M M
implicit none
double precision::    x(51), grad(51), ereal
double precision::    r12, rk
logical::             gradt
integer::             i,j,k, totdof

!Convert these from crappy units to atomic units
beta(:,:)= betaang(:,:)*0.529177d0
delta6(:,:)= deltaang6(:,:)*0.529177d0
delta8(:,:)= deltaang8(:,:)*0.529177d0
A0(:,:)= Aang0(:,:)*1.59362d-3
AM(:,:)= AangM(:,:)*1.59362d-3/0.529177d0
A1(:,:)= Aang1(:,:)*1.59362d-3*0.529177d0
C6(:,:)= Cang6(:,:)*1.59362d-3/(0.529177d0**6)
C8(:,:)= Cang8(:,:)*1.59362d-3/(0.529177d0**8)
C10(:,:)= Cang10(:,:)*1.59362d-3/(0.529177d0**10)

totdof=waterdof+methanedof
! do i=1, totdof
!    write(*,*) (x(3*(i-1)+j), j=1,3)
! end do
ereal=0.0d0
do i=1,waterdof
   do j=1, methanedof
      r12= calcr(x(3*(i-1)+1:3*(i-1)+3),x(3*(j-1)+25:3*(j-1)+27))
      ereal= ereal+ tangtoennies(r12, i,j)
   end do
end do

grad(:)=0.0d0
if (gradt) then
   do i=1,waterdof
      do j=1, methanedof
         r12= calcr(x(3*(i-1)+1:3*(i-1)+3),x(3*(j-1)+25:3*(j-1)+27))
         do k=1,3
            rk= x(3*(i-1) + k) - x(3*(j-1) + k+ 24)
            grad(3*(i-1) + k) = grad(3*(i-1) + k) + (rk*gradtangtoennies(r12, i, j)/r12)
            grad(3*(j-1) + k + 24) = grad(3*(j-1) + k+24) - (rk*gradtangtoennies(r12, i, j)/r12)
         end do
      end do
   end do
end if
return
end subroutine wmrb
!-------------------------------------------------------------------------------------------
function calcr(x1, x2)
implicit none
double precision::   x1(:), x2(:), calcr
integer::            i

calcr= 0.0d0
do i=1,3
   calcr= calcr+(x1(i)-x2(i))**2
end do
calcr= sqrt(calcr)
return
end function calcr
!-------------------------------------------------------------------------------------------
function tangtoennies(r, a, b)
!Function returning the diatomic potential of Tang-Toennies, taking indeces
!a and b as a label to the fitted parameters defined at the start.
implicit none
double precision::    r, eint, tangtoennies
integer::             a, b

eint= exp(-beta(a,b)*r)*(A0(a,b) + A1(a,b)*r + AM(a,b)/r)
eint= eint + watercharge(a)*methanecharge(b)/(r)
eint= eint + (C6(a,b)/r**6.0d0)*gammp(7.0d0,delta6(a,b)*r)
eint= eint + (C8(a,b)/r**8.0d0)*gammp(9.0d0,delta8(a,b)*r)
eint= eint + (C10(a,b)/r**10.0d0)*gammp(11.0d0,delta8(a,b)*r)

tangtoennies= eint
return
end function tangtoennies
!-------------------------------------------------------------------------------------------
function gradtangtoennies(r, a, b)
!Function returning the diatomic potential of Tang-Toennies, taking indeces
!a and b as a label to the fitted parameters defined at the start.
implicit none
double precision::    r, grad,gradtangtoennies
integer::             a, b, k

grad= -beta(a,b)*exp(-beta(a,b)*r)*(A0(a,b) + A1(a,b)*r + AM(a,b)/r)
grad= grad + exp(-beta(a,b)*r)*(A1(a,b) - AM(a,b)/r**2)
grad= grad - watercharge(a)*methanecharge(b)/(r**2)
grad= grad - 6.0d0*(C6(a,b)/r**7.0d0)*gammp(7.0d0,delta6(a,b)*r)
grad= grad - 8.0d0*(C8(a,b)/r**9.0d0)*gammp(9.0d0,delta8(a,b)*r)
grad= grad - 10.0d0*(C10(a,b)/r**11.0d0)*gammp(11.0d0,delta8(a,b)*r)
grad= grad + C6(a,b)*(delta6(a,b)**7)*exp(-delta6(a,b)*r)/exp(gammln(7.0d0))
grad= grad + C8(a,b)*(delta8(a,b)**9)*exp(-delta8(a,b)*r)/exp(gammln(9.0d0))
grad= grad + C10(a,b)*(delta8(a,b)**11)*exp(-delta8(a,b)*r)/exp(gammln(11.0d0))

gradtangtoennies= grad
return
end function gradtangtoennies

!-----------------------------------------------------------------------------------------
!Numerical recipes functions to calculate normalized incomplete gamma function P
! (the Tang-Toennies damping function)

function gammp(a,x)
  double precision:: a,gammp,x
  !Returns the incomplete gamma function P (a, x).
  double precision:: gammcf,gamser,gln
  if(x.lt.0.d0.or.a.le.0.d0) then
     write(*,*) 'gammp'
     stop
  end if
  if(x.lt.a+1.)then
     ! Use the series representation.
     call gser(gamser,a,x,gln)
     gammp=gamser
  else
     ! Use the continued fraction representation
     call gcf(gammcf,a,x,gln)
     gammp=1.d0-gammcf
     ! and take its complement.
  endif
  return
end function gammp

subroutine gser(gamser,a,x,gln)
  integer:: ITMAX
  double precision:: a,gamser,gln,x,EPS
  PARAMETER (ITMAX=100,EPS=3.e-7)
  ! Returns the incomplete gamma function P (a, x) evaluated by its series representation as
  ! gamser. Also returns ln Γ(a) as gln.
  integer:: n
  double precision:: ap,del,sum
  gln=gammln(a)
  if(x.le.0.0d0) then
     if(x.lt.0.0d0) stop 'gser'
     gamser=0.0d0
     return
  endif
  ap=a
  sum=1.0d0/a
  del=sum
  do n=1,ITMAX
     ap=ap+1.
     del=del*x/ap
     sum=sum+del
     if(abs(del).lt.abs(sum)*EPS) goto 1
  enddo
  write(*,*) 'gser'
  stop
1  gamser=sum*exp(-x+a*log(x)-gln)
  return
END subroutine gser

subroutine gcf(gammcf,a,x,gln)
  integer:: itmax
  double precision:: a,gammcf,gln,x,eps,fpmin
  parameter (itmax=100,eps=3.e-7,fpmin=1.e-30)
  ! Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction repre-
  ! sentation as gammcf. Also returns ln Γ(a) as gln .
  ! Parameters: ITMAX is the maximum allowed number of iterations; EPS is the relative accu-
  ! racy; FPMIN is a number near the smallest representable floating-point number.
  integer:: i
  double precision:: an,b,c,d,del,h
  gln=gammln(a)
  b=x+1.-a
  ! Set up for evaluating continued fraction by modified Lentz’s method (§5.2) with b0 = 0.
  c=1.0d0/fpmin
  d=1.0d0/b
  h=d
  do i=1,itmax
     ! Iterate to convergence.
     an=-i*(i-a)
     b=b+2.
     d=an*d+b
     if(abs(d).lt.fpmin) d=fpmin
     c=b+an/c
     if(abs(c).lt.fpmin) c=fpmin
     d=1.0d0/d
     del=d*c
     h=h*del
     if(abs(del-1.d0).lt.EPS) goto 1
  enddo
  stop 'gcf'
1 gammcf=exp(-x+a*log(x)-gln)*h
  ! Put factors in front.
  return
end subroutine gcf

function gammln(xx)
  double precision:: gammln,xx
  ! Returns the value ln[Gamma(xx)] for xx > 0.
  integer j
  double precision:: ser,stp,tmp,x,y,cof(6)
  ! Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
  ! accuracy is good enough.
  data cof,stp/76.18009172947146d0,-86.50532032941677d0,&
       24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,&
       -.5395239384953d-5,2.5066282746310005d0/
  x=xx
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
  do j=1,6
     y=y+1.d0
     ser=ser+cof(j)/y
  enddo
  gammln=tmp+log(stp*ser/x)
  return
end function gammln


end module watermethane_mod
