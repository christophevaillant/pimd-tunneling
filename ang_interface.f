c------------------------------------------------------------------------------
c
      subroutine Potential_ang(Etot,omegA,omegB,R_AB,valA,valB)
      implicit real*8 (a-h,o-z)
      dimension Oa(3),Ha1(3),Ha2(3),Ob(3),Hb1(3),Hb2(3)
      dimension omegA(3),omegB(3),valA(3),valB(3)

      common /ddaattaa/ param(18,6),parab(84,6,6),
     .       c(1000),cc(2000),params(1000),
     .       chrg(25),sites(3,25),ind_data(6250),
     .       nsitea,nsiteb,nparsall,iembed,ipotparts,icc,
     .       iemonomer,iembedinterf,ixz

      call interface(omegA,omegB,R_AB,valA,valB,
     .               Oa,Ha1,Ha2,Ob,Hb1,Hb2,
     .               iembedinterf,ixz)

      call CCpol_xyz
     .      (nnr,Oa,Ha1,Ha2,Ob,Hb1,Hb2,iembed,ipotparts,icc,Etot,
     .      param,parab,nsitea,nsiteb,c,
     .      cc,params,nparsall,
     .      ind_data,chrg,sites)

      if (iemonomer.eq.1) then

c        pi=dacos(-1.d0)
c        rad2deg=180.d0/pi
c        a0=0.529177249d0
        h2kcal=627.510d0

c        rA1=valA(1)*a0
c        rA2=valA(2)*a0
c        thA=valA(3)*rad2deg
c        rB1=valB(1)*a0
c        rB2=valB(2)*a0
c        thB=valB(3)*rad2deg

c        write(6,'(a,3f13.8)') "rA1, rA2, thA:  ",rA1, rA2, thA
c        write(6,'(a,3f13.8)') "rB1, rB2, thB:  ",rB1, rB2, thB

        call pots(vA,valA(1),valA(2),valA(3))
        call pots(vB,valB(1),valB(2),valB(3))
c        write(6,'(a,f13.8)') "vA= ",vA*h2kcal
c        write(6,'(a,f13.8)') "vB= ",vB*h2kcal

c        write(6,'(a,f10.5)') "Etot(dimer only): ",Etot
        Etot=Etot+(vA+vB)*h2kcal

      endif

      return
      end
c
c -----------------------------------------------------------------------------
c The Cartesian coordinates of the atoms in the complex are calculated from
c the geometry given using angles and distances.
c To perform this procedure one has to define the type of the embedding which
c is to be used.
c
      subroutine interface(omegA,omegB,R_AB,valA,valB,
     .                     OOOa,HHHa1,HHHa2,OOOb,HHHb1,HHHb2,
     .                     iembedang,ixz)
      implicit real*8 (a-h,o-z)

      dimension O1(3),H1(3),H2(3)

      dimension Rot(3,3),Rota(3,3),Rotb(3,3)
      dimension omegA(3),omegB(3),valA(3),valB(3)
      dimension Oa(3),Ha1(3),Ha2(3)
      dimension Ob(3),Hb1(3),Hb2(3)
      dimension OOOa(3),HHHa1(3),HHHa2(3)
      dimension OOOb(3),HHHb1(3),HHHb2(3)
      parameter (bohr2A=0.52917724924d0)

      A2bohr=1.0d0/bohr2A

      R_AB_ang=R_AB*bohr2A

c calculate the rotation matrix which puts the Oz axis at the position 
c indicated by the Euler angles
      phi  =omegA(1)
      theta=omegA(2)
      chi  =omegA(3)
      call eulerrot(Rota,phi,theta,chi,1)

c put the molecule A at the initial position in the xz (yz) surface (embedding)
      call put_flex(iembedang,valA,Oa,Ha1,Ha2,ixz)
c rotate the molecule
      do i=1,3
        OOOa(i)=0.0d0
        HHHa1(i)=0.0d0
        HHHa2(i)=0.0d0
        do k=1,3
         OOOa(i)=OOOa(i)+Rota(i,k)*Oa(k)
         HHHa1(i)=HHHa1(i)+Rota(i,k)*Ha1(k)
         HHHa2(i)=HHHa2(i)+Rota(i,k)*Ha2(k)
        enddo
      enddo

      phi  =omegB(1)
      theta=omegB(2)
      chi  =omegB(3)
      call eulerrot(Rotb,phi,theta,chi,1)
      call put_flex(iembedang,valB,Ob,Hb1,Hb2,ixz)
      do i=1,3
        OOOb(i)=0.0d0
        HHHb1(i)=0.0d0
        HHHb2(i)=0.0d0
        do k=1,3
         OOOb(i)=OOOb(i)+Rotb(i,k)*Ob(k)
         HHHb1(i)=HHHb1(i)+Rotb(i,k)*Hb1(k)
         HHHb2(i)=HHHb2(i)+Rotb(i,k)*Hb2(k)
        enddo
      enddo
c - for the molecule B shift the coordinates by R_AB
      OOOb(3)=OOOb(3)+R_AB_ang
      HHHb1(3)=HHHb1(3)+R_AB_ang
      HHHb2(3)=HHHb2(3)+R_AB_ang

c check the coordinates of the nuclei
c      write(6,*) "positions of the nuclei:"
c      write(6,'(3d18.10)') OOOa
c      write(6,'(3d18.10)') HHHa1
c      write(6,'(3d18.10)') HHHa2
c      write(6,'(3d18.10)') OOOb
c      write(6,'(3d18.10)') HHHb1
c      write(6,'(3d18.10)') HHHb2

      return
      end
c ---------------------------------------------------------

!     ===+=========+=========+=========+=========+=========+=========+==
      Subroutine Eulerrot(rot,alpha,beta,gamma,transpose)
!     ===+=========+=========+=========+=========+=========+=========+==
*     calculate rotation matrix given the euler angles (in radians)
*       alpha  corresponds to polar coordinate phi
*       beta   corresponds to polar coordinate theta
*       gamma  corresponds to the "spin" angle chi
*
      Implicit real*8 (a-h,o-z)
      double precision rot(3,3)
      integer transpose
c
      cosa=dcos(alpha)
      sina=dsin(alpha)
      cosb=dcos(beta)
      sinb=dsin(beta)
      cosc=dcos(gamma)
      sinc=dsin(gamma)
c
      if( transpose .eq. 0 ) then 
        rot(1,1) = cosb*cosa*cosc-sina*sinc
        rot(1,2) = cosb*sina*cosc+cosa*sinc
        rot(1,3) = -sinb*cosc
        rot(2,1) = -cosb*cosa*sinc-sina*cosc
        rot(2,2) = -cosb*sina*sinc+cosa*cosc
        rot(2,3) = sinb*sinc
        rot(3,1) = sinb*cosa
        rot(3,2) = sinb*sina
        rot(3,3) = cosb
      else
        rot(1,1) = cosb*cosa*cosc-sina*sinc
        rot(2,1) = cosb*sina*cosc+cosa*sinc
        rot(3,1) = -sinb*cosc
        rot(1,2) = -cosb*cosa*sinc-sina*cosc
        rot(2,2) = -cosb*sina*sinc+cosa*cosc
        rot(3,2) = sinb*sinc
        rot(1,3) = sinb*cosa
        rot(2,3) = sinb*sina
        rot(3,3) = cosb
      endif
      return
      End 
c ---------------------------------------------------------
c------------------------------------------------------------------
c - znieksztalcona czasteczke umiesc na plaszczyznie YZ (XZ) w taki sposob,
c   aby srodek masy znalazl sie w poczatku ukladu wspolrzednych,
c   a wektor dwusieczny skierowany byl wzdluz osi Z, 
c - zawolaj procedure eck_rad_tst lub radau_f1_tst ze wspolrzednymi kartezjanskimi
c   uzyskanymi w poprzednim punkcie
c - czasteczke obroc na plaszczyznie wokol poczatku ukladu wspolrzednych 
c   w taki sposob, aby wektor dwusieczny pokryl sie z vecIa 
!          valA(:)  : valence coordinates r1,r2,theta of monomer A

      subroutine put_flex(iembedang,val,O,H1,H2,ixz)
      implicit real*8 (a-h,o-z)
      dimension vi1(3),vi2(3)
      dimension O(3),H1(3),H2(3)
      dimension Or(3),H1r(3),H2r(3)
      dimension w1(3),w2(3),vshift(3)
      dimension Oposition(3),COM(3)

      dimension vecI(3),vecJ(3)
      dimension vz(3)
      dimension val(3),val_bohr(3)
      parameter (bohr2A=0.52917724924d0)
c ds=sin(th_ref_half), dc=cos(th_ref_half)
      a0=0.529177249d0

      pi=dacos(-1.d0)
      rad2deg=180.d0/pi
c distorted molecules put in the plane xz or yz in such a way that COM
c is in the origin and the bisection vector is oriented along the 0z axis.
      r1=val(1)*bohr2A
      r2=val(2)*bohr2A
      th=val(3)
      th_deg=th*rad2deg
      th2=th/2.0d0
      th2_deg=th_deg/2.0d0
      do i=1,3
        Or(i)=0.0d0
        H1r(i)=0.0d0
        H2r(i)=0.0d0
      enddo
      ds=dsin(th2)
      dc=dcos(th2)
      if (ixz.eq.1) then
       H1r(1)= r1*ds   
       H1r(3)=-r1*dc
       H2r(1)=-r2*ds
       H2r(3)=-r2*dc
      else
       H1r(2)=-r1*ds    ! Claude's convention
       H1r(3)=-r1*dc
       H2r(2)= r2*ds
       H2r(3)=-r2*dc
      endif

      call COMcalc(Or,H1r,H2r,COM)
      do i=1,3
        Or(i)=Or(i)-COM(i)
        H1r(i)=H1r(i)-COM(i)
        H2r(i)=H2r(i)-COM(i)
      enddo
c - call embedding procedure eck_rad_tst or radau_f1_tst 
      if (iembedang.eq.1) then
        call eck_rad_tst(Or,H1r,H2r,vecI,vecJ)
      endif
      if (iembedang.eq.2) then
        call radau_f1_tst(Or,H1r,H2r,vecI,vecJ)
      endif
c rotate the molecule around the origin to make the bisection vector
c aligned with the vector vecIa
      vz(1)=0.0d0
      vz(2)=0.0d0
      vz(3)=-1.0d0
c vecI is normalized
      ss=0.0d0
      vzn=0.0d0
      vecIn=0.0d0
      do i=1,3
        ss=ss+vz(i)*vecI(i)
        vzn=vzn+vz(i)*vz(i)
        vecIn=vecIn+vecI(i)*vecI(i)
      enddo
      vznn=dsqrt(vzn)
      vecInn=dsqrt(vecIn)
      angle=dacos(ss/(vznn*vecInn)) 

      call cross(w1,vz,vecI)
      zero=0.0d0
      if (ixz.eq.1) then
        if (w1(2).lt.zero) angle=-angle
      else
        if (w1(1).lt.zero) angle=-angle
      endif

      if (ixz.eq.1) then
c rotate the molecule in the XZ plane by the angle
       O(1)= dcos(angle)*Or(1)-dsin(angle)*Or(3)
       O(2)= 0.0d0    
       O(3)= dsin(angle)*Or(1)+dcos(angle)*Or(3)
       H1(1)= dcos(angle)*H1r(1)-dsin(angle)*H1r(3)
       H1(2)= 0.0d0
       H1(3)= dsin(angle)*H1r(1)+dcos(angle)*H1r(3)
       H2(1)= dcos(angle)*H2r(1)-dsin(angle)*H2r(3)
       H2(2)= 0.0d0
       H2(3)= dsin(angle)*H2r(1)+dcos(angle)*H2r(3)

c       O(1)= dcos(angle)*Or(1)+dsin(angle)*Or(3)
c       O(2)= 0.0d0    
c       O(3)=-dsin(angle)*Or(1)+dcos(angle)*Or(3)
c       H1(1)= dcos(angle)*H1r(1)+dsin(angle)*H1r(3)
c       H1(2)= 0.0d0
c       H1(3)=-dsin(angle)*H1r(1)+dcos(angle)*H1r(3)
c       H2(1)= dcos(angle)*H2r(1)+dsin(angle)*H2r(3)
c       H2(2)= 0.0d0
c       H2(3)=-dsin(angle)*H2r(1)+dcos(angle)*H2r(3)
      else
c rotate the molecule in the YZ plane by the angle
       O(1)= 0.0d0     ! for Claude's convention
       O(2)= dcos(angle)*Or(2)+dsin(angle)*Or(3)
       O(3)=-dsin(angle)*Or(2)+dcos(angle)*Or(3)
       H1(1)= 0.0d0
       H1(2)= dcos(angle)*H1r(2)+dsin(angle)*H1r(3)
       H1(3)=-dsin(angle)*H1r(2)+dcos(angle)*H1r(3)
       H2(1)= 0.0d0
       H2(2)= dcos(angle)*H2r(2)+dsin(angle)*H2r(3)
       H2(3)=-dsin(angle)*H2r(2)+dcos(angle)*H2r(3)
      endif
c sprawdzmy liczac vecI - powinien byc wzdluz osi OZ
c      if (iembedang.eq.1) then
c        call eck_rad_tst(O,H1,H2,vecI,vecJ)
c      endif
c      if (iembedang.eq.2) then
c        call radau_f1_tst(O,H1,H2,vecI,vecJ)
c      endif
c      write(37,*) ' check:',(vecI(i),i=1,3)
c sprawdzam geometrie czasteczki
c      rr1=0.0d0
c      rr2=0.0d0
c      sss=0.0d0
c      do j=1,3
c        rr1=rr1+(H1(j)-O(j))**2
c        rr2=rr2+(H2(j)-O(j))**2
c        sss=sss+(H1(j)-O(j))*(H2(j)-O(j))
c      enddo
c      rr1=dsqrt(rr1)
c      rr2=dsqrt(rr2)
c      ddd=sss/(rr1*rr2)
c      aaa=dacos(ddd)*rad2deg
c      write(6,'(a,3f13.8)') "rr1,rr2,angle:  ",rr1,rr2,aaa

      return
      end
c------------------------------------------------------------------
