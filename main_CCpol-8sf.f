      subroutine init_ccpol(isurf,iemon,iembedang,ixyz)
      implicit real*8 (a-h,o-z)
      dimension Oaa(3),Haa1(3),Haa2(3),Obb(3),Hbb1(3),Hbb2(3)
      dimension carta(3,3),cartb(3,3)
      character*50 nameofsaptdatafile

      common /ddaattaa/ param(18,6),parab(84,6,6),
     .       c(1000),cc(2000),params(1000),
     .       chrg(25),sites(3,25),ind_data(6250),
     .       nsitea,nsiteb,nparsall,iembed,ipotparts,icc,
     .       iemonomer,iembedinterf,ixz

      write(6,'(a)') "-----------------------------------------"
c one can choose from the following list of full-dimensional surfaces
c (the definitions are given in the paper):
c     isurf          short description:            
c       1            CCpol-8sf[2014] (Eckart embedding)
c       2            CCpol-8sfIR[2014] (Eckart embedding)
c       3            CCpol-8sfIR[2012] (Radau f=1 embedding)
c       4            CCpol-8sfIR[2012] (Eckart embedding)
c       5            SAPT-5s'f[2014]         ! 
c       6            SAPT-5s'fIR[2014]       !
c       7            SAPT-5s'fIR[2012]       !
c       8            SAPT-5s'f[2006]         !
c       9            SAPT-5s'fIR[2006]       !
c      10            CCpol-8sfIR[2014] (Radau f=1 embedding)

      if (isurf.eq.1) then
        write(6,*) 
     .      "generating CCpol-8sf[2014] (Eckart embedding) energies"
        iembed=1
        ipotparts=1
        nameofsaptdatafile='./data_SAPT5spf_2014'
        icc=1

      elseif (isurf.eq.2)then
        write(6,*) 
     .      "generating CCpol-8sfIR[2014] (Eckart embedding) energies"
        iembed=1
        ipotparts=1
        nameofsaptdatafile='./data_SAPT5spfIR_2014'
        icc=1

      elseif (isurf.eq.3)then
        write(6,*) 
     .    "generating CCpol-8sfIR[2012] (Radau f=1 embedding) energies"
        iembed=2
        ipotparts=1
        nameofsaptdatafile='./data_SAPT5spfIR_2006'
        icc=1

      elseif (isurf.eq.4)then
        write(6,*) 
     .    "generating CCpol-8sfIR[2012] (Eckart embedding) energies"
        iembed=1
        ipotparts=1
        nameofsaptdatafile='./data_SAPT5spfIR_2006'
        icc=1

      elseif (isurf.eq.5)then
        write(6,*) "generating SAPT-5s'f[2014] energies"
        iembed=1
        ipotparts=1
        nameofsaptdatafile='./data_SAPT5spf_2014'
        icc=0

      elseif (isurf.eq.6)then
        write(6,*) "generating SAPT-5s'fIR[2014] energies"
        iembed=1
        ipotparts=1
        nameofsaptdatafile='./data_SAPT5spfIR_2014'
        icc=0

      elseif (isurf.eq.7)then
        write(6,*) "generating SAPT-5s'fIR[2012] energies"
        iembed=1
        ipotparts=1
        nameofsaptdatafile='./data_SAPT5spfIR_2006'
        icc=0

      elseif (isurf.eq.8)then
        write(6,*) "generating SAPT-5s'f[2006] energies"
        iembed=1
        ipotparts=0
        nameofsaptdatafile='./data_SAPT5spf_2006'
        icc=0

      elseif (isurf.eq.9)then
        write(6,*) "generating SAPT-5s'fIR[2006] energies"
        iembed=1
        ipotparts=0
        nameofsaptdatafile='./data_SAPT5spfIR_2006'
        icc=0

      elseif (isurf.eq.10)then
        write(6,*) 
     .     "generating CCpol-8sfIR[2014] (Radau f=1 embedding) energies"
        iembed=2
        ipotparts=1
        nameofsaptdatafile='./data_SAPT5spfIR_2014'
        icc=1

      else
        write(6,*) "wrong value of isurf: ",isurf
        stop
      endif

c      iembed=2       ! 1 for Eckart     <--------!!!!!!!!!!!
c                     ! 2 for Radau f=1  <--------!!!!!!!!!!!

      pi=dacos(-1.d0)
      rad2deg=180.d0/pi

c..... initialize the interaction energy procedures .....
      open(55,file=nameofsaptdatafile,status='unknown')
      imode = -1
      call poten(r1,r2,theta1,r3,r4,theta2,
     .      rin,alphaAv,betaAv,gammaAv,alphaBv,betaBv,gammaBv,val,imode,
     .      icart,carta,cartb,ipotparts,
     .      param,parab,nsitea,nsiteb,c)
      close(55)

      imode = -1
      call ccpol8s_dimer(imode,Oaa,Haa1,Haa2,Obb,Hbb1,Hbb2,Erigid,
     .                   cc,params,nparsall,
     .                   ind_data,chrg,sites)
c........................................................

      iemonomer=iemon

      if (iembedang.ne.0) then
      iembedinterf=iembedang
      ixz=ixyz
      if (iembedang.eq.1) then
        write(6,"(a,a)") " Eckart embedding is applied to translate",
     .             " (angles + bonds) to the Cartesian coordinates"
      else
       if (iembedang.eq.2) then
        write(6,"(a,a)") " Radau f=1 embedding is applied to translate",
     .                 " (angles + bonds) to the Cartesian coordinates"
       else
          write(6,*) "wrong value of iembedang"
          stop
        endif
      endif

      if (ixyz.eq.1) then
        write(6,*) "initially the molecule is put in the xz plane"
      else
        if (ixyz.eq.0) then
          write(6,*) "initially the molecule is put in the yz plane"
        else
          write(6,*) "wrong value of ixz"
          stop
        endif
      endif
      endif ! <---- if (iembedang.eq.0) then

      if (iemon.eq.1) then
        write(6,*) 
     .    "energy of the monomers is added to the interaction energy"
      else
        if (iemon.ne.0) then
          write(6,*) "wrong value of iemonomer"
          stop
        endif
      endif
      write(6,'(a)') "-----------------------------------------"

      call test_parameters(isurf,iemon)

      return
      end
c
c------------------------------------------------------------------------------
      subroutine test_parameters(isurf,iemon)
      implicit real*8 (a-h,o-z)
      dimension val(10),valm(10)
      dimension Oa(3),Ha1(3),Ha2(3),Ob(3),Hb1(3),Hb2(3)
      data val / -1.45754, -1.36660, -1.34676, -1.34688, -0.76978,
     .           -0.67884, -0.65912, -0.75268, -0.65451, -1.36648/
      data valm/ 25.04582, 25.13676, 25.15660, 25.15648, 25.73358,
     .           25.82452, 25.84424, 25.75068, 25.84885, 25.13688/
      data Oa  /  0.6458557220D-01, 0.3399054992D-02,-0.1782922818D-01/
      data Ha1 / -0.5505396894D+00, 0.6383283738D-01,-0.6241648475D+00/
      data Ha2 / -0.4744802670D+00,-0.1177783097D+00, 0.9071276528D+00/
      data Ob  / -0.5658499752D-01,-0.1827211353D-03, 0.2523332441D+01/
      data Hb1 /  0.6334615998D+00, 0.2269642816D+00, 0.2055825014D+01/
      data Hb2 /  0.2645834235D+00,-0.2240643644D+00, 0.3293430563D+01/

      call ccpol(Etot,Oa,Ha1,Ha2,Ob,Hb1,Hb2)

      write(6,'(a)') "to check if the parameters are read properly:"
      if (iemon.eq.0) then
        write(6,'(2(a,f10.5))') 
     .      "Calculated value: ",Etot,"  and should be:",val(isurf)
      else
        write(6,'(2(a,f10.5))') 
     .      "Calculated value: ",Etot,"  and should be:",valm(isurf)
      endif
      write(6,'(a)') "-----------------------------------------"

      return
      end
c
c------------------------------------------------------------------------------
c
c The input Cartesian coordinates should be given in Angstroms.
c
      subroutine ccpol(Etot,Oa,Ha1,Ha2,Ob,HB1,Hb2)
      implicit real*8 (a-h,o-z)
      dimension Oa(3),Ha1(3),Ha2(3),Ob(3),Hb1(3),Hb2(3)

      common /ddaattaa/ param(18,6),parab(84,6,6),
     .       c(1000),cc(2000),params(1000),
     .       chrg(25),sites(3,25),ind_data(6250),
     .       nsitea,nsiteb,nparsall,iembed,ipotparts,icc,
     .       iemonomer,iembedinterf,ixz

      call CCpol_xyz
     .      (nnr,Oa,Ha1,Ha2,Ob,Hb1,Hb2,iembed,ipotparts,icc,Etot,
     .      param,parab,nsitea,nsiteb,c,
     .      cc,params,nparsall,
     .      ind_data,chrg,sites)

      if (iemonomer.eq.1) then
        rA1=0.0d0
        rA2=0.0d0
        rB1=0.0d0
        rB2=0.0d0
        ssA=0.0d0
        ssB=0.0d0
        do j=1,3
          rA1=rA1+(Ha1(j)-Oa(j))**2
          rA2=rA2+(Ha2(j)-Oa(j))**2
          rB1=rB1+(Hb1(j)-Ob(j))**2
          rB2=rB2+(Hb2(j)-Ob(j))**2
          ssA=ssA+(Ha1(j)-Oa(j))*(Ha2(j)-Oa(j))
          ssB=ssB+(Hb1(j)-Ob(j))*(Hb2(j)-Ob(j))
        enddo
        rA1=dsqrt(rA1)
        rA2=dsqrt(rA2)
        rB1=dsqrt(rB1)
        rB2=dsqrt(rB2)
        thA=dacos(ssA/(rA1*rA2))
        thB=dacos(ssB/(rB1*rB2))

c        pi=dacos(-1.d0)
c        rad2deg=180.d0/pi
        a0=0.529177249d0
        h2kcal=627.510d0

c        write(6,'(a,3f13.8)') "rA1, rA2, thA:  ",rA1, rA2, thA*rad2deg
c        write(6,'(a,3f13.8)') "rB1, rB2, thB:  ",rB1, rB2, thB*rad2deg
        rA1=rA1/a0
        rA2=rA2/a0
        rB1=rB1/a0
        rB2=rB2/a0
        call pots(vA,rA1,rA2,thA)
        call pots(vB,rB1,rB2,thB)
c        write(6,'(a,f13.8)') "vA= ",vA*h2kcal
c        write(6,'(a,f13.8)') "vB= ",vB*h2kcal

c        write(6,'(a,f13.8)') "Etot(dimer only): ",Etot
        Etot=Etot+(vA+vB)*h2kcal

      endif

      return
      end
c
c -----------------------------------------------------------------------------
      Subroutine CCpol_xyz(nnr,Oa,Ha1,Ha2,Ob,HB1,Hb2,
     .                 iembed,ipotparts,icc,Etot,
     .                 param,parab,nsitea,nsiteb,c,
     .                 cc,params,nparsall,
     .                 ind_data,chrg,sites)
      implicit real*8 (a-h,o-z)
      dimension Oa(3),Ha1(3),Ha2(3),COMa(3)
      dimension Ob(3),Hb1(3),Hb2(3),COMb(3)
      dimension Oaa(3),Haa1(3),Haa2(3)
      dimension Obb(3),Hbb1(3),Hbb2(3)
      character*1 c1
      dimension vecI(3),vecIa(3),vecIb(3),vecJa(3),vecJb(3)
      dimension carta(3,3),cartb(3,3)
      dimension cartaa(3,3),cartbb(3,3)
      character*50 nameofsaptdatafile
      parameter (nsitemax=8,ntypemax=6,ntypemax2=ntypemax*ntypemax)
      parameter (maxpar1=18,maxpar2=84)
c
      dimension param(maxpar1,ntypemax),parab(maxpar2,ntypemax,ntypemax)
      dimension c(1000)
      PARAMETER (maxlin=2000,mmaxpar=1000)
      dimension cc(maxlin),params(mmaxpar)

      dimension chrg(25)
      dimension sites(3,25)
      dimension ind_data(6250)
      integer omp_get_num_threads,omp_get_thread_num

c      ns=omp_get_num_threads()
c      id=omp_get_thread_num()

      a0=0.529177249d0
      vecI(1)=0.0d0
      vecI(2)=0.0d0
      vecI(3)=-1.0d0

      imode = 0

        call align_on_z_axis(Oa,Ha1,Ha2,Ob,Hb1,Hb2,Rcom)

        do jj=1,3      
          carta(1,jj)=Oa(jj)           
          carta(2,jj)=Ha1(jj)           
          carta(3,jj)=Ha2(jj)           
          cartb(1,jj)=Ob(jj)           
          cartb(2,jj)=Hb1(jj)           
          cartb(3,jj)=Hb2(jj)           
        enddo

c make embedding for molecule A
        if (iembed.eq.1) then
          call eck_rad_tst(Oa,Ha1,Ha2,vecIa,vecJa)
        endif
        if (iembed.eq.2) then
          call radau_f1_tst(Oa,Ha1,Ha2,vecIa,vecJa)
        endif

        call put_rigid(vecIa,vecJa,Oaa,Haa1,Haa2)

c make embedding for molecule B
        Ob(3)=Ob(3)-Rcom
        Hb1(3)=Hb1(3)-Rcom
        Hb2(3)=Hb2(3)-Rcom

        if (iembed.eq.1) then
          call eck_rad_tst(Ob,Hb1,Hb2,vecIb,vecJb)
        endif
        if (iembed.eq.2) then
          call radau_f1_tst(Ob,Hb1,Hb2,vecIb,vecJb)
        endif

        Ob(3)=Ob(3)+Rcom
        Hb1(3)=Hb1(3)+Rcom
        Hb2(3)=Hb2(3)+Rcom

        call put_rigid(vecIb,vecJb,Obb,Hbb1,Hbb2)
        Obb(3)=Obb(3)+Rcom
        Hbb1(3)=Hbb1(3)+Rcom
        Hbb2(3)=Hbb2(3)+Rcom

        do jj=1,3          
          cartaa(1,jj)=Oaa(jj)           
          cartaa(2,jj)=Haa1(jj)           
          cartaa(3,jj)=Haa2(jj)           
          cartbb(1,jj)=Obb(jj)           
          cartbb(2,jj)=Hbb1(jj)           
          cartbb(3,jj)=Hbb2(jj)           
        enddo

c        icc=1
c        icc=0
        if (icc.eq.1) then
         call driver_potss_sapt5sf(carta,cartb,val,
     .                     ipotparts,param,parab,nsitea,nsiteb,c)
         call driver_potss_sapt5sf(cartaa,cartbb,vall,
     .                     ipotparts,param,parab,nsitea,nsiteb,c)
         call ccpol8s_dimer(imode,Oaa,Haa1,Haa2,Obb,Hbb1,Hbb2,Erigid,
     .                     cc,params,nparsall,
     .                     ind_data,chrg,sites)
         Etot= Erigid + (val-vall)
        else 
         call driver_potss_sapt5sf(carta,cartb,val,
     .                     ipotparts,param,parab,nsitea,nsiteb,c)
         Etot= val
        endif

      return
      end
c------------------------------------------------------------------
c The rigid (reference) molecule is embedded into the plane defined by the original
c (distorted) molecule and by the axis I related to the embedding method. The bisection
c vector of the rigid molecule is put along the I vector and COM of the rigid molecule
c is put at the position of COM of the distorted one.
c
c Information about two normalized vectors is passed to the procedure.
c vi1 - the bisection vector of the reference water is put in this direction
c vi2 - it is perpendicular to vi1 and lies in the plane of the molecule
c
      subroutine put_rigid (vi1,vi2,O,H1,H2)
      implicit real*8 (a-h,o-z)
      dimension vi1(3),vi2(3)
      dimension O(3),H1(3),H2(3)
      dimension w1(3),w2(3),vshift(3)
      dimension Oposition(3),COM(3)
c ds=sin(th_ref_half), dc=cos(th_ref_half)
      ds=0.79170358110560535d0
      dc=0.61090542612139243d0
      rOHref=0.97162570027717354d0
      com_shift=0.66429466101803d-01
      a0=0.529177249d0

      pi=dacos(-1.d0)
      rad2deg=180.d0/pi

c Initially O is put in the origin, but next the whole rigid molecule should be shifted
c to put its COM in the origin. This shif is in the vi1 direction by com_shift.
      Oposition(1)=0.0d0
      Oposition(2)=0.0d0
      Oposition(3)=0.0d0

      do j=1,3
        w1(j)=dc*vi1(j)+ds*vi2(j)
        w2(j)=dc*vi1(j)-ds*vi2(j)
        vshift(j)=-com_shift*vi1(j)
      enddo
      do j=1,3
        w1(j)=rOHref*w1(j)
        w2(j)=rOHref*w2(j)
      enddo

      do j=1,3
        w1(j)=w1(j)+vshift(j)
        w2(j)=w2(j)+vshift(j)
        Oposition(j)=Oposition(j)+vshift(j)
      enddo
      do j=1,3
        O(j)=Oposition(j)
        H1(j)=w1(j)
        H2(j)=w2(j)
      enddo

      return
      end
c------------------------------------------------------------------
c This procedure is based on the Konrad's antigetgeo.py program.
c It transforms the cartesian coordinates of two water molecules given in the following
c order (OA, HA1, HA2, OB, HB1, HB2), into cartesian coordinates but in the system, where
c 1. COM of the first molecule is at the (0,0,0) position
c 2. COM of the second one is at (0,0,R), where R is the distance between COMs
c
      subroutine align_on_z_axis(O1A,H1A,H2A,O1B,H1B,H2B,Rcom)
      implicit real*8 (a-h,o-z)
      dimension xyzA(3,3), xyzB(3,3)
      dimension xyzAA(3,3), xyzBB(3,3)
      dimension xyzAAA(3,3), xyzBBB(3,3)
      dimension ama(3),comA(3),comB(3),s(3),s1(3)
      dimension O1A(3),H1A(3),H2A(3),O1B(3),H1B(3),H2B(3)
      character*1 c1
c      data ama/15.994915d0, 1.007825d0, 1.007825d0/
c      save ama

      bohr2a=0.529177d0
      ama(1)=15.994915d0
      ama(2)=1.007825d0 
      ama(3)=1.007825d0 

c read the cartesian coordiantes
      thr=1.0d-9

      do j=1,3
        xyzA(1,j)=O1A(j)
        xyzA(2,j)=H1A(j)
        xyzA(3,j)=H2A(j)
        xyzB(1,j)=O1B(j)
        xyzB(2,j)=H1B(j)
        xyzB(3,j)=H2B(j)
      enddo

c test a procedure
      call COMcalc3 (xyzA,comA)
      call COMcalc3 (xyzB,comB)

      sss=0.0d0
      do j=1,3
        sss=sss+(comB(j)-comA(j))**2
      enddo
      Rcom=dsqrt(sss)

c shift the whole system in such a way that comA coincides with the origin of the 
c system of coordinates
      do i=1,3
        do j=1,3
          xyzA(i,j)=xyzA(i,j)-comA(j)
          xyzB(i,j)=xyzB(i,j)-comA(j)
        enddo
        comB(i)=comB(i)-comA(i)
      enddo
      do i=1,3
        comA(i)=comA(i)-comA(i)
      enddo

      ss=dsqrt(comB(1)*comB(1)+comB(2)*comB(2))
      if (ss.lt.thr) then
c        write(66,*) "    both center of masses are on the Z axis, ss=",ss
c "Konrad problem"
c Even if both are on the Z axis, take care about position of comB with respect to comA: 
c the value of the z-coordinate of B should be larger than of A (means 0 in our case)
c Make an inversion when needed
        if (comB(3).lt.0.0d0) then
c          write(66,*) "I have to flip the complex"
          do i=1,3
            do j=1,3
              xyzA(i,j)=-xyzA(i,j)
              xyzB(i,j)=-xyzB(i,j)
            enddo
            comA(i)=-comA(i)
            comB(i)=-comB(i)
          enddo
        endif
        do i=1,3
          do j=1,3
            xyzAAA(i,j)=xyzA(i,j)
            xyzBBB(i,j)=xyzB(i,j)
          enddo
        enddo
      else
c        write(66,*) "not both center of masses are on Z axis, ss=",ss
c find a vector s (rotation axis) which is perpendicular both to the vector connecting
c comA and comB, and to z-axis.
c s = comB x (0,0,1) = (comB(2),-comB(1),0)
c + normalization
        xnorm=dsqrt(comB(1)*comB(1)+comB(2)*comB(2))
        s(1)=comB(2)/xnorm
        s(2)=-comB(1)/xnorm
        s(3)=0.0d0
c find an angle between comB and (0,0,1)
        rr=dsqrt(comB(1)*comB(1)+comB(2)*comB(2)+comB(3)*comB(3))
        ccos=comB(3)/rr
        ssin=dsqrt(1.0d0-ccos*ccos)
c now we have a new orthogonal basis {s,s',(0,0,01)}, where
c s'=(-comB(1)/sqrt(comB(1)**2+comB(2)**2),-comB(2)/sqrt(comB(1)**2+comB(2)**2),0)
        s1(1)=-comB(1)/xnorm
        s1(2)=-comB(2)/xnorm
        s1(3)=0.0d0

c let's express the positions of atoms in the new basis:
        do i=1,3
          xyzAA(i,1)=xyzA(i,1)*s(1)+xyzA(i,2)*s(2)
          xyzAA(i,2)=xyzA(i,1)*s1(1)+xyzA(i,2)*s1(2)
          xyzAA(i,3)=xyzA(i,3)
        enddo
        do i=1,3
          xyzBB(i,1)=xyzB(i,1)*s(1)+xyzB(i,2)*s(2)
          xyzBB(i,2)=xyzB(i,1)*s1(1)+xyzB(i,2)*s1(2)
          xyzBB(i,3)=xyzB(i,3)
        enddo

c rotate everything:
        do i=1,3
          xyzAAA(i,1)= xyzAA(i,1)
          xyzAAA(i,2)= xyzAA(i,2)*ccos+xyzAA(i,3)*ssin
          xyzAAA(i,3)=-xyzAA(i,2)*ssin+xyzAA(i,3)*ccos
        enddo
        do i=1,3
          xyzBBB(i,1)= xyzBB(i,1)
          xyzBBB(i,2)= xyzBB(i,2)*ccos+xyzBB(i,3)*ssin
          xyzBBB(i,3)=-xyzBB(i,2)*ssin+xyzBB(i,3)*ccos
        enddo
      endif

      do j=1,3
        O1A(j)=xyzAAA(1,j)
        H1A(j)=xyzAAA(2,j)
        H2A(j)=xyzAAA(3,j)
        O1B(j)=xyzBBB(1,j)
        H1B(j)=xyzBBB(2,j)
        H2B(j)=xyzBBB(3,j)
      enddo

      return
      end
c-------------------------------------------------------------------------------------
      subroutine COMcalc3 (xyz,COM)
      implicit real*8 (a-h,o-z)
      dimension xyz(3,3), xyzB(3,3)
      dimension O1(3),H1(3),H2(3),COM(3)
      real*8 mO,mH,M
      mO=15.994 914 6221d0
      mH= 1.007 825 032 1d0

      do j=1,3
        O1(j)=xyz(1,j)
        H1(j)=xyz(2,j)
        H2(j)=xyz(3,j)
      enddo
      M=mO+mH+mH
      do i=1,3
        COM(i)= (mO*O1(i) +mH*H1(i) +mH*H2(i))/M
      enddo 
      end
c-------------------------------------------------------------------------------------
c
c I assume that the molecule is already shifted to fit COM to (0,0,0) 
c 
      subroutine eck_rad_tst(r0,r1,r2,vecI,vecJ)
      implicit real*8 (a-h,o-z)
      dimension rot(3,3)
      dimension q1(3),q2(3),bv(3),bbv(3)
      dimension rcom(3),r0(3),r1(3),r2(3)
      dimension br0(3)
      dimension pom1(3),pom2(3)
      dimension temp1(3),temp2(3),vecI(3),vecJ(3)
      character*1 c1
c      data xmO /15.994915d0/, xmH /1.007825d0/
      xmO=15.9949146221d0
      xmH=1.0078250321d0
      xq1e=0.95111822d0
      xq2e=0.95111822d0
      theta_r_e_deg=107.952612d0
      theta_r_e=1.88412851d0
      a0=0.529177249d0

      deg2rad=dacos(-1.d0)/180.d0

c To perform the Eckart embedding we need info about the Radau coordinates of 
c the molecule in its equilibrium geometry. It was precalculated with the 
c data given above and hard-coded as xq1e, xq2e, and theta_r_e

      xm12=2*xmH
      xm=xm12+xmO
      xm0=xmO
      xm1=xmH
      xm2=xmH
      alpha=dsqrt(xm0/xm)               ! Eq. (26)    can be precalculated
      b=(alpha-alpha*alpha)*xm/xm12     ! Eq. (27)

c move the COM of the molecule to the origin of the coordinate system
      do j=1,3
        q1(j)=r1(j)-b*r0(j)    ! Eq. (28)
        q2(j)=r2(j)-b*r0(j)    ! Eq. (29)
        br0(j)=b*r0(j)
      enddo

c calculate length of q1 and q2 vectors, and the angle between them
        xq1=0.0d0
        xq2=0.0d0
        sss=0.0d0
        do j=1,3
          xq1=xq1+q1(j)*q1(j)
          xq2=xq2+q2(j)*q2(j)
          sss=sss+q1(j)*q2(j)
        enddo
        xq1=dsqrt(xq1)
        xq2=dsqrt(xq2)
        theta_r=dacos( sss/(xq1*xq2) )

        sss=0.0d0
        do j=1,3
          pom1(j)=q1(j)/xq1
          pom2(j)=q2(j)/xq2
          bv(j)=pom1(j)+pom2(j)
          sss=sss+bv(j)*bv(j)
        enddo
        sss=dsqrt(sss)
        do j=1,3
          bv(j)=bv(j)/sss
        enddo

c formula (33) from Wei & Carrington, JCP 107, 2813 (1997)
c masses are the same and can be skipped, in fact xq1e and xq2e could be also
        eta_e=0.5d0*theta_r_e
        ang=theta_r-theta_r_e+eta_e
        sss=(xq2e*xq2*dsin(ang) + xq1e*xq1*dsin(eta_e))/
     .      (xq2e*xq2*dcos(ang) + xq1e*xq1*dcos(eta_e))
        eta=datan(sss)

c To obtain the I vector which bisects the angle made by q1e and q2e I do not
c need an explicit information about these vectors, but I can get it by the
c rotation  of the q1 vector in direction of q2 by the angle equal to eta.
c Temporary axes temp1=norm(q1) and temp2=norm[q2-((q1*q2)/(q1*q1)) q1]
c                                or temp2=norm[q2-(temp1*q2) temp1]

        sss=0.0d0
        do j=1,3
          temp1(j)=q1(j)/xq1
          sss=sss+temp1(j)*q2(j)
        enddo
        ttt=0.0d0
        do j=1,3
          temp2(j)=q2(j)-sss*temp1(j)  
          ttt=ttt+temp2(j)**2
        enddo
        ttt=dsqrt(ttt)
        do j=1,3
          temp2(j)=temp2(j)/ttt
        enddo
        
c        write(67,666) (temp1(j),j=1,3)," temp1"
c        write(67,666) (temp2(j),j=1,3)," temp2"
        s1=0.0d0
        s2=0.0d0
        s3=0.0d0
        do j=1,3
          s1=s1+temp1(j)*temp1(j)
          s2=s2+temp2(j)*temp2(j)
          s3=s3+temp1(j)*temp2(j)
        enddo
c        write(67,'(a,3f12.6)') "s1, s2, s3:", s1, s2, s3
c        write(67,*)
c vecI=cos(eta)*temp1 + sin(eta)*temp2
        do j=1,3
          vecI(j)= dcos(eta)*temp1(j)+dsin(eta)*temp2(j)
          vecJ(j)=-dsin(eta)*temp1(j)+dcos(eta)*temp2(j)
        enddo

      return

  666 format(3f12.6,3x,a)
  777 format(a,2f12.8,f12.6,f12.8)
  778 format(a,f12.6,f12.8)

      end
c-----------------------------------------------------------------------------------
c
c I assume that the molecule is already shifted to fit COM to (0,0,0) 
c 
      subroutine radau_f1_tst(r0,r1,r2,vecI,vecJ)
      implicit real*8 (a-h,o-z)
      dimension rot(3,3)
      dimension q1(3),q2(3),bv(3),bbv(3)
      dimension rcom(3),r0(3),r1(3),r2(3)
      dimension br0(3)
      dimension pom1(3),pom2(3)
      dimension temp1(3),temp2(3),vecI(3),vecJ(3)
      character*1 c1
c      data xmO /15.994915d0/, xmH /1.007825d0/
      xmO=15.9949146221d0
      xmH=1.0078250321d0
      xq1e=0.95111822d0
      xq2e=0.95111822d0
      theta_r_e_deg=107.952612d0
      theta_r_e=1.88412851d0
      a0=0.529177249d0

      deg2rad=dacos(-1.d0)/180.d0

c data given above and hard-coded as xq1e, xq2e, and theta_r_e


      xm12=2*xmH
      xm=xm12+xmO
      xm0=xmO
      xm1=xmH
      xm2=xmH
      alpha=dsqrt(xm0/xm)               ! Eq. (26)    can be precalculated
      b=(alpha-alpha*alpha)*xm/xm12     ! Eq. (27)

c move the COM of the molecule to the origin of the coordinate system
      do j=1,3
        q1(j)=r1(j)-b*r0(j)    ! Eq. (28)
        q2(j)=r2(j)-b*r0(j)    ! Eq. (29)
        br0(j)=b*r0(j)
      enddo

c calculate length of q1 and q2 vectors, and the angle between them
        xq1=0.0d0
        xq2=0.0d0
        sss=0.0d0
        do j=1,3
          xq1=xq1+q1(j)*q1(j)
          xq2=xq2+q2(j)*q2(j)
          sss=sss+q1(j)*q2(j)
        enddo
        xq1=dsqrt(xq1)
        xq2=dsqrt(xq2)
        theta_r=dacos( sss/(xq1*xq2) )

        sss=0.0d0
        do j=1,3
          pom1(j)=q1(j)/xq1
          pom2(j)=q2(j)/xq2
          bv(j)=pom1(j)+pom2(j)
          sss=sss+bv(j)*bv(j)
        enddo
        sss=dsqrt(sss)
        do j=1,3
          bv(j)=bv(j)/sss
          vecI(j)=bv(j)
        enddo

c I need also a vector which is perpendicular to vecI and lies in the plane
c of the molecule.
c I can obtain it from the following expression:
c temp2=norm[q2-(vecI*q2) vecI]     (q2 or q1 - does not matter)

        sss=0.0d0
        do j=1,3
          sss=sss+vecI(j)*q2(j)
        enddo
        ttt=0.0d0
        do j=1,3
          temp2(j)=q2(j)-sss*vecI(j)
          ttt=ttt+temp2(j)**2
        enddo
        ttt=dsqrt(ttt)
        do j=1,3
          temp2(j)=temp2(j)/ttt
c          vecJ(j)=temp2(j)          mozna wybrac zwrot
          vecJ(j)=-temp2(j)
        enddo

      return

  666 format(3f12.6,3x,a)
  777 format(a,2f12.8,f12.6,f12.8)
  778 format(a,f12.6,f12.8)

      end
c-----------------------------------------------------------------------------------
      subroutine vecprod(v1,v2,v3)
      implicit real*8 (a-h,o-z)
      dimension v1(3), v2(3), v3(3)

      v3(1)=v1(2)*v2(3)-v1(3)*v2(2)
      v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
      v3(3)=v1(1)*v2(2)-v1(2)*v2(1)
      return
      end
c-----------------------------------------------------------------------
      subroutine read_cc_data(ind_data,chrg,sites)
      implicit real*8 (a-h,o-z)

      character*1 cnull
      dimension ind_charge(25)
      dimension ind_beta(25,25)
      dimension ind_d1(25,25)
      dimension ind_d6(25,25)
      dimension ind_d8(25,25)
      dimension ind_d10(25,25)
      dimension ind_C6(25,25)
      dimension ind_C8(25,25)
      dimension ind_C10(25,25)
      dimension chrg(25)
      dimension sites(3,25)

      dimension ind_data(6250)
      dimension ind_data1(6250)
      equivalence (ind_data1(1),ind_charge)
      equivalence (ind_data1(26),ind_beta)
      equivalence (ind_data1(651),ind_d1)
      equivalence (ind_data1(1276),ind_d6)
      equivalence (ind_data1(1901),ind_d8)
      equivalence (ind_data1(2526),ind_d10)
      equivalence (ind_data1(3151),ind_c6)
      equivalence (ind_data1(3776),ind_c8)
      equivalence (ind_data1(4401),ind_c10)

      do i=1,6250
        ind_data1(i)=0
      enddo

      do i=1,25
       ind_charge(i)=0
       chrg(i)=0.0d0
       do j=1,25
        ind_beta(i,j)=0
        ind_d1(i,j)=0
        ind_d6(i,j)=0
        ind_d8(i,j)=0
        ind_d10(i,j)=0
        ind_C6(i,j)=0
        ind_C8(i,j)=0
        ind_C10(i,j)=0
       enddo
       do j=1,3
        sites(j,i)=0.0d0
       enddo
      enddo

      open (8,file='data_ccdata')

      read(8,*) cnull
      do i=1,25
        read(8,*) (sites(j,i),j=1,3)
      enddo
      read(8,*) cnull
      read(8,*) (chrg(j),j=1,5)

      read(8,*) cnull
      read(8,*) (ind_charge(j),j=1,5)

      read(8,*) cnull
      do i=1,25
        read(8,*) (ind_beta(i,j),j=1,25)
      enddo
      read(8,*) cnull
      do i=1,5
        read(8,*) (ind_d1(i,j),j=1,5)
      enddo
      read(8,*) cnull
      do i=1,3
        read(8,*) (ind_d6(i,j),j=1,3)
      enddo
      read(8,*) cnull
      do i=1,3
        read(8,*) (ind_d8(i,j),j=1,3)
      enddo
      read(8,*) cnull
      do i=1,3
        read(8,*) (ind_d10(i,j),j=1,3)
      enddo
      read(8,*) cnull
      do i=1,3
        read(8,*) (ind_c6(i,j),j=1,3)
      enddo
      read(8,*) cnull
      do i=1,3
        read(8,*) (ind_c8(i,j),j=1,3)
      enddo
      read(8,*) cnull
      do i=1,3
        read(8,*) (ind_c10(i,j),j=1,3)
      enddo

      do i=1,6250
        ind_data(i)=ind_data1(i)
      enddo

      close(8)
      return
      end
c-------------------------------------------------------
