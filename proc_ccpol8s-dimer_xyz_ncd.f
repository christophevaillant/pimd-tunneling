c----------------------------------------------------------
      subroutine ccpol8s_dimer(imode,Oa,Ha1,Ha2,Ob,Hb1,Hb2,Erigid,
     .                         c,params,nparsall,
     .                         ind_data,chrg,sites)
      implicit real*8 (a-h,o-z) 
      PARAMETER (nsite=25,maxdat=3000,maxlin=2000)
      PARAMETER (maxpar=1000)
      parameter (nsites=2*nsite*maxdat)
      dimension params(maxpar)
c      save params,nparsall

      dimension c(maxlin)
      dimension ElA(3),ElB(3)
      dimension RA(3),RB(3)
      dimension T(3,3),vec(3)
      dimension sites(3,nsite)     ! site coordinates for standard molecule
      dimension sitesAB(3,nsites) ! site coords for each of translated and rotated 2*ndat molecules
      dimension aj(maxlin)
      dimension rsites(3,nsite)     ! site coordinates for standard molecule - robocza
      dimension Oa(3),Ha1(3),Ha2(3)
      dimension Ob(3),Hb1(3),Hb2(3)
      dimension COMa(3),COMb(3)

      integer omp_get_num_threads,omp_get_thread_num
      character*1 cnull
      dimension ind_data(6250)
      dimension chrg_omp(25)
      dimension sites_omp(3,25)

c      save sites 

      zero=0.d0
      bohr2a=0.529177249d0
      h2kcal=627.510d0

      pi180=dacos(-1.d0)/180.d0

c        id=omp_get_thread_num()

      if (imode.eq.-1) then
      open (7,file='data_CCpol8s')
      read (7,*) nparsall
      if (nparsall.gt.maxpar) stop 010
      do 10 i=1,nparsall
         read (7,*) ii,params(i)
         if (ii.ne.i) stop 020
 10   continue

      read (7,*) nlin0
      if (nlin0.gt.maxlin) stop 030
      do 20 i=1,nlin0
         read (7,*) ii,c(i)
         if (ii.ne.i) stop 040
 20   continue
      call read_cc_data(ind_data,chrg,sites)
      close(7)
      return
      endif

c      write (*,'(a)') 'dimer energies (hartree,kcal/mol)'
c      read (*,*) ndat
c      if (ndat.gt.maxdat) stop 111

        do j=1,3
          Oa(j)=Oa(j)/bohr2a
          Ha1(j)=Ha1(j)/bohr2a
          Ha2(j)=Ha2(j)/bohr2a
          Ob(j)=Ob(j)/bohr2a
          Hb1(j)=Hb1(j)/bohr2a
          Hb2(j)=Hb2(j)/bohr2a
        enddo

c If our program is used independently, use this subroutine to check whether our
c molecule is indeed the water in the reference geometry
c      call check_molecule(Oa,Ha1,Ha2)   
c      call check_molecule(Ob,Hb1,Hb2)   

        call fill_sites(Oa,Ha1,Ha2,sites,rsites)
        indA = 0
        do ns=1,nsite
          do jj=1,3
            sitesAB(jj,indA+ns) = rsites(jj,ns)
          enddo
        enddo

        call fill_sites(Ob,Hb1,Hb2,sites,rsites)
        indB = nsite
        do ns=1,nsite
          do jj=1,3
            sitesAB(jj,indB+ns) = rsites(jj,ns)
          enddo
        enddo

      itwo=2
      call indN_iter (itwo, sitesAB(1,indA+1), nsite, Eind, chrg)

      Eind = Eind           ! hartree
c      Eind = Eind*h2kcal    ! kcal/mol

      indA = 1
      indB = nsite +1
      call U0(nsite,sitesAB(1,indA),sitesAB(1,indB),dummy,a0,aj,nlin,
     .        params,nparsall,
     .        ind_data,chrg,sites)

      E=Eind
      do 210 nl=1,nlin
        E=E +c(nl)*aj(nl)
 210  continue
      E=E+a0

c      write (*,'(i5,2f20.8)') nd,E,E*h2kcal
      Erigid=E*h2kcal

      return
      end
c----------------------------------------------------------------------
      subroutine U0 (nsite,sitesA,sitesB,energy, a0,aj,nlin,
     .               params,nparsall,
     .               ind_data,chrg,sites)
      implicit real*8 (a-h,o-z)
      PARAMETER (maxpar=1000)

c      dimension sitesA(3,*),sitesB(3,*), aj(*)
      dimension sitesA(3,25),sitesB(3,25), aj(144)
      dimension params(maxpar)

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
      integer omp_get_num_threads,omp_get_thread_num

      zero=0.d0
      h2kcal=627.510d0
      do i=1,6250
        ind_data1(i)=ind_data(i)
      enddo
      
      nlin=144

      do 5 i=1,nlin
      aj(i)=zero
 5    continue

      E_exp=zero
      E_ele=zero
      E_ind=zero

      do 10 nsA=1,nsite
      do 10 nsB=1,nsite

      call distan (sitesA(1,nsA),sitesB(1,nsB),R)

      if (ind_beta(nsA,nsB).ne.0) then     ! exponential
         beta = params(ind_beta(nsA,nsB))
         eks = dexp(-beta*R)
         indlin=ind_beta(nsA,nsB)-98
         if (indlin.lt.0) indlin=indlin+65
         ind0= indlin
         ind1= ind0+36
         ind2= ind1+36
         ind3= ind2+36
         aj(ind0) = aj(ind0) +eks
         aj(ind1) = aj(ind1) +eks*R
         aj(ind2) = aj(ind2) +eks*R*R
         aj(ind3) = aj(ind3) +eks*R*R*R
c         E_exp=...
      end if

      if (ind_charge(nsA)*ind_charge(nsB).ne.0) then     ! elst
         qA = params(ind_charge(nsA))
         qB = params(ind_charge(nsB))
         d1 = params(ind_d1(nsA,nsB))
         f1 = damp(1,d1,R)
         E_ele = E_ele+ f1* qA*qB/R
      end if

      if (ind_d6(nsA,nsB).ne.0) then     ! ind-disp
         d6 = params(ind_d6(nsA,nsB))
         d8 = params(ind_d8(nsA,nsB))
         d10 = params(ind_d10(nsA,nsB))
         C6 = params(ind_C6(nsA,nsB))
         C8 = params(ind_C8(nsA,nsB))
         C10 = params(ind_C10(nsA,nsB))
         f6 = damp(6,d6,R)
         f8 = damp(8,d8,R)
         f10 = damp(10,d10,R)
         R2 = R*R
         R6 = R2*R2*R2
         R8 = R6*R2
         R10 = R8*R2
         E_ind = E_ind -f6*C6/R6 -f8*C8/R8 -f10*C10/R10
      end if

 10   continue

c      energy = E_exp +E_ele +E_ind
      a0 = E_ele +E_ind

c!test1
c      write (8,*) 'E_exp=',E_exp*h2kcal
c      write (8,*) 'E_ele=',E_ele*h2kcal
c      write (8,*) 'E_ind=',E_ind*h2kcal
c      write (8,*) 'energy=',energy*h2kcal
c      sum=zero
c      do 200 i=1,40
c      sum=sum + aj(i)*params(i+43)
c 200  continue
c      write (8,*) 'test:',sum*h2kcal
c!test0      


      end
c---------------------------------------------------------------------
      subroutine indN_iter(N,sites,nsite,energy,chrg)
c Compute the iterative induction between N water molecules.
c The electrostatic field supplied by routine efield.
c
c Assume that the dipole center is at NOT the center
c of mass, but in special point.
      implicit real*8 (a-h,o-z)
      PARAMETER (Nmax=1024)
      PARAMETER (maxsite=30)
      parameter (maxNsite=Nmax*maxsite)
      parameter (maxit=200)

      dimension sites(3,maxNsite) ! site coords for all molecules
      dimension R(3,Nmax)         ! positions of polarizable sites
      dimension G2(3,Nmax)        ! induced dipoles
      dimension E0x(Nmax),E0y(Nmax),E0z(Nmax) ! fields from static charges
      dimension epom(3)
      dimension dist(Nmax,Nmax)   ! distances between pol.centers (**(-3))

      dimension chrg(25)
      integer omp_get_num_threads,omp_get_thread_num

c        ns=omp_get_num_threads()
c        id=omp_get_thread_num()
c        nn=30+id
c        write(nn,*) "in indN_iter  0  id=",id
c        call flush(nn)
      
      pol=9.922d0                ! exprmt average polarizability
      sig=0.367911875040999981d0 ! position of polarizable center 
                                 ! same as in 3B fits
      plen=1.1216873242d0
      dmpfct=1.d0
      zero=0.d0

      do j=1,N
       do i=1,3
        G2(i,j) = zero
       enddo
      enddo

c Compute positions of polarizable sites and the distance between them
c      dist(i,j)= 1.d0/rij**3

      do 10 i=1,N
      ind1= (i-1)*nsite +1
      ind2= (i-1)*nsite +2
      ind3= (i-1)*nsite +3
      do 11 ii=1,3
       pom = 0.5d0* (sites(ii,ind2) + sites(ii,ind3))
       R(ii,i) = sites(ii,ind1) + sig*(pom-sites(ii,ind1))/plen
 11   continue
 10   continue

      do 20 i=1,N
      do 21 j=i+1,N
      dist(i,j) = 0.d0
      do 25 ii=1,3
 25   dist(i,j) = dist(i,j) +(R(ii,i)-R(ii,j))**2.d0
      dist(i,j)=dist(i,j)**(-1.5d0)
      dist(j,i)=dist(i,j)
 21   continue
 20   continue

c
c Permanent charges fields
c
      do 30 i=1,N
      E0x(i) = zero
      E0y(i) = zero
      E0z(i) = zero
      do 35 j=1,N
      if (j.eq.i) goto 35
      index= (j-1)*nsite +1
      call efield_bohr(R(1,i),sites(1,index),nsite,epom,chrg) ! field of j on i
      E0x(i) = E0x(i) +epom(1)
      E0y(i) = E0y(i) +epom(2)
      E0z(i) = E0z(i) +epom(3)
 35   continue
 30   continue


C   LETS DO IT (the iterations, that is...)

      thr_iter = 1.d-20
      change = 10.d0
      isteps = 0

      DO WHILE (change.gt.thr_iter.and.isteps.lt.maxit)   
c     -----------------------------------------------

      energy= 0.0d0
      change= 0.0d0

      do 100 i=1,N
       E1X = E0x(i)
       E1Y = E0y(i) 
       E1Z = E0z(i) 
       do 200 j=1,N
          if(j.eq.i) go to 200
c calculate the total electric field at G1 of molecule i
         call TTTprod(R(1,i),R(1,j),G2(1,j),dist(i,j),epom)
          E1X = E1X + dmpfct*epom(1)
          E1Y = E1Y + dmpfct*epom(2)
          E1Z = E1Z + dmpfct*epom(3)
 200   continue

c total field at ii ready. Place induced dipole....
      polE1x = pol*E1X
      polE1y = pol*E1Y
      polE1z = pol*E1Z

      change=(G2(1,i)-polE1X)**2
     >      +(G2(2,i)-polE1Y)**2
     >      +(G2(3,i)-polE1Z)**2 + change
      G2(1,i)=polE1X
      G2(2,i)=polE1Y
      G2(3,i)=polE1Z
      energy =-0.5d0*pol*(E1X*E0X(i)
     >                   +E1Y*E0Y(i)
     >                   +E1Z*E0Z(i)) + energy

 100  continue

c      write (*,*) isteps,energy
      isteps = isteps +1
      END DO     !  end while
c     -----------------------

      if (isteps.ge.maxit) then
         write (*,*) 'No convergence in indN_iter'
         write (*,'(a,g12.3)') 'energy change=',change
         write (*,'(a,g12.3)') 'thr_iter=',thr_iter
c$$$         stop
         energy=0.0d0
      end if

c      energy = 627.510d0*energy   ! hartree--->kcal/mol
      end
c----------------------------------------------------------
c     
c Calculate the contribution of molecule b to the field on a
c Use the CC-pol site charges.
c Units: input positions of bohr, 
c output fields in au
c
      subroutine efield_bohr(veci,sitebt,nsiteb,e,chrg)
      implicit real*8 (a-h,o-z)
      PARAMETER (nsitemax=30)

      dimension veci(3),sitebt(3,1),e(3),chrg(25)
      dimension sep(3,nsitemax),sepl(nsitemax)

      integer omp_get_num_threads,omp_get_thread_num
      dimension sites(3,25)
      character*1 cnull
c       save chrg
c      data a0 /0.529177249d0/
      a0=1.d0
      crgfct=18.22262373d0

c
c Position of polarizable center of A is in veci(*)
      do isite=1,nsiteb
       if(chrg(isite).ne.0.d0) then   ! if charged site
       sepl(isite) = 0.d0
       do k=1,3
        sep(k,isite) = veci(k) - sitebt(k,isite)
        sepl(isite) = sepl(isite) + sep(k,isite)*sep(k,isite)
       end do
       sepl(isite) = sepl(isite)**(-1.5d0)
       endif                            ! if charged site
      end do

      do k=1,3
       e(k) = 0.d0
      end do

      do isite=1,nsiteb
       if(chrg(isite).ne.0.d0) then   ! if charged site
        do k=1,3
         e(k) = e(k) + a0*a0*chrg(isite)*sep(k,isite)*sepl(isite)
        end do
       endif      ! if charged site
      end do
c
      return
      end
c----------------------------------------------------------
      subroutine Av (A,v,b)
c calculates A*v=b
      implicit real*8 (a-h,o-z)
      dimension A(3,3),v(3),b(3)
      zero=0.d0

      do 10 i=1,3
      sum=zero
      do 20 j=1,3
      sum=sum +a(i,j)*v(j)
 20   continue
      b(i)=sum
 10   continue

      end
c----------------------------------------------------------
      subroutine distan (R1,R2,d)
      implicit real*8 (a-h,o-z)
      dimension R1(3),R2(3)

      d=0.d0
      do 10 i=1,3
      R12i=R1(i)-R2(i)
      d=d +R12i*R12i
 10   continue
      d=sqrt(d)
      end
c----------------------------------------------------------
      function damp (n,beta,r)
c
c     calculate the damping factor (small R correct)
c
      implicit real*8 (a-h,o-z)
      br=beta*r
c The following line added by RB, Sept. 18 1997
      if(br.eq.0.d0) then
       damp = 0.d0
       return
      endif
      sum=1.0d0
      term=1.0d0
      ncn=n
      do i=1,ncn
        term=term*br/i
        sum=sum+term
      enddo
      damp=1.0d0 - dexp(-br)*sum
c     in case of d --> 0 use
c     d=1.0d0 - dexp(-br)*sum = sum_m=ncn+1^\infty br^m/m!
      if(dabs(damp).lt.1.0d-8) then
        damp=0.0d0
        do i=ncn+1,1000
          term=term*br/i
          damp=damp+term
          if(term/damp .lt. 1.0d-8) go to 111
        enddo
        write(6,*) 'No convergence in damp'
  111 continue
      damp=damp*dexp(-br)
      endif
c     write(6,'(i4,2f10.5,e20.10)') n,beta,r,d
      return
      end
c----------------------------------------------------------
      subroutine fill_sites(O,H1,H2,sites,rsites)
      implicit real*8 (a-h,o-z)
      PARAMETER (nsite=25,maxdat=3000,maxlin=2000)
      dimension sites(3,nsite)     ! site coordinates for standard molecule
      dimension rsites(3,nsite)  
      dimension O(3),H1(3),H2(3)
      dimension COM(3)
      dimension ex(3),ey(3),ez(3)
      dimension v1(3),v2(3),v3(3)

c      save sites 

      bohr2a=0.529177249d0

c ss1 - bond length
c      dv1pv2=2.0d0*ss1*dcos(theta/2.0d0)
c      dv1mv2=2.0d0*ss1*dsin(theta/2.0d0)
      dv1pv2=1.99230765895d0
      dv1mv2=2.907303924565d0

      pi180=dacos(-1.d0)/180.d0

c First prepare two vectors, which span the plane defined by the molecule.
c Let's call them ez and ex, and require that they are oriented with respect to 
c the molecule in a standard way.
c ez is a bisector vector  

      call COMcalc(O,H1,H2,COM)

      do j=1,3
        v1(j)=H1(j)-COM(j)
        v2(j)=H2(j)-COM(j)
      enddo
      do j=1,3
        ez(j)=-(v1(j)+v2(j))
        ex(j)=v2(j)-v1(j)
      enddo
      do j=1,3
        ez(j)=ez(j)/dv1pv2
        ex(j)=ex(j)/dv1mv2
      enddo
      call cross(ey,ez,ex)

      do kk=1,nsite
        rsites(1,kk)=
     .       ex(1)*sites(1,kk) + ey(1)*sites(2,kk) + ez(1)*sites(3,kk)
        rsites(2,kk)=
     .       ex(2)*sites(1,kk) + ey(2)*sites(2,kk) + ez(2)*sites(3,kk)
        rsites(3,kk)=
     .       ex(3)*sites(1,kk) + ey(3)*sites(2,kk) + ez(3)*sites(3,kk)
c        write(6,666) (rsites(j,kk),j=1,3)," rsites(*,kk)"
      enddo

      do kk=1,nsite
        do j=1,3
          rsites(j,kk)=rsites(j,kk)+COM(j)
        enddo
      enddo

  666 format(3f18.12,a) 
      return
      end
c----------------------------------------------------------
      subroutine cross (z,x,y)
c calculates z = vector cross poduct of x an y
      implicit real*8 (a-h,o-z)
      dimension x(3),y(3),z(3)

      z(1)=x(2)*y(3)-x(3)*y(2)
      z(2)=x(3)*y(1)-x(1)*y(3)
      z(3)=x(1)*y(2)-x(2)*y(1)
      end
c------------------------------------------------------------------
      subroutine COMcalc (O1,H1,H2,COM)
      implicit real*8 (a-h,o-z)
      dimension O1(3),H1(3),H2(3),COM(3)
      real*8 mO,mH,M
      mO=15.994 914 6221d0
      mH=1.007 825 032 1d0

      M=mO+mH+mH
      do 10 i=1,3
      COM(i)= (mO*O1(i) +mH*H1(i) +mH*H2(i))/M
 10   continue

      end
c------------------------------------------------------------------
