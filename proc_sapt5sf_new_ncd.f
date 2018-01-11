       subroutine driver_potss_sapt5sf(carta,cartb,val,ipotparts,
     .                            param,parab,nsitea,nsiteb,c)
       implicit real*8 (a-h,o-z)
      parameter (naamax=14,nbbmax=14)
      parameter (nsitemax=8,ntypemax=6,ntypemax2=ntypemax*ntypemax)
      parameter (ntypemax_9=ntypemax*9,ntypemax2_9=ntypemax2*9)
      parameter (maxpar1=18,maxpar2=84)
c
       dimension carta(3,3),cartb(3,3)
      dimension param(maxpar1,ntypemax),parab(maxpar2,ntypemax,ntypemax)
      dimension c(1000)

      a0=0.529177249d0

       icart=1  ! choose this flag if the cartesian coordinates are provided
                ! in the following format 
                ! (for (O1,H1,H2) interacting with (O2,H3,H4)):
                ! xO1  yO1  zO1
                ! xH1  yH1  zH1
                ! xH2  yH2  zH2
                ! xO2  yO2  zO2
                ! xH3  yH3  zH3
                ! xH4  yH4  zH4
                ! Angstrom is a default unit

        pi = dacos(-1.d0)
        d2rad = pi/180.d0
        rad2d = 1.d0/d2rad

        r0=0.9716257d0
        theta0=104.69d0

      imode = 0
      alphaAv=0.0d0

        do ii=1,3
          do jj=1,3
            carta(ii,jj)=carta(ii,jj)/a0
            cartb(ii,jj)=cartb(ii,jj)/a0
          enddo
        enddo

        call poten(r1,r2,theta1,r3,r4,theta2,
     .      rin,alphaAv,betaAv,gammaAv,alphaBv,betaBv,gammaBv,val,imode,
     .      icart,carta,cartb,ipotparts,
     .        param,parab,nsitea,nsiteb,c)

201   format(2(1x,3(1x,f9.4)),7(1x,f9.4),f15.8,f5.1)
202   format(1x,3f15.8/1x,3f15.8/1x,3f15.8/
     .       1x,3f15.8/1x,3f15.8/1x,3f15.8,f15.8)

      return
      end
c========================================================================
      subroutine poten(r1_ang,r2_ang,theta1_deg,
     .                 r3_ang,r4_ang,theta2_deg,
     .        R,alphaAv,betaAv,gammaAv,alphaBv,betaBv,gammaBv,val,mode,
     .        icart,carta,cartb,ipotparts,
     .        param,parab,nsitea,nsiteb,c)
      implicit real*8 (a-h,o-z)
c
      parameter (naamax=14,nbbmax=14)
      parameter (nsitemax=8,ntypemax=6,ntypemax2=ntypemax*ntypemax)
      parameter (ntypemax_9=ntypemax*9,ntypemax2_9=ntypemax2*9)
      parameter (maxpar1=18,maxpar2=84)
c
      dimension oa(3), ob(3), values(100)  !     ???
c
      dimension siteat(3,nsitemax,-1:naamax),
     .          sitebt(3,nsitemax,-1:nbbmax)
      dimension sa(3,-1:naamax),sb(3,-1:nbbmax)
      dimension itypea(nsitemax),itypeb(nsitemax)
      dimension itypus(ntypemax,ntypemax,2)
      dimension numtm(2)   ! to be set in potparts....
      dimension sitex(3,nsitemax),sx(3)
      dimension carta(3,3),cartb(3,3)
c      save pi,d2rad,rad2d,iaa,ibb

      dimension param(maxpar1,ntypemax),parab(maxpar2,ntypemax,ntypemax)
      dimension c(1000)
      integer omp_get_num_threads,omp_get_thread_num
      character*50 nameofsaptdatafile
c      save param,parab,nsitea,nsiteb,c

c        ns=omp_get_num_threads()
c        id=omp_get_thread_num()
c
c If mode = -1 -- read input
c
      if(mode.eq.-1) then
        call data1(param,parab,nsitea,nsiteb,c)
       return
      endif

        pi = dacos(-1.d0)
        d2rad = pi/180.d0
        rad2d = 1.d0/d2rad
        iaa=-1
        ibb=-1

       call set_sites(carta,sitex,sx,itypea,nsitea)

       do j=1,3
        sa(j,iaa)=sx(j)
       end do
       do ia=1,nsitea
        do j=1,3
         siteat(j,ia,iaa)=sitex(j,ia)
        end do
       end do

       call set_sites(cartb,sitex,sx,itypeb,nsiteb)

       do j=1,3
        sb(j,ibb)=sx(j)
       end do
       do ib=1,nsiteb
        do j=1,3
         sitebt(j,ib,ibb)=sitex(j,ib)
        end do
       end do

c
c Zero out the used-type matrix
c
      do k1=1,ntypemax
       do k2=1,ntypemax
        itypus(k1,k2,1) = 0
        itypus(k1,k2,2) = 0
       end do
      end do
c
c Loop over sites in A and B
c
      iii = 1
      val = 0.d0
c
      do ia=1,nsitea
       do ib=1,nsiteb
        valp = 0.d0
c
c Compute the distance between the two sites
c
        rij = dist(siteat(1,ia,iaa),sitebt(1,ib,ibb))
c        write(44,*) "r(",ia,",",ib,")=",rij
c
c Compute the "basis functions" values pertaining to the
c required form of the potential. numtm(1) is the number of terms 
c returned for the "symmetric" part, numtm(2) - number
c of antisymmetric terms, numt = total number of terms, including constant part 
c
       if (ipotparts.eq.1) then
        call potparts(iaa,ibb,rij,ia,ib,numt,numtm,values,
     .                c,sa,sb,param,parab,itypea,itypeb)
        elseif(ipotparts.eq.0) then
        call potparts_old(iaa,ibb,rij,ia,ib,numt,numtm,values,
     .                c,sa,sb,param,parab,itypea,itypeb)
       else
        write(6,*) "wrong value of ipotparts:",ipotparts
        stop
       endif
c
c If abs(ntpot)>10, the last values is the fixed part.
c
        ntpot=124161
        numt0 = numt
        if(abs(ntpot).gt.10) then
         numt0 = numt-1
         valp = valp + values(numt)
        endif
c
c decide which basis function it is (basis function is determined
c by the pair of site types and the position in values vector).
c iii is the first basis function for this pair of types. It is
c updated into the next available position for the next pair of types.
c
c RB !!!!!!! change the following if numtm matrix added to potparts !!!
c
c GM+RB: additional loop over permutational symmetry
c
        itsmax = 1
        if(itypea(ia).ne.itypeb(ib)) itsmax = 3
        do its=1,itsmax,2
c
         itsm = min(its,2)
         itu = itypus(itypea(ia),itypeb(ib),itsm)
         if(itu.eq.0) then
          itypus(itypea(ia),itypeb(ib),itsm) = iii
          itypus(itypeb(ib),itypea(ia),itsm) = iii
          itu = iii
          iii = iii + numtm(itsm)
         endif
c
         do i=1,numtm(itsm)
          itu1 = itu+i-1
          if(itsm.eq.1) then
           valp = valp + c(itu1)*values(i)
          else
           valp = valp + c(itu1)*values(i+numtm(1))
          endif
         end do
c
        end do   ! loop over its
        val = val + valp

       end do
      end do
c
c Don't forget to add asymptotics when the fitting part tested...
c
      call dipind(iaa,ibb,R,fcind,
     .     sa,sb,param,parab,siteat,sitebt,nsitea,nsiteb,itypea,itypeb)
      val = val + fcind

c         write(94,200)
c     .    iaa,r1_ang,r2_ang,theta1_deg,(sa(j,iaa),j=1,3),
c     .    ibb,r3_ang,r4_ang,theta2_deg,(sb(j,ibb),j=1,3),
c     .    id,R,(oa(j)*rad2d,j=2,3),(ob(j)*rad2d,j=1,3),val,enfitpj
200   format(2(1(1x,i4),6(1x,f9.4)),1(1x,i6),6(1x,f9.4),f15.8,f5.1)

      return
      end
c============================================================================
c
      function dist(a,b)
      implicit real*8 (a-h,o-z)
      dimension a(3),b(3)
      ttt = 0.d0
      do i=1,3
       diff = a(i)-b(i)
       ttt = ttt + diff*diff
      end do
      dist = dsqrt(ttt)
      return
      end
c============================================================================
c
      subroutine potparts(iaa,ibb,rij,ia,ib,numt,numtm,values,
     .                    c,sa,sb,param,parab,itypea,itypeb)
      implicit real*8 (a-h,o-z)
c
      parameter (naamax=14,nbbmax=14)
      parameter (nsitemax=8,ntypemax=6,ntypemax2=ntypemax*ntypemax)
      parameter (ntypemax_9=ntypemax*9,ntypemax2_9=ntypemax2*9)
      parameter (maxpar1=18,maxpar2=84)
c
      dimension values(100),numtm(2)
      dimension param(maxpar1,ntypemax),parab(maxpar2,ntypemax,ntypemax)
      dimension itypea(nsitemax),itypeb(nsitemax)
      dimension c(1000)
      dimension sa(3,-1:naamax),sb(3,-1:nbbmax)
       
      a0=0.529177249d0
      au2cal=627.510d0
      
c
c ***********************************************************
c GM-modified: more linear basis functions inc squares (but now with couplings
c between monomers--differs from 24161 in that accidently neglected a3 
c coplings put in)
c RB-modified
c ntpot = 124161: exponential*(1+sum r^n, n=1...3) + r^-(6-8-10) + elst, totally
c non-linear      FLEXIBILIZED
c ***********************************************************
      ntpot=124161
      if(ntpot.eq.124161) then
c
c Previously fitted rigid contributions
c hen dealing with type=6, assume the parameters of H.
c
      ita = itypea(ia)
      itb = itypeb(ib)
c      if(ita.eq.6) ita = 2
c      if(itb.eq.6) itb = 2
      beta = parab(1,ita,itb)
      alpha = parab(2,ita,itb)
      a = dexp(alpha)
c Rigid C6,C8,C10 understood as symmetric in s_A, s_B
      c6 = parab(3,itypea(ia),itypeb(ib))
      c8 = parab(4,itypea(ia),itypeb(ib))
      c10 = parab(5,itypea(ia),itypeb(ib))
      dmp1 = parab(6,itypea(ia),itypeb(ib))
      dmp6 = parab(7,itypea(ia),itypeb(ib))
      dmp8 = parab(8,itypea(ia),itypeb(ib)) 
      dmp10= parab(9,itypea(ia),itypeb(ib))
c parameter index numbers updated
      a1 = parab(38,itypea(ia),itypeb(ib))
      a2 = parab(39,itypea(ia),itypeb(ib))
      a3 = parab(40,itypea(ia),itypeb(ib))
c ccccc
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
c
c and previously fitted flexible contributions
c
      s1=sa(1,iaa)
      s2=sa(2,iaa)
      s3=sa(3,iaa)
      s4=sb(1,ibb)
      s5=sb(2,ibb)
      s6=sb(3,ibb)
      signa=1.0d0
      signb=1.0d0
      if (ia .eq. 3) signa=-1.0d0
      if (ib .eq. 3) signb=-1.0d0
      s3 = signa*s3
      s6 = signb*s6
      qa
     & =qa
     & +param(2,itypea(ia))*s1
     & +param(3,itypea(ia))*s2
     & +param(4,itypea(ia))*s3
     & +param(5,itypea(ia))*s1*s2
     & +param(6,itypea(ia))*s2*s3
     & +param(7,itypea(ia))*s1*s1
     & +param(8,itypea(ia))*s2*s2
     & +param(9,itypea(ia))*s3*s3
      qb
     & =qb
     & +param(2,itypeb(ib))*s4
     & +param(3,itypeb(ib))*s5
     & +param(4,itypeb(ib))*s6
     & +param(5,itypeb(ib))*s4*s5
     & +param(6,itypeb(ib))*s5*s6
     & +param(7,itypeb(ib))*s4*s4
     & +param(8,itypeb(ib))*s5*s5
     & +param(9,itypeb(ib))*s6*s6
c Make basis functions even in s3, s6 if site is NOT H but not for 
c charges--thus charges have been moved above this--or in dipind 
c since charges and dipole in dipind were not fited with this definition
c      if(itypea(ia).ne.2) s3 = dabs(s3)
c      if(itypeb(ib).ne.2) s6 = dabs(s6)
      if(itypea(ia).ne.2) s3 = s3*s3
      if(itypeb(ib).ne.2) s6 = s6*s6
c Add the flexible part of beta here....
      if(itypea(ia).eq.itypeb(ib)) then
       beta = beta + parab(41,itypea(ia),itypeb(ib))*(s3+s6)
       beta = beta + parab(46,itypea(ia),itypeb(ib))*(s3*s3+s6*s6)
      elseif(itypea(ia).lt.itypeb(ib)) then
       beta = beta + parab(41,itypea(ia),itypeb(ib))*s3
       beta = beta + parab(42,itypea(ia),itypeb(ib))*s6
       beta = beta + parab(46,itypea(ia),itypeb(ib))*s3*s3
       beta = beta + parab(47,itypea(ia),itypeb(ib))*s6*s6
      elseif(itypea(ia).gt.itypeb(ib)) then
       beta = beta + parab(41,itypea(ia),itypeb(ib))*s6
       beta = beta + parab(42,itypea(ia),itypeb(ib))*s3
       beta = beta + parab(47,itypea(ia),itypeb(ib))*s3*s3
       beta = beta + parab(46,itypea(ia),itypeb(ib))*s6*s6
      endif
      beta = dabs(beta)   ! make it positive
c end of flex beta part....
c Add the flexible part of alpha here....
      if(itypea(ia).eq.itypeb(ib)) then
       alpha = alpha + parab(43,itypea(ia),itypeb(ib))*(s3+s6)
       alpha = alpha + parab(48,itypea(ia),itypeb(ib))*(s3*s3+s6*s6)
      elseif(itypea(ia).lt.itypeb(ib)) then
       alpha = alpha + parab(43,itypea(ia),itypeb(ib))*s3
       alpha = alpha + parab(44,itypea(ia),itypeb(ib))*s6
       alpha = alpha + parab(48,itypea(ia),itypeb(ib))*s3*s3
       alpha = alpha + parab(49,itypea(ia),itypeb(ib))*s6*s6
      elseif(itypea(ia).gt.itypeb(ib)) then
       alpha = alpha + parab(43,itypea(ia),itypeb(ib))*s6
       alpha = alpha + parab(44,itypea(ia),itypeb(ib))*s3
       alpha = alpha + parab(48,itypea(ia),itypeb(ib))*s6*s6
       alpha = alpha + parab(49,itypea(ia),itypeb(ib))*s3*s3
      endif
      a = dexp(alpha)   ! make it positive
c end of flex alpha part....
      d1 = d(1,dmp1,rij)
      d6 = d(6,dmp6,rij)
      d8 = d(8,dmp8,rij)
      d10 = d(10,dmp10,rij)
      c6
     & =c6
     & +parab(11,itypea(ia),itypeb(ib))*(s3+s6)
     & +parab(14,itypea(ia),itypeb(ib))*(s1+s4)
     & +parab(17,itypea(ia),itypeb(ib))*(s2+s5)
     & +parab(20,itypea(ia),itypeb(ib))*(s3*s6)
     & +parab(23,itypea(ia),itypeb(ib))*(s1*s4)
     & +parab(26,itypea(ia),itypeb(ib))*(s2*s5)
      c8
     & =c8
     & +parab(12,itypea(ia),itypeb(ib))*(s3+s6)
     & +parab(15,itypea(ia),itypeb(ib))*(s1+s4)
     & +parab(18,itypea(ia),itypeb(ib))*(s2+s5)
     & +parab(21,itypea(ia),itypeb(ib))*(s3*s6)
     & +parab(24,itypea(ia),itypeb(ib))*(s1*s4)
     & +parab(27,itypea(ia),itypeb(ib))*(s2*s5)
      c10
     & =c10
     & +parab(13,itypea(ia),itypeb(ib))*(s3+s6)
     & +parab(16,itypea(ia),itypeb(ib))*(s1+s4)
     & +parab(19,itypea(ia),itypeb(ib))*(s2+s5)
     & +parab(22,itypea(ia),itypeb(ib))*(s3*s6)
     & +parab(25,itypea(ia),itypeb(ib))*(s1*s4)
     & +parab(28,itypea(ia),itypeb(ib))*(s2*s5)
c Update c6, c8, c10 with asymptotic terms asymmetric in s_A, s_B:
      c6as = 0.d0
      c8as = 0.d0
      c10as = 0.d0
      if(itypea(ia).ne.itypeb(ib)) then
      c6as
     & =c6as
     & +parab(29,itypea(ia),itypeb(ib))*(s3-s6)
     & +parab(32,itypea(ia),itypeb(ib))*(s1-s4)
     & +parab(35,itypea(ia),itypeb(ib))*(s2-s5)
      c8as
     & =c8as
     & +parab(30,itypea(ia),itypeb(ib))*(s3-s6)
     & +parab(33,itypea(ia),itypeb(ib))*(s1-s4)
     & +parab(36,itypea(ia),itypeb(ib))*(s2-s5)
      c10as
     & =c10as
     & +parab(31,itypea(ia),itypeb(ib))*(s3-s6)
     & +parab(34,itypea(ia),itypeb(ib))*(s1-s4)
     & +parab(37,itypea(ia),itypeb(ib))*(s2-s5)
       if(itypea(ia).gt.itypeb(ib)) then
        c6as = -c6as
        c8as = -c8as
        c10as = -c10as
       endif
      endif
      c6 = c6 + c6as
      c8 = c8 + c8as
      c10 = c10 + c10as
c
      if (beta.gt.0.d0) then 
      if(itypea(ia).ne.2.and.itypeb(ib).ne.2) then ! no H involved
       if(itypea(ia).eq.itypeb(ib)) then
        numtm(1) = 40
        numtm(2) =  0 !28
       else   ! test: make it the same as symmetric
        numtm(1) = 40
        numtm(2) = 28
       endif
      else
       if(itypea(ia).eq.itypeb(ib)) then
        numtm(1) = 40
        numtm(2) =  0 !28  ! no s-expansion if none of sites is H
       else
        numtm(1) = 40
        numtm(2) = 28
       endif
      endif
      numt = numtm(1) + numtm(2) + 1
c
      val0=a*dexp(-beta*rij)
      val1=val0*rij
      val2=val1*rij
      val3=val2*rij
       values(numt) = val0 + a1*val1 + a2*val2 + a3*val3
     & + d1*qa*qb/rij
     & - d6*c6/rij**6 - d8*c8/rij**8 - d10*c10/rij**10
c
c Construct planning matrix by filling values vector
c s-Symmetric terms
      if(itypea(ia).ne.2.and.itypeb(ib).ne.2) then ! no H involved
c delta a0,1,2,3 contributions to planning matrix
      values(1)  = (s1+s4)*val0
      values(2)  = (s1+s4)*val1
      values(3)  = (s1+s4)*val2
      values(4)  = (s1+s4)*val3

      values(5)  = (s2+s5)*val0
      values(6)  = (s2+s5)*val1
      values(7)  = (s2+s5)*val2
      values(8)  = (s2+s5)*val3

      values(9)  = (s3+s6)*val0
      values(10) = (s3+s6)*val1
      values(11) = (s3+s6)*val2
      values(12) = (s3+s6)*val3


      values(13) = (s1*s2+s4*s5)*val0
      values(14) = (s1*s2+s4*s5)*val1
      values(15) = (s1*s2+s4*s5)*val2     !x
      values(16) = (s1*s2+s4*s5)*val3     !x

      values(17) = (s2*s3+s5*s6)*val0
      values(18) = (s2*s3+s5*s6)*val1
      values(19) = (s2*s3+s5*s6)*val2
      values(20) = (s2*s3+s5*s6)*val3

      values(21) = (s1*s1+s4*s4)*val0
      values(22) = (s1*s1+s4*s4)*val1
      values(23) = (s1*s1+s4*s4)*val2
      values(24) = (s1*s1+s4*s4)*val3

      values(25) = (s2*s2+s5*s5)*val0
      values(26) = (s2*s2+s5*s5)*val1
      values(27) = (s2*s2+s5*s5)*val2
      values(28) = (s2*s2+s5*s5)*val3


      values(29) = (s1*s4)*val0
      values(30) = (s1*s4)*val1
      values(31) = (s1*s4)*val2
      values(32) = (s1*s4)*val3

      values(33) = (s2*s5)*val0
      values(34) = (s2*s5)*val1
      values(35) = (s2*s5)*val2
      values(36) = (s2*s5)*val3

      values(37) = (s3*s6)*val0
      values(38) = (s3*s6)*val1
      values(39) = (s3*s6)*val2
      values(40) = (s3*s6)*val3
c For uneqal types, add s-antisymmetric terms
      if(itypea(ia).lt.itypeb(ib)) then
      values(41) = (s1-s4)*val0
      values(42) = (s1-s4)*val1
      values(43) = (s1-s4)*val2
      values(44) = (s1-s4)*val3

      values(45) = (s2-s5)*val0
      values(46) = (s2-s5)*val1
      values(47) = (s2-s5)*val2
      values(48) = (s2-s5)*val3

      values(49) = (s3-s6)*val0
      values(50) = (s3-s6)*val1
      values(51) = (s3-s6)*val2
      values(52) = (s3-s6)*val3


      values(53) = (s1*s2-s4*s5)*val0
      values(54) = (s1*s2-s4*s5)*val1
      values(55) = (s1*s2-s4*s5)*val2     !x
      values(56) = (s1*s2-s4*s5)*val3     !x

      values(57) = (s2*s3-s5*s6)*val0
      values(58) = (s2*s3-s5*s6)*val1
      values(59) = (s2*s3-s5*s6)*val2
      values(60) = (s2*s3-s5*s6)*val3

      values(61) = (s1*s1-s4*s4)*val0
      values(62) = (s1*s1-s4*s4)*val1
      values(63) = (s1*s1-s4*s4)*val2
      values(64) = (s1*s1-s4*s4)*val3

      values(65) = (s2*s2-s5*s5)*val0
      values(66) = (s2*s2-s5*s5)*val1
      values(67) = (s2*s2-s5*s5)*val2
      values(68) = (s2*s2-s5*s5)*val3
      endif
      if(itypea(ia).gt.itypeb(ib)) then
      values(41) =-(s1-s4)*val0
      values(42) =-(s1-s4)*val1
      values(43) =-(s1-s4)*val2
      values(44) =-(s1-s4)*val3

      values(45) =-(s2-s5)*val0
      values(46) =-(s2-s5)*val1
      values(47) =-(s2-s5)*val2
      values(48) =-(s2-s5)*val3

      values(49) =-(s3-s6)*val0
      values(50) =-(s3-s6)*val1
      values(51) =-(s3-s6)*val2
      values(52) =-(s3-s6)*val3


      values(53) =-(s1*s2-s4*s5)*val0
      values(54) =-(s1*s2-s4*s5)*val1
      values(55) =-(s1*s2-s4*s5)*val2     !x
      values(56) =-(s1*s2-s4*s5)*val3     !x

      values(57) =-(s2*s3-s5*s6)*val0
      values(58) =-(s2*s3-s5*s6)*val1
      values(59) =-(s2*s3-s5*s6)*val2
      values(60) =-(s2*s3-s5*s6)*val3

      values(61) =-(s1*s1-s4*s4)*val0
      values(62) =-(s1*s1-s4*s4)*val1
      values(63) =-(s1*s1-s4*s4)*val2
      values(64) =-(s1*s1-s4*s4)*val3

      values(65) =-(s2*s2-s5*s5)*val0
      values(66) =-(s2*s2-s5*s5)*val1
      values(67) =-(s2*s2-s5*s5)*val2
      values(68) =-(s2*s2-s5*s5)*val3
      endif
      else    ! at least 1 H involved
c delta a0,1,2,3 contributions to planning matrix
      values(1)  = (s1+s4)*val0
      values(2)  = (s1+s4)*val1
      values(3)  = (s1+s4)*val2
      values(4)  = (s1+s4)*val3

      values(5)  = (s2+s5)*val0
      values(6)  = (s2+s5)*val1
      values(7)  = (s2+s5)*val2
      values(8)  = (s2+s5)*val3

      values(9)  = (s3+s6)*val0
      values(10) = (s3+s6)*val1
      values(11) = (s3+s6)*val2
      values(12) = (s3+s6)*val3


      values(13) = (s1*s2+s4*s5)*val0
      values(14) = (s1*s2+s4*s5)*val1
      values(15) = (s1*s2+s4*s5)*val2     !x
      values(16) = (s1*s2+s4*s5)*val3     !x

      values(17) = (s2*s3+s5*s6)*val0
      values(18) = (s2*s3+s5*s6)*val1
      values(19) = (s2*s3+s5*s6)*val2
      values(20) = (s2*s3+s5*s6)*val3

      values(21) = (s1*s1+s4*s4)*val0
      values(22) = (s1*s1+s4*s4)*val1
      values(23) = (s1*s1+s4*s4)*val2
      values(24) = (s1*s1+s4*s4)*val3

      values(25) = (s2*s2+s5*s5)*val0
      values(26) = (s2*s2+s5*s5)*val1
      values(27) = (s2*s2+s5*s5)*val2
      values(28) = (s2*s2+s5*s5)*val3


      values(29) = (s1*s4)*val0
      values(30) = (s1*s4)*val1
      values(31) = (s1*s4)*val2
      values(32) = (s1*s4)*val3

      values(33) = (s2*s5)*val0
      values(34) = (s2*s5)*val1
      values(35) = (s2*s5)*val2
      values(36) = (s2*s5)*val3

      values(37) = (s3*s6)*val0
      values(38) = (s3*s6)*val1
      values(39) = (s3*s6)*val2
      values(40) = (s3*s6)*val3
c For uneqal types, add s-antisymmetric terms
      if(itypea(ia).lt.itypeb(ib)) then
      values(41) = (s1-s4)*val0
      values(42) = (s1-s4)*val1
      values(43) = (s1-s4)*val2
      values(44) = (s1-s4)*val3

      values(45) = (s2-s5)*val0
      values(46) = (s2-s5)*val1
      values(47) = (s2-s5)*val2
      values(48) = (s2-s5)*val3

      values(49) = (s3-s6)*val0
      values(50) = (s3-s6)*val1
      values(51) = (s3-s6)*val2
      values(52) = (s3-s6)*val3


      values(53) = (s1*s2-s4*s5)*val0
      values(54) = (s1*s2-s4*s5)*val1
      values(55) = (s1*s2-s4*s5)*val2     !x
      values(56) = (s1*s2-s4*s5)*val3     !x

      values(57) = (s2*s3-s5*s6)*val0
      values(58) = (s2*s3-s5*s6)*val1
      values(59) = (s2*s3-s5*s6)*val2
      values(60) = (s2*s3-s5*s6)*val3

      values(61) = (s1*s1-s4*s4)*val0
      values(62) = (s1*s1-s4*s4)*val1
      values(63) = (s1*s1-s4*s4)*val2
      values(64) = (s1*s1-s4*s4)*val3

      values(65) = (s2*s2-s5*s5)*val0
      values(66) = (s2*s2-s5*s5)*val1
      values(67) = (s2*s2-s5*s5)*val2
      values(68) = (s2*s2-s5*s5)*val3
      endif
      if(itypea(ia).gt.itypeb(ib)) then
      values(41) =-(s1-s4)*val0
      values(42) =-(s1-s4)*val1
      values(43) =-(s1-s4)*val2
      values(44) =-(s1-s4)*val3

      values(45) =-(s2-s5)*val0
      values(46) =-(s2-s5)*val1
      values(47) =-(s2-s5)*val2
      values(48) =-(s2-s5)*val3

      values(49) =-(s3-s6)*val0
      values(50) =-(s3-s6)*val1
      values(51) =-(s3-s6)*val2
      values(52) =-(s3-s6)*val3


      values(53) =-(s1*s2-s4*s5)*val0
      values(54) =-(s1*s2-s4*s5)*val1
      values(55) =-(s1*s2-s4*s5)*val2     !x
      values(56) =-(s1*s2-s4*s5)*val3     !x

      values(57) =-(s2*s3-s5*s6)*val0
      values(58) =-(s2*s3-s5*s6)*val1
      values(59) =-(s2*s3-s5*s6)*val2
      values(60) =-(s2*s3-s5*s6)*val3

      values(61) =-(s1*s1-s4*s4)*val0
      values(62) =-(s1*s1-s4*s4)*val1
      values(63) =-(s1*s1-s4*s4)*val2
      values(64) =-(s1*s1-s4*s4)*val3

      values(65) =-(s2*s2-s5*s5)*val0
      values(66) =-(s2*s2-s5*s5)*val1
      values(67) =-(s2*s2-s5*s5)*val2
      values(68) =-(s2*s2-s5*s5)*val3
      endif
c
      endif   ! if H involved
 
      else
c
      numt=1
      numtm(1) = 0
      numtm(2) = 0
      values(numt) =
     & + d1*qa*qb/rij
     & - d6*c6/rij**6 - d8*c8/rij**8 - d10*c10/rij**10
c
      end if 
c
      return
      endif
      end
c
c============================================================================
c
      subroutine potparts_old(iaa,ibb,rij,ia,ib,numt,numtm,values,
     .                        c,sa,sb,param,parab,itypea,itypeb)
      implicit real*8 (a-h,o-z)
c
      parameter (naamax=14,nbbmax=14)
      parameter (nsitemax=8,ntypemax=6,ntypemax2=ntypemax*ntypemax)
      parameter (ntypemax_9=ntypemax*9,ntypemax2_9=ntypemax2*9)
      parameter (maxpar1=18,maxpar2=84)
c
      dimension values(100),numtm(2)
      dimension param(maxpar1,ntypemax),parab(maxpar2,ntypemax,ntypemax)
      dimension itypea(nsitemax),itypeb(nsitemax)
      dimension c(1000)
      dimension sa(3,-1:naamax),sb(3,-1:nbbmax)
       
      a0=0.529177249d0
      au2cal=627.510d0
      
c
c ***********************************************************
c GM-modified: more linear basis functions inc squares (but now with couplings
c between monomers--differs from 24161 in that accidently neglected a3 
c coplings put in)
c RB-modified
c ntpot = 124161: exponential*(1+sum r^n, n=1...3) + r^-(6-8-10) + elst, totally
c non-linear      FLEXIBILIZED
c ***********************************************************
      ntpot=124161
      if(ntpot.eq.124161) then
c
c Previously fitted rigid contributions
c hen dealing with type=6, assume the parameters of H.
c
      ita = itypea(ia)
      itb = itypeb(ib)
c      if(ita.eq.6) ita = 2
c      if(itb.eq.6) itb = 2
      beta = parab(1,ita,itb)
      alpha = parab(2,ita,itb)
      a = dexp(alpha)
c Rigid C6,C8,C10 understood as symmetric in s_A, s_B
      c6 = parab(3,itypea(ia),itypeb(ib))
      c8 = parab(4,itypea(ia),itypeb(ib))
      c10 = parab(5,itypea(ia),itypeb(ib))
      dmp1 = parab(6,itypea(ia),itypeb(ib))
      dmp6 = parab(7,itypea(ia),itypeb(ib))
      dmp8 = parab(8,itypea(ia),itypeb(ib)) 
      dmp10= parab(9,itypea(ia),itypeb(ib))
c parameter index numbers updated
      a1 = parab(38,itypea(ia),itypeb(ib))
      a2 = parab(39,itypea(ia),itypeb(ib))
      a3 = parab(40,itypea(ia),itypeb(ib))
c ccccc
      qa = param(1,itypea(ia))
      qb = param(1,itypeb(ib))
c
c and previously fitted flexible contributions
c
      s1=sa(1,iaa)
      s2=sa(2,iaa)
      s3=sa(3,iaa)
      s4=sb(1,ibb)
      s5=sb(2,ibb)
      s6=sb(3,ibb)
      signa=1.0d0
      signb=1.0d0
      if (ia .eq. 3) signa=-1.0d0
      if (ib .eq. 3) signb=-1.0d0
      s3 = signa*s3
      s6 = signb*s6
      qa
     & =qa
     & +param(2,itypea(ia))*s1
     & +param(3,itypea(ia))*s2
     & +param(4,itypea(ia))*s3
     & +param(5,itypea(ia))*s1*s2
     & +param(6,itypea(ia))*s2*s3
     & +param(7,itypea(ia))*s1*s1
     & +param(8,itypea(ia))*s2*s2
     & +param(9,itypea(ia))*s3*s3
      qb
     & =qb
     & +param(2,itypeb(ib))*s4
     & +param(3,itypeb(ib))*s5
     & +param(4,itypeb(ib))*s6
     & +param(5,itypeb(ib))*s4*s5
     & +param(6,itypeb(ib))*s5*s6
     & +param(7,itypeb(ib))*s4*s4
     & +param(8,itypeb(ib))*s5*s5
     & +param(9,itypeb(ib))*s6*s6
c Make basis functions even in s3, s6 if site is NOT H but not for 
c charges--thus charges have been moved above this--or in dipind 
c since charges and dipole in dipind were not fited with this definition
c      if(itypea(ia).ne.2) s3 = dabs(s3)
c      if(itypeb(ib).ne.2) s6 = dabs(s6)
      if(itypea(ia).ne.2) s3 = s3*s3
      if(itypeb(ib).ne.2) s6 = s6*s6
c Add the flexible part of beta here....
      if(itypea(ia).eq.itypeb(ib)) then
       beta = beta + parab(41,itypea(ia),itypeb(ib))*(s3+s6)
       beta = beta + parab(46,itypea(ia),itypeb(ib))*(s3*s3+s6*s6)
      elseif(itypea(ia).lt.itypeb(ib)) then
       beta = beta + parab(41,itypea(ia),itypeb(ib))*s3
       beta = beta + parab(42,itypea(ia),itypeb(ib))*s6
       beta = beta + parab(46,itypea(ia),itypeb(ib))*s3*s3
       beta = beta + parab(47,itypea(ia),itypeb(ib))*s6*s6
      elseif(itypea(ia).gt.itypeb(ib)) then
       beta = beta + parab(41,itypea(ia),itypeb(ib))*s6
       beta = beta + parab(42,itypea(ia),itypeb(ib))*s3
       beta = beta + parab(47,itypea(ia),itypeb(ib))*s3*s3
       beta = beta + parab(46,itypea(ia),itypeb(ib))*s6*s6
      endif
      beta = dabs(beta)   ! make it positive
c end of flex beta part....
c Add the flexible part of alpha here....
      if(itypea(ia).eq.itypeb(ib)) then
       alpha = alpha + parab(43,itypea(ia),itypeb(ib))*(s3+s6)
       alpha = alpha + parab(48,itypea(ia),itypeb(ib))*(s3*s3+s6*s6)
      elseif(itypea(ia).lt.itypeb(ib)) then
       alpha = alpha + parab(43,itypea(ia),itypeb(ib))*s3
       alpha = alpha + parab(44,itypea(ia),itypeb(ib))*s6
       alpha = alpha + parab(48,itypea(ia),itypeb(ib))*s3*s3
       alpha = alpha + parab(49,itypea(ia),itypeb(ib))*s6*s6
      elseif(itypea(ia).gt.itypeb(ib)) then
       alpha = alpha + parab(43,itypea(ia),itypeb(ib))*s6
       alpha = alpha + parab(44,itypea(ia),itypeb(ib))*s3
       alpha = alpha + parab(48,itypea(ia),itypeb(ib))*s6*s6
       alpha = alpha + parab(49,itypea(ia),itypeb(ib))*s3*s3
      endif
      a = dexp(alpha)   ! make it positive
c end of flex alpha part....
      d1 = d(1,dmp1,rij)
      d6 = d(6,dmp6,rij)
      d8 = d(8,dmp8,rij)
      d10 = d(10,dmp10,rij)
      c6
     & =c6
     & +parab(11,itypea(ia),itypeb(ib))*(s3+s6)
     & +parab(14,itypea(ia),itypeb(ib))*(s1+s4)
     & +parab(17,itypea(ia),itypeb(ib))*(s2+s5)
     & +parab(20,itypea(ia),itypeb(ib))*(s3*s6)
     & +parab(23,itypea(ia),itypeb(ib))*(s1*s4)
     & +parab(26,itypea(ia),itypeb(ib))*(s2*s5)
      c8
     & =c8
     & +parab(12,itypea(ia),itypeb(ib))*(s3+s6)
     & +parab(15,itypea(ia),itypeb(ib))*(s1+s4)
     & +parab(18,itypea(ia),itypeb(ib))*(s2+s5)
     & +parab(21,itypea(ia),itypeb(ib))*(s3*s6)
     & +parab(24,itypea(ia),itypeb(ib))*(s1*s4)
     & +parab(27,itypea(ia),itypeb(ib))*(s2*s5)
      c10
     & =c10
     & +parab(13,itypea(ia),itypeb(ib))*(s3+s6)
     & +parab(16,itypea(ia),itypeb(ib))*(s1+s4)
     & +parab(19,itypea(ia),itypeb(ib))*(s2+s5)
     & +parab(22,itypea(ia),itypeb(ib))*(s3*s6)
     & +parab(25,itypea(ia),itypeb(ib))*(s1*s4)
     & +parab(28,itypea(ia),itypeb(ib))*(s2*s5)
c Update c6, c8, c10 with asymptotic terms asymmetric in s_A, s_B:
      c6as = 0.d0
      c8as = 0.d0
      c10as = 0.d0
      if(itypea(ia).ne.itypeb(ib)) then
      c6as
     & =c6as
     & +parab(29,itypea(ia),itypeb(ib))*(s3-s6)
     & +parab(32,itypea(ia),itypeb(ib))*(s1-s4)
     & +parab(35,itypea(ia),itypeb(ib))*(s2-s5)
      c8as
     & =c8as
     & +parab(30,itypea(ia),itypeb(ib))*(s3-s6)
     & +parab(33,itypea(ia),itypeb(ib))*(s1-s4)
     & +parab(36,itypea(ia),itypeb(ib))*(s2-s5)
      c10as
     & =c10as
     & +parab(31,itypea(ia),itypeb(ib))*(s3-s6)
     & +parab(34,itypea(ia),itypeb(ib))*(s1-s4)
     & +parab(37,itypea(ia),itypeb(ib))*(s2-s5)
       if(itypea(ia).gt.itypeb(ib)) then
        c6as = -c6as
        c8as = -c8as
        c10as = -c10as
       endif
      endif
      c6 = c6 + c6as
      c8 = c8 + c8as
      c10 = c10 + c10as
c
      if (beta.gt.0.d0) then 
      if(itypea(ia).ne.2.and.itypeb(ib).ne.2) then ! no H involved
       if(itypea(ia).eq.itypeb(ib)) then
        numtm(1) = 40
        numtm(2) =  0 !28
       else   ! test: make it the same as symmetric
        numtm(1) = 40
        numtm(2) = 28
       endif
      else
       if(itypea(ia).eq.itypeb(ib)) then
        numtm(1) = 40
        numtm(2) =  0 !28  ! no s-expansion if none of sites is H
       else
        numtm(1) = 40
        numtm(2) = 28
       endif
      endif
      numt = numtm(1) + numtm(2) + 1
c
      val0=a*dexp(-beta*rij)
      val1=val0*rij
      val2=val1*rij
      val3=val2*rij
       values(numt) = val0 + a1*val1 + a2*val2 + a3*val3
     & + d1*qa*qb/rij
     & - d6*c6/rij**6 - d8*c8/rij**8 - d10*c10/rij**10
c
c Construct planning matrix by filling values vector
c s-Symmetric terms
      if(itypea(ia).ne.2.and.itypeb(ib).ne.2) then ! no H involved
c delta a0,1,2,3 contributions to planning matrix
      values(1)  = (s1+s4)*val0
      values(2)  = (s1+s4)*val1
      values(3)  = (s1+s4)*val2
      values(4)  = (s1+s4)*val3

      values(5)  = (s2+s5)*val0
      values(6)  = (s2+s5)*val1
      values(7)  = (s2+s5)*val2
      values(8)  = (s2+s5)*val3

      values(9)  = (s3+s6)*val0
      values(10) = (s3+s6)*val1
      values(11) = (s3+s6)*val2
      values(12) = (s3+s6)*val3


      values(13) = (s1*s2+s4*s5)*val0
      values(14) = (s1*s2+s4*s5)*val1
      values(15) = (s1*s2+s4*s4)*val2
      values(16) = (s1*s2+s4*s4)*val3

      values(17) = (s2*s3+s5*s6)*val0
      values(18) = (s2*s3+s5*s6)*val1
      values(19) = (s2*s3+s5*s6)*val2
      values(20) = (s2*s3+s5*s6)*val3

      values(21) = (s1*s1+s4*s4)*val0
      values(22) = (s1*s1+s4*s4)*val1
      values(23) = (s1*s1+s4*s4)*val2
      values(24) = (s1*s1+s4*s4)*val3

      values(25) = (s2*s2+s5*s5)*val0
      values(26) = (s2*s2+s5*s5)*val1
      values(27) = (s2*s2+s5*s5)*val2
      values(28) = (s2*s2+s5*s5)*val3


      values(29) = (s1*s4)*val0
      values(30) = (s1*s4)*val1
      values(31) = (s1*s4)*val2
      values(32) = (s1*s4)*val3

      values(33) = (s2*s5)*val0
      values(34) = (s2*s5)*val1
      values(35) = (s2*s5)*val2
      values(36) = (s2*s5)*val3

      values(37) = (s3*s6)*val0
      values(38) = (s3*s6)*val1
      values(39) = (s3*s6)*val2
      values(40) = (s3*s6)*val3
c For uneqal types, add s-antisymmetric terms
      if(itypea(ia).lt.itypeb(ib)) then
      values(41) = (s1-s4)*val0
      values(42) = (s1-s4)*val1
      values(43) = (s1-s4)*val2
      values(44) = (s1-s4)*val3

      values(45) = (s2-s5)*val0
      values(46) = (s2-s5)*val1
      values(47) = (s2-s5)*val2
      values(48) = (s2-s5)*val3

      values(49) = (s3-s6)*val0
      values(50) = (s3-s6)*val1
      values(51) = (s3-s6)*val2
      values(52) = (s3-s6)*val3


      values(53) = (s1*s2-s4*s5)*val0
      values(54) = (s1*s2-s4*s5)*val1
      values(55) = (s1*s2-s4*s4)*val2
      values(56) = (s1*s2-s4*s4)*val3

      values(57) = (s2*s3-s5*s6)*val0
      values(58) = (s2*s3-s5*s6)*val1
      values(59) = (s2*s3-s5*s6)*val2
      values(60) = (s2*s3-s5*s6)*val3

      values(61) = (s1*s1-s4*s4)*val0
      values(62) = (s1*s1-s4*s4)*val1
      values(63) = (s1*s1-s4*s4)*val2
      values(64) = (s1*s1-s4*s4)*val3

      values(65) = (s2*s2-s5*s5)*val0
      values(66) = (s2*s2-s5*s5)*val1
      values(67) = (s2*s2-s5*s5)*val2
      values(68) = (s2*s2-s5*s5)*val3
      endif
      if(itypea(ia).gt.itypeb(ib)) then
      values(41) =-(s1-s4)*val0
      values(42) =-(s1-s4)*val1
      values(43) =-(s1-s4)*val2
      values(44) =-(s1-s4)*val3

      values(45) =-(s2-s5)*val0
      values(46) =-(s2-s5)*val1
      values(47) =-(s2-s5)*val2
      values(48) =-(s2-s5)*val3

      values(49) =-(s3-s6)*val0
      values(50) =-(s3-s6)*val1
      values(51) =-(s3-s6)*val2
      values(52) =-(s3-s6)*val3


      values(53) =-(s1*s2-s4*s5)*val0
      values(54) =-(s1*s2-s4*s5)*val1
      values(55) =-(s1*s2-s4*s4)*val2
      values(56) =-(s1*s2-s4*s4)*val3

      values(57) =-(s2*s3-s5*s6)*val0
      values(58) =-(s2*s3-s5*s6)*val1
      values(59) =-(s2*s3-s5*s6)*val2
      values(60) =-(s2*s3-s5*s6)*val3

      values(61) =-(s1*s1-s4*s4)*val0
      values(62) =-(s1*s1-s4*s4)*val1
      values(63) =-(s1*s1-s4*s4)*val2
      values(64) =-(s1*s1-s4*s4)*val3

      values(65) =-(s2*s2-s5*s5)*val0
      values(66) =-(s2*s2-s5*s5)*val1
      values(67) =-(s2*s2-s5*s5)*val2
      values(68) =-(s2*s2-s5*s5)*val3
      endif
      else    ! at least 1 H involved
c delta a0,1,2,3 contributions to planning matrix
      values(1)  = (s1+s4)*val0
      values(2)  = (s1+s4)*val1
      values(3)  = (s1+s4)*val2
      values(4)  = (s1+s4)*val3

      values(5)  = (s2+s5)*val0
      values(6)  = (s2+s5)*val1
      values(7)  = (s2+s5)*val2
      values(8)  = (s2+s5)*val3

      values(9)  = (s3+s6)*val0
      values(10) = (s3+s6)*val1
      values(11) = (s3+s6)*val2
      values(12) = (s3+s6)*val3


      values(13) = (s1*s2+s4*s5)*val0
      values(14) = (s1*s2+s4*s5)*val1
      values(15) = (s1*s2+s4*s4)*val2
      values(16) = (s1*s2+s4*s4)*val3

      values(17) = (s2*s3+s5*s6)*val0
      values(18) = (s2*s3+s5*s6)*val1
      values(19) = (s2*s3+s5*s6)*val2
      values(20) = (s2*s3+s5*s6)*val3

      values(21) = (s1*s1+s4*s4)*val0
      values(22) = (s1*s1+s4*s4)*val1
      values(23) = (s1*s1+s4*s4)*val2
      values(24) = (s1*s1+s4*s4)*val3

      values(25) = (s2*s2+s5*s5)*val0
      values(26) = (s2*s2+s5*s5)*val1
      values(27) = (s2*s2+s5*s5)*val2
      values(28) = (s2*s2+s5*s5)*val3


      values(29) = (s1*s4)*val0
      values(30) = (s1*s4)*val1
      values(31) = (s1*s4)*val2
      values(32) = (s1*s4)*val3

      values(33) = (s2*s5)*val0
      values(34) = (s2*s5)*val1
      values(35) = (s2*s5)*val2
      values(36) = (s2*s5)*val3

      values(37) = (s3*s6)*val0
      values(38) = (s3*s6)*val1
      values(39) = (s3*s6)*val2
      values(40) = (s3*s6)*val3
c For uneqal types, add s-antisymmetric terms
      if(itypea(ia).lt.itypeb(ib)) then
      values(41) = (s1-s4)*val0
      values(42) = (s1-s4)*val1
      values(43) = (s1-s4)*val2
      values(44) = (s1-s4)*val3

      values(45) = (s2-s5)*val0
      values(46) = (s2-s5)*val1
      values(47) = (s2-s5)*val2
      values(48) = (s2-s5)*val3

      values(49) = (s3-s6)*val0
      values(50) = (s3-s6)*val1
      values(51) = (s3-s6)*val2
      values(52) = (s3-s6)*val3


      values(53) = (s1*s2-s4*s5)*val0
      values(54) = (s1*s2-s4*s5)*val1
      values(55) = (s1*s2-s4*s4)*val2
      values(56) = (s1*s2-s4*s4)*val3

      values(57) = (s2*s3-s5*s6)*val0
      values(58) = (s2*s3-s5*s6)*val1
      values(59) = (s2*s3-s5*s6)*val2
      values(60) = (s2*s3-s5*s6)*val3

      values(61) = (s1*s1-s4*s4)*val0
      values(62) = (s1*s1-s4*s4)*val1
      values(63) = (s1*s1-s4*s4)*val2
      values(64) = (s1*s1-s4*s4)*val3

      values(65) = (s2*s2-s5*s5)*val0
      values(66) = (s2*s2-s5*s5)*val1
      values(67) = (s2*s2-s5*s5)*val2
      values(68) = (s2*s2-s5*s5)*val3
      endif
      if(itypea(ia).gt.itypeb(ib)) then
      values(41) =-(s1-s4)*val0
      values(42) =-(s1-s4)*val1
      values(43) =-(s1-s4)*val2
      values(44) =-(s1-s4)*val3

      values(45) =-(s2-s5)*val0
      values(46) =-(s2-s5)*val1
      values(47) =-(s2-s5)*val2
      values(48) =-(s2-s5)*val3

      values(49) =-(s3-s6)*val0
      values(50) =-(s3-s6)*val1
      values(51) =-(s3-s6)*val2
      values(52) =-(s3-s6)*val3


      values(53) =-(s1*s2-s4*s5)*val0
      values(54) =-(s1*s2-s4*s5)*val1
      values(55) =-(s1*s2-s4*s4)*val2
      values(56) =-(s1*s2-s4*s4)*val3

      values(57) =-(s2*s3-s5*s6)*val0
      values(58) =-(s2*s3-s5*s6)*val1
      values(59) =-(s2*s3-s5*s6)*val2
      values(60) =-(s2*s3-s5*s6)*val3

      values(61) =-(s1*s1-s4*s4)*val0
      values(62) =-(s1*s1-s4*s4)*val1
      values(63) =-(s1*s1-s4*s4)*val2
      values(64) =-(s1*s1-s4*s4)*val3

      values(65) =-(s2*s2-s5*s5)*val0
      values(66) =-(s2*s2-s5*s5)*val1
      values(67) =-(s2*s2-s5*s5)*val2
      values(68) =-(s2*s2-s5*s5)*val3
      endif
c
      endif   ! if H involved
 
      else
c
      numt=1
      numtm(1) = 0
      numtm(2) = 0
      values(numt) =
     & + d1*qa*qb/rij
     & - d6*c6/rij**6 - d8*c8/rij**8 - d10*c10/rij**10
c
      end if 
c
      return
      endif
      end
c
c============================================================================
c
c     calculate the damping factor (small R correct)
c
      function d(n,beta,r)
c
      implicit real*8 (a-h,o-z)
      br=beta*r
c The following line added by RB, Sept. 18 1997
      if(br.eq.0.d0) then
       d = 0.d0
       return
      endif
      sum=1.0d0
      term=1.0d0
      ncn=n
      do i=1,ncn
        term=term*br/i
        sum=sum+term
      enddo
      d=1.0d0 - dexp(-br)*sum
c     in case of d --> 0 use
c     d=1.0d0 - dexp(-br)*sum = sum_m=ncn+1^\infty br^m/m!
      if(dabs(d).lt.1.0d-8) then
        d=0.0d0
        do i=ncn+1,1000
          term=term*br/i
          d=d+term
          if(term/d .lt. 1.0d-8) go to 111
        enddo
        write(66,*) 'No convergence in d'
  111 continue
      d=d*dexp(-br)
      endif
      return
      end
c============================================================================
c
c The parameters of the fit are read in this subroutine.
c
      subroutine data1(param,parab,nsitea,nsiteb,c)
      implicit real*8 (a-h,o-z)
c
      parameter (naamax=14,nbbmax=14)
      parameter (nsitemax=8,ntypemax=6,ntypemax2=ntypemax*ntypemax)
      parameter (ntypemax_9=ntypemax*9,ntypemax2_9=ntypemax2*9)
      parameter (maxpar1=18,maxpar2=84)
c
      character*8 label(9)
      dimension param(maxpar1,ntypemax),parab(maxpar2,ntypemax,ntypemax)
      dimension c(1000)
      integer omp_get_num_threads,omp_get_thread_num

      a0=0.529177249d0

c      ns=omp_get_num_threads()
c      id=omp_get_thread_num()
      
c
c Zero out the nonlinear parameters matrices
c
      do i=1,ntypemax
       do k=1,maxpar1
        param(k,i) = 0.d0
       end do
      end do
      do i=1,ntypemax
       do j=1,ntypemax
        do k=1,maxpar2
         parab(k,i,j) = 0.d0
        end do
       end do
      end do
c
      nsitea=8
      nsiteb=8

c
c Read in the one-site nonlinear parameters and opt indicators..
c 
c nparm used only here
      read(55,*) nparm
      do i=1,nparm
       read(55,*) ityp, inumpar, val, indopt
       param(inumpar,ityp) = val
c Convert charges to the proper units 
c and DELTA charge!!! GM 3/29/05
c Switch off for the total charge density variant.... (everything
c in atomic units)
       if(inumpar.le.9) then
        param(inumpar,ityp) = 18.22262373d0*param(inumpar,ityp)
       endif
      end do
c
c Read in the two-site nonlinear parameters and opt indicators..
c
c nparab used only here
      read(55,*) nparab
      do i=1,nparab
       read(55,*) ityp1,ityp2,inumpar, val, indopt
       parab(inumpar,ityp1,ityp2) = val
       parab(inumpar,ityp2,ityp1) = val
      end do
c
c Read the type of the potential form to be used, inear coeff. reading
c indicator, and the optimization parameters....
c
c w wersji ostatecznej ustawic:
c ntpot=124161
c idonl=1 (w procedurze wolajacej)
c wartosci iopt, iweight, iasdone niewazne - mozna pominac
c oraz dodac do commonu dataoffit
      read(55,*) ntpot, idonl, iopt, iweight, iasdone
c      ntpot=124161
c the following parameters are not necessary for evaluation of the interaction
c energy value, but are left to use the standard input - can be removed
c in the final "paper" version.
      read(55,*) TOLF,TOLR,ANOISE,TIMLIM,LINP,LOUT,ALFST,IPR,SAFETL
      read(55,*) R_0,isyst,npowers       ! to be removed in the final version
      read(55,*) RCOND                   ! to be removed in the final version
c      ntpot=124161
c      iasdone=1

c
        read(55,*)numlin
        do i=1,numlin
         read(55,*) c(i)
        end do

       return
       end
c============================================================================
c
c A simple routine to calculate the dipole-dipole induction,
c now in the flexible case
c units: charges and polarizability calculated in au; energy returned in kcal/mol
c
      subroutine dipind(iaa,ibb,R,energy,
     .    sa,sb,param,parab,siteat,sitebt,nsitea,nsiteb,itypea,itypeb)
      implicit real*8 (a-h,o-z)
cccccc Parameters taken from potparts
      parameter (naamax=14,nbbmax=14)
      parameter (nsitemax=8,ntypemax=6,ntypemax2=ntypemax*ntypemax)
      parameter (ntypemax_9=ntypemax*9,ntypemax2_9=ntypemax2*9)
      parameter (maxpar1=18,maxpar2=84)
c
      dimension param(maxpar1,ntypemax),parab(maxpar2,ntypemax,ntypemax)
      dimension sa(3,-1:naamax),sb(3,-1:nbbmax)
      dimension siteat(3,nsitemax,-1:naamax),
     .          sitebt(3,nsitemax,-1:nbbmax)
      dimension itypea(nsitemax),itypeb(nsitemax)
cccccc end of parameters taken from potparts
c      dimension dma_unrot(3), dmb_unrot(3), dma(3), dmb(3), u(3)
      dimension dma(3), dmb(3), u(3)
      a0=0.529177249d0
      har2kcal=627.510d0
c
c Calculate static dipole moments of A and B in these orientations and
c Let's place the polarizability on Oxygen, just
c for simplicity
c
c monomer A
c
      s1=sa(1,iaa)
      s2=sa(2,iaa)
      s3=sa(3,iaa)
c      write(31,*) iaa,s1,s2,s3
      signa=1.0d0
      do i=1,3
c       dma_unrot(i)=0.0d0  !for checking
       dma(i)=0.0d0
      end do
      polisa=0.0d0
c
      do 11 ia=1,nsitea
      if (ia .eq. 3) signa=-1.0d0
      s3=signa*s3
c Make basis functions even in s6 if site is NOT H
c     if(itypea(ia).ne.2) s3=dabs(s3)
      qa
     & =param(1,itypea(ia))         !already in units of [kcal/mol & Ang]/Ang^n ; multiply by Ang^n to get bohr^3
     & +param(2,itypea(ia))*s1
     & +param(3,itypea(ia))*s2
     & +param(4,itypea(ia))*s3
     & +param(5,itypea(ia))*s1*s2
     & +param(6,itypea(ia))*s2*s3
     & +param(7,itypea(ia))*s1*s1
     & +param(8,itypea(ia))*s2*s2
     & +param(9,itypea(ia))*s3*s3
      qa=qa/18.22262373d0           !convert back to atomic units
c      write(31,101)
c     &iaa,ia,param(1,itypea(ia))/18.22262373d0,qa
      do i=1,3
c       dma_unrot(i)=dma_unrot(i)+qa*(sitea(i,ia,iaa))/a0
       dma(i)=dma(i)+qa*(siteat(i,ia,iaa))/a0
      end do
      if (ia.eq.1) then !Oxygen only
      polisa
     & =param(10,itypea(ia))
     & +param(11,itypea(ia))*s1
     & +param(12,itypea(ia))*s2
     & +param(13,itypea(ia))*s3
     & +param(14,itypea(ia))*s1*s2
     & +param(15,itypea(ia))*s2*s3
     & +param(16,itypea(ia))*s1*s1
     & +param(17,itypea(ia))*s2*s2
     & +param(18,itypea(ia))*s3*s3  !already in units of bohr^3/Ang^n ; multiply by Ang^n to get bohr^3
      end if
11    continue
c      write(31,201) iaa, (dma_unrot(i),i=1,3)
c      write(31,301) iaa, (dma(i),i=1,3)
c      write(31,*) 'iaa, polisa:', iaa, polisa
c
101   format('iaa,ia,qa_0,qa:',
     &2(1x,i6),2(1x,f12.7))
201   format('iaa,dma_unrot(i):',
     &1(1x,i6),3(1x,f12.7))
301   format('iaa,dma(i):',
     &1(1x,i6),3(1x,f12.7))
c
c monomer B
c
      s4=sb(1,ibb)
      s5=sb(2,ibb)
      s6=sb(3,ibb)
c      write(32,*) ibb,s4,s5,s6
      signb=1.0d0
      do i=1,3
c       dmb_unrot(i)=0.0d0  !for checking
       dmb(i)=0.0d0
      end do
      polisb=0.0d0
c
      do 12 ib=1,nsiteb
      if (ib .eq. 3) signb=-1.0d0
      s6=signb*s6
c Make basis functions even in s6 if site is NOT H
c     if(itypeb(ib).ne.2) s6=dabs(s6)
      qb
     & =param(1,itypeb(ib))         !already in units of [kcal/mol & Ang]/Ang^n ; multiply by Ang^n to get bohr^3
     & +param(2,itypeb(ib))*s4
     & +param(3,itypeb(ib))*s5
     & +param(4,itypeb(ib))*s6
     & +param(5,itypeb(ib))*s4*s5
     & +param(6,itypeb(ib))*s5*s6
     & +param(7,itypeb(ib))*s4*s4
     & +param(8,itypeb(ib))*s5*s5
     & +param(9,itypeb(ib))*s6*s6
      qb=qb/18.22262373d0           !convert back to atomic units
c      write(32,102)
c     &ibb,ib,param(1,itypeb(ib))/18.22262373d0,qb
      Rtemp=0.0
      do i=1,3
       if(i.eq.3)Rtemp=R
c       dmb_unrot(i)=dmb_unrot(i)+qb*(siteb(i,ib,ibb))/a0
       dmb(i)=dmb(i)+qb*(sitebt(i,ib,ibb)-Rtemp)/a0
      end do
      if (ib.eq.1) then !Oxygen only
      polisb
     & =param(10,itypeb(ib))
     & +param(11,itypeb(ib))*s4
     & +param(12,itypeb(ib))*s5
     & +param(13,itypeb(ib))*s6
     & +param(14,itypeb(ib))*s4*s5
     & +param(15,itypeb(ib))*s5*s6
     & +param(16,itypeb(ib))*s4*s4
     & +param(17,itypeb(ib))*s5*s5
     & +param(18,itypeb(ib))*s6*s6  !already in units of bohr^3/Ang^n ; multiply by Ang^n to get bohr^3
      end if
12    continue
c      write(32,202) ibb, (dmb_unrot(i),i=1,3)
c      write(32,302) ibb, (dmb(i),i=1,3)
c      write(32,*) 'ibb, polisb:', ibb, polisb
c
102   format('ibb,ib,qb_0,qb:',2(1x,i6),2(1x,f12.7))
202   format('ibb,dmb_unrot(i):',1(1x,i6),3(1x,f12.7))
302   format('ibb,dmb(i):',1(1x,i6),3(1x,f12.7))
c
c Calculate the distance between the O's (pol centers)
c
      dlen = 0.d0
      do i=1,3
       pom = sitebt(i,1,ibb) - siteat(i,1,iaa)
       dlen = dlen + pom*pom
      end do
      dlen = dsqrt(dlen)
c
c Compute the damping factor for induction
c
      dmpind = d(6,parab(10,1,1),dlen)
c
      dlen = dlen**(-3.d0)   ! inverse cube of the distance
c
c u(*): field of A on B and the induction energy
c
      call TTTprod(siteat(1,1,iaa),sitebt(1,1,ibb),dma,dlen,u)
      energy_a_on_b = polisa*scalp(u,u)
c
c u(*): field of B on A and the induction energy
c
      call TTTprod(siteat(1,1,iaa),sitebt(1,1,ibb),dmb,dlen,u)
      energy_b_on_a = polisb*scalp(u,u)
c
      energy = energy_a_on_b + energy_b_on_a
      energy = -0.5d0*(a0**6)*har2kcal*energy*dmpind
c
      return
      end
c=============================================================================
c
c Calculate the vector V resulting from the action
c of the dipole propagator tensor T_ij on U
c Ri, Rj - position vectors of molecules i and j, respectively
c rij = |Ri - Rj|^(-3)
c
        subroutine TTTprod(Ri,Rj,u,rij,v)
        implicit real*8 (a-h,o-z)
        dimension Ri(3),Rj(3),u(3),v(3)
c
c      goto 777
        ddd = rij**(0.66666666666666666d0)
        scal = 0.d0
        do i=1,3
         v(i) = Ri(i) - Rj(i)
         scal = scal + v(i)*u(i)
        end do
        do i=1,3
         v(i) = (3.d0*v(i)*scal*ddd - u(i))*rij
        end do
c
 777   continue
        return
        end
c=============================================================================
c
        function scalp(a,b)
        implicit real*8 (a-h,o-z)
        dimension a(3),b(3)
        scalp = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
        return
        end 
c=============================================================================
c In this procedure positions of sites are established on the base of 
c the cartesian coordinates of atoms. 
c Written on the base of data from the subroutine getgeo
c
c Also the sa and itypea matrices are filled in this routine

      subroutine set_sites(carta,sitea,sa,itypea,nsitea)
      implicit real*8 (a-h,o-z)

      parameter (lay_claude=0)
      parameter (nsitemax=8)
c
      dimension sitea(3,nsitemax),sa(3)
      dimension itypea(nsitemax)
      dimension carta(3,3),cartb(3,3)
      dimension v1(3),vn1(3),v2(3),vn2(3),v(3),vn(3),vv(3),vsm(3)
      dimension vb(3),vp(3),vd1a(3),vd1b(3),vd2a(3),vd2b(3)

      a0=0.529177249d0
      zero=0.d0
      r0_ang=0.9716257d0
      theta0_deg=104.69d0

      iaa=-1
      ibb=-1
c           
c Info regarding site positions for an undistorted monomer with approximation 
c that amO=16 amu, amH=1 amu 
c  sig   is the usual O site-COM site separation
c  sig2  is the O sep along the z coordinate of Bunny 1
c  sig3  is the y coordinate of Bunny 1
c  sig4  is the O sep along the z coordinate of Bunny 2
c  sig5  is the y coordinate of Bunny 2
      sig=  1.24631924913843d-01                 ! in bohr
      sig2= 0.371792435d0                        !
      sig3= 0.2067213d0                          !
      sig4= 0.125368076d0                        !
      sig5= 0.2d0                                !

      shift=9.01563628739252d-4  !accounts for shift to precise masses (in bohr)

      pi=dacos(-1.d0)
      rad2d=180.d0/pi

c Distances in carta are given in bohr, so we will use bohr as the default unit
c They should be changed to Angstroms at the end
c--- O1
      sitea(1,1)=carta(1,1)
      sitea(2,1)=carta(1,2)
      sitea(3,1)=carta(1,3)
c--- H1
      sitea(1,2)=carta(2,1)
      sitea(2,2)=carta(2,2)
      sitea(3,2)=carta(2,3)
c--- H2
      sitea(1,3)=carta(3,1)
      sitea(2,3)=carta(3,2)
      sitea(3,3)=carta(3,3)

c vn1 is a normalized vector pointing from O to H1
      v1(1)=carta(2,1)-carta(1,1)
      v1(2)=carta(2,2)-carta(1,2)
      v1(3)=carta(2,3)-carta(1,3)
      xnv1=dsqrt(v1(1)**2+v1(2)**2+v1(3)**2)
      vn1(1)=v1(1)/xnv1
      vn1(2)=v1(2)/xnv1
      vn1(3)=v1(3)/xnv1
      xnv1_ang=xnv1*a0
c      write(44,*) "xnv1_ang=",xnv1_ang," AA"

c vn2 is a normalized vector pointing from O to H2
      v2(1)=carta(3,1)-carta(1,1)
      v2(2)=carta(3,2)-carta(1,2)
      v2(3)=carta(3,3)-carta(1,3)
      xnv2=dsqrt(v2(1)**2+v2(2)**2+v2(3)**2)
      vn2(1)=v2(1)/xnv2
      vn2(2)=v2(2)/xnv2
      vn2(3)=v2(3)/xnv2
      xnv2_ang=xnv2*a0
c      write(44,*) "xnv2_ang=",xnv2_ang," AA"

c components of the bisector (vb) vector (normalized):
      v(1)=vn1(1)+vn2(1)
      v(2)=vn1(2)+vn2(2)
      v(3)=vn1(3)+vn2(3)
      xnv=dsqrt(v(1)**2+v(2)**2+v(3)**2)
      vb(1)=v(1)/xnv
      vb(2)=v(2)/xnv
      vb(3)=v(3)/xnv

c a normalized vector, perpendicular (vp) to the surface defined by 
c the v1 and v2 vectors (or vn1 and vn2)
      v(1)=v1(2)*v2(3)-v1(3)*v2(2)
      v(2)=v1(3)*v2(1)-v1(1)*v2(3)
      v(3)=v1(1)*v2(2)-v1(2)*v2(1)
      xn=dsqrt(v(1)**2+v(2)**2+v(3)**2)
      vp(1)=v(1)/xn
      vp(2)=v(2)/xn
      vp(3)=v(3)/xn

c set bunny_ratio_A
      r0=r0_ang/a0
      theta0=theta0_deg/rad2d
      cta=dcos(0.5d0*theta0)
      prodv1vb=v1(1)*vb(1)+v1(2)*vb(2)+v1(3)*vb(3)
      prodv2vb=v2(1)*vb(1)+v2(2)*vb(2)+v2(3)*vb(3)
      IF (LAY_CLAUDE .EQ. 0) THEN
        bunny_ratio_A=(0.5d0*(prodv1vb+prodv2vb))/(r0*cta)
c        write(44,*) "bunny_ratio=",bunny_ratio_A
      END IF
      IF (LAY_CLAUDE .EQ. 1) THEN
        bunny_ratio_A=1.0d0
      END IF

c--- Bunny1 1 -- charged one
      vd1a(1)=sig3*vp(1)+sig2*vb(1)*bunny_ratio_A
      vd1a(2)=sig3*vp(2)+sig2*vb(2)*bunny_ratio_A
      vd1a(3)=sig3*vp(3)+sig2*vb(3)*bunny_ratio_A
      sitea(1,4)=carta(1,1)+vd1a(1)
      sitea(2,4)=carta(1,2)+vd1a(2)
      sitea(3,4)=carta(1,3)+vd1a(3)
c--- Bunny1 2 -- charged one
      vd1b(1)=-sig3*vp(1)+sig2*vb(1)*bunny_ratio_A
      vd1b(2)=-sig3*vp(2)+sig2*vb(2)*bunny_ratio_A
      vd1b(3)=-sig3*vp(3)+sig2*vb(3)*bunny_ratio_A
      sitea(1,5)=carta(1,1)+vd1b(1)
      sitea(2,5)=carta(1,2)+vd1b(2)
      sitea(3,5)=carta(1,3)+vd1b(3)
c--- Bunny2 1 -- exp-type one
      vd2a(1)=sig5*vp(1)-sig4*vb(1)*bunny_ratio_A
      vd2a(2)=sig5*vp(2)-sig4*vb(2)*bunny_ratio_A
      vd2a(3)=sig5*vp(3)-sig4*vb(3)*bunny_ratio_A
      sitea(1,6)=carta(1,1)+vd2a(1)
      sitea(2,6)=carta(1,2)+vd2a(2)
      sitea(3,6)=carta(1,3)+vd2a(3)
c--- Bunny2 2 -- exp-type one
      vd2b(1)=-sig5*vp(1)-sig4*vb(1)*bunny_ratio_A
      vd2b(2)=-sig5*vp(2)-sig4*vb(2)*bunny_ratio_A
      vd2b(3)=-sig5*vp(3)-sig4*vb(3)*bunny_ratio_A
      sitea(1,7)=carta(1,1)+vd2b(1)
      sitea(2,7)=carta(1,2)+vd2b(2)
      sitea(3,7)=carta(1,3)+vd2b(3)

c calculate COM position with Robert''s masses
      xm16=15.994915d0
      xm1=1.007825d0
      sm=xm16+2.0d0*xm1
      vsm(1)=(xm16*carta(1,1)+xm1*carta(2,1)+xm1*carta(3,1))/sm
      vsm(2)=(xm16*carta(1,2)+xm1*carta(2,2)+xm1*carta(3,2))/sm
      vsm(3)=(xm16*carta(1,3)+xm1*carta(2,3)+xm1*carta(3,3))/sm
c and then shift the position of the com-site to the position of COM calculated
c with the approximated masses.
c--- com
      sitea(1,8)= vsm(1) - shift*vb(1)
      sitea(2,8)= vsm(2) - shift*vb(2)
      sitea(3,8)= vsm(3) - shift*vb(3)
c      sss=sqrt(vsm(1)**2+vsm(2)**2+vsm(3)**2)
c      write(44,*) "|vsm|=",sss

c transform to Angstroms...
      do ia=1,nsitea
       do i=1,3
        sitea(i,ia) = sitea(i,ia)*a0
       end do
      end do

c As in the subroutine "monomer":
c distorts monomer A (or B) and writes the distorted symmetry coordinates sa
c (passed back as sa or sb) 

      sprod=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
      ccos=sprod/(xnv1*xnv2)
      theta1=dacos(ccos)
      theta1_deg=theta1*rad2d
c      write(44,*) "theta1_deg=",theta1_deg
      dsqrt2=dsqrt(2.0d0)
      sa(1)=( (xnv1_ang-r0_ang) + (xnv2_ang-r0_ang) )/dsqrt2
      sa(2)=dsqrt(xnv1_ang*xnv2_ang)*(theta1_deg-theta0_deg)/rad2d
      sa(3)=( (xnv1_ang-r0_ang) - (xnv2_ang-r0_ang) )/dsqrt2
     
      itypea(1)=1
      itypea(2)=2
      itypea(3)=2
      itypea(4)=3
      itypea(5)=3
      itypea(6)=4
      itypea(7)=4
      itypea(8)=5

      return
      end
c==============================================================================
