ccc   randomizes order of VEGAS unweighted events and
ccc   prints nev events to record
      subroutine unwprint
      implicit double precision(a-y)
      integer i,j,k,l,m
      integer evfill(2000000)

      include 'pdg.f'
      include 'unweighted.f'
      include 'mom.f'
      include 'proc.f'
      include 'decay.f'
      include 'record.f'
      include 'hepevt.f'
      include 'leshouches.f'
      include 'inparticle.f'
      include 'ewpars.f'
      include 'mp.f'
      include 'vars.f'

      do i=1,evnum
         evfill(i)=1
      enddo

      range=dble(evnum)

      do i=1,nev

 555     call r2455(r)

         j=nint(r*range)
         if(dble(j).lt.r*range)j=j+1
         if(evfill(j).eq.0)goto 555
         evfill(j)=0

c         if(i.eq.20)then
c            print*,j
c         endif

ccccccccccccccccccccccccccccccccccccccccccccccc
ccccc Les Houches
ccccccccccccccccccccccccccccccccccccccccccccccc

         if(erec.eq.'lhe')then
            
            do k=1,nup+2
               idup(k)=pdgid(k)
            enddo
            
           do k=3,nup+2
               do l=1,4
                  pup(l,k)=evrec(j,k,l)
               enddo
               pup(5,k)=dsqrt(dabs(pup(4,k)**2-pup(3,k)**2
     &              -pup(2,k)**2-pup(1,k)**2))
            enddo

            if(inparticle.eq.'el')then
               pup(5,3)=me
               pup(5,4)=me
            elseif(inparticle.eq.'prot')then
               pup(5,3)=mp
               pup(5,4)=mp
            endif

            if(i.eq.1)then
               write(45,*)'<LesHouchesEvents version="1.0">'
               write(45,*)'<header>'
               call headerlhe
               write(45,*)'</header>'
               write(45,*)'<init>'
               write(45,302)idup(1),idup(2),pup(4,1),pup(4,2),pdfgup(1)
     &              ,pdfgup(2),pdfsup(1),pdfsup(2),idwtup,nprup
               write(45,312)xsecup(1),xerrup(1),xmaxup(1),1
               write(45,*)'</init>'
            endif
c
 302        format(5x,i4,2x,i4,2x,E16.9,2x,E16.9,2x,i5,2x,i5,2x,i6,2x,i6
     &           ,2x,i1,2x,i1)
            
 312        format(E16.9,2x,E16.9,2x,E16.9,2x,i1)
           

            scalup=mx
            aqcdup=alphas(mx**2)
            
            write(45,*)'<event>'
            write(45,304)nup,idprup,xwgtup,scalup,aqedup,aqcdup
            do m=3,nup+2
               write(45,303)idup(m),istup(m),mothup(1,m),
     &              mothup(2,m),icolup(1,m),icolup(2,m),pup(1,m)
     &              ,pup(2,m),pup(3,m),pup(4,m),pup(5,m),vtimup(m)
     &              ,spinup(m)
            enddo
            write(45,*)'</event>'
            
c            write(45,*)''

 304        format(i2,4x,i1,3x,F2.0,3x,E16.9,3x,E16.9,3x,E16.9)
            
         endif
         
ccccccccccccccccccccccccccccccccccccccccccccccc
cccc  HEPEVT
ccccccccccccccccccccccccccccccccccccccccccccccc
         
         if(erec.eq.'hepevt')then
            
            nevhep=nev
            
            do k=1,nhep
               idhep(k)=pdgid(k)
            enddo
            
           do k=3,nhep
               do l=1,4
                  phep(l,k)=evrec(j,k,l)
               enddo
               phep(5,k)=dsqrt(dabs(phep(4,k)**2-phep(3,k)**2
     &              -phep(2,k)**2-phep(1,k)**2))
            enddo

            if(inparticle.eq.'el')then
               phep(5,3)=me
               phep(5,4)=me
            elseif(inparticle.eq.'prot')then
               phep(5,3)=mp
               phep(5,4)=mp
            endif

            do k=1,2
               do m=5,nhep
                  jmohep(k,m)=mothup(k,m)
               enddo
            enddo
               
            write(45,*)i

            
            do m=1,nhep
               write(45,300)m,idhep(m),isthep(m),jmohep(1,m),
     &              jmohep(2,m),jdahep(1,m),jdahep(2,m),
     &              phep(1,m),phep(2,m),phep(3,m),phep(4,m)
     &              ,phep(5,m),vhep(1,m),vhep(2,m),vhep(3,m),vhep(4,m)
            enddo
            
            write(45,*)''
            
         endif
               
      enddo

      if(erec.eq.'lhe')then
         write(45,*)'</LesHouchesEvents>'
      endif
         

 300  format(i4,1x,i8,1x,i4,1x,i4,1x,i4,1x,i4,1x,i4,1x,E16.9,1x,E16.9
     &,1x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9,1x
     &,E16.9)

 301  format(i5,1x,i4,1x,i8,1x,i4,1x,i4,1x,i4,1x,i4,1x,E16.9,1x,
     &E16.9,1x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9)

 303  format(7x,i9,1x,i8,1x,i4,1x,i4,1x,i4,1x,i4,1x,E16.9,1x,
     &E16.9,1x,E16.9,1x,E16.9,1x,E16.9,1x,F2.0,1x,F2.0)

      return
      end
