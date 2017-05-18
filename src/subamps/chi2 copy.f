ccc   gg --> chi_0 subprocess amplitude
      subroutine chi2(p,mqq,mxx,q1,q2,echi2,out)
      implicit double precision (a-z)
      double precision q1(2),q2(2)
      double precision qt1(4),qt2(4)
      complex*16 out,echi2(5,4,4),cpp
      integer p,i,j
      double precision pcm(4),pboo(4),plb(4)
      double precision q1b(4),q2b(4)

      include 'pi.f'
      include 'zi.f'
      include 'vars.f'
      include 'mom.f'
      include 'quarkonia.f'

      qt1(4)=0d0
      qt2(4)=0d0
      qt1(3)=0d0
      qt2(3)=0d0

      do i=1,4
         pboo(i)=-q(i,5)
      enddo
         pboo(4)=q(4,5)

c$$$      call boost(mx,pboo,q1,plb)
c$$$      do i=1,2
c$$$         qt1(i)=plb(i)
c$$$      enddo
c$$$
c$$$      call boost(mx,pboo,q2,plb)
c$$$      do i=1,2
c$$$         qt2(i)=plb(i)
c$$$      enddo

      do i=1,4
         pcm(i)=q(i,1)
      enddo
      call boost(mx,pboo,pcm,plb)
      do i=1,4
         q1b(i)=plb(i)
      enddo

      do i=1,4
         pcm(i)=q(i,2)
      enddo
      call boost(mx,pboo,pcm,plb)
      do i=1,4
         q2b(i)=plb(i)
      enddo
      
      do i=1,2
         qt1(i)=q1(i)
         qt2(i)=q2(i)
      enddo

      call boost(mx,pboo,qt1,plb)
      do i=1,4
         qt1(i)=plb(i)
      enddo

      call boost(mx,pboo,qt2,plb)
      do i=1,4
         qt2(i)=plb(i)
      enddo
  
      qt1sq=-(qt1(4)**2-qt1(3)**2-qt1(2)**2-qt1(1)**2)
      qt2sq=-(qt2(4)**2-qt2(3)**2-qt2(2)**2-qt2(1)**2)

c      qt1sq=qt1(1)**2+qt1(2)**2
c      qt2sq=qt2(1)**2+qt2(2)**2
      q1q2=(mx**2+qt1sq+qt2sq)/2d0

      cchi=dsqrt(pi*mx**3*gamchi0/3d0)

c$$$      cpp=s*(qt1(1)*qt2(1)*echi2(p,1,1)+qt1(1)*qt2(2)*echi2(p,1,2)
c$$$     &     +qt1(2)*qt2(1)*echi2(p,2,1)+qt1(2)*qt2(2)*echi2(p,2,2))
c$$$      cpp=cpp-2d0*(qt1(1)*qt2(1)+qt1(2)*qt2(2))*
c$$$     &     (q(3,1)*q(3,2)*echi2(p,3,3)+q(4,1)*q(4,2)*echi2(p,4,4)
c$$$     &     -q(3,1)*q(4,2)*echi2(p,3,4)-q(4,1)*q(3,2)*echi2(p,4,3))

  
            
      cpp=s*(qt1(1)*qt2(1)*echi2(p,1,1)+qt1(1)*qt2(2)*echi2(p,1,2)
     &     +qt1(2)*qt2(1)*echi2(p,2,1)+qt1(2)*qt2(2)*echi2(p,2,2)
     &     +qt1(3)*qt2(1)*echi2(p,3,1)+qt1(3)*qt2(2)*echi2(p,3,2)
     &     +qt1(3)*qt2(3)*echi2(p,3,3)-qt1(3)*qt2(4)*echi2(p,3,4)
     &     -qt1(4)*qt2(1)*echi2(p,4,1)-qt1(4)*qt2(2)*echi2(p,4,2)
     &     -qt1(4)*qt2(3)*echi2(p,4,3)+qt1(4)*qt2(4)*echi2(p,4,4)
     &     +qt1(1)*qt2(3)*echi2(p,1,3)-qt1(1)*qt2(4)*echi2(p,1,4)
     &     +qt1(2)*qt2(3)*echi2(p,2,3)-qt1(2)*qt2(4)*echi2(p,2,4))
c$$$      cpp=cpp-2d0*(qt1(1)*qt2(1)+qt1(2)*qt2(2)+qt1(3)*qt2(3)
c$$$     &     -qt1(4)*qt2(4))*
c$$$     &     (q1b(3)*q2b(3)*echi2(p,3,3)+q1b(4)*q2b(4)*echi2(p,4,4)
c$$$     &     -q1b(3)*q2b(4)*echi2(p,3,4)-q1b(4)*q2b(3)*echi2(p,4,3)
c$$$     &     +q1b(1)*q2b(1)*echi2(p,1,1)+q1b(1)*q2b(2)*echi2(p,1,2)
c$$$     &     +q1b(1)*q2b(3)*echi2(p,1,3)-q1b(1)*q2b(4)*echi2(p,1,4)
c$$$     &     +q1b(2)*q2b(1)*echi2(p,2,1)+q1b(2)*q2b(2)*echi2(p,2,2)
c$$$     &     +q1b(2)*q2b(3)*echi2(p,2,3)-q1b(2)*q2b(4)*echi2(p,2,4))

      do i=1,3
         qt1(i)=-qt1(i)
         qt2(i)=-qt2(i)
         q1b(i)=-q1b(i)
         q2b(i)=-q2b(i)
      enddo

      print*,cpp

      cpp=(0d0,0d0)
      
      do i=1,4
         do j=1,4
c            cpp=cpp+s*qt1(i)*qt2(j)*echi2(p,i,j)-2d0*(qt1(1)*qt2(1)+
c     &           qt1(2)*qt2(2)+qt1(3)*qt2(3)-qt1(4)*qt2(4))*
c     &           q1b(i)*q2b(j)*echi2(p,i,j)

            cpp=cpp+s*qt1(i)*qt2(j)*echi2(p,i,j)

            
         enddo
      enddo

      print*,cpp
      print*,''
      
      
c$$$      cpp=s*(qt1(1)*qt2(1)*echi2(p,1,1)+qt1(1)*qt2(2)*echi2(p,1,2)
c$$$     &     +qt1(2)*qt2(1)*echi2(p,2,1)+qt1(2)*qt2(2)*echi2(p,2,2))
c$$$      cpp=cpp-2d0*(qt1(1)*qt2(1)+qt1(2)*qt2(2))*
c$$$     &     (q1b(3)*q2b(3)*echi2(p,3,3)+q1b(4)*q2b(4)*echi2(p,4,4)
c$$$     &     -q1b(3)*q2b(4)*echi2(p,3,4)-q1b(4)*q2b(3)*echi2(p,4,3))
 
      cpp=cpp*cchi*dsqrt(2d0)*mx/s
      cpp=cpp/(2d0*mqq*mx+qt1sq+qt2sq)**2*4d0
      cpp=cpp*dsqrt(mx/mxx)
      cpp=cpp*mx**2/2d0

      out=cpp

      return
      end
