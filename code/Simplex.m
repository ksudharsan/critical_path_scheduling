C=[0;0.50;-0.40;0.60;-0.40;-0.10;-0.50;-0.20;9.90];
C=[C;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
A=zeros(16,25);
A(1,2)=-1;
A(1,10)=1;
A(2,3)=-1;A(2,11)=1;
A(3,2)=1;A(3,3)=-1;A(3,12)=-5;
A(4,2)=1;A(4,4)=-1;A(4,13)=1;
A(5,2)=1;A(5,5)=-1;A(5,14)=1;
A(6,2)=1;A(6,6)=-1;A(6,15)=1;
A(7,3)=1;A(7,4)=-1;A(7,16)=1;
A(8,3)=1;A(8,8)=-1;A(8,17)=1;
A(9,4)=1;A(9,5)=-1;A(9,18)=1;
A(10,4)=1;A(10,8)=-1;A(10,19)=1;
A(11,5)=1;A(11,6)=-1;A(11,20)=1;
A(12,5)=1;A(12,7)=-1;A(12,21)=1;
A(13,5)=1;A(13,8)=-1;A(13,22)=1;
A(14,6)=1;A(14,7)=-1;A(14,23)=1;
A(15,7)=1;A(15,8)=-1;A(15,24)=1;
A(16,8)=1;A(16,9)=-1;A(16,25)=1;
B=[-4;-9;-5;-4;-2;-8;-1;-4;-7;-11;-6;-12;-7;-3;-6;-5];
#[xn,sn,ln,mn,kn,nxs,rp,rd] = pclp(A,B,C);
[xn,zmin,p,pp]=glpk(C,A,B);
disp('optimal solution')
xn
disp('optimal value')
zmin = C'*xn