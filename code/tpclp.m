%    TEST FUNCTION FOR PCLP MACRO : FULL AND SPARSE MATRICES
% 
%  Scilab syntax:   exec tpclp.sci

clear 

% getf('pclp.sci');
% getf('concats.sci');
n=500;
%  p=  round(n/2);
p=5;

e = 0.01

%  Dense matrices
  A = rand(p,n);

  x = rand(n,1);
  s = rand(n,1);
  xi=x;
  si=s;

  b = A*x;
  c = s ;
  l = zeros(p,1); 

   minf = 1e-10;
   km=300;

  x=ones(n,1);
  s=ones(n,1);

  m=x'*s/n;
  verbose=3;
  redlin = 1;

disp('STARTING PCLP: TESTING THE FULL MATRICES CASE');


[xn,sn,ln,mn,kn,nxs,r1,r2] = pclp(A,b,c,x,s,km,minf,e,redlin,verbose);

%  Check optimality:

minxn = min(xn)

minsn = min(sn)

xnsnoverm = xn'*sn / n

primalfeas = norm(A*xn-b)

dualfeas = norm(c + A'*ln  - sn)
% 

disp('STARTING PCLP: TESTING THE SPARSE MATRICES CASE');

A =sparse(A); 

[xn,sn,ln,mn,kn,nxs,r1,r2] = pclp(A,b,c,x,s,km,minf,e,redlin,verbose);

%  Check optimality:

minxn = min(xn)

minsn = min(sn)

xnsnoverm = xn'*sn / n

primalfeas = norm(A*xn-b)

dualfeas = norm(c + A'*ln  - sn)
% 
