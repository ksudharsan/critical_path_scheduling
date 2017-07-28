%    TEST FUNCTION FOR PCLPM MACRO : FULL MATRICES AND SPARSE MATRICES
%               VERSION WITH MULTIPLE CORRECTOR STEPS
% 
%  Scilab syntax:   exec tpclpm.sci

clear 

% getf('pclpm.sci');
% getf('chsolvem.sci');
% getf('concats.sci');

n=50;
%  p=  round(n/2);
p=3;

e = 0.01

spar=-1
if spar > 0 % then
%  Sparse matrices
  den=0.01; A = sprand(p,n,den);
  a=rand(p,1);
  for i=1:p, A(i,i) = a(i); end; 
else
%  Dense matrices
  A = rand(p,n);
end;

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
x = normest(b) * x;
  s=ones(n,1);

  m=x'*s/n;
  redlin = 1;
  verbose= -1;

disp('STARTING PCLPM');

for redlin=[1,-1],
  for ifull=[-1,1], 
    AA = A;
    if ifull > 0 % then 
      AA = full(A);
    end;
%   
%  With factorization
factor = 1;
[xn,sn,ln,mn,kiter,kfact,klin]=pclpm(AA,b,c,x,s,km,minf,e,redlin,factor,verbose);

%  Check optimality:
minxn = min(xn)
minsn = min(sn)
xnsnoverm = xn'*sn / n
primalfeas = normest(A*xn-b)
%  normest(c + A'*ln + H*xn - sn)
dualfeas = normest(c + A'*ln  - sn)

%  Without factorization
factor = -1;

[xn,sn,ln,mn,kiterf,kfactf,klinf]=pclpm(AA,b,c,x,s,km,minf,e,redlin,factor,verbose);

%  Check optimality:
minxn = min(xn)
minsn = min(sn)
xnsnoverm = xn'*sn / n
primalfeas = normest(A*xn-b)
%  normest(c + A'*ln + H*xn - sn)
dualfeas = normest(c + A'*ln  - sn)

disp('Number of major iterations, factorizations, linear systems solved');

dispac = 'With refactorizations:   ';
dispac = concats(dispac, concats(num2str(kiter),'  '));
dispac = concats(dispac, concats(num2str(kfact),'  '));
dispac = concats(dispac, num2str(klin));
disp(dispac);
% 
dispac = 'Without refactorizations:   ';
dispac = concats(dispac, concats(num2str(kiterf),'  '));
dispac = concats(dispac, concats(num2str(kfactf), '  '));
dispac = concats(dispac,num2str(klinf));
disp(dispac);


  end;
end;



%  EOF



