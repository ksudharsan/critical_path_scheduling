function [xn,sn,ln,mn,kiter,kfact,klin,nxs,rp,rd]=pclpm(A,b,c,x,s,km,minf,e,redlin,factor,verbose)
%                            PCLPM.SCI
%               A SCILAB/MATLAB IMPLEMENTATION OF THE  
%    PREDICTOR-CORRECTOR ALGORITHM WITH MULTIPLE CORRECTOR STEPS
%                 FOR LINEAR PROGRAMMING PROBLEMS
%                  VERSION OF SEPTEMBER 17, 2000
% 
%                     
%  ORIGIN
%  Frederic BONNANS, INRIA-Rocquencourt
%  Mounir HADDOU,  INRIA-Rocquencourt and Universite d'Orleans
% 
%  REFERENCE
%  F. Bonnans, J.-Ch. Gilbert, C. Lemarechal, C. Sagastizabal,
%  Optimisation numerique : aspects th\'eoriques et pratiques
%  Springer-Verlag, Paris, 1997. 
% 
%  DESCRIPTION
%   This function  solves the LP problem 
% 
%            min c'x    s.t Ax=b;  x>=0  
% 
%   by a large neihghborhood infeasible predictor_corrector algorithm
%   making Newton steps on the perturbed optimality system
%            x.*s         =  m*un
%            Ax           =  b
%            - s + A'l    = -c
%            x>= 0 ,    s>=0
%            m--->0
%  The matrix A may be either full or sparse; computations are made
%  accordingly. 
% 
%  CALLING SEQUENCE
%     [xn,sn,ln,mn,kiter,kfact,klin,nxs,rp,rd] = pclpm(A,b,c)
%     [xn,sn,ln,mn,kiter,kfact,klin,nxs,rp,rd] = 
%             pclpm(A,b,c,x,s,km,minf,e,redlin,factor,verbose)
% 
%  CALLING PARAMETERS
%  (Only A,b,c are compulsory; other parameters are optional)
% 
%  A (p x n)            constraint matrix (sparse or full) 
%  b (p x 1), c (n x 1) rhs of constraints and cost (full vectors)
%  x > 0, s > 0 (n x 1) starting points; default to ones  (full vectors)
%  km (integer)         maximum number of major iterations; default value 100
%  minf > 0             precision; default value  1.e-10
%  e, 0 < e < 1/2       size of the large neighborhood; default value 0.01
%  redlin               reduced linear system solved if redlin > 0 (default)
%  factor               single factorization in corrections if > 0 (default)
%  verbose              if positive, information displayed at each iteration
%                       if greater than 1: additional displays; 
%                       default value -1
% 
%  RETURN PARAMETERS
%  xn, sn (n x 1)       primal-dual optimal solutions
%  ln     (p x 1)       column vector of Lagrange multipliers
%  mn                   penalty parameter
%  kiter                number of major iterations
%  kfact                number of factorizations
%  klin                 number of linear systems solved
%  nxs                  complementarity: nxs = x'*s
%  rp  (resp. rd)       norm of primal (resp. dual) constraints

%  The optimality conditions are xn and sn nonnegative, and  
%  nxs, rp, rd close to 0. 

% 
%  EXAMPLE 
%  % getf('pclpm.sci');
%  n=100; p=20; c=rand(n,1); x=rand(n,1); A=rand(p,n); b=A*x;
%  [xn,sn,ln,mn,kn,nxs,rp,rd] = pclpm(A,b,c);
%  min(xn), min(sn), xn'*sn / n 
%  normest(A*xn-b), normest(c + A'*ln  - sn)

software = 'matlab';

[p,n] = size(A);
un = ones(n,1);

%  ALGORITHM PARAMETERS
%  Upper bound for m
minfdef=1.e-10;
%  Size of neighborhood (default value)
edefault=0.01;
%  Reduction wrt maximum step
fr = 0.999;
%  Maximum number of inner iterations in linesearch for correction step
kiterc = 100; 
%  Relative reduction of step in linesearch for correction step
redc=.5;
%  epsilon value in computing maximal affine step (see function affstep)
epsaff = 1.e-10;
%  Threshold value of max(x.*s) / min(x.*s) for additional 
%  correction step
ratmax = 1.e2
%  Maximum number of correction steps in a major iteration
ksolveinm = 3; 

%  INITIAL VALUES

%  
ra = nargin;

%  SETTING DEFAULT VALUES
%  Algorithm parameters
if ra<11 % then  
  verbose = -1; 
end;
if ra<10 % then  
  factor  =  1; 
end;
if ra<9  % then  
  redlin  =  1; 
end;
if ra<8  % then  
  e = edefault; 
end;
if ra<7  % then
  minf = minfdef; 
end;
if ra<6  % then
  km = 100; 
end;

%  Starting point
if ra<5  % then
  s = (normest(c) + normest(A) ) * un;  
end;
if ra<4  % then
  x = (1 + normest(b) / normest(A)) * un;  
end;

%  INCONSISTANT INPUT DATA
if ra<3  % then  
  error('pclpm: minimum of 3 input argument needed');  
end;


%  CONSTANTS
einv = 1 / e ;
l = zeros(p,1);
m = x'*s/n;
gp = (b-A*x)/m;
%  gd = -(c+A'*l-s)/m;
gd = -(c+ (l' * A)' -s)/m;

%  Arrays needed when solving non reduced system
if redlin <= 0 % then 
  if issparse(A) % then 
%  Sparse matrices
    Opp=sparse(p,p);
    Onp=sparse(n,p);
    Opn=sparse(p,n);
    Onn=sparse(n,n);
    In= speye(n,n);
  else
%  Full matrices
   dispac = 'PCLPM Warning: full data and non reduced linear systems option';
    disp(dispac);
    Opp=zeros(p,p);
    Onp=zeros(n,p);
    Opn=zeros(p,n);
    Onn=zeros(n,n);
    In= eye(n,n);
  end;
  else
    Opp=[];
    Onp=[];
    Opn=[];
    Onn=[];
    In =[];
end;

%  Working arrays, used when simplified steps based on past
%  factorizations are used
Bfact = [];
B = [];
shat = [];
d = [];
sqrtxs = [];

%  Initialization of counters (see return parameters) 
kiter = 0;
kfact = 0;
klin = 0;

while ( (m > minf) & (kiter < km ) ) % then
kiter = kiter + 1;
%  STARTING MAJOR ITERATION
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if verbose > 0 % then 
  dispmes = concats('PCLPM: ITERATION  ',num2str(kiter) );
  disp(dispmes)
  end;
% %%%%%%%%%%%%%%%%%%
%   CORRECTION STEP     
% %%%%%%%%%%%%%%%%%%
  kfact=kfact+1;
  klin=klin+1;
  dofactor=1; 
  [x,l,s,B,Bfact,shat,sqrtxs,d] = solcm(A,b,c,gp,gd,m,x,l,s,e,factor,dofactor,B,Bfact,shat,sqrtxs,d,p,n,redlin,un,In,Onn,Onp,Opn,Opp,fr,kiterc,redc,verbose,software);
% %%%%%%%%%%%%%%%%%%
%   ADDITIONAL CORRECTION STEP    
% %%%%%%%%%%%%%%%%%%
  y = x.*s; 
  raty=max(y)/min(y);
  if verbose > 0 % then 
      dispac = concats('max(x.*s)/mu = ', num2str(max(y)/m));
      dispac = concats(dispac,concats('  min(x.*s)/mu = ',num2str(min(y)/m)));
      dispac = concats(dispac,concats('  ratio = ',num2str(max(y)/min(y))));
    disp(dispac); 
  end;
  ksolvein=0;
% 
  while ( raty > ratmax & ksolvein < ksolveinm) % then
  ksolvein=ksolvein+1;    
  if verbose > 0 % then 
      dispac = concats('major iteration = ',  num2str(kiter));
      dispac = concats(dispac,concats('  additional correction step ',num2str(ksolvein))); 
      disp(dispac); 
    end;
    if factor <=0 % then
      kfact=kfact+1;
    end
    klin=klin+1;
    dofactor=-1;
    [x,l,s,B,Bfact,shat,sqrtxs,d] = solcm(A,b,c,gp,gd,m,x,l,s,e,factor,dofactor,B,Bfact,shat,sqrtxs,d,p,n,redlin,un,In,Onn,Onp,Opn,Opp,fr,kiterc,redc,verbose,software);
    y = x.*s; 
    raty=max(y)/min(y);
    if verbose > 0 % then 
      dispac = concats('max(x.*s)/mu = ', num2str(max(y)/m));
      dispac = concats(dispac,concats('  min(x.*s)/mu = ', num2str(min(y)/m)));
      dispac = concats(dispac,concats('  ratio = ', num2str(max(y)/min(y))));
      disp(dispac); 
    end;
 end;
  % %%%%%%%%%%%%%%%%%%%%%
  % AFFINE (PREDICTION) STEP      sola_1.sci
  % %%%%%%%%%%%%%%%%%%%%%
  if factor <=0 % then
    kfact=kfact+1;
  end
  klin=klin+1;
  [x,l,s,m] = solam(A,b,c,gp,gd,m,x,l,s,e,factor,B,Bfact,shat,sqrtxs,d,p,n,redlin,un,In,Onn,Onp,Opn,Opp,epsaff,fr,verbose,software);
%  END OF MAJOR ITERATION
end

% %%%%%%%%%%%%%%%%%%
%  Results
% %%%%%%%%%%%%%%%%%%
xn = x;
sn = s;
ln = l;
mn = m;
nxs = x' * s;
rp = normest(A*x-b);
% rd = normest(c+A'*l+H*x-s);
rd = normest(c + (l' * A)' - s);
% 

function [x,l,s,B,Bfact,shat,sqrtxs,d]=solcm(A,b,c,gp,gd,m,x,l,s,e,factor,dofactor,B,Bfact,shat,sqrtxs,d,p,n,redlin,un,In,Onn,Onp,Opn,Opp,fr,kiterc,redc,verbose,software)
% 
%  MACRO COMPUTING THE CORRECTION DISPLACEMENT
% 
einv = 1/e;
%  1) COMPUTATION OF WEIGHTED CONSTRAINT MATRIX AND SCALES IF NECESSARY
if (factor <= 0 | dofactor > 0)
  d = sqrt(x./s)';
  shat=s;
  sqrtxs =  sqrt(x .* s);
  if issparse(A) % then 
    if software == 'scilab' % then 
      B = A * sparse([(1:n)',(1:n)'],d);
    else
      B = A * sparse((1:n)',(1:n)',d);
    end;
  else
    B(p,n)=0;
    for i=1:p,
      B(i,:) = A(i,:) .* d;
    end;
  end;
  d = d';
end;
% 
%  2) COMPUTATION OF NEWTON DIRECTION
% 
if redlin > 0 % then 
%  2a) REDUCED SYSTEM
%   B*B' dl = rhs
  r = B*((m*un - x.*s ) ./ sqrtxs );
  if factor > 0 % then
%   REDUCED SYSTEM: FACTORIZE
    if dofactor > 0   % then
%   Factorize sparse matrix of reduced system
      if issparse(B) % then 
        Bfact = chfact(B*B');
      else 
%     Factorize full matrix of reduced system
        Bfact = chol(B*B'); 
      end;
    end;
%   SOLVE
%   Sparse matrices
      if issparse(B) % then 
        dl = chsolvem(Bfact,r);
      else 
%   Full matrices
        dl = chsolvem(Bfact,r); 
      end
%   REDUCED SYSTEM: NO FACTORIZATION
    else 
      dl = (B*B')\r;
    end
    V = (dl' * A)';
    U = (m*un - x.*s)./shat - V.*(d.*d);
  else
%  2b) NON REDUCED SYSTEM
%  Matrix of original system (H=0 if LP)
%       M=[S    X    Onp
%          A    Opn  Opp
%          H    In    A' ];
%  Scaled matrix 
      M=[In    In   Onp;
         B     Opn  Opp;
         Onn  -In    B' ];
% 
      dir=M\[ (m*un - x.*s ) ./ sqrtxs; zeros(n+p,1) ];
      U  = dir(1:n) .* d;
      V  = dir((n+1):(2*n)) ./ d;
      dl = dir((2*n+1):(2*n+p));
  end   
%  
%  3) COMPUTATION OF STEP-SIZE
%  Initializations
y = x.*s/m; 
z = U.*V/m; 
q = (x.*V+s.*U)/m; 
tx = 1; 
ts = 1; 
%  x must be >0 : 
if min(U)<0 % then 
  tx = 1/max(-U./x); 
end 
%  s must be >0 : 
if min(V)<0 % then 
  ts = 1/max(-V./s); 
end 
t = min([1,fr * tx, fr * ts]); 
% 
%  Find the first feasible step-size in large neigborhood: 
%  Armijo-like linear search 
kv = 0; 
vois1 = y+t*q+(t^2)*z- e *un; 
vois1 = min(vois1); 
while (vois1 < 0) & (kv < kiterc) % then 
  t = redc*t; 
  vois1 = y+t*q+(t^2)*z- e *un; 
  vois1 = min(vois1); 
  kv = kv+1; 
end 
%  kv = 0; 
vois2 = y+t*q+(t^2)*z - einv * un; 
vois2 = max(vois2); 
while (vois2 > 0) & (kv < kiterc) % then 
  t = redc*t; 
  vois2 = y+t*q+(t^2)*z - einv * un; 
  vois2 = max(vois2); 
  kv = kv+1; 
end 

%  4) NEXT POINT (x,l,s) 

if verbose > 0 % then 
  disptc = concats('correction step: tc = ', num2str(t));
%   disptc = concats(disptc,concats(' reduction steps   = ' , num2str(kv)));
  disp( disptc );
end;
x = x+t*U; 
l = l+t*dl; 
s = s+t*V;


function [x,l,s,m]=solam(A,b,c,gp,gd,m,x,l,s,e,factor,B,Bfact,shat,sqrtxs,d,p,n,redlin,un,In,Onn,Onp,Opn,Opp,epsaff,fr,verbose,software);
% 
%  MACRO COMPUTING THE AFFINE DISPLACEMENT
% 
%  1) COMPUTATION OF WEIGHTED CONSTRAINT MATRIX AND SCALES IF NECESSARY
%  We could force refactorization in the affine step; 
%  For that set factora = -1; 
factora = factor;
% 
%  1) COMPUTATION OF WEIGHTED CONSTRAINT MATRIX AND SCALES IF NECESSARY
if factora <= 0 % then
  shat = s; 
  sqrtxs=sqrt(x.*s);
  d = sqrt(x./s)';
  if issparse(A) % then 
    if software == 'scilab' % then 
      B = A * sparse([(1:n)',(1:n)'],d);
    else
      B = A * sparse((1:n)',(1:n)',d);
    end;
  else
    B(p,n)=0;
    for i=1:p,
      B(i,:) = A(i,:) .* d;
    end;
  end;  
  d = d';
end;
% 
%  2) COMPUTATION OF NEWTON DIRECTION
% 
if redlin > 0 % then 
%  2a) REDUCED SYSTEM
%   B*B' dl = rhs
  r = B*(m * d .* gd +  (m*un - x.*s ) ./ sqrtxs ) - m*gp;
%   REDUCED SYSTEM: FACTORIZATION ALREADY DONE
  if factora > 0 % then
%   Sparse matrices
    if issparse(B) % then 
      dl = chsolvem(Bfact,r);
    else 
%   Full matrices
        dl = chsolvem(Bfact,r); 
    end
  else
%   REDUCED SYSTEM: NO FACTORIZATION
    dl = (B*B')\r;
  end
  V = (dl' * A)' - m*gd;
  U = (m*un - x.*s)./shat - V.*(d.*d);
  else
%  2B) NON REDUCED SYSTEM
%  Scaled matrix 
      M=[In    In   Onp;
         B     Opn  Opp;
         Onn  -In    B' ];
% 
      dir=M\[ (m*un - x.*s ) ./ sqrtxs; m * gp; m * d .* gd ];
      U  = dir(1:n) .* d;
      V  = dir((n+1):(2*n)) ./ d;
      dl = dir((2*n+1):(2*n+p));
  end; 
%  3) COMPUTING STEP-SIZE
y = x.*s/m; 
z = U.*V/m; 
q = (x.*V+s.*U)/m; 
tax = 1; 
tas = 1; 
%   x must be >0: 
if min(U)<0 % then 
  tax = 1/max(-U./x); 
end 
%   pause
%   s must be >0 : 
if min(V)<0 % then 
  tas = 1/max(-V./s); 
end 
ta = fr * min([1,tax,tas]);  

%  (x,s,m) must stay in Ne
t = aff(y,z,q,e,epsaff,verbose);

%  4) NEX POINT (x,l,s,m) 
x = x+t*U; 
l = l+t*dl; 
s = s+t*V; 
m = (1-t)*m;
if verbose > 0 % then 
  dispta = concats('affine     step: ta = ',num2str(t));
  dispta = concats(dispta,concats('   mu = ',num2str(m)));
  disp( dispta ); 
end;
% 
% 
function[t]=aff(y,z,q,e,epsaff,verbose);
%  Computation of maximum affine step 
%     (x,s,m) must stay in Ne
%     compute the largest t satisfying : t^2*z + t*(q+e*un) + y-e*un>=0 
      z1 = z./(y-e);
      q1 =(q+e)./(y-e);
%     t must satisfy t^2*z1 + t*q1 +un >= 0  
      t1 = affstep(z1,q1,epsaff,verbose);
%     t1<0:   error
% 
%     t  must satisfy t^2*z + t*(q+e^(-1)*un) + y-e^(-1)*un<=0
     z2 = z./(y-e^(-1));
     q2 =(q+e^(-1))./(y-e^(-1));
 %    t such that  t^2*z2 + t*q2 +un >= 0 
     t2 = affstep(z2,q2,epsaff,verbose);
%      t1<0: error
     t = min (t1,t2);
% 
%  
function [t]=affstep(z1,q1,epsaff,verbose)
t1 = 1; 
t2 = 1; 
t3 = 1; 
%  t must satisfy t^2*z1 + t*q1 +un >=0
dis1 = q1.^2 - 4*z1;
%    if dis1<0 there is no problem
%    if some component of of  dis1 are >=0
I1 = find(dis1>=0);
z1 = z1(I1);   
q1 = q1(I1);
  dis1 = dis1(I1);
%   the case  z1=0 is particular  
Io = find(z1<=epsaff & z1>=-epsaff);
qo = q1(Io);
%    solve the  inequality t*q1 +un >=0
%   si q1>=0 ==>t*q1 +un >=0 pour tout t>=0 
Jo = find(qo<0);
if ~isempty(Io) & ~isempty(Jo) 
% if Io ~= [] & Jo ~= [] % then
  t1=-1/min(qo(Jo));
end
%  
Ip=find(z1 > epsaff & q1 < 0);
 if ~isempty(Ip)
% if Ip ~= [] % then
  zp = z1(Ip);
  qp = q1(Ip);
  disip = dis1(Ip);
  solp =min( (-qp-sqrt(disip))./(2*zp) , (-qp+sqrt(disip))./(2*zp));
  t2 = min (solp);
end;
% 
In=find(z1 < -epsaff);
 if ~isempty(In)
% if In ~= [] % then 
  zn = z1(In);
  qn = q1(In);
  disn = dis1(In);
  soln =max( (-qn-sqrt(disn))./(2*zn) , (-qn+sqrt(disn))./(2*zn));
  t3 = min (soln);
end
if verbose > 2 % then  
  disp123 = concats(num2str(t1),concats(num2str(t2),num2str(t3)));
%   disp(disp123); 
end;
t = min([t1,t2,t3]);
% 
