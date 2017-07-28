function [xn,sn,ln,mn,kn,nxs,rp,rd]=pclp(A,b,c,x,s,km,minf,e,redlin,verbose)
% 
%                            PCLPM.SCI
%               A SCILAB/MATLAB IMPLEMENTATION OF THE  
%                   PREDICTOR-CORRECTOR ALGORITHM 
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
%     [xn,sn,ln,mn,kn,nxs,rp,rd] = pclp(A,b,c)
%     [xn,sn,ln,mn,kn,nxs,rp,rd] = pclp(A,b,c,x,s,km,minf,e,redlin,verbose)
% 
%  CALLING PARAMETERS
%  (Only A,b,c are compulsory; other parameters are optional)
% 
%  A (p x n)            sparse or full matrix
%  b (p x 1), c (n x 1) rhs of constraints and cost
%  x > 0, s > 0 (n x 1) starting points; default to ones 
%  km (integer)         max number of iterations; default value 100
%  minf > 0             precision; default value sqrt(%eps) ~ 1.e-10
%  e, 0 < e < 1/2       size of the large neighborhood; default value 0.01
%  redlin               reduced linear system solved if redlin > 0 (default)
%  verbose              if positive, information displayed at each iteration
%                       if greater than 1: additional displays 
% 
%  RETURN PARAMETERS
%  xn, sn (n x 1)       primal-dual optimal solutions
%  ln     (p x 1)       column vector of Lagrange multipliers
%  mn                   penalty parameter
%  kn                   number of iterations
%  nxs                  complementarity: nxs = x'*s
%  rp  (resp. rd)       norm of primal (resp. dual) constraints

%  rd                   if unbounded rd = -1

%  The optimality conditions are xn and sn nonnegative, and  
%  nxs, rp, rd close to 0. 

% 
%  EXAMPLE 
%  % getf('pclp.sci');
%  n=100; p=20; c=-rand(n,1); x=rand(n,1); A=rand(p,n); b=A*x;
%  [xn,sn,ln,mn,kn,nxs,rp,rd] = pclp(A,b,c);
%  min(xn), min(sn), xn'*sn / n 
%  norm(A*xn-b), norm(c + A'*ln  - sn)


software = 'matlab';
unbounded = 0;
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
%  Maximum number of inner iterations in linesearch for affine step
kitera = 100; 
%  Relative reduction of step in linesearch for affine step
reda = .95;
%  epsilon value in computing maximal affine step (see function affstep)
epsaff = 1.e-10;
%  INIT VALUES

 ra = nargin;
if ra<10 % then
  verbose = -1; 
end
if ra<9  % then
  redlin  =  1; 
end
if ra<8  % then
  e = edefault; 
end
if ra<7  % then
  minf = minfdef; 
end
if ra<6  % then
  km = 100; 
end  
if ra<5  % then
  s =  (norm(c) + norm(A) ) * un; 
end  
if ra<4  % then
  x=  (1 + norm(b) / norm(A)) * un; 
end


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

k=0;

while ((m>=minf)&(k<km)) % then
  k = k+1;
%  STARTING MAJOR ITERATION
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose > 0 % then 
  dispmes = concats(' PCLP: ITERATION  ',num2str(k) );
  disp(dispmes);
%   fprintf(' PCLP: major iteration %6.0f\n ',k);
end;
  % %%%%%%%%%%%%%%%%%%
  % CORRECTION STEP     solc_1.sci
  % %%%%%%%%%%%%%%%%%%
  [x,l,s] = solc(A,b,c,gp,gd,m,x,l,s,e,p,n,redlin,un,In,Onn,Onp,Opn,Opp,fr,kiterc,redc,verbose,software);
  if norm(x) > 1e10
	  disp('unbounded')
	  unbounded = 1;
	  break;
  end
  % %%%%%%%%%%%%%%%%%%%%%
  % PREDICTION STEP      sola_1.sci
  % %%%%%%%%%%%%%%%%%%%%%

  [x,l,s,m] = sola(A,b,c,gp,gd,m,x,l,s,e,p,n,redlin,un,In,Onn,Onp,Opn,Opp,epsaff,fr,verbose,software);
%  END OF MAJOR ITERATION
end;

% %%%%%%%%%%%%%%%%%%
%  Results
% %%%%%%%%%%%%%%%%%%
xn = x;
sn = s;
ln = l;
mn = m;
nxs = x' * s;
rp = norm(A*x-b);
% rd = norm(c+A'*l+H*x-s);
rd = norm(c + (l' * A)' - s);
kn = k;
if unbounded == 1
	rd = -1;
end
% 
% 
% 
% 
%  called functions
% 
function [x,l,s]=solc(A,b,c,gp,gd,m,x,l,s,e,p,n,redlin,un,In,Onn,Onp,Opn,Opp,fr,kiterc,redc,verbose,software)
%  CORRECTION STEP
einv = 1/e;
%  Newton direction 
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
% 
  if redlin > 0 % then 
%   Reduced system B*B' dl = rhs
    r = B*((m*un - x.*s ) ./ sqrt(x .* s) );
    dl = (B*B')\r;
%   V = A'*dl;
    V = (dl' * A)';
    U = (m*un - x.*V)./s -x;
  else
%  Non reduced system
%  Matrix of original system (H=0 if LP)
%       M=[S    X    Onp
%          A    Opn  Opp
%          H    In    A' ];
%  Scaled matrix 
      M=[In    In   Onp;
         B     Opn  Opp;
         Onn  -In    B' ];
% 
      dir=M\[ (m*un - x.*s ) ./ sqrt(x .* s); zeros(n+p,1) ];
      U  = dir(1:n) .* d;
      V  = dir((n+1):(2*n)) ./ d;
      dl = dir((2*n+1):(2*n+p));
  end   
%  step-size 
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
%  The first feasible step-size in Ne : 
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

%  new (x,l,s) 

if verbose > 0 % then 
  disptc = concats('correction step: tc = ', num2str(t));
  disptc = concats(disptc,concats('   inner iterations = ', num2str(kv)));
  disp( disptc );
end;
x = x+t*U; 
l = l+t*dl; 
s = s+t*V;


function [x,l,s,m]=sola(A,b,c,gp,gd,m,x,l,s,e,p,n,redlin,un,In,Onn,Onp,Opn,Opp,epsaff,fr,verbose,software) 
%  AFFINE STEP
    [p,n] = size(A);
%  Newton direction 
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
  if redlin > 0 % then 
%   Reduced system B*B' dl = rhs
    r = B*(m * d .* gd +  (m*un - x.*s ) ./ sqrt(x .* s) ) - m*gp;
    dl = (B*B')\r;
    V = (dl' * A)' - m*gd;
    U = (m*un - x.*V)./s -x; 
  else
%  Non reduced system
%  Scaled matrix 
      M=[In    In   Onp;
         B     Opn  Opp;
         Onn  -In    B' ];
% 
      dir=M\[ (m*un - x.*s ) ./ sqrt(x .* s); m * gp; m * d .* gd ];
      U  = dir(1:n) .* d;
      V  = dir((n+1):(2*n)) ./ d;
      dl = dir((2*n+1):(2*n+p));
  end; 
%  step-size 
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

%  (x,s,m) must stay in Ne
tne = aff(y,z,q,e,epsaff,verbose);
ta = fr * min([1,tax,tas]);  
t = min(ta,tne);

%  new (x,l,s,m) 
x = x+t*U; 
l = l+t*dl; 
s = s+t*V; 
m = (1-t)*m;
if verbose > 0 % then 
  dispta = concats('affine     step: ta = ', num2str(t));
  dispta = concats(dispta, concats('   mu = ', num2str(m)));
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
%     if t1<0       keyboard     end
% 
%     t  must satisfy t^2*z + t*(q+e^(-1)*un) + y-e^(-1)*un<=0
     z2 = z./(y-e^(-1));
     q2 =(q+e^(-1))./(y-e^(-1));
 %    t such that  t^2*z2 + t*q2 +un >= 0 
     t2 = affstep(z2,q2,epsaff,verbose);
%      if t1<0       keyboard     end
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
if length(I1) > 0
z1 = z1(I1);   
q1 = q1(I1);
  dis1 = dis1(I1);
%   the case  z1=0 is particular  
Io = find(z1<=epsaff & z1>=-epsaff);
qo = q1(Io);
%    solve the  inequality t*q1 +un >=0
%   si q1>=0 ==>t*q1 +un >=0 pour tout t>=0 
Jo = find(qo<0);
if ~isempty(Io) & ~isempty(Jo) % then
% if Io ~= [] & Jo ~= [] % then
  t1=-1/min(qo(Jo));
end
%  
Ip=find(z1 > epsaff & q1 < 0);
 if ~isempty(Ip) % then
% if Ip ~= [] % then
  zp = z1(Ip);
  qp = q1(Ip);
  disip = dis1(Ip);
  solp =min( (-qp-sqrt(disip))./(2*zp) , (-qp+sqrt(disip))./(2*zp));
  t2 = min (solp);
end;
% 
In=find(z1 < -epsaff);
 if ~isempty(In) % then
% if In ~= [] % then 
  zn = z1(In);
  qn = q1(In);
  disn = dis1(In);
  soln =max( (-qn-sqrt(disn))./(2*zn) , (-qn+sqrt(disn))./(2*zn));
  t3 = min (soln);
end
if verbose > 1 % then  
  disp123 = concats(num2str(t1),concats(num2str(t2),num2str(t3)));
%   disp(disp123); 
end;
  t = min([t1,t2,t3]);
else
	t = 1e10;
end
% 
