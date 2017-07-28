function x = chsolvem(r,b)
%  SOLVER OF LINEAR SYSTEM BASED ON THE 
%  (ALREADY COMPUTED) CHOLESKY FACTOR
%  This is a naive implementation, useful
%  only in absence of a 'chsolve' function. 
% 
%  INPUT
%  r Cholesky factor
%  b right-hand side
% 
%  OUTPUT
%  x solution of linear system
% 
x = r \ ( r' \ b ); 
%  EOF
