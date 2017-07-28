function sortie = concats(stra,strb)
%  CONCATENATION OF STRINGS 
software = 'matlab';
if software == 'scilab' % then 
  sortie = stra + strb;
else
 sortie = [stra strb];
end;
