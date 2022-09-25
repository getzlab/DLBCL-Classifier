function [X1]=trimStruct(X,k) 
f=fieldnames(X);
kmax=max(k);
if ((kmax<2)&(length(k)>2))|any(k==0)  % k is array flag not list of indices
    k=find(k);
    kmax=max(k);
end

X1=X;
if (isfield(X,'N'))
    if (X.N>0)
        N=X.N;
        if (length(N)>1)
            error(' field N should be element count ');
        end
    end
end
if ~exist('N','var')
    N=0;
    for n=1:length(f)
        f1=char(f(n));
        x=X.(f1);
        m=size(x);
        N=max(N,m(1));
    end
end

N1=N;

for n=1:length(f)
    f1=char(f(n));
    x=X.(f1);
    m=size(x); 
    if ( (m(1)==N) && ~(strcmp(f1,'vlab')) && ~(strcmp(f1,'N'))  && ~(strcmp(f1,'head')) && ~(strcmp(f1,'hd')) )
       if (m(1)<kmax) 
            display(sprintf('ooops  trimStruct: %s\t%dx%d\n',f1,m))
            %continue;
            X1=[];
            return;
      end
      x1=x(k,:);
      X1.(f1)=x1;
      N1=length(x1);
    else
      X1.(f1)=x;    
    end
end

X1.N=N1;
X1.N=lengthStruct(X1); 
 
     
function N=lengthStruct(X) 
%function N=lengthStruct(X) 
% N=length of stuctural elements in X
f=fieldnames(X);
N=0;
for n=1:length(f)
    f1=f{n};
    m=size(X.(f1));
    % skip magic structure names
    if (~(strcmp(f1,'vlab')) && ~(strcmp(f1,'N')) && ~(strncmp(f1,'head',4)) && ~(strcmp(f1,'hd')))
       if (m(1)>N) & (~isstruct(X.(f1)))
            N=m(1);
       end
    end
end

function test()
 
 area='~/Projects/MobileElement/test/sim2/Spanner/build'
 f=[area '/sim2_fused.20.retro.span']
 what='sim2'

 % load all retro fragments 
A=loadRetroSpan(f);
k=find(bitget(A.element,1)>0);
A1=trimSPanStruct(A,k);