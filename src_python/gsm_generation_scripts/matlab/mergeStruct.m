function [X]=mergeStruct(X1,X2) 
% function [X]=mergeStruct(X1,X2) 
% X= X2 tacked on to X1 with fields matching X1
X2=fieldStruct(X2,X1);
X1.N=lengthStruct(X1);
X2.N=lengthStruct(X2);
if (X2.N<1)
    X=X1;
    return;
end
f=fieldnames(X1);
f2=fieldnames(X2);
headlines=[];
if any(ismember(f,{'headline'}))
    headlines=X1.headline(:)';
end
if any(ismember(f2,{'headline'}))
    X2.headline = X2.headline(:)';
    headlines=[headlines(1:(end-1)) X2.headline];
end
headlines=unique(sort(headlines));
q=ismember(f,{'vlab','headline','header'});
X1=rmfield(X1,f(q));
f=f(~q);
if (X1.N<1)
    XT=X1;
    X1=X2;
    X2=XT;
end
X=X1;
N=X1.N+X2.N;
%pastN=false;
k=[];
for n=1:length(f)
    if ~(strcmp(f,'vlab'))
        m=size(X1.(char(f(n))));
        k=[k  m(1)];
    end
end    
kmax=max(k);
fu='';
for n=1:length(f)
    f1=char(f(n));
    m=size(X1.(f1));
    %if ( (m(1)>1) && ~(strcmp(f1,'vlab'))  )
    if ( (m(1)>=kmax) && ~(strcmp(f1,'vlab')) && ~(strcmp(f1,'N')) && ~(strncmp(f1,'head',4)) && ~(strcmp(f1,'hd')))
       if (m(1)<kmax) 
            X=[];
            display('ooops')
            return;
       end
    %if ((m(1)>1) && ~strcmp(f1,'vlab'))
    %if (isnumeric(X1.(f1)) && pastN)
    %if ( (~strcmp(f1,'vlab')) && pastN)
      L1=size(X1.(f1),2);
      L2=size(X2.(f1),2);
      fu=f1;
      if (L1==L2)
          if (iscellstr(X1.(f1))&~iscellstr(X2.(f1)))
              X2.(f1)=cellstr(num2str(X2.(f1)));
          end
          %if (iscellstr(X2.(f1))&~iscellstr(X1.(f1)))
          if (iscellstr(X2.(f1))&isnumeric(X1.(f1)))
              X1.(f1)=cellstr(num2str(X1.(f1)));
          end          
          X.(f1)=[X1.(f1); X2.(f1)];
          N=size(X.(f1),1);
      else
          X.(f1)=strvcat(X1.(f1), X2.(f1));
      end
    end
    %pastN=pastN||strcmp(f1,'N');
end
if length(headlines)>0
    X.headline=headlines;
end
X.N=N;

function [X]=fieldStruct(X1,XF) 
% function [X]=fieldStruct(X1,XF) 
% X = revised X1 such that X has same fields as XF
%     adds empty fields missing in X1 present in XF
%     removes extra fields not in XF

f1=fieldnames(X1);
ff=fieldnames(XF);
X=X1;
af=setdiff(ff, f1);
af(strcmp(af,'N'))=[];
X.N=lengthStruct(X);
for i=1:length(af)
    x=XF.(af{i});
    if (iscellstr(x))
       X.(af{i})=repmat({''},X.N,1);
    else
       X.(af{i})=NaN*zeros(X.N,1);
    end   
end
f=fieldnames(X);
mf=setdiff(f,ff);
for i=1:length(mf)
    X=rmfield(X,mf{i});
end


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
A1=trimSpanStruct(A,k);
k=find(bitget(A.element,4)>0);
A2=trimSpanStruct(A,k);
A12=mergeSpanStruct(A1,A2);
