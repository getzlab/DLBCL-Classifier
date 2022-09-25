function k = strfindk(str,pat,opt)
inverse=0;
if (nargin>2)
    if (opt(1)=='e')
        k=strmatch(pat, str, 'exact');
        return
    end
    if (opt(1)=='v')
        inverse=1;
    end
end   
k1=strfind(str,pat);
L=length(k1); k=zeros(L,1);
if (iscell(k1))
   k=cellfun(@length,k1)>0;
else
   k=k1>0;
end

% for i=1:L 
%     if (iscell(k1))
%         if (~isempty(k1{i}))
%             k(i)=1;
%         end
%     else
%         if (~isempty(k1(i)))
%             k(i)=1;
%         end
%     end       
% end 
k=find(k);
if (inverse)
    k1=k;
    k=1:numel(str);
    k(k1)=[];
end
     