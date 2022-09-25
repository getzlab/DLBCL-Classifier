function X = load_tsv_v2(tsv_file, comment,delimeter,autoformat)
%LOAD_TSV_V2 load table into struct in UTF-8 coding
% inputs: 
%     tsv_file : input delimited table with header line of field names
%     comment  : character used at start of line to mark comment lines (default '#')
%     delimeter: delimiter character (',',' '.':'... default '\t' tab)
%     autoformat: 0: cellstr output fields, 1=number when possible (default 1) 
%
X=[];
if (nargin<2)
    comment='#';
end
if (nargin<3)
    delimeter='\t';
end
if (nargin<4)
    autoformat=1;
end

n=0;
fid = fopen(tsv_file, 'rt','n','UTF-8');
while true
    thisline = fgetl(fid);
    if ~ischar(thisline); break; end  %end of file
    if thisline(1)=='#'; continue;  end % comment
    n=n+1;
    fields=regexp(thisline,delimeter,'split');
    if (n==1)
        headers=fields;
    else
        val(n-1,:)=fields;
    end
end
fclose(fid);
nv=length(fields);
for v=1:nv
    f1=genvarname(headers{v});
    X.(f1)=val(:,v);
    if autoformat & ~ismember(lower(f1),{'chromosome','chr','chrom'})
        v=X.(f1);
        x=str2double(v);
        k=find(cellfun(@length,v)>0);
        if any(isnan(x(k))), continue;end
        X.(f1)=x;
    end
end

end

function test
maf1='/Users/stewart/Downloads/16403_COH-006_W1D1_PL_T_v1-16403_COH-006_W1D1_PBMC_N_v1.MuTect1_SNV.MuTect2_INDEL.Strelka_INDEL.annotated.maf'
X=load_tsv_v2(maf1)
X=load_tsv_v2(maf1,'#','\t',0)
end

