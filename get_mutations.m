% Code to determine mutation rates from OR NGS data
% EAJ 10/31/2023

clear
cutoff = 0.000;
names = 'triple_mismatch_figure_3'; %Name of experiment

fignum = 1;
expt = 'nucsko'; %Datasets labeled with 'wt' for wild-type or 'nucsko' for NucS Knock-out strain
directory = ['C:\Users\eajoseph\Box\NucS manuscript\data\' names '\' expt '\']; %Directory with data

mlist = dir([directory '*_R1*.fastq']);
s = 'CCGCAGACCCTGATCAACATCCGTCCCGTCGTGGCGGCGATCAAGGAGTTCTTCGGCACCAGCCAGCTGTCGCAGTTCATGGACCAGAACAACCCGCTGTCGGGTCTGACCCACAAGCGTCGTCTTTCGGCGCTGGGCCCCGGCGGTCTGTCCCGTGAGCGCGCCGGCCTCGAGGTCCGCGACGTGCACCCCAGCCACTACGGCCGCATGTGCCCGATCGAGACCCCTGAGGGTCCCAACATCGGTCTG';

three = ('GAAAGACGACGCTTGTGGGT')
five = ('CTCCTTGATCGCCGCCACGA')    
oligos = ['CAGACCCGACAGCGGGTTGTTCTGGTCCATGAACTGCGACAGCTGGCTGGTGCCGAAGAA'];
rifr = 35;
rifrs = 'A';

target = oligos(1,:);
fid = fopen(['C:\Users\eajoseph\Box\NucS manuscript\data\' names '.txt'],'r')
while ~feof(fid)
    oligos = [oligos;fgetl(fid)]
end
fclose(fid)
target = oligos(1,:);
mm = find(sum(oligos~=target)~=0);
target(mm)='X';
target = strrep(target,'X','[ATCG]');
%
c=0;

mutG = zeros(numel(mlist),length(oligos(1,:)));
mutA = zeros(numel(mlist),length(oligos(1,:)));
mutT = zeros(numel(mlist),length(oligos(1,:)));
mutC = zeros(numel(mlist),length(oligos(1,:)));

for x = 1:numel(mlist);
    fid = fopen([directory 'oligo_region' num2str(x) '.txt'],'w+');
    n1o = [directory mlist(x).name]
    n2o = strrep(n1o,'_R1_','_R2_');

    f1 = fastqread([n1o]);
    f2 = fastqread([n2o]);
    
    for y = 1:numel(f1)
       forw = char(cellstr(regexp(f2(y).Sequence,target,'match')));
       rev = char(cellstr(regexp(seqrcomplement(f1(y).Sequence),target,'match')));

    if (length(forw)>0).*(length(rev)>0)
    if sum(forw==rev)==length(forw)
        if forw(rifr)==rifrs;
           fprintf(fid,'%s\n',forw);
        end
    end
    end
    end
    fclose(fid);
    fid = fopen([directory 'oligo_region' num2str(x) '.txt'],'r');
    fid2 = fopen([directory 'mutations_at_hotspots' num2str(x) '.txt'],'w+');

    while ~feof(fid)
        r = fgetl(fid);
        for m = 1:numel(mm);
            if r(1,mm(m))~=oligos(1,mm(m))
                if r(1,mm(m))=='A'
                    mutA(x,mm(m))=mutA(x,mm(m))+1;
                    fprintf(fid2,'%s','A');
                elseif r(1,mm(m))=='T'
                    mutT(x,mm(m))=mutT(x,mm(m))+1;
                    fprintf(fid2,'%s','T');
                elseif r(1,mm(m))=='G'
                    mutG(x,mm(m))=mutG(x,mm(m))+1;
                    fprintf(fid2,'%s','G');
                elseif r(1,mm(m))=='C'
                    mutC(x,mm(m))=mutC(x,mm(m))+1;
                    fprintf(fid2,'%s','C');
                end
            else
                    fprintf(fid2,'%s','n');
            end
        end
         fprintf(fid2,'\n');
    end
    fclose(fid) ;
    fclose(fid2) ;
    
   end
mut = [mutG+mutC+mutA+mutT]';
mA = mutA(:,sum(mut')~=0)./(max(mut))'*100;

mG = mutG(:,sum(mut')~=0)./(max(mut))'*100;
mC = mutC(:,sum(mut')~=0)./(max(mut))'*100;
mT = mutT(:,sum(mut')~=0)./(max(mut))'*100;
mA(mA==0)=nan;
mC(mC==0)=nan;
mG(mG==0)=nan;
mT(mT==0)=nan;
L = 1:length(oligos(1,:));
L = L(sum(mut')~=0);
figure(fignum);
hold off
%stem(L-.75/2,mA','g.')
scatter(L-.75/2,mA',20,'g','filled','MarkerEdgeColor',[0 0 0],'LineWidth',1)
hold on
scatter(L-0.25/2,mT',20,'r','filled','MarkerEdgeColor',[0 0 0],'LineWidth',1)
scatter(L+.25/2,mG',20,'k','filled','MarkerEdgeColor',[0 0 0],'LineWidth',1)
scatter(L+.75/2,mC',20,'c','filled','MarkerEdgeColor',[0 0 0],'LineWidth',1)
%stem(L-0.25/2,mT','r.')
%stem(L+.25/2,mG','k.')
%stem(L+.75/2,mC','b.')


xticks([1 mm length(oligos(1,:))])
axis([min(mm)-1,max(mm)+1,-0.5 20.5])
gca.FontSize = 6;
drawnow
variant_list='';
for x = 1:numel(mlist);
    variants='';
fid = fopen([directory 'mutations_at_hotspots' num2str(x) '.txt'],'r');
while ~feof(fid)
    variants = [variants; fgetl(fid)];
end
eval(['v' num2str(x) ' = variants;']); 
variant_list = unique([variants;variant_list],'rows');
end
fclose(fid);

v = zeros(size(variant_list,1), 3);
for x = 1:size(variant_list,1)
for y = 1:numel(mlist)
    eval(['vv = v' num2str(y) ';']);
    v(x,y)=sum(sum(vv==variant_list(x,:),2)==length(variant_list(1,:)));
end
end
meanv = (mean(v./sum(v),2)')';
goodv = find(meanv>=0.003);
goodvar = variant_list(goodv,:);
vorder = sortrows([meanv(goodv) (1:numel(goodv))' ],'descend');
mordered = vorder(:,1);
vorder = vorder(:,2);
sortgoodvar = goodvar(vorder,:);
varmat = nt2int(sortgoodvar);
varmat(varmat==15) = 0;
hold off
figure(fignum+1)
imagesc(varmat)
 xticks([1:length(mm)])
 xticklabels(mm)
 yticks(1:numel(goodv));
 yticklabels(mordered);
