clear
cutoff = 0.000;

directory = 'C:\Users\eajoseph\Box\NucS manuscript\data\single_mismatch_figure_2\wt\';
mlist = dir([directory '*_R1*.fastq']);
s = 'CCGCAGACCCTGATCAACATCCGTCCCGTCGTGGCGGCGATCAAGGAGTTCTTCGGCACCAGCCAGCTGTCGCAGTTCATGGACCAGAACAACCCGCTGTCGGGTCTGACCCACAAGCGTCGTCTTTCGGCGCTGGGCCCCGGCGGTCTGTCCCGTGAGCGCGCCGGCCTCGAGGTCCGCGACGTGCACCCCAGCCACTACGGCCGCATGTGCCCGATCGAGACCCCTGAGGGTCCCAACATCGGTCTG';
t = 'GACCGCCGGGGCCCAGCGCCGAAAGACGACGCTTGCGGGTCAGACCCGACAGTGGGTTGTTCTGGTCCATG';
three = seqrcomplement('GACCGCCGGGGC');
five =  seqrcomplement('CTGGTCCATG');
oligos=['GACCGCCGGGGCCCAGCGCCGAAAGACGACGCTTGTGGGTCAGACCCGACAGCGGGTTGTTCTGGTCCATG'];
target = oligos(1,:);
% target(mm)='X';
% target = strrep(target,'X','[ATCG]');
%regexp(oligos,target,'match')
c=0;
lengths=zeros(numel(mlist),300);
for x = 1:numel(mlist);
    n1o = [directory mlist(x).name]
    n2o = strrep(n1o,'_R1_','_R2_');
    plot(lengths')
    drawnow
    f1 = fastqread([n1o]);
    f2 = fastqread([n2o]);
    
    
    
    for y = 1:numel(f1)
 
   f51= strfind(f2(y).Sequence,five);
   f31= strfind(f2(y).Sequence,three);
   f52= strfind(seqrcomplement(f1(y).Sequence),five);
   f32= strfind(seqrcomplement(f1(y).Sequence),three); 
   if (numel(f51)==1).*(numel(f31)==1).*(numel(f52)==1).*(numel(f32)==1)
   if length(f1(y).Sequence)==length(f2(y).Sequence)
       if sum(seqrcomplement(f2(y).Sequence)==(f1(y).Sequence))==length(f1(y).Sequence)
           L = length((f2(y).Sequence));
           lengths(x,L)=lengths(x,L)+1;
       end
   end
   end
   end
end