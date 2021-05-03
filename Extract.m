clearvars
Name = ["atk1_count.txt","atk3_count.txt","braf_count.txt","kras12_count.txt","kras61_count.txt","map2k1_57_count.txt","map2k1_121_count.txt","map2k1_203_count.txt","map2k2_57_count.txt","map2k2_125_count.txt","map2k2_207_count.txt","nras12_count.txt","nras61_count.txt","pik542_count.txt","pik1047_count.txt"];
Blocker = ["CGCTGCCTGCAGTGGACCACT", "GCCTCCAGTTTTTTATATATTCTCCTACATGAGG", "CCATCGAGATTTCACTGTAGCTAGACCAAAA", "TGCCTACGCCACCAGCTCCA", "CACTGTACTCCTCTTGACCTGCTGTG", "TCTTACCCAGAAGCAGAAGGTGGGA", "CTGCATGAGTGCAACTCTCCGTACA", "AAAGTCACAGAGCTTGATCTCCCCAC", "TTGGCTTTCTGGGTGAGAAAGGCTT", "CTGCACGAATGCAACTCGCCGTA", "TCACACAGCTTGATCTCCCCTCTAGA", "CCCAACACCACCTGCTCCAACC", "CACTGTACTCTTCTTGTCCAGCTGTATC", "ATCCTCTCTCTGAAATCACTGAGCAGG", "CAAATGAATGATGCACATCATGGTGGC"];
Start = [16, 19, 12, 14, 116, 13, 12, 11, 11, 11, 80, 13, 88, 17, 19];
%mp1_57 shifted 1 position down since text file missing the first nucleotide position, enrichment region star and end
End = [36, 52, 42, 33, 141, 37, 36, 36, 35, 33, 103, 34, 115, 43, 45];
pv_all = cell(1,1);
%read each text file
for j = 1:15
T = readtable(Name(1,j));
Total = T.Var2 + T.Var3 + T.Var4 + T.Var5;
%convert blocker seq to numeric
B = convertStringsToChars(Blocker(1,j));
Bl= length(B);
Bp = fromGCAT(B);
n = 1;
A = [];
%store number of variant nucleotides within the blocker region in A
for i = Start(1,j):End(1,j)
       if Bp(1,n) == 1
           A(n,1) = T.Var3(i,1);
           A(n,2) = T.Var4(i,1);
           A(n,3) = T.Var5(i,1);
           n = n+1;
       elseif Bp(1,n) == 2
           A(n,1) = T.Var2(i,1);
           A(n,2) = T.Var4(i,1);
           A(n,3) = T.Var5(i,1);
           n = n+1;
       elseif Bp(1,n) == 3
           A(n,1) = T.Var2(i,1);
           A(n,2) = T.Var3(i,1);
           A(n,3) = T.Var5(i,1);
           n = n+1;
       elseif Bp(1,n) == 4
           A(n,1) = T.Var2(i,1);
           A(n,2) = T.Var4(i,1);
           A(n,3) = T.Var3(i,1);
           n = n+1;
       end 

    
end
%find the max variant nucleotide and calculate var percent
v = max(A,[],2);
tn = Total (Start(1,j):End(1,j));
pv = (v./tn)*100;
%add percent var matrix into array
pv_all = [pv_all pv];


end
%delete the first empty cell
pv_all(1) = [];  
%transpose array and convert to table with rownames
Var_tab = cell2table(pv_all.', 'RowNames',{'Atk1' 'Atk3' 'Braf' 'kras12' 'kras61' 'map1_57' 'map1_121' 'map1_203' 'map2_57' 'map2_125' 'map2_207' 'nras12' 'nras61' 'pik542' 'pik1047'});
writetable(Var_tab, 'Var_tab.csv','WriteRowNames',true);
function [output]=fromGCAT(inputseq)
for i = 1:length(inputseq)
    if ((inputseq(i) == 'A')||(inputseq(i) == 'a'))
        output(i) = 1;
    elseif ((inputseq(i) == 'T')||(inputseq(i) == 't')||(inputseq(i) == 'U')||(inputseq(i) == 'u'))
        output(i) = 4;
    elseif ((inputseq(i) == 'C')||(inputseq(i) == 'c'))
        output(i) = 2;
    elseif ((inputseq(i) == 'G')||(inputseq(i) == 'g'))
        output(i) = 3;
	elseif (inputseq(i) == ' ')
		output(i) = 0;
	else
        fprintf('Unexpected nucleotide!\n');
	end
       

end
end


