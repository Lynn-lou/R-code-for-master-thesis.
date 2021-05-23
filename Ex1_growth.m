% Example 1:  Berkeley Growth Data
clear all;
file_in = 'growth.dat';          
file_out = 'growth.out';     
nrow = 2883;  
[idsubj,idclass,data] = datap(file_in,nrow);   
nc = 2;                 
pcopt.ops   =  0;  
clustopt.PP = 0.9;
clustopt.MM = 0.2; 
[nc_kcfc,Pout,Mout,idfpca,idkcfc] = kcfc(nc,pcopt,clustopt,data);
[aRandfpca,cRatefpca,orderfpca,tablefpca] = cmp2p(idfpca,idclass)
[aRandkcfc,cRatekcfc,orderkcfc,tablekcfc] = cmp2p(idkcfc,idclass)
% idclass : gender classification (1 = boy, 2 = girl)   
fidout = fopen(file_out,'w');
n = length(idsubj);
for i = 1:n;
  fprintf(fidout,'%u \t', idsubj(i),idclass(i),idfpca(i),idkcfc(i));
  fprintf(fidout,'\n');
end;
fclose(fidout);
