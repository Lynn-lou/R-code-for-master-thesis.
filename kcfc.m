function [nc_kcfc,Pout,Mout,idfpca,idkcfc] = kcfc(nc_in,pcopt,clustopt,data);

% Cluster n curves into nc_in clusters by k-centers FC algorithm.
% Input
%  nc_in   - total number of clusters for initial clustering
%  pcopt - parameters for FPCA, including:
%    1. If pcopt.ops = 0, then use pc1.m without smoothing. 
%       If pcopt.ops = 1, then use pc2.m with smoothing steps.
%    2. When pcopt.ops = 1, then user must set the input arguments:
%       pcopt.bw1d and pcopt.bw2d (bandwidths for 1d and 2d smoothing)
%  clustdopt - parameters of clustering method, including: 
%    clustopt.PP - threshold value (\tau_{\lambda}) for selecting  
%                        the number FPC scores for initial clustering.  
%    clustopt.MM - threshold value (\tau_D) for selecting the number 
%                        functional principal components for iterative
%                        reclassificaiton. 
%  data    - includeing:
%          data.isobs - n x m matrix of indicator variables.
%                       (m is the max number of time points among subjects)
%                       isobs(i,j) = 1: the ith curve is observed 
%                       isobs(i,j) = 0, the ith curve is missing or not  
%                                       observed at time point Tin(i,j).
%            data.Tin - n x m matrix of observed time points of n subjects. 
%            data.Yin - n x m matrix of observations corresponding to Tin
% Output
%   nc_kcfc - total number of clusters
%   Pout   - the number of random components used in initial clustering
%   Mout   - (1 x nc_in) vector for the number of components used in 
%             the truncated K-L expansion for each group.
%   idfpca  - n x 1 vector of initial clustering results (cluster labels)
%   idkcfc  - n x 1 vector of final clustering results (cluster labels)
  
      ops = pcopt.ops;
      if ops == 0;
         bw1d = 0;   bw2d = 0;        
      elseif ops == 1;
         bw1d = pcopt.bw1d;  bw2d = pcopt.bw2d;
      end;
      PP = clustopt.PP;                  
      MM = clustopt.MM;      
      iter = 100;    
      isobs = data.isobs; Tin = data.Tin; Yin = data.Yin;
      m = size(Yin,2); n = size(Yin,1);
      nobs = sum(isobs,2);
       
      if ops == 0;
         tin = Tin(1,:)';  t_unique = tin;   nt = m;
         [meanfcn,covfcn,eigfcn,eigval,Pcs,varprop] = ...
              pc1(t_unique,Yin);
      elseif ops == 1;                                          
         [t_unique] = unique(Tin(isobs == 1));
         nt = length(t_unique);
         [meanfcn,covfcn,eigfcn,eigval,Pcs,varprop,evar] = ...
            pc2(nt,bw1d,bw2d,isobs,Tin,Yin,t_unique);
      end;
      
      for i = 1:m; peval(i) = sum(varprop(1:i)); end; 
      idp = find(peval >= PP);  np = min(idp);  prop = peval(np);
      Pout = np; 
      iseed = 232316; rand('state',iseed);
      [idc,C,sumd,D] = kmeans(Pcs(:,[1:np]),nc_in,'distance','sqEuclidean');
      idfpca = idc;
      
      nclust = nc_in;
      btin = zeros(n,m);   btin(:,[2:m]) = Tin(:,[1:m-1]);   
      Delta = Tin - btin;  Delta(isobs == 0) = 0;
      indx_tin = zeros(n,m);
      for j = 1:nt;
          id = Tin == t_unique(j);  indx_tin(id) = j;
      end; 
      
      reCluster = zeros(n,iter);
      ISE = zeros(iter+1,nclust);
      groupM = ISE;
      [uniqclust] = unique(idfpca);
      for kk = 1:iter;     
          idold = idc;
          [groupM(kk,:),groupMufcn,groupEigfcn,groupEigval,groupevar,...
            ISE(kk,:)] = allsse(nt,nclust,idold,uniqclust,ops,...
            MM,bw1d,bw2d,isobs,Tin,Yin,Delta,t_unique,indx_tin);  
          diff = 10000000 * ones(n,nclust);
          for i = 1:n;  
              Ypred = zeros(nclust,m);   
              for k = 1:nclust;   
               if length(find(idold == uniqclust(k))) <= 1; continue; end;
               npcs = groupM(kk,k);  
               if idold(i) == k; 
                  id = find(idold == uniqclust(k)); id(id == i) = [];  
                  ncurve = length(id);
                  if ncurve <= 1; continue; end;
                  tempisobs = isobs(id,:); tempTin = Tin(id,:);
                  tempYin = Yin(id,:);     tempD = Delta(id,:);
                  if ops == 0;
                     [mufcn,YCov,tempEigfcn,tempEigval,PCs] = ...
                         pc1(t_unique,tempYin);
                  elseif ops == 1; 
                     [mufcn,covfcn,tempEigfcn,tempEigval,tempFpcs,varprop,evar] = ...
                         pc2(nt,bw1d,bw2d,tempisobs,tempTin,tempYin,t_unique);
                  end;                 
               else;    
                  id = find(idold == uniqclust(k));
                  mufcn = groupMufcn(:,k);
                  tempEigfcn = groupEigfcn(:,:,k);
                  tempEigval = groupEigval(:,k);  
              end; 
              if ops == 0;
                 auxmufcn = mufcn';
                 if npcs > 0; auxefcn = tempEigfcn(:,[1:npcs]); end;
               elseif ops == 1;
                 auxmufcn = zeros(1,m);
                 auxmufcn(1,[1:nobs(i)]) = mufcn(indx_tin(i,[1:nobs(i)]))';
                 if npcs > 0;
                    auxefcn = zeros(m,npcs);
                    auxefcn([1:nobs(i)],:) = ...
                    tempEigfcn(indx_tin(i,[1:nobs(i)]),[1:npcs]);
                 end;
              end;
              Ypred(k,:) = auxmufcn;
              if npcs > 0;
                 Y_i = Yin(i,:);  T_i = Tin(i,:);  isobs_i = isobs(i,:);
                 index_i = indx_tin(i,:);  temp_i = [1:nobs(i)];
                 pcs_i = ((Y_i(temp_i) - auxmufcn(temp_i)).*...
                             Delta(i,temp_i))*auxefcn(temp_i,:);
                 if ops == 1;   
                    auxeigval = tempEigval(1:npcs)';
                    range_t = range(t_unique);
                    rho_0 = range_t*groupevar(k);
                    temp = range_t*groupevar(k)/sum(isobs_i);
                    sk = temp(:,ones(npcs,1));
                    temp1 = (auxeigval + sk);
                    id1 = find(temp1 ~= 0); id2 = find(temp1 == 0);
                    sk(id1) = auxeigval(id1)./temp1(id1);
                    sk(id2) = 1;
                    pcs_i = sk.*pcs_i;
                 end;
                 Ypred(k,:) = Ypred(k,:) + pcs_i*auxefcn(:,[1:npcs])';
              end;
              diff(i,k) = sum(((Yin(i,:) - Ypred(k,:)).^2).*Delta(i,:));
           end;   
           tempdiff = diff(i,:);
           idtemp = find( tempdiff == min(tempdiff) ); 
           idc(i) = idtemp(1);
          end;   
          reCluster(:,kk) = idc;   
          idkcfc = idc; 
          true = 0;
          if idfpca == idkcfc;  %^^
             true = 1;
          else;
             for j = 1:kk-2;
                 if reCluster(:,j)==idkcfc; true = 2; break; end;
             end;
          end;
          if true == 1; idselect = 1; break;
          elseif true == 2;
             idnan = find(ISE == -99);
             if length(idnan) > 0;
                tempSSE = ISE; tempSSE(idnan) = 0;
                sumSSE = sum(tempSSE,2);
             else;
                sumSSE = sum(ISE,2); 
             end;   
             idselect = find(sumSSE(1:kk-1)==min(sumSSE(1:kk-1)));
             idselect = idselect(1);
             if idselect==1; idkcfc = idfpca;
             else; idkcfc = reCluster(:,idselect-1);
             end; break;
          end;                 %==
             idmove = logical(idc == idold);
          nmove = length(find(idmove == 0));
          if nmove < 1; idselect = kk;  break;   end;
          if kk == iter;    
             idselect = kk;
             disp(['Warning- kcfc no covergence after ', ...
                 num2str(iter), ' iterations.']);
          end;
      end; 
 
      Mout = groupM(idselect,:);
      [unique_c] = unique(idkcfc);
      nc_kcfc = length(unique_c);
     

function [groupM,groupMufcn,groupEigfcn,groupEigval,groupevar,groupISE] = ...
          allsse(nt,nclust,idclust,uniqclust,ops,MM,bw1d,bw2d,...
          isobs,Tin,Yin,Delta,t_unique,indx_tin);
  nobs = sum(isobs,2);            groupM = zeros(1,nclust);
  groupMufcn = zeros(nt,nclust);  groupEigfcn = zeros(nt,nt,nclust);
  groupEigval = zeros(nt,nclust); groupevar = zeros(nclust,1);
  groupISE = zeros(1,nclust);     
  for k = 1:nclust;     
      id = find(idclust == uniqclust(k)); 
      ncurve = length(id);
      if ncurve <= 1;      
         groupISE(k) = -99;  gruopM(k) = -99; 
      elseif ncurve > 1;
         tempisobs = isobs(id,:);   tempTin = Tin(id,:);
         tempYin = Yin(id,:);       tempD = Delta(id,:);
         if ops == 0;
            [ymean,ycov,eigvec,eigval,PCs] = pc1(t_unique,tempYin);
            groupMufcn(:,k) = ymean;   
            groupEigfcn(:,:,k) = eigvec;
            groupEigval(:,k) = eigval; 
            auxmufcn =  ymean(:,ones(ncurve,1))';
            auxefcn = eigvec(:,:,ones(ncurve,1));
            auxeval = eigval;
            m=nt; %--
         elseif ops == 1;                
            [mufcn,covfcn,tempEigfcn,eigval,PCs,varprop,evar] = ...
            pc2(nt,bw1d,bw2d,tempisobs,tempTin,tempYin,t_unique); 
            groupMufcn(:,k) = mufcn;     
            groupEigfcn(:,:,k) = tempEigfcn;
            groupEigval(:,k) = eigval;
            groupevar(k) = evar;
            m = size(Tin,2);
            auxmufcn = zeros(ncurve,m);
            auxefcn = zeros(m,nt,ncurve); 
            for i = 1:ncurve;
              j = id(i);
              auxmufcn(i,[1:nobs(j)]) = mufcn(indx_tin(j,[1:nobs(j)]))'; 
              auxefcn([1:nobs(j)],:,i) = ...
              tempEigfcn(indx_tin(j,[1:nobs(j)]),[1:nt]); 
            end;
            auxeval = eigval;
        end; 
        ncurve = size(tempYin,1); 
        for nk = 1:nt+1;  
          PY = auxmufcn;
          if nk-1 > 0;
            for j = 1:nk-1;
              temp1 = PCs(:,j);
              temp2 = zeros(m,ncurve); %--
              temp2(:,:) = auxefcn(:,j,:);
              PY = PY + temp1(:,ones(m,1)).*temp2'; %--
            end;
          end;
          tempsse(nk) = sum(sum(((tempYin - PY).^2).*tempD,2));      
          if nk == 1;
            Ratio_sse(nk) = 1;
          else;
            Ratio_sse(nk) = (tempsse(nk-1)-tempsse(nk))/tempsse(1);
          end;
          if Ratio_sse(nk) >= MM;
            groupM(k) = nk-1; groupISE(k) = tempsse(nk);
          else;
            break;
          end;
        end;  
     end;  
  end;  
    
  
function [YMean,YCov,eigvec,eigval,PCs,varprop] = pc1(tin,Yin);
  n = size(Yin,1);  m = size(Yin,2); 
  Tin = tin(:,ones(n,1))';
  YMean = sum(Yin,1)/n; 
  YMean = YMean';
  YCov = (Yin-YMean(:,ones(n,1))')'*(Yin-YMean(:,ones(n,1))')/n;
  eigval = zeros(m,1);   
  eigvec = zeros(m,m);
  bt = zeros(m,1);    
  bt(2:m) = tin(1:m-1);   
  delta = tin-bt; 
  id = find(delta == 0);
  delta(id) = 0.0001;
  quadwt = diag(delta);
  tempA = sqrt(quadwt)*YCov*sqrt(quadwt);
  [eigvecA,eval] = eig(tempA);
  temp = 1./sqrt(delta);
  invW = diag(temp);
  evec = invW*eigvecA;
  eval = diag(eval);
  [seval,index] = sort(eval);
  for i = 1:m;
      eigval(i) = eval(index(m-i+1));  
      eigvec(:,i) = evec(:,index(m-i+1));
  end;
  eigval(find(eigval< 0)) = 0;
  varprop = eigval/sum(eigval); 
  PCs = zeros(n,m);
  btin = zeros(n,m);   
  btin(:,[2:m]) = Tin(:,[1:m-1]);   
  Delta = Tin - btin;
  for k = 1:m;
      PCs(:,k) = ((Yin-YMean(:,ones(n,1))').*Delta)*eigvec(:,k); 
  end;


function [Meanfcn,Covfcn,Eigfcn,Eigval,Pcs,varprop,evar] = ...
    pc2(npcs,bw1d,bw2d,isobs,Tin,Yin,Tout);
  n = size(Yin,1);  m = size(Yin,2);  mp = size(Tout,1);
  nobs = sum(isobs,2);
  tempTin = reshape(Tin',[],1);
  tempisobs = reshape(isobs',[],1);
  tempTin(tempisobs == 0) = [];
  t_all = cat(1,tempTin,Tout);
  [t_unique] = unique(t_all);
  nt = length(t_unique);               
  indx_tin = ones(n,m);  indx_tout = zeros(mp,1);
  for k = 1:nt;
      id = find(Tin == t_unique(k));  indx_tin(id) = k;
      id = find(Tout == t_unique(k)); indx_tout(id) = k;
  end;           
  indx_tin(find(isobs == 0)) = -1; 
  delta = zeros(n,nt);  auxy = zeros(n,nt); win = zeros(nt,1);
  for j = 1:nt;
      [idr,idc] = find(Tin == t_unique(j)); 
      delta(idr,j) = 1;
      for i = 1:length(idr); auxy(idr(i),j) = Yin(idr(i),idc(i)); end;
      if length(idr) > 0; win(j) = 1; end;
  end;    
  sumauxy = sum(auxy,1); sumdelta = sum(delta,1); %^^
  id1 = find(sumdelta > 0);
  auxmfc = zeros(1,nt); 
  auxmfc(id1) = sumauxy(id1)./sumdelta(id1); %==
  npoly = 1; 
  [tempmu] = sm1(nt,nt,npoly,bw1d,t_unique,auxmfc',win,t_unique);
  Meanfcn = tempmu(indx_tout);
  Mufcn = zeros(n,m);
  for i = 1:n;
      for j = 1:nobs(i);
          if indx_tin(i,j) > 0;  Mufcn(i,j) = tempmu(indx_tin(i,j))'; end
      end;
  end;
  if (npcs == 0);  
    Covfcn = 0; Eigfcn = 0; Eigval = 0; Pcs = 0; varprop = 0; evar = 0;
  end;
  if (npcs > 0);  
    efcn = zeros(nt,npcs); Eigval = zeros(npcs,1);  
    auxr = zeros(n,nt); aumfc = zeros(n,nt);
    delta = zeros(n,nt); auxy = zeros(n,nt); auxmufcn = zeros(n,nt);
    for j = 1:nt;
      [idr,idc] = find(Tin == t_unique(j));    
      delta(idr,j) = 1;                        
      for i = 1:length(idr);
          auxy(idr(i),j) = Yin(idr(i),idc(i));      
          auxmufcn(idr(i),j) = Mufcn(idr(i),idc(i));  %--
      end;
    end; 
    auxr = (auxy-auxmufcn).*delta;  %^^
    Wsum = delta'*delta; 
    [idwr,idwc] = find(Wsum == 0); 
    id0 = find(Wsum == 0);  Wsum(id0) = 1;
    covest_raw = auxr'*auxr./Wsum;  %==
    count = nt*(nt-1);
    temp1 = t_unique(:,ones(nt,1))'; temp1 = reshape(temp1,[],1);
    temp2 = t_unique(:,ones(nt,1));  temp2 = reshape(temp2,[],1);
    id = find(temp1-temp2 == 0); 
    temp1(id) = []; temp2(id) = [];
    xin = cat(2,temp1,temp2); 
    win = ones(count,1);  %^^
    for i = 1:length(idwr); 
        temp3 = find(xin(:,1) == t_unique(idwr(i)) & xin(:,2) == t_unique(idwc(i))); 
        win(temp3) = 0; 
    end;                  %==
    tempid = diag(ones(nt,1));
    id = find(tempid == 0);
    yin = covest_raw(id);            
    mm = nt*(nt+1)/2;
    tempid = ones(nt,nt); tempid = triu(tempid);
    [idr,idc] = find(tempid' == 1);
    xou = cat(2,t_unique(idc),t_unique(idr)); 
    npoly = [1 1]'; bw = [bw2d bw2d]';    
    [you] = sm2(count,mm,npoly,bw,xin,yin,win,xou);
    i1 = 0;
    for i = 1:nt;
      j = nt-i+1;
      covest(i,[i:nt]) = you(i1+1:i1+j)';  
      covest([i:nt],i) = you(i1+1:i1+j);  
      if (covest(i,i) < 0); covest(i,i) = 0; end;
      i1 = i1+j;
    end;  
    trmat = [1/sqrt(2) 1/sqrt(2); -1/sqrt(2)  1/sqrt(2)]; 
    trxin = trmat*xin';  trxin = trxin';                
    trxou = zeros(nt,2); trxou(:,1) = sqrt(2)*t_unique;    
    npoly = [1 2]'; 
    [tildeG] = sm2(count,nt,npoly,bw,trxin,yin,win,trxou);
    id = find(tildeG < 0); tildeG(id) = 0;
    yin = diag(covest_raw);
    win = ones(nt,1);
    npoly = 1;  
    [Vfcn] = sm1(nt,nt,npoly,bw1d,t_unique,yin,win,t_unique);
    id = find(Vfcn < 0); Vfcn(id) = 0;  
    range_t = max(t_unique)-min(t_unique); 
    i1 = range_t/4;  i2 = 3*range_t/4;   
    id = find(t_unique >= i1 & t_unique <= i2);
    Delta = zeros(nt,1); Delta(1) = t_unique(1);
    Delta(2:nt) = t_unique(2:nt) - t_unique(1:nt-1);
    diff = (Vfcn(id)-tildeG(id)).*Delta(id);   
    evar = 2*sum(diff)/range_t;
    if evar < 0; evar = 0; end; % 

    bt = zeros(nt,1); bt(2:nt) = t_unique(1:nt-1); Delta = t_unique-bt;  
    id = find(Delta == 0); Delta(id) = 0.0001;
    quadwt = diag(Delta);
    tempA = sqrt(quadwt)*covest*sqrt(quadwt);
    [eigvecA,eigval] = eig(tempA);
    temp = 1./sqrt(Delta);
    invW = diag(temp);
    eigvec = invW*eigvecA;
    eval = diag(eigval);
    [seval,index] = sort(eval);
    eval(find(eval<0)) = 0; 
    for i = 1:npcs;
      Eigval(i) = eval(index(nt-i+1));  
      efcn(:,i) = eigvec(:,index(nt-i+1));
    end;
    varprop = Eigval/sum(eval);
    Covfcn = covest(indx_tout,indx_tout);  
    Eigfcn = zeros(mp,npcs);   
    Eigfcn = efcn(indx_tout,:);   
  
    auxefcn = zeros(npcs,m,n); 
    for i = 1:n;
      auxefcn(:,[1:nobs(i)],i) = efcn(indx_tin(i,[1:nobs(i)]),[1:npcs])'; 
    end;   
    btin = zeros(n,m);    
    btin(:,[2:m]) = Tin(:,[1:m-1]);
    Delta = Tin - btin;  Delta(isobs == 0) = 0;      
    temp = zeros(m,n);   Pcs = zeros(n,npcs);
    for k = 1:npcs;
      temp(:,:) = auxefcn(k,:,:);
      Pcs(:,k) =  diag(((Yin - Mufcn).*Delta)*temp);   
    end;
    auxeigval = Eigval(1:npcs);
    auxeigval = auxeigval(:,ones(n,1))';
    range_t = ones(n,1)*(max(t_unique) - min(t_unique)); 
    auxevar = evar*ones(n,1);            
    rho_0 = range_t(1)*evar;
    temp = range_t.*auxevar./nobs;
    sk = temp(:,ones(npcs,1));
    temp1 = (auxeigval + sk);
    id1 = find(temp1 ~= 0);
    id2 = find(temp1 == 0);
    sk(id1) = auxeigval(id1)./temp1(id1);     
    sk(id2) = 1;
    Pcs = sk.*Pcs;
  end;  
      

function [you] = sm1(n,m,npoly,bw,xin,yin,win,xou)
  aa = 1;
  npoly2 = max(npoly,1);
  xin_sort = zeros(n,1);
  you = zeros(m,1);
  xin_sort(1:n-1) = xin(2:n);
  if any( xin_sort([1:n-1]) < xin([1:n-1]) ) == 1;
     [b1,ix1] = sort(xin);
     xin(1:n) = b1;
     yin(1:n) = yin(ix1(1:n));
     win(1:n) = win(ix1(1:n));
  end;
  aux1 = []; aux2 = []; aux3 = [];
  id = win > 0;
  tempx = xin(id);
  tempy = yin(id);
  tempw = win(id);
  for j = 1:m;
      uu = xou(j);
      ul = uu - aa*bw;
      uh = uu + aa*bw;
      you(j) = -99;
      count = 0;
      id = tempx >= ul;
      aux1 = tempx(id);
      aux2 = tempy(id);
      aux3 = tempw(id);
      id =  aux1 <= uh;
      aux1 = (aux1(id)-uu)./bw;
      aux2 = aux2(id);
      aux3 = aux3(id);
      count = length(aux2);
      if count < (2 + npoly2);
          if ( npoly2 == 1);
              if count == 1;  you(j) = aux2(1);  end;
              if count == 2;
                 if aux1(1) >= 0;  you(j) = aux2(1); end;
                 if aux1(2) <= 0;  you(j) = aux2(2); end;
                 if (aux1(1) < 0 & aux1(2) > 0);
                     xh = aux1(2) - aux1(1);
                     if aux2(2) > aux2(1);
                        you(j) = ( -aux1(1) / xh * ( aux2(2) - aux2(1) ) ) + aux2(1);
                     elseif aux2(2) < aux2(1);
                        you(j) = ( aux1(2) / xh * ( aux2(1) - aux2(2) ) ) + aux2(2);
                     else;
                        you(j) = (aux2(1) + aux2(2))/2;
                     end;
                 end;
              end;
          end;
          continue;
      end;
      aux3= aux3.*(1-aux1.*aux1);  
      aux3 = max(0,aux3);
      aux8 = zeros(count,npoly2+1);
      aux8(:,1) = 1;
      aux1 = aux1*bw;
      for k = 1:npoly2;
          aux8(:,k+1) = aux1.^k;
      end;     
     nxmat = npoly2 + 1;
     aux7 = zeros(nxmat,1);  
     swt = sqrt(aux3([1:count]));
     b = aux2([1:count]).*swt;
     SW = swt(:,ones(1,nxmat));
     A = aux8.*SW;
     aux7 = pinv(A)*b;
     you(j) = aux7(1);
end; 

  
function [you] = sm2(n,m,npoly,bw,xin,yin,win,xou)
  spoly = 0;
  for k =1:2;
      npoly(k) = max(npoly(k),1);
      spoly = spoly + npoly(k);
  end;  
  you = zeros(m,1);
  id = win > 0;
  tempx = xin(id,:); tempy = yin(id); tempw = win(id);
  for j = 1:m;
      you(j) = -99;
      count = 0;  
      x = tempx;
      y = tempy;
      w = tempw;
      aa = 1;
      for k = 1:2;
          id1 = []; id2 = [];
          id1 = (x(:,k) >= xou(j,k) - aa*bw(k));
          x = x(id1,:);  y = y(id1,:);     w = w(id1,:); 
          id2 = (x(:,k) <= xou(j,k) + aa*bw(k));
          x = x(id2,:);  y = y(id2,:);     w = w(id2,:); 
      end;
      count = length(y);
      temp = xou(j,:);       
      temp = temp(ones(count,1),:);
      tempd = bw';           
      tempd = tempd(ones(count,1),:);
      x = (x - temp)./tempd;
      if count < (2 + spoly);  continue;   end;
      we = ones(count,1);
      we = we.*prod(1-x.^2,2);
      we = w.*we;
      we = max(zeros(count,1),we);
      temp = bw';
      x = x.*temp(ones(count,1),:);
      xmat = ones(count,1);
      for k = 1:2;
          temp1 = x(:,k);
          temp2 = [];
          for i1 = 1:npoly(k);
              temp2 = cat(2,temp2,temp1.^i1);
          end;
         xmat = cat(2,xmat,temp2);
     end;
     nxmat = size(xmat,2);
     auy = zeros(nxmat,1); 
     swt = sqrt(we([1:count]));
     b = y([1:count]).*swt;
     SW = swt(:,ones(1,nxmat));
     A = xmat.*SW;
     auy = pinv(A)*b;
     you(j) = auy(1);
end; 
