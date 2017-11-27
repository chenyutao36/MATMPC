function [Hc,gc,Cc,cc,tcond] = Full_Condensing(Q,S,R,A,B,gx,gu,Cx,Cu,c,a,ds0,input,sizes)

N=sizes.N;
nx=sizes.nx;
nu=sizes.nu;
nc=sizes.nc*2;
ncN=sizes.ncN*2;

tic;

%% G aknd L

G=input.data.G;
L=input.data.L;

L(:,1)=ds0;

for i=1:N   
    G{i,i}=B{i};                 
    for k=i+1:N
          G{k,i}=A{k}*G{k-1,i};
    end   
    L(:,i+1)=A{i}*L(:,i)+a(:,i);
end


%% Hc

Hc=input.data.Hc;

for i=1:N
        
    W=Q{N+1}*G{N,i};
    
    for k=N:-1:i+1
                
        tmp=G{k-1,i};
                             
        Hc{k,i}=S{k}'*tmp+B{k}'*W;
        Hc{i,k}=Hc{k,i}';
        
        W=Q{k}*tmp+A{k}'*W;
             
    end
    
    Hc{i,i} = R{i}+B{i}'*W;
end

%% gc

gc=input.data.gc;

w=gx(:,N+1)+Q{N+1}*L(:,N+1);
for i=N:-1:2
    gc(:,i)=gu(:,i)+S{i}'*L(:,i)+B{i}'*w;
    w=gx(:,i)+Q{i}*L(:,i)+A{i}'*w;
end
gc(:,1)=gu(:,1)+S{1}'*L(:,1)+B{1}'*w;

%% Cc and cc

Cc=input.data.Cc;
cc=input.data.cc;

for i=1:N
    Cc{i,i}=Cu{i};
    cc((i-1)*nc+1:i*nc)=c((i-1)*nc+1:i*nc)+Cx{i}*L(:,i);
    for k=i+1:N
        Cc{k,i}=Cx{k}*G{k-1,i};        
    end
    Cc{N+1,i}=Cx{N+1}*G{N,i};
end

cc(N*nc+1:N*nc+ncN)=c(N*nc+1:N*nc+ncN)+Cx{N+1}*L(:,N+1);
%%

gc=reshape(gc,[N*nu,1]);
Hc=cell2mat(Hc);
Cc=cell2mat(Cc);

%%
tcond=toc*1e3;

end