function [Hc,gc,Cc,cc,tCOND] = Full_Condensing_Mex(Q,S,R,A,B,gx,gu,Cx,Cu,c,a,ds0,sizes)

tcond=tic;

[Hc,gc, Cc,cc]=Condensing_mex_new(A,B,Q,S,R,Cx,Cu,ds0,a,c,gx,gu,sizes);

tCOND=toc(tcond)*1e3;

end