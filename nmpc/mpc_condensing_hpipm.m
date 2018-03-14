function [mem,tCOND] = mpc_condensing_hpipm(mem,settings)
    tcond=tic;
    condensing_hpipm(mem,settings);
    tCOND=toc(tcond)*1e3;
   
    mem.Hc= tril(mem.Hc,-1)'+mem.Hc;
    
    mem.Hc=flipud(mem.Hc);
    mem.Hc=fliplr(mem.Hc);
    
    mem.gc=flipud(mem.gc);
    
    mem.Cc=tril(mem.Cc);
    
    mem.lcc=flipud(mem.lcc);
    
    mem.ucc=flipud(mem.ucc);
    
end

