function [ OBJ ] = GetObj( z,xN,od,y,yN,Q,QN,input,sizes )

    N=sizes.N;   
        
    obj_vec=input.data.obj_vec;
    
    for i=1:N
        zi = z(i,:)';
        refi = y(i,:)';
        
        if ~isempty(od)
            parai = od(i,:)';
        else
            parai = 0;
        end
        
        obj_vec{i}=obji_vec_fun('obji_vec_fun',zi,refi,parai,Q);        
    end
    
    i=N+1;
    zN = xN;
    refN = yN;
    if ~isempty(od)
         paraN = od(i,:)';
     else
         paraN = 0;
    end
    
    obj_vec{i}=objN_vec_fun('objN_vec_fun',zN,refN,paraN,QN);
    
    obj_vec=cell2mat(obj_vec);
    
    OBJ=0.5*norm(full(obj_vec),2)^2;
end

