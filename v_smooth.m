function [v_smooth]=v_smooth(v,s1,s2,nit)

s = 1./v ;

s0 = s ;
s(:) = 0 ;

    hh1=(s1-1)/2; hh2=(s2-1)/2;
    [n1,n2]=size(s0);
    for it=1:nit
        s0 = [repmat(s0(:,1),1,hh2) , s0 , repmat(s0(:,end),1,hh2)] ;
        s0 = [repmat(s0(1,:),hh1,1) ; s0 ; repmat(s0(end,:),hh1,1)] ; 
        
        for i2=hh2+1:hh2+n2
            for i1=hh1+1:hh1+n1
                i11=i1-hh1; i12=i1+hh1 ;
                i21=i2-hh2; i22=i2+hh2 ;
                ss = s0(i11:i12,i21:i22) ;
                s(i1,i2) = sum(sum(ss))/(s1*s2) ;
            end
        end
        s0 = s(hh1+1:hh1+n1,hh2+1:hh2+n2) ;
    end
    s = s0 ;
    
v_smooth = 1./s ;

end




    