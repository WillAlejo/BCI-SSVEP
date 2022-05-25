function [prec,error,exactitud]=eficacia(X,Y,k)
vp=0;
fn=0;
fp=0;
vn=0;
for i=1:length(X)
    if (X(i)==k) && (Y(i)==k)
        vp=vp+1;
    end
    if (X(i)==k) && (Y(i)~=k)
        fn=fn+1; 
    end
    if (X(i)~=k) && (Y(i)==k)
        fp=fp+1; 
    end
    if (X(i)~=k) && (Y(i)~=k)
        vn=vn+1; 
    end
end
prec=vp/(vp+fp);
error=(fp+fn)/(vp+vn+fn+fp);
exactitud=(vp+vn)/(vp+fp+vn+fn);
end