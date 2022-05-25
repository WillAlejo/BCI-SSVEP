
function [n,v_minimo,v_maximo]=normalizacion(x)
    v_minimo=min(x);
    r=x-v_minimo;
    v_maximo=max(r);
    n=r/v_maximo;
end