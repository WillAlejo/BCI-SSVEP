function [repetidos, z] = Valores_repetidos(vector)

%{
clear;
clc;
close all;


vector =[7.5 6.0 7.5 8.0 7.5 3.0 2.0 7.5 4];
%}

valores_unicos = unique(vector);

N = length(valores_unicos);

vector_str = num2str(vector);

if N > 1
    for i = 1:N
        
        caracter = num2str(valores_unicos(1,i));
        
        S1 = count(vector_str,caracter);
        
        repetidos(1,i) = S1;
        
    end
    
    [x, y] = max(repetidos);
    
    z = valores_unicos(1,y);
else
    repetidos = N;
    z = valores_unicos;
end

end