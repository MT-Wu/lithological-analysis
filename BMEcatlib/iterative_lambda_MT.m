function [Lfit]=iterative_lambda_MT(Lfit,ptheor,pbeta);
%% iterative�G���NLambda���즬��

delta=[];
PLfit=exp(Lfit);

if ptheor>pbeta,
    delta=log(((ptheor-pbeta)*(1-PLfit))/((1-ptheor+pbeta)*PLfit));
    
else
    delta=log(((ptheor+pbeta)*(1-PLfit))/((1+ptheor-pbeta)*PLfit));
    
end

Lfit=Lfit+delta;


end