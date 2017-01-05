function [Pfit,niter]=iterativerescaling_MT(Pbiv,ndim,nc,Pbeta,tol);

% iterativerescaling     - iterative rescaling algorithm for maximum entropy fitting
%                          (December 1, 2003)
%
% Iterative rescaling algorithm for the fitting of a multidimensional
% probability table. This function is used by maxentropytable.m
%
% SYNTAX :
%
% [Pfit]=iterativerescaling(Pbiv,ndim,nc,tol);
%
% INPUT :
%
% Pbiv     ndim by ndim     square cell array where each upper diagonal cell (i,j) is
%                           a symmetric nc by nc matrix of bivariate probabilities
%                           for categories considered at locations separated by a
%                           given distance, where ndim is the number of locations.
% ndim     scalar           dimensionality of Pbiv. ��nmax+1(������m�I���Сϭ�l�������I�^
% nc       scalar           number of elements along each dimension of Pbiv, that also
%                           corresponds to the number of categories for the variable.
% tol      scalar           value of the stopping criterion for the iterative scaling,
%                           that corresponds to the maximum of the absolute differences
%                           between joint probabilities estimates for two successive
%                           iterations. Default value is equal to 1e-3.
%
% OUTPUT :
%
% Pfit    nc by ... by nc   ndim-dimensional table with nc elements along each
%                           dimension, that contains joint probability estimates
%                           after iterative rescaling is completed.

test=0;
Pfit=createones(ndim,nc)/(nc^ndim); %�ؤ@��table�������O1
Lfit=createzeros(ndim,nc);  %�ؤ@��table�������O0
Pfitold=Pfit;

niter=0;
while test==0,
  niter=niter+1;
  for i=1:ndim, %������m�I�ƶq
    for j=i+1:ndim,
      Pmar=sumoverallexcepttwo(Pfit,i,j);
      Pmar=squeeze(Pmar);
      index=ones(ndim,1);
      for k=1:nc, %�������O1?3
        for l=1:nc, %�������O1?3
          index(i)=k;
          index(j)=l;
          if ~isnan(Pbiv{i,j}(k,l))
            [Pfit,niter]=multiplysubtable_MT(Pfit,Lfit,[i j],[index(i) index(j)],Pbiv{i,j}(k,l),Pmar(k,l),Pbeta{i,j}(k,l),tol);
          end;
        end;
      end;
    end;
  end;
  if max(abs(Pfit(:)-Pfitold(:)))<tol,  %?��Lamda�N���ħ����ҥH�o�q������
    test=1;
  else
    Pfitold=Pfit;
  end;
end;