function [sdmax,sdmin,tdmax,tdmin]=maxmind(cMS,date1)
%% calculate spatial-temporal max and mix distance
%calculate spatial max distance
c1=cMS;
c1=combinedupli(c1,c1(:,1));
D=kron(c1,ones(length(c1),1))-kron(ones(length(c1),1),c1);
D=sqrt(sum((D.^2)')');
D=reshape(D,length(c1),length(c1));
index=find(D~=0);
D=D(index);
D=reshape(D,length(c1)-1,length(c1));
D=max(D);
sdmax=max(D');
%calculate spatial min distance
c1=combinedupli(c1,c1(:,1));
D=kron(c1,ones(length(c1),1))-kron(ones(length(c1),1),c1);
D=sqrt(sum((D.^2)')');
D=reshape(D,length(c1),length(c1));
index=find(D~=0);
D=D(index);
D=reshape(D,length(c1)-1,length(c1));
D=min(D);
sdmin=max(D');
%calculate temporal max distance
% % % % if ~isempty(ck)
% % % % t=sort([date1;ck(:,3)]);
% % % % else
 t=sort(date1); 
% % % % end
tdmax=t(end,1);
%calculate temporal min distance
tdmin=t(1,1); 