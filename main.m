clear
% clc

load('100leaves.mat')
% if size(Y,1)~=size(X{1},2)
%     for i=1:length(X)
%         X{i} = double(X{i}');
%         X{i} = X{i}./repmat(sqrt(sum(X{i}.^2)),[size(X{i},1) 1]);% Ïû³ýÁ¿¸Ù
%         X{i}=X{i}';
%     end
% else
%     for i=1:length(X)
%         X{i} = double(X{i});
%         X{i} = X{i}./repmat(sqrt(sum(X{i}.^2)),[size(X{i},1) 1]);
%         X{i}=X{i}';
%     end
% end
% for i = 1 :size(X,2)
%     for  j = 1:size(X{1},1)
%         X{i}(j,:) = ( X{i}(j,:) - mean( X{i}(j,:) ) ) / std( X{i}(j,:) );
%     end
% end
data = X;
labels = Y;

% parameter setting
alpha = [1e-5,1e-4,1e-3,1e-2,1e-1,1,10,100,1e3,1e4,1e5];
beta = [1e-5,1e-4,1e-3,1e-2,1e-1,1,10,100,1e3,1e4,1e5];
% for i = 1:length(alpha)
%     for ii = 1:length(beta)

                    [res ,Obj,E,C,A,mu,obj1,obj3,y] = FAMvC(data, labels,10,  1e4,10);

%     end
% end
% 