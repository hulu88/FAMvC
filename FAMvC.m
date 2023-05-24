function [res,Obj,E,A,S,w,obj1,obj3,y] = FAMvC(data, labels, alpha, beta, knn)
%% The performance of different computer may be different

k = max(labels);
n = length(labels);
V = length(data);

for v = 1:V
    S{v} = constructW_PKN(data{v}',knn);
    S{v} = (S{v}+S{v}')/2;
    D10 = diag(1./sqrt(sum(S{v}, 2)));
    S{v} = D10*S{v}*D10;
    [Ft,~,~] = svds(S{v},k);
    Fs{v} = Ft;
end
Qs = Fs;
w = ones(V,1) / V;

T = zeros(n,k);
for idx = 1:v
    Rt = Qs{idx};
    St = S{idx};
    temp = St*Rt;
    w_v = w(idx);
    T = T + w_v*temp;
end
[Uy, ~, Vy] = svds(T,k);
F = Uy*Vy';
F = max(F,0);


% ... update ... %

for iter = 1:200
%     if alpha~=0 && beta 
    B = alpha*ones(V) - diag(alpha*ones(1,V)) + diag(beta*ones(1,V));
    M = B.*(w*w') + 2*diag(w);
    %
    commom_baA = zeros(n);
    for v = 1:V
        baA{v} = alpha*w(v)*S{v};
        special_baA{v} = beta*w(v)*S{v};
        commom_baA = commom_baA + baA{v};
    end
    for v = 1:V
        true_baA{v} = commom_baA - baA{v} + special_baA{v};
        temp_1 = full(w(v)*(2*F*Qs{v}' + true_baA{v}));
        P{v} = temp_1(:);
    end
    right_t = cat(2, P{:})';
    if det(M) == 0
        solution = (pinv(M) * right_t)';
        fprintf('------------,%d',iter)
    else
        solution = (M \ right_t)';
    end
    solution(solution<0) = 0;
    for v = 1:V
        temp_2 = solution(:,v);

        A{v} = zeros(n);
        A{v} = reshape(temp_2, [n n]);
        A{v} = max(A{v}, A{v}');
        A{v} = min(A{v}, S{v});

    end
    
    %% Z
    T = zeros(n,k);
    for idx = 1:v
        Rt = Qs{idx};
        St = A{idx};
        temp_3 = St*Rt;
        w_v = w(idx);
        T = T + w_v*temp_3;
    end
    [Uy, ~, Vy] = svds(T,k);
    F = Uy*Vy';
    F = max(F,0);
    obj = 0;
    %% Qs
    for idx = 1:v
        St = A{idx};
        Rt = St'*F;   % without
        Qs{idx} = Rt*pinv(F'*F);
       
        temp_4 = St - F*Rt';
        E{idx} = S{idx} - A{idx};
        % update weightes
        w(idx) = 0.5/norm(temp_4,'fro');
        obj = obj + norm(temp_4,'fro');
    end

    coef2 = zeros(V);
    for i = 1:V
        for j = i:V
            coef2(i,j) = sum(sum(E{i}.*E{j}));
            coef2(j,i) = coef2(i,j);
        end
    end
    coef2 = coef2.*B;
    obj2 = sum(sum(coef2.*(w*w')));
    obj1(iter)=obj;
    obj3(iter)=obj2;
    Obj(iter) = obj + obj2;
    % convergence checking
    if iter>1
        temp_obj = Obj(iter -1);
    else
        temp_obj = 0;
    end

    if abs(Obj(iter) - temp_obj)/temp_obj <1e-6 %1e-6
        fprintf('------------ ÊÕÁ²\n')
        break;
    end
    
% end
end
[~,y] = max(F,[],2);
res = EvaluationMetrics(labels, y);



