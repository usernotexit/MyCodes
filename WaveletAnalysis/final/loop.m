function [x2, t2, T] = loop(x, t)
    [n_v,~] = size(x);
    T = [];
    t_e = [];
%% making tables
% for each triangle: k1,k2,k3
% T:
%   k1,k2 -> k3,(k4), e_{1,2}
%   k1,k3 -> k2,(k5), e_{1,3}
%   k2,k3 -> k1,(k6), e_{2,3}
%   if k4==0, then edge (k1,k2) is on the boundary
% t_e:
%   [x1,e_{1,2},e_{1,3}]; [x2,e_{2,3},e_{1,2}]; [x3,e_{1,3},e_{2,3}];
%   [e_{2,3},e_{1,3},e_{1,2}]
    new = n_v+1;
    %tris = [t;t(:,[2,3,1]);t(:,[3,1,2])];
    for tris = transpose(t)       
        t_e_i = [];
        for tri = {tris([2,3,1]),tris([3,1,2]),tris([1,2,3])}
            tri = tri{1};
            k1 = tri(1); k2=tri(2); k3=tri(3);
        
            % take k1,k2 -> k3 as example          
            if isempty(T)
                to_insert = false;
            else
                to_insert = Insert(T, [k1,k2]); %to_insert = sum(ismember(T(:,1:2),edge),2)==2;
            end
            if any(to_insert)
                T(to_insert, 4) = k3;% no more than 1
                t_e_i = [t_e_i, T(to_insert, 5)];
            else
                T = [T; [k1,k2, k3,0, new]];
                t_e_i = [t_e_i, new];
                new = new+1;
            end
        end
        x1 = tris(1); x2=tris(2); x3=tris(3);
        e1 = t_e_i(1); e2=t_e_i(2); e3=t_e_i(3);
        t_e = [t_e; [x1,e3,e2]; [x2,e1,e3]; [x3,e2,e1]; [e1,e2,e3]];
    end

%% fill t2
t2 = t_e;

%% calculate x2
    x2 = zeros(new-1, 3);
    %x2(1:n_v,:) = x;
    
    Boun_e = T(:,4)==0; % boundary mask
    x2(T(~Boun_e,5),1:3) = 3/8*(x(T(~Boun_e,1),1:3) + x(T(~Boun_e,2),1:3)) ...
        + 1/8*(x(T(~Boun_e,3),1:3) + x(T(~Boun_e,4),1:3)); % 右边 x和x2 等效
    x2(T(Boun_e,5),1:3) = 1/2*(x(T(Boun_e,1),1:3) + x(T(Boun_e,2),1:3));
    
    n_neighbor = zeros(n_v, 1);
    for Ti=transpose(T)
        n_neighbor(Ti(1)) = n_neighbor(Ti(1))+1;
        n_neighbor(Ti(2)) = n_neighbor(Ti(2))+1;
    end
    Boun_x = false(n_v,1);
    Boun_x(T(Boun_e,1)) = true;    Boun_x(T(Boun_e,2)) = true;
    % fixed the coefficients gamma_x at 2023/12/7
    gamma_x = 3/8+(1/8+cos(2*pi./n_neighbor)/4).^2; %8/5*(3/8+1/4*cos(2*pi./n_neighbor).^2);
    delta_x = (1 - gamma_x)./n_neighbor;
    x_to_e_mat = sparse([T(:,1) T(:,2)], [T(:,5) T(:,5)], ...
        1, n_v, new-1);
        
    x2(Boun_x,:) = x(Boun_x,:);        
    x2(~Boun_x,:) = gamma_x(~Boun_x) .* x(~Boun_x,:) + delta_x(~Boun_x) .* (x_to_e_mat(~Boun_x,:) * x2);
end

%%
function to_insert = Insert(T, edge)
[n,~] = size(T);
to_insert = false(n,1);
to_insert((T(:,1)==edge(1)&T(:,2)==edge(2))|(T(:,1)==edge(2)&T(:,2)==edge(1))) = true;
end