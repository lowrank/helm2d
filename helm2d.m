function [x] = helm2d()
%HELM2D computes solution of Helmholtz equation in 2D domain [0,1]x[0,1]
% boundary condition is taken either absorbing or Dirichlet.
% Finite element method based.

%% paramaters
k = 15; % freq
sig = 0; % absorption
margin = 2;

opt = struct('deg', 3, 'qdeg', 6, 'min_area', 1e-3, 'edge', [0 1 1 0; 0 0 1 1], 'hull',...
    [0 0 1 0 1 1 0 1 ...
    -margin -margin 1+margin -margin ...
    1+margin 1+margin -margin 1+margin]',...
    'idx', [0 1 1 2 2 3 3 0 4 5 5 6 6 7 7 4]');
V = femm(opt);


S = V.build('s', 1);
M = V.build('m', 1);

f =@(x)(exp(-((x(1,:) - 0.5).^2 + (x(2,:) - 0.5).^2)/100)); 
load_f = f(V.quad2d);
L = V.build('l', load_f);


bcs = BC('dirichlet');
bcs.set_constraint('x-1.25');
bcs.set_constraint('x+0.25');
bcs.set_constraint('y-1.25');
bcs.set_constraint('y+0.25');

[e1, e2, e3, e4] = bcs.get_boundary(V.space.edges, V.space.nodes, 4);

E = unique([e1 e2 e3 e4]);
N = size(S, 1);
dofs = setdiff(1:N, E);

A = -S+ (k^2 + sqrt(-1)*k * sig) * M;
x = zeros(N, 1);
x(E) = 0;
x(dofs)  = A(dofs, dofs)\L(dofs);


trisurf(V.space.elems(1:3,:)', V.space.nodes(1,:), V.space.nodes(2,:), real(x), 'edgeColor', 'none');
colormap jet; view(2); shading interp;
end

