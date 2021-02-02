%% Sierpinski_triangle
%
% Function to compute, display, and save a Sierpinski triangle defined by
% three given points of the 2D or 3D space, at any iteration number / depth level.
%
% Author & support : nicolas.douillet (at) free.fr, 2020-2021.
%
%% Syntax
%
% Sierpinski_triangle;
%
% Sierpinski_triangle(nb_it);
%
% Sierpinski_triangle(nb_it, M1, M2, M3);
%
% Sierpinski_triangle(nb_it, M1, M2, M3, option_display);
%
% [V,T] = Sierpinski_triangle(nb_it, M1, M2, M3, option_display);
%
%% Description
%
% Sierpinski_triangle computes and displays the
% 3-iterations / depth level Sierpinski triangle.
%
% Sierpinski_triangle(nb_it) computes and displays the
% nb_it iterations / depth level Sierpinski triangle.
%
% Sierpinski_triangle(nb_it, M1, M2, M3) uses M1, M2, and M3 as the three
% triangle summits.
%
% Sierpinski_triangle(nb_it, M1, M2, M3, option_display) displays it when
% option_display is set to logical *true/1 (default), and doesn't
% when it is set to  logical false/0.
%
% [V,T] = Sierpinski_triangle(nb_it, M1, M2, M3, option_display) stores the resulting
% vertices coordinates in the array V, and the corresponding triplet indices list in the array T.
%
%% See also
%
% <https://fr.mathworks.com/matlabcentral/fileexchange/73178-n-level-sierpinski-tetrahedron Sierpinski_tetrahedron> |
%
%% Input arguments
%
% - nb_it : positive integer scalar double, the number of iterations / depth level.
%
% - Mi = [Mix Miy Miz], real row vector double, the coordinates of one of the three triangle summits. Size(Mi) = [1,3].
%
% - option_display : either logical *true / false or numeric *1/0.
%
%% Output arguments
%
%        [ |  |  |]
% - V = [Vx Vy Vz], real matrix double, the vertex coordinates. Size(V) = [nb_vertices,3].
%        [ |  |  |]
%
%        [ |  |  |]
% - T = [T1 T2 T3], positive integer matrix double, the triangulation. Size(T) = [nb_triangles,3].
%        [ |  |  |]
%
%% Example #1
% Computes and displays the default 2D Sierpinski equilateral triangle at iteration 3.

Sierpinski_triangle;
axis square;
view(2);

%% Example #2
% Computes and displays the default 2D Sierpinski equilateral triangle at iteration 5.

Sierpinski_triangle(5);
axis square;
view(2);

%% Example #3
% Computes and displays a defined non equilateral 2D Sierpinski triangle at iteration 3.

M1 = [-0.5 0 0];
M2 = [0 0 0];
M3 = [0 2 0];
Sierpinski_triangle(3,M1,M2,M3);
axis equal;
view(2);

%% Example #4
% Computes and displays a random 3D double Sierpinski triangle at iteration 2.

M1 = cat(2,rand(1,2),0);
M2 = cat(2,rand(1,2),0);
M3 = cat(2,rand(1,2),0);
Sierpinski_triangle(2,M1,M2,M3,'double');
axis equal;
view(2);

%% Example #5
% Computes, displays, and saves the 3D (M1,M2,M3) Sierpinski triangle at iteration 4.

M1 = [0 1 1];
M2 = [0 -1 1];
M3 = [0.2 0 0.5];
[V,T] = Sierpinski_triangle(4,M1,M2,M3);
view(28,15);