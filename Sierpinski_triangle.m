function [V, T] = Sierpinski_triangle(nb_it, M1, M2, M3, type, option_display)
%% Sierpinski_triangle : function to compute, display, and save a Sierpinski triangle
% defined by three given points of the 2D or 3D space, at any iteration number / depth level.
%
% Author & support : nicolas.douillet (at) free.fr, 2020-2021.
%
%
% Syntax
%
% Sierpinski_triangle;
% Sierpinski_triangle(nb_it);
% Sierpinski_triangle(nb_it, M1, M2, M3);
% Sierpinski_triangle(nb_it, M1, M2, M3, type);
% Sierpinski_triangle(nb_it, M1, M2, M3, type, option_display);
% [V,T] = Sierpinski_triangle(nb_it, M1, M2, M3, type, option_display);
%
%
% Description
%
% Sierpinski_triangle computes and display the
% 3-iterations / depth level Sierpinski triangle.
%
% Sierpinski_triangle(nb_it) computes and display the
% nb_it iterations / depth level Sierpinski triangle.
%
% Sierpinski_triangle(nb_it, M1, M2, M3) uses M1, M2, and M3 as the three
% triangle summits. 
%
% Sierpinski_triangle(nb_it, M1, M2, M3, type) computes the simple
% Sierpinski triangle by default, or when type is set to 'simple', and
% computes the double Sierpinski fractal when set to double.
%
% Sierpinski_triangle(nb_it, M1, M2, M3, type, option_display) displays it when
% option_display is set to logical *true/1 (default), and doesn't
% when it is set to  logical false/0.
%
% [V,T] = Sierpinski_triangle(nb_it, M1, M2, M3, type, option_display) stores the resulting
% vertices coordinates in the array V, and the corresponding triplet indices list in the array T.
%
%
% Input arguments
%
% - nb_it : positive integer scalar double, the number of iterations / depth level.
%
% - Mi = [Mix Miy Miz], real row vector double, the coordinates of one of the three triangle summits. Size(Mi) = [1,3].
%
% - type : characters string in the set {*'simple','double'}, the type of Sierpinski fractal. Case insensitive.
%
% - option_display : either logical *true / false or numeric *1/0.
%
%
% Output arguments
%
%       [ |  |  |]
% - V = [Vx Vy Vz], real matrix double, the vertex coordinates. Size(V) = [nb_vertices,3].
%       [ |  |  |]
%
%       [ |  |  |]
% - T = [T1 T2 T3], positive integer matrix double, the triangulation. Size(T) = [nb_triangles,3].
%       [ |  |  |]
%
%
% Example #1 : computes and displays the default 2D Sierpinski equilateral triangle at iteration 3.
% 
% Sierpinski_triangle;
% axis square;
% view(2);
%
%
% Example #2 : computes and displays the default 2D Sierpinski equilateral triangle at iteration 5.
% 
% Sierpinski_triangle(5);
% axis square;
% view(2);
%
%
% Example #3 : computes and displays a defined non equilateral 2D Sierpinski triangle at iteration 3.
% 
% M1 = [-0.5 0 0];
% M2 = [0 0 0];
% M3 = [0 2 0];
% Sierpinski_triangle(3,M1,M2,M3);
% axis equal;
% view(2);
%
%
% Example #4 : computes and displays a random 3D double Sierpinski triangle at iteration 2.
% 
% M1 = cat(2,rand(1,2),0);
% M2 = cat(2,rand(1,2),0);
% M3 = cat(2,rand(1,2),0);
% Sierpinski_triangle(2,M1,M2,M3,'double');
% axis equal;
% view(2);
%
%
% Example #5 : computes, displays, and saves the 3D (M1,M2,M3) Sierpinski triangle at iteration 4.
% 
% M1 = [0 1 1];
% M2 = [0 -1 1];
% M3 = [0.2 0 0.5];
% [V,T] = Sierpinski_triangle(4,M1,M2,M3);
% view(28,15);


%% Input parsing
assert(nargin < 6,'Too many input arguments.');

if ~nargin
    nb_it = 3;
    type = 'simple';
    option_display = true;
    M1 = [0.5*sqrt(3) 0 0];
    M2 = [0 1 0];
    M3 = [-0.5*sqrt(3) 0 0];
elseif nargin > 0
    assert(isnumeric(nb_it) && nb_it == floor(nb_it) && nb_it >= 0,'Error : nb_it parameter value must be numeric positive or null integer.');
    if nargin > 1
        if nargin < 4
            error('The three input arguments M1, M2 and M3 must defined together.');
        else % check M1, M2, M3 same size, dim
            assert(isequal(size(M1),size(M2),size(M3)),     'Error : all inputs points must have the same size.');
            assert(isequal(numel(M1),numel(M2),numel(M3)),  'Error : all inputs points must have the same number of elements (2 or 3).');
            assert(isequal(ndims(M1),ndims(M2),ndims(M3),2),'Error : all inputs points must have the same number of dimensions (2).');
            assert(isreal(M1) && isreal(M2) && isreal(M3),  'Error : all inputs points must contain real numbers only.');
            if nargin > 4 % check type char
                assert(ischar(type),'Error : type must be a character string.');                
            else
                type = 'simple';
            end
            if nargin > 5
                assert(islogical(option_display) || isnumeric(option_display),'Error : option_display parameter type must be either logical or numeric.');
            else
                option_display = true;
            end
        end
    else
        M1 = [0.5*sqrt(3) 0 0];
        M2 = [0 1 0];
        M3 = [-0.5*sqrt(3) 0 0];
        type = 'simple';
        option_display = true;
    end
end

warning('on');
if option_display && nb_it > 7
    warning('%s triangles to display ! Make sure your graphic card has enough memory.',num2str(3^(nb_it+1)))
end
warning('off');


%% Body
nb_max_it = 9;
sample_step = 2^(nb_max_it-nb_it);

% Create a meshed triangle
[V1,T_array_1] = sample_triangle(M1',M2',M3',sample_step);

% Middle edge vertices
edg_idx1 = 1 + sample_step/2;
edg_idx3_vect = cumsum(edg_idx1:1+sample_step);
edg_idx3 = edg_idx3_vect(end);
edg_idx2 = edg_idx3 - sample_step/2;

V_array_1 = V1;

% Iterate over nb_it
p = 0;

while p ~= nb_it
    
    if strcmpi(type,'simple')
        
        new_V_array_1 = repmat(V_array_1,[1 1 3]);
        
        for j = 1:size(V_array_1,3) % loop on current nb Sierpinski triangles
            
            new_V_array_1(:,:,3*(j-1)+1) = sample_triangle(V_array_1(1,:,j)',V_array_1(edg_idx1,:,j)',V_array_1(edg_idx2,:,j)',sample_step);               % bottom left triangle (#1)
            new_V_array_1(:,:,3*(j-1)+2) = sample_triangle(V_array_1(1 + sample_step,:,j)',V_array_1(edg_idx3,:,j)',V_array_1(edg_idx1,:,j)',sample_step); % bottom right triangle (#2)
            new_V_array_1(:,:,3*(j-1)+3) = sample_triangle(V_array_1(end,:,j)',V_array_1(edg_idx2,:,j)',V_array_1(edg_idx3,:,j)',sample_step);             % top triangle (#3)
            
        end
        
    elseif strcmpi(type,'double')
        
        new_V_array_1 = repmat(V_array_1,[1 1 4]);
        
        for j = 1:size(V_array_1,3) % loop on current nb Sierpinski triangles
            
            % bottom left triangle (#1)
            new_V_array_1(:,:,4*(j-1)+1) = sample_triangle(V_array_1(1,:,j)',V_array_1(edg_idx1,:,j)',V_array_1(edg_idx2,:,j)',sample_step);
            
            % bottom right triangle (#2)
            new_V_array_1(:,:,4*(j-1)+2) = sample_triangle(V_array_1(edg_idx1,:,j)',V_array_1(1 + sample_step,:,j)',V_array_1(edg_idx3,:,j)',sample_step);
            
            % top triangle (#3)
            new_V_array_1(:,:,4*(j-1)+3) = sample_triangle(V_array_1(edg_idx2,:,j)',V_array_1(edg_idx3,:,j)',V_array_1(end,:,j)',sample_step);
            
            % centre triangle (#4)
            new_V_array_1(:,:,4*(j-1)+4) = sample_triangle(new_V_array_1(edg_idx3,:,4*(j-1)+1)',new_V_array_1(edg_idx1,:,4*(j-1)+2)',new_V_array_1(edg_idx2,:,4*(j-1)+3)',sample_step);
            
        end
        
    end
    
    V_array_1 = new_V_array_1;    
    p = p+1;
    
end

V = V_array_1(:,:,1);
T = T_array_1(:,:,1);

for k = 1:size(V_array_1,3)
    
    T = cat(1,T,T_array_1(:,:,1)+size(V,1));
    V = cat(1,V,V_array_1(:,:,k));  
    
end


%% Display
if option_display
    
    figure;
    set(gcf,'Color',[0 0 0]);
    t = trisurf(T,V(:,1),V(:,2),V(:,3)); shading interp, hold on;
    t.EdgeColor = 'g';
    set(gca,'Color',[0 0 0],'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1]);
    colormap([0 1 0]);
    axis equal, axis tight;
    camlight right;
    camlight head;
    
end


end % Sierpinski_triangle


%% sample_triangle subfunction
function [V, T] = sample_triangle(V1, V2, V3, nbstep)


% Create sampling grid
global Ndim;

Ndim = size(V1,1);

% (V1V2, V1V3) base
u = (V2 - V1);
v = (V3 - V1);

V = zeros(sum(1:nbstep+1),Ndim);

nu = u / norm(u);
nv = v / norm(v);
stepu = norm(u) / nbstep;
stepv = norm(v) / nbstep;
k = 1;

% Sampling & vertices generation
for m = 0:nbstep
    
    for n = 0:nbstep
        
        if m+n <= nbstep % in (V1,V2,V3) triangle conditions ; indices # nb segments
            
            % translation vector
            tv = m*stepu*nu + n*stepv*nv;
            V(k,:) = (V1 + tv)';
            k = k+1;
            
        end
        
    end
    
end


% Index triplets list construction
T = zeros(nbstep^2,3);
row_length = 1 + nbstep;
cum_row_length = row_length;
row_idx = 1;
p = 1;

while p <= nbstep^2 && row_length > 1
    
     i = p;
    
    if p < 2 % "right" triangle serie only
        
        while i < cum_row_length
            
            T(row_idx,:) = [i i+1 i+row_length];
            row_idx = row_idx + 1;
            i = i +1;
            
        end
        
        row_length = row_length - 1;
        cum_row_length = cum_row_length + row_length;
        p = p + row_length+1;
        
    else
        
        % Since p >= 2
        while i < cum_row_length % both triangle series
            
            T(row_idx,:) = [i i+1 i+row_length];
            row_idx = row_idx + 1;            
            T(row_idx,:) = [i i+1 i-row_length]; % + upside-down triangles serie
            row_idx = row_idx + 1;            
            i = i +1;
            
        end
        
        row_length = row_length - 1;
        cum_row_length = cum_row_length + row_length;
        p = p + row_length+1;
        
    end
    
end

T = sort(T,2);
T = unique(T,'rows','stable');


end % sample_triangle