clc; clear; close all;
addpath('C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Calculation codes\ndSparse_G4_2021_03_16')
%% INPUTS
% FILES

% fn1 = "C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Backdiagonal_Mises.txt";
% fn2 = "C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\model_Midlle9.txt";
% fn3 = "C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\mises_Diagonal.txt";
% fn4 = "C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Beam_StiffSoft9.txt";
% fn5 = "C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\eams9_softStiff_x.txt";
% fn6 = "C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\new_BeamModel_Middle9.txt";
% fn9 = "C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\New_MiddleVertical_Beam9.txt";
fn1 = "C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\New_Rising_Beam9.txt";
fn2 = "C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\New_middle_Beam9.txt";
fn3 = "C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\New_Falling_Beam9.txt";
fn4 = "C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\New_outer_Beam9.txt";

fni={fn1, fn2, fn3, fn4};

% FLAGS
isGaussianSmooth = true;
smoothing_length_scale = [0.0125 0.0125]; % Give values as fraction of domain size
% smoothing_length_scale = [0.1 0.1];
N_gaussian_smooth = 1;

%% BODY
for jj = 1: length(fni)
    fn = fni{jj};
    % Load the data
    b=importdata(fn,' ',9);
    data=b.data;
    
    % Get the spatial limits of the data
    for i = 1:3
        lims(i,:) = [min(data(:,i)) max(data(:,i))];
    end
    
    % Bin coords into voxel indices
    % Define N_grid based on how many unique coordinates there are in the data
    N_grid = [];
    for i = 1:3
        tol = 1e-8; % Finding unique values within this tolerance
        N_grid(i) = numel(uniquetol(data(:,i),tol));
    end
    
    disp(['N_grid = ' num2str(N_grid)])
    P = [data(:,1) data(:,2) data(:,3)];
    V = data(:,4);
    [x_idx, y_idx, z_idx] = bin_points(P, N_grid);
    coords = [x_idx,y_idx,z_idx];
    
    % Construct a 3-D sparse matrix containing stress values where they exist, and zero where no stress values exist.
    S = ndSparse.build(coords,V); % Download ndSparse from matlab file exchange (author is Matt J)
    S_2D = mean(S,3);
    S_2D = S_2D';
    S_2D = full(S_2D);
    
    % Smooth the data with a Gaussian filter
    gauss_filt_sigma = N_grid(1:2).*smoothing_length_scale;
    if isGaussianSmooth
        for i = 1:N_gaussian_smooth
            S_2D = imgaussfilt(S_2D,gauss_filt_sigma);
        end
    end
    
%     % Plot the result
%     fig = figure();
%     ax = axes(fig);
%     imagesc(ax,lims(1,:),lims(2,:),S_2D)
%     set(ax,'YDir','normal')
%     xlabel(ax,'x')
%     ylabel(ax,'y')
%     colorbar(ax)
%     colormap(ax,'hot')
%     daspect(ax,[1 1 1])

% Plot the result
fig = figure();
ax = axes(fig);
imagesc(ax,lims(1,:),lims(2,:),S_2D./10^6)
set(ax,'YDir','normal')
colorbar_handle = colorbar(ax);
colormap(ax,'hot')

% Add grid lines
grid(ax, 'on')
ax.GridColor = [0.5,1,0.5]; % Light gray color
ax.GridAlpha = 0.5; % Transparency of the grid lines

% Set the grid to be 9 by 9
x_ticks = linspace(lims(1,1), lims(1,2), 10); % 9 intervals create 10 tick marks
y_ticks = linspace(lims(2,1), lims(2,2), 10); % 9 intervals create 10 tick marks
set(ax, 'XTick', x_ticks, 'YTick', y_ticks)

% Remove axis labels, numbers, and ticks
set(ax, 'XTickLabel', [], 'YTickLabel', [], 'XColor', 'none', 'YColor', 'none')

% Customize the colorbar
set(colorbar_handle, 'Ticks', [min(S_2D(:)), max(S_2D(:))], 'TickLabels', {'Min', 'Max'})

% Add text annotations for "Min" and "Max"
% ylabel(colorbar_handle, 'Stress')
% yticks(colorbar_handle, [min(S_2D(:)), max(S_2D(:))])
% yticklabels(colorbar_handle, {'Min', 'Max'})

set(colorbar_handle, 'FontSize', 14)

end
    % set(ax,'colorscale','log')
    
    % This function is written by ChatGPT so hopefully it's right!
    function [x_idx, y_idx, z_idx] = bin_points(P, N_grid)
    % bin_points - Bins input coordinates into voxel indices in a 3D grid.
    %
    % Syntax:
    %   [x_idx, y_idx, z_idx] = bin_points(P, N_grid)
    %
    % Inputs:
    %   P      - N x 3 matrix of points (coordinates).
    %   N_grid - 3 x 1 vector specifying the number of bins along x, y, z.
    %
    % Outputs:
    %   x_idx  - N x 1 vector of voxel indices along the x-dimension.
    %   y_idx  - N x 1 vector of voxel indices along the y-dimension.
    %   z_idx  - N x 1 vector of voxel indices along the z-dimension.

    % Ensure N_grid is a column vector
    N_grid = N_grid(:);

    % Compute the minimum and maximum coordinates in each dimension
    P_min = min(P, [], 1);  % 1 x 3 vector of minimum x, y, z
    P_max = max(P, [], 1);  % 1 x 3 vector of maximum x, y, z

    % Avoid division by zero if P_max equals P_min
    range = P_max - P_min;
    range(range == 0) = 1;

    % Map the points P to grid indices
    x_idx = floor((P(:,1) - P_min(1)) ./ range(1) * N_grid(1)) + 1;
    y_idx = floor((P(:,2) - P_min(2)) ./ range(2) * N_grid(2)) + 1;
    z_idx = floor((P(:,3) - P_min(3)) ./ range(3) * N_grid(3)) + 1;

    % Ensure indices are within the valid range [1, N_grid]
    x_idx = min(max(x_idx, 1), N_grid(1));
    y_idx = min(max(y_idx, 1), N_grid(2));
    z_idx = min(max(z_idx, 1), N_grid(3));
    end

