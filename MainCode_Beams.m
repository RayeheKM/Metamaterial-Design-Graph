clc
clear all

% addpath('/shared/Apps/COMSOL61/multiphysics/mli') %On the server
% addpath('/usr/local/comsol61/multiphysics/mli') %On the server
% mphstart(20001); %On the server

import com.comsol.model.*
import com.comsol.model.util.*

rownumber=9;
columnnumber=9;

Base_Element = 'Rod'; % 'Rod', 'Triangle', 'Pyramid', 'Cube', 'PeriodicCube'

% NameData = 'RodCheck';
Grid_Number = 6;
cube_size = 1; %mm
rho = 1190; %kg/m^3
E = 3.3; %GPa
nu = 0.39;
u_disp = 0.01*cube_size;

% randi([60,100])/100
p1 = 0.6.*ones(rownumber,columnnumber);
p2= ones(rownumber,columnnumber);
bx = 0.1.*ones(rownumber,columnnumber);
by = 0.1.*ones(rownumber,columnnumber);
bz = 0.1.*ones(rownumber,columnnumber);

% bx(1:3, :) = 1;       % First rows
% bx(4:6, :) = 0.1;     % Middle rows
% bx(7:9, :) = 1;       % Last rows
% 
% 
% by(1:3, :) = 0.1;     % First rows
% by(4:6, :) = 1;       % Middle rows
% by(7:9, :) = 0.1;     % Last rows

% by(1:3, :) = 1;       % First rows
% by(4:6, :) = 0.1;     % Middle rows
% by(7:9, :) = 1;       % Last rows
% 
% 
% bx(1:3, :) = 0.1;     % First rows
% bx(4:6, :) = 1;       % Middle rows
% bx(7:9, :) = 0.1;     % Last rows

% 
% bx(1:3, 7:9) = 1;       % First rows
% bx(4:6, 4:6) = 1;     % Middle rows
% bx(7:9, 1:3) = 1;       % Last rows
% 
% 
% by(1:3, 1:6) = 1;     % First rows
% by(4:6, [1:3,7:9]) = 1;       % Middle rows
% by(7:9, 4:9) = 1;     % Last rows

% bx(4:6, 4:6) = 1;     % Middle rows
% 
% by([1:3,7:9], [1:3,7:9]) = 1;     % First rows
% by(4:6, [1:3,7:9]) = 1;       % Middle rows


% bx(1:3, 1:3) = 1;       % First rows
% bx(4:6, 4:6) = 1;       % Middle rows
% bx(7:9, 7:9) = 1;       % Last rows
% 
% 
% by(1:3, 4:9) = 1;             % First rows
% by(4:6, [1:3,7:9]) = 1;       % Middle rows
% by(7:9, 1:6) = 1;             % Last rows

% bx(1:3, 9) = 1;
% bx(1:4, 8) = 1;
% bx(1:5, 7) = 1;
% bx(2:6, 6) = 1;
% bx(3:7, 5) = 1;
% bx(4:8, 4) = 1;
% bx(5:9, 3) = 1;
% bx(6:9, 2) = 1;
% bx(7:9, 1) = 1;

% bx(1:3, 1) = 1;
% bx(1:4, 2) = 1;
% bx(1:5, 3) = 1;
% bx(2:6, 4) = 1;
% bx(3:7, 5) = 1;
% bx(4:8, 6) = 1;
% bx(5:9, 7) = 1;
% bx(6:9, 8) = 1;
% bx(7:9, 9) = 1;

% bx(4:6, 1:9) = 1;

bx([1:3,7:9], 1:9) = 1;

for i=1:rownumber
    for j=1:columnnumber

        disp([num2str(i), ',', num2str(j)])
        NameData = [Base_Element,num2str((i-1)*rownumber+j)];

        % system(['python3 /home/rkm41/GraphDataZ/Main_Graph.py --NameData ',NameData, ' --p1 ', num2str(p1), ' --p2 ', num2str(p2), ' --bx ', num2str(bx), ' --by ', num2str(by), ' --bz ', num2str(bz), ' --Grid_Number ', num2str(Grid_Number), ' --cube_size ', num2str(cube_size), ' --Base_Element ', Base_Element]);
        system(['python C:\Users\BrinsonLab\PycharmProjects\model3d\Main_Graph.py --NameData ',NameData, ' --p1 ', num2str(p1(i,j)), ' --p2 ', num2str(p2(i,j)), ' --bx ', num2str(bx(i,j)), ' --by ', num2str(by(i,j)), ' --bz ', num2str(bz(i,j)), ' --Grid_Number ', num2str(Grid_Number), ' --cube_size ', num2str(cube_size), ' --Base_Element ', Base_Element]);

    end
end

Beams(rownumber, columnnumber, cube_size, u_disp, rho, E, nu);