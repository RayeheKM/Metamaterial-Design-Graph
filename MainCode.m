clc
clear all

% addpath('/shared/Apps/COMSOL61/multiphysics/mli') %On the server
% addpath('/usr/local/comsol61/multiphysics/mli') %On the server
% mphstart(12344); %On the server

import com.comsol.model.*
import com.comsol.model.util.*

for i=1:66

    disp(i)
    
    
    p1 = randi([70,100])/100;
    p2 = randi([0,100])/100;
    Base_Element = 'Rod'; % 'Rod', 'Triangle', 'Pyramid', 'Cube', 'PeriodicCube'
    NameData = [Base_Element,num2str(i)];
    % NameData = 'RodCheck';
    Grid_Number = 6;
    cube_size = 5; %mm
    rho = 1190; %kg/m^3
    E = 3.3; %GPa
    nu = 0.39;

    % system(['python /home/rkm41/GraphData/Main_Graph.py --NameData ',NameData, ' --p1 ', num2str(p1), ' --p2 ', num2str(p2), ' --Grid_Number ', num2str(Grid_Number), ' --cube_size ', num2str(cube_size), ' --Base_Element ', Base_Element]);
    system(['python C:\Users\BrinsonLab\PycharmProjects\model3d\Main_Graph.py --NameData ',NameData, ' --p1 ', num2str(p1), ' --p2 ', num2str(p2), ' --Grid_Number ', num2str(Grid_Number), ' --cube_size ', num2str(cube_size), ' --Base_Element ', Base_Element]);
    
    u_disp = 0.01*cube_size;

    if strcmp(Base_Element, 'Rod')
        disp('Rod selected.')
        [Stiffness, Impedance, VolumeFraction] = Rod(NameData, Grid_Number, cube_size, u_disp, rho, E, nu);
        % disp(['Stiffness: '  num2str(Stiffness)  ' (Pa), Impedance: '  num2str(Impedance)  ' (Rayl), Volume Fraction: '  num2str(VolumeFraction*100)  ' (%)']);
    elseif strcmp(Base_Element, 'Triangle')
        disp('Triangle selected.')
        [Stiffness, Impedance, VolumeFraction] = Triangle(NameData, Grid_Number, cube_size, u_disp, rho, E, nu);
        disp(['Stiffness: '  num2str(Stiffness)  ' (Pa), Impedance: '  num2str(Impedance)  ' (Rayl), Volume Fraction: '  num2str(VolumeFraction*100)  ' (%)']);
    elseif strcmp(Base_Element, 'Pyramid')
        disp('Pyramid selected.')
        [Stiffness, Impedance, VolumeFraction] = Pyramid(NameData, Grid_Number, cube_size, u_disp, rho, E, nu);
        disp(['Stiffness: '  num2str(Stiffness)  ' (Pa), Impedance: '  num2str(Impedance)  ' (Rayl), Volume Fraction: '  num2str(VolumeFraction*100)  ' (%)']);
    elseif strcmp(Base_Element, 'Cube')
        disp('Cube selected.')
        [Stiffness, Impedance, VolumeFraction] = Cube(NameData, Grid_Number, cube_size, u_disp, rho, E, nu);
        disp(['Stiffness: '  num2str(Stiffness)  ' (Pa), Impedance: '  num2str(Impedance)  ' (Rayl), Volume Fraction: '  num2str(VolumeFraction*100)  ' (%)']);
    elseif strcmp(Base_Element, 'PeriodicCube')
        disp('Periodic Cube selected.')
        [Stiffness, Impedance, VolumeFraction] = PeriodicCube(NameData, Grid_Number, cube_size, u_disp, rho, E, nu);
        disp(['Stiffness: '  num2str(Stiffness)  ' (Pa), Impedance: '  num2str(Impedance)  ' (Rayl), Volume Fraction: '  num2str(VolumeFraction*100)  ' (%)']);
    else
        error(['Invalid value for Base_Element: ' Base_Element]);
    end

end