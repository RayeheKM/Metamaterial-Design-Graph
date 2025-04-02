clc
clear all
% close all

volumeFraction_1 = [];
impedance1_1 = [];
impedance2_1 = [];
impedance3_1 = [];
stiffness1_1 = [];
stiffness2_1 = [];
stiffness3_1 = [];

volumeFraction_11 = [];
impedance1_11 = [];
impedance2_11 = [];
impedance3_11 = [];
stiffness1_11 = [];
stiffness2_11 = [];
stiffness3_11 = [];

volumeFraction_12 = [];
impedance1_12 = [];
impedance2_12 = [];
impedance3_12 = [];
stiffness1_12 = [];
stiffness2_12 = [];
stiffness3_12 = [];

volumeFraction_13 = [];
impedance1_13 = [];
impedance2_13 = [];
impedance3_13 = [];
stiffness1_13 = [];
stiffness2_13 = [];
stiffness3_13 = [];

volumeFraction_2 = [];
impedance1_2 = [];
impedance2_2 = [];
impedance3_2 = [];
stiffness1_2 = [];
stiffness2_2 = [];
stiffness3_2 = [];

volumeFraction_21 = [];
impedance1_21 = [];
impedance2_21 = [];
impedance3_21 = [];
stiffness1_21 = [];
stiffness2_21 = [];
stiffness3_21 = [];


volumeFraction_3 = [];
impedance1_3 = [];
impedance2_3 = [];
impedance3_3 = [];
stiffness1_3 = [];
stiffness2_3 = [];
stiffness3_3 = [];

volumeFraction_31 = [];
impedance1_31 = [];
impedance2_31 = [];
impedance3_31 = [];
stiffness1_31 = [];
stiffness2_31 = [];
stiffness3_31 = [];

volumeFraction_32 = [];
impedance1_32 = [];
impedance2_32 = [];
impedance3_32 = [];
stiffness1_32 = [];
stiffness2_32 = [];
stiffness3_32 = [];

volumeFraction_4 = [];
impedance1_4 = [];
impedance2_4 = [];
impedance3_4 = [];
stiffness1_4 = [];
stiffness2_4 = [];
stiffness3_4 = [];

volumeFraction_41 = [];
impedance1_41 = [];
impedance2_41 = [];
impedance3_41 = [];
stiffness1_41 = [];
stiffness2_41 = [];
stiffness3_41 = [];

volumeFraction_43 = [];
impedance1_43 = [];
impedance2_43 = [];
impedance3_43 = [];
stiffness1_43 = [];
stiffness2_43 = [];
stiffness3_43 = [];

volumeFraction_42 = [];
impedance1_42 = [];
impedance2_42 = [];
impedance3_42 = [];
stiffness1_42 = [];
stiffness2_42 = [];
stiffness3_42 = [];

volumeFraction_S = [];
impedance1_S = [];
impedance2_S = [];
impedance3_S = [];
stiffness1_S = [];
stiffness2_S = [];
stiffness3_S = [];

volumeFraction_N = [];
impedance1_N = [];
impedance2_N = [];
impedance3_N = [];
stiffness1_N = [];
stiffness2_N = [];
stiffness3_N = [];

%This is the one used for 4 together
for i = 1: 100
    % data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\Rod_Data\Rod',num2str(i),'.mat']);
    data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\GraphData2\Rod',num2str(i),'.mat']);
    volumeFraction_1 = [volumeFraction_1; data.VolumeFraction];
    impedance1_1 = [impedance1_1; data.Impedance(1)];
    impedance2_1 = [impedance2_1; data.Impedance(2)];
    impedance3_1 = [impedance3_1; data.Impedance(3)];
    stiffness1_1 = [stiffness1_1; data.Stiffness(1)];
    stiffness2_1 = [stiffness2_1; data.Stiffness(2)];
    stiffness3_1 = [stiffness3_1; data.Stiffness(3)];
end

%This is the one used for new bias
for i = 1: 200
    % data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\Rod_Data\Rod',num2str(i),'.mat']);
    data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\NewDataParametrization\Rod_1x01y01z_',num2str(i),'.mat']);
    volumeFraction_11 = [volumeFraction_11; data.VolumeFraction];
    impedance1_11 = [impedance1_11; data.Impedance(1)];
    impedance2_11 = [impedance2_11; data.Impedance(2)];
    impedance3_11 = [impedance3_11; data.Impedance(3)];
    stiffness1_11 = [stiffness1_11; data.Stiffness(1)];
    stiffness2_11 = [stiffness2_11; data.Stiffness(2)];
    stiffness3_11 = [stiffness3_11; data.Stiffness(3)];
end

% for i = 1: 100
%     data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\GraphData2\Rod',num2str(i),'.mat']);
%     volumeFraction_12 = [volumeFraction_12; data.VolumeFraction];
%     impedance1_12 = [impedance1_12; data.Impedance(1)];
%     impedance2_12 = [impedance2_12; data.Impedance(2)];
%     impedance3_12 = [impedance3_12; data.Impedance(3)];
%     stiffness1_12 = [stiffness1_12; data.Stiffness(1)];
%     stiffness2_12 = [stiffness2_12; data.Stiffness(2)];
%     stiffness3_12 = [stiffness3_12; data.Stiffness(3)];
% end

% This is the one used for rod thickness
for i = 1: 100
    data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\GraphData3\Rod',num2str(i),'.mat']);
    volumeFraction_13 = [volumeFraction_13; data.VolumeFraction];
    impedance1_13 = [impedance1_13; data.Impedance(1)];
    impedance2_13 = [impedance2_13; data.Impedance(2)];
    impedance3_13 = [impedance3_13; data.Impedance(3)];
    stiffness1_13 = [stiffness1_13; data.Stiffness(1)];
    stiffness2_13 = [stiffness2_13; data.Stiffness(2)];
    stiffness3_13 = [stiffness3_13; data.Stiffness(3)];
end

for i = 1: 100
    % data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\Pyramid_Data\Pyramid',num2str(i),'.mat']);
    data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\GraphData2\Pyramid',num2str(i),'.mat']);
    volumeFraction_2 = [volumeFraction_2; data.VolumeFraction];
    impedance1_2 = [impedance1_2; data.Impedance(1)];
    impedance2_2 = [impedance2_2; data.Impedance(2)];
    impedance3_2 = [impedance3_2; data.Impedance(3)];
    stiffness1_2 = [stiffness1_2; data.Stiffness(1)];
    stiffness2_2 = [stiffness2_2; data.Stiffness(2)];
    stiffness3_2 = [stiffness3_2; data.Stiffness(3)];
end

%This is the one used for new parametric study
for i = 1: 200
    % data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\Pyramid_Data\Pyramid',num2str(i),'.mat']);
    data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\NewDataParametrization\Pyramid',num2str(i),'.mat']);
    volumeFraction_21 = [volumeFraction_21; data.VolumeFraction];
    impedance1_21 = [impedance1_21; data.Impedance(1)];
    impedance2_21 = [impedance2_21; data.Impedance(2)];
    impedance3_21 = [impedance3_21; data.Impedance(3)];
    stiffness1_21 = [stiffness1_21; data.Stiffness(1)];
    stiffness2_21 = [stiffness2_21; data.Stiffness(2)];
    stiffness3_21 = [stiffness3_21; data.Stiffness(3)];
end

% for i = 1: 49
%     data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\GraphData2\Pyramid',num2str(i),'.mat']);
%     volumeFraction_2 = [volumeFraction_2; data.VolumeFraction];
%     impedance1_2 = [impedance1_2; data.Impedance(1)];
%     impedance2_2 = [impedance2_2; data.Impedance(2)];
%     impedance3_2 = [impedance3_2; data.Impedance(3)];
%     stiffness1_2 = [stiffness1_2; data.Stiffness(1)];
%     stiffness2_2 = [stiffness2_2; data.Stiffness(2)];
%     stiffness3_2 = [stiffness3_2; data.Stiffness(3)];
% end

% This is the one used for all together
for i = 1: 100
    % data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\Cube_Data\Cube',num2str(i),'.mat']);
    data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\GraphDataZ\Cube',num2str(i),'.mat']);
    volumeFraction_3 = [volumeFraction_3; data.VolumeFraction];
    impedance1_3 = [impedance1_3; data.Impedance(1)];
    impedance2_3 = [impedance2_3; data.Impedance(2)];
    impedance3_3 = [impedance3_3; data.Impedance(3)];
    stiffness1_3 = [stiffness1_3; data.Stiffness(1)];
    stiffness2_3 = [stiffness2_3; data.Stiffness(2)];
    stiffness3_3 = [stiffness3_3; data.Stiffness(3)];
end

% This is the one used for new parametric study
for i = 1: 100
    % data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\Cube_Data\Cube',num2str(i),'.mat']);
    data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\NewDataParametrization\Cube_6_',num2str(i),'.mat']);
    volumeFraction_31 = [volumeFraction_31; data.VolumeFraction];
    impedance1_31 = [impedance1_31; data.Impedance(1)];
    impedance2_31 = [impedance2_31; data.Impedance(2)];
    impedance3_31 = [impedance3_31; data.Impedance(3)];
    stiffness1_31 = [stiffness1_31; data.Stiffness(1)];
    stiffness2_31 = [stiffness2_31; data.Stiffness(2)];
    stiffness3_31 = [stiffness3_31; data.Stiffness(3)];
end

% for i = 1: 100
%     data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\GraphData2\Cube',num2str(i),'.mat']);
%     volumeFraction_32 = [volumeFraction_32; data.VolumeFraction];
%     impedance1_32 = [impedance1_32; data.Impedance(1)];
%     impedance2_32 = [impedance2_32; data.Impedance(2)];
%     impedance3_32 = [impedance3_32; data.Impedance(3)];
%     stiffness1_32 = [stiffness1_32; data.Stiffness(1)];
%     stiffness2_32 = [stiffness2_32; data.Stiffness(2)];
%     stiffness3_32 = [stiffness3_32; data.Stiffness(3)];
% end

% for i = 1: 100
%     data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\Triangle_Data\Triangle',num2str(i),'.mat']);
%     volumeFraction_4 = [volumeFraction_4; data.VolumeFraction];
%     impedance1_4 = [impedance1_4; data.Impedance(1)];
%     impedance2_4 = [impedance2_4; data.Impedance(2)];
%     impedance3_4 = [impedance3_4; data.Impedance(3)];
%     stiffness1_4 = [stiffness1_4; data.Stiffness(1)];
%     stiffness2_4 = [stiffness2_4; data.Stiffness(2)];
%     stiffness3_4 = [stiffness3_4; data.Stiffness(3)];
% end

% This is the one used for all together
for i = 1: 100
    data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\GraphData2\Triangle',num2str(i),'.mat']);
    volumeFraction_4 = [volumeFraction_4; data.VolumeFraction];
    impedance1_4 = [impedance1_4; data.Impedance(1)];
    impedance2_4 = [impedance2_4; data.Impedance(2)];
    impedance3_4 = [impedance3_4; data.Impedance(3)];
    stiffness1_4 = [stiffness1_4; data.Stiffness(1)];
    stiffness2_4 = [stiffness2_4; data.Stiffness(2)];
    stiffness3_4 = [stiffness3_4; data.Stiffness(3)];
end

% This is the one used for new parametric study
for i = 1: 200
    data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\NewDataParametrization\Triangle',num2str(i),'.mat']);
    volumeFraction_41 = [volumeFraction_41; data.VolumeFraction];
    impedance1_41 = [impedance1_41; data.Impedance(1)];
    impedance2_41 = [impedance2_41; data.Impedance(2)];
    impedance3_41 = [impedance3_41; data.Impedance(3)];
    stiffness1_41 = [stiffness1_41; data.Stiffness(1)];
    stiffness2_41 = [stiffness2_41; data.Stiffness(2)];
    stiffness3_41 = [stiffness3_41; data.Stiffness(3)];
end

% for i = 1: 100
%     data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\GraphData2\Triangle',num2str(i),'.mat']);
%     volumeFraction_42 = [volumeFraction_42; data.VolumeFraction];
%     impedance1_42 = [impedance1_42; data.Impedance(1)];
%     impedance2_42 = [impedance2_42; data.Impedance(2)];
%     impedance3_42 = [impedance3_42; data.Impedance(3)];
%     stiffness1_42 = [stiffness1_42; data.Stiffness(1)];
%     stiffness2_42 = [stiffness2_42; data.Stiffness(2)];
%     stiffness3_42 = [stiffness3_42; data.Stiffness(3)];
% end

% plot(volumeFraction, impedance1./10^6, 'o','MarkerFaceColor', 'b')
% hold on
% plot(volumeFraction, impedance2./10^6, 'o','MarkerFaceColor', 'r')
% plot(volumeFraction, impedance3./10^6, 'o','MarkerFaceColor', 'g')
% legend('Z_x', 'Z_y', 'Z_z', 'Location','northwest')
% xlabel('Volume fraction')
% ylabel('Impedance (MRayl)')
% 
% figure
% plot(volumeFraction, stiffness1./10^9, 'o','MarkerFaceColor', 'b')
% hold on
% plot(volumeFraction, stiffness2./10^9, 'o','MarkerFaceColor', 'r')
% plot(volumeFraction, stiffness3./10^9, 'o','MarkerFaceColor', 'g')
% legend('C_{11}', 'C_{22}', 'C_{33}', 'Location','northwest')
% xlabel('Volume fraction')
% ylabel('Stiffness (GPa)')

% figure
% plot3(impedance1_1./10^6./1.5, impedance2_1./10^6./1.5, impedance3_1./10^6./1.5, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none')
% hold on
% plot3(impedance1_4./10^6./1.5, impedance2_4./10^6./1.5, impedance3_4./10^6./1.5, 'o', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'none')
% plot3(impedance1_2./10^6./1.5, impedance2_2./10^6./1.5, impedance3_2./10^6./1.5, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none')
% plot3(impedance1_3./10^6./1.5, impedance2_3./10^6./1.5, impedance3_3./10^6./1.5, 'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'none')
% grid on
% box on
% legend('Rod','Triangle','Pyramid','Cube','Location','Best')
% xlabel('Z_x (MRayl)')
% ylabel('Z_y (MRayl)')
% zlabel('Z_z (MRayl)')
% 
% figure
% plot3(stiffness1_1./10^9./2, stiffness2_1./10^9./2, stiffness3_1./10^9./2, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none')
% hold on
% plot3(stiffness1_4./10^9./2, stiffness2_4./10^9./2, stiffness3_4./10^9./2, 'o', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'none')
% plot3(stiffness1_2./10^9./2, stiffness2_2./10^9./2, stiffness3_2./10^9./2, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none')
% plot3(stiffness1_3./10^9./2, stiffness2_3./10^9./2, stiffness3_3./10^9./2, 'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'none')
% grid on
% box on
% xlim([0 3.5])
% ylim([0 3.5])
% zlim([0 3.5])
% legend('Rod','Triangle','Pyramid','Cube','Location','Best')
% xlabel('C_{11} (GPa)')
% ylabel('C_{22} (GPa)')
% zlabel('C_{33} (GPa)')

% figure
% plot3(impedance1_1./10^6, impedance2_1./10^6, impedance3_1./10^6, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none')
% hold on
% plot3(impedance1_4./10^6, impedance2_4./10^6, impedance3_4./10^6, 'o', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'none')
% plot3(impedance1_2./10^6, impedance2_2./10^6, impedance3_2./10^6, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none')
% plot3(impedance1_3./10^6, impedance2_3./10^6, impedance3_3./10^6, 'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'none')
% grid on
% box on
% legend('Rod','Triangle','Pyramid','Cube','Location','Best')
% xlabel('Z_x (MRayl)')
% ylabel('Z_y (MRayl)')
% zlabel('Z_z (MRayl)')
% 
% figure
% plot3(stiffness1_1./10^9, stiffness2_1./10^9, stiffness3_1./10^9, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none')
% hold on
% plot3(stiffness1_4./10^9, stiffness2_4./10^9, stiffness3_4./10^9, 'o', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'none')
% plot3(stiffness1_2./10^9, stiffness2_2./10^9, stiffness3_2./10^9, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none')
% plot3(stiffness1_3./10^9, stiffness2_3./10^9, stiffness3_3./10^9, 'o', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'none')
% grid on
% box on
% % xlim([0 3.5])
% % ylim([0 3.5])
% % zlim([0 3.5])
% legend('Rod','Triangle','Pyramid','Cube','Location','Best')
% xlabel('C_{11} (GPa)')
% ylabel('C_{22} (GPa)')
% zlabel('C_{33} (GPa)')
% 
% figure
% plot3(impedance1_1./10^6, impedance2_1./10^6, impedance3_1./10^6, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none')
% hold on
% plot3(impedance1_13./10^6, impedance2_13./10^6, impedance3_13./10^6, 'o', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'none')
% grid on
% box on
% legend('Rod','Triangle','Pyramid','Cube','Location','Best')
% xlabel('Z_x (MRayl)')
% ylabel('Z_y (MRayl)')
% zlabel('Z_z (MRayl)')
% 
% figure
% plot3(stiffness1_1./10^9, stiffness2_1./10^9, stiffness3_1./10^9, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none')
% hold on
% plot3(stiffness1_13./10^9, stiffness2_13./10^9, stiffness3_13./10^9, 'o', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'none')
% grid on
% box on
% legend('Rod','Triangle','Pyramid','Cube','Location','Best')
% xlabel('C_{11} (GPa)')
% ylabel('C_{22} (GPa)')
% zlabel('C_{33} (GPa)')

%%
figure
scatter3(impedance1_1./10^6, impedance2_1./10^6, impedance3_1./10^6, 50, 'filled', 'MarkerFaceColor', [0 0 0.5], 'MarkerEdgeColor', [0.5 0.5 1], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
hold on
scatter3(impedance1_4./10^6, impedance2_4./10^6, impedance3_4./10^6, 50, 'filled', 'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', [1 0.75 0.5], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
scatter3(impedance1_2./10^6, impedance2_2./10^6, impedance3_2./10^6, 50, 'filled', 'MarkerFaceColor', [0.5 0 0], 'MarkerEdgeColor', [1 0.5 0.5], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
scatter3(impedance1_3./10^6, impedance2_3./10^6, impedance3_3./10^6, 50, 'filled', 'MarkerFaceColor', [0 0.5 0], 'MarkerEdgeColor', [0.5 1 0.5], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
grid on
box on
legend('Rod','Triangle','Pyramid','Cube','Location','Best')
xlabel('Z_x (MRayl)')
ylabel('Z_y (MRayl)')
zlabel('Z_z (MRayl)')

figure
scatter3(stiffness1_1./10^9, stiffness2_1./10^9, stiffness3_1./10^9, 50, 'filled', 'MarkerFaceColor', [0 0 0.5], 'MarkerEdgeColor', [0.5 0.5 1], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
hold on
scatter3(stiffness1_4./10^9, stiffness2_4./10^9, stiffness3_4./10^9, 50, 'filled', 'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', [1 0.75 0.5], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
scatter3(stiffness1_2./10^9, stiffness2_2./10^9, stiffness3_2./10^9, 50, 'filled', 'MarkerFaceColor', [0.5 0 0], 'MarkerEdgeColor', [1 0.5 0.5], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
scatter3(stiffness1_3./10^9, stiffness2_3./10^9, stiffness3_3./10^9, 50, 'filled', 'MarkerFaceColor', [0 0.5 0], 'MarkerEdgeColor', [0.5 1 0.5], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
grid on
box on
xlim([0 8])
ylim([0 8])
zlim([0 8])
legend('Rod','Triangle','Pyramid','Cube','Location','Best')
xlabel('C_{11} (GPa)')
ylabel('C_{22} (GPa)')
zlabel('C_{33} (GPa)')

figure
scatter3(impedance1_1./10^6, impedance2_1./10^6, impedance3_1./10^6, 50, 'filled', 'MarkerFaceColor', [0 0 0.5], 'MarkerEdgeColor', [0.5 0.5 1], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
hold on
scatter3(impedance1_13./10^6, impedance2_13./10^6, impedance3_13./10^6, 50, 'filled', 'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', [1 0.75 0.5], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
grid on
box on
legend('r = 0.05 mm','r = 0.1 mm','Location','Best')
xlabel('Z_x (MRayl)')
ylabel('Z_y (MRayl)')
zlabel('Z_z (MRayl)')

figure
scatter3(stiffness1_1./10^9, stiffness2_1./10^9, stiffness3_1./10^9, 50, 'filled', 'MarkerFaceColor', [0 0 0.5], 'MarkerEdgeColor', [0.5 0.5 1], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
hold on
scatter3(stiffness1_13./10^9, stiffness2_13./10^9, stiffness3_13./10^9, 50, 'filled', 'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', [1 0.75 0.5], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
grid on
box on
legend('r = 0.05 mm','r = 0.1 mm','Location','Best')
xlabel('C_{11} (GPa)')
ylabel('C_{22} (GPa)')
zlabel('C_{33} (GPa)')



%% New parametric study
% figure
% scatter3(impedance1_1./10^6, impedance2_1./10^6, impedance3_1./10^6, 50, 'filled', 'MarkerFaceColor', [0 0 0.5], 'MarkerEdgeColor', [0.5 0.5 1], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
% hold on
% scatter3(impedance1_11./10^6, impedance2_11./10^6, impedance3_11./10^6, 50, 'filled', 'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', [1 0.75 0.5], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
% grid on
% box on
% legend('r=0.05 mm','r=0.1 mm','Location','Best')
% xlabel('Z_x (MRayl)')
% ylabel('Z_y (MRayl)')
% zlabel('Z_z (MRayl)')
% 
% figure
% scatter3(stiffness1_1./10^9, stiffness2_1./10^9, stiffness3_1./10^9, 50, 'filled', 'MarkerFaceColor', [0 0 0.5], 'MarkerEdgeColor', [0.5 0.5 1], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
% hold on
% scatter3(stiffness1_11./10^9, stiffness2_11./10^9, stiffness3_11./10^9, 50, 'filled', 'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', [1 0.75 0.5], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
% grid on
% box on
% legend('r=0.05 mm','r=0.1 mm','Location','Best') % should be adjusted
% xlabel('C_{11} (GPa)')
% ylabel('C_{22} (GPa)')
% zlabel('C_{33} (GPa)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
scatter3(impedance1_2./10^6, impedance2_2./10^6, impedance3_2./10^6, 50, 'filled', 'MarkerFaceColor', [0 0 0.5], 'MarkerEdgeColor', [0.5 0.5 1], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
hold on
scatter3(impedance1_21./10^6, impedance2_21./10^6, impedance3_21./10^6, 50, 'filled', 'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', [1 0.75 0.5], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
grid on
box on
legend('Domain size = 5 mm','Domain size = 2 mm','Location','Best')
xlabel('Z_x (MRayl)')
ylabel('Z_y (MRayl)')
zlabel('Z_z (MRayl)')

figure
scatter3(stiffness1_2./10^9, stiffness2_2./10^9, stiffness3_2./10^9, 50, 'filled', 'MarkerFaceColor', [0 0 0.5], 'MarkerEdgeColor', [0.5 0.5 1], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
hold on
scatter3(stiffness1_21./10^9, stiffness2_21./10^9, stiffness3_21./10^9, 50, 'filled', 'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', [1 0.75 0.5], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
grid on
box on
xlim([0 8])
ylim([0 8])
zlim([0 8])
legend('Domain size = 5 mm','Domain size = 2 mm','Location','Best')
xlabel('C_{11} (GPa)')
ylabel('C_{22} (GPa)')
zlabel('C_{33} (GPa)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
scatter3(impedance1_3./10^6, impedance2_3./10^6, impedance3_3./10^6, 50, 'filled', 'MarkerFaceColor', [0 0 0.5], 'MarkerEdgeColor', [0.5 0.5 1], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
hold on
scatter3(impedance1_31./10^6, impedance2_31./10^6, impedance3_31./10^6, 50, 'filled', 'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', [1 0.75 0.5], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
grid on
box on
legend('Number of grids = 10','Number of grids = 5','Location','Best')
xlabel('Z_x (MRayl)')
ylabel('Z_y (MRayl)')
zlabel('Z_z (MRayl)')

figure
scatter3(stiffness1_3./10^9, stiffness2_3./10^9, stiffness3_3./10^9, 50, 'filled', 'MarkerFaceColor', [0 0 0.5], 'MarkerEdgeColor', [0.5 0.5 1], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
hold on
scatter3(stiffness1_31./10^9, stiffness2_31./10^9, stiffness3_31./10^9, 50, 'filled', 'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', [1 0.75 0.5], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
grid on
box on
xlim([2,7])
ylim([2,7])
zlim([2,7])

legend('Number of grids = 10','Number of grids = 5','Location','Best')
xlabel('C_{11} (GPa)')
ylabel('C_{22} (GPa)')
zlabel('C_{33} (GPa)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
scatter3(impedance1_4./10^6, impedance2_4./10^6, impedance3_4./10^6, 50, 'filled', 'MarkerFaceColor', [0 0 0.5], 'MarkerEdgeColor', [0.5 0.5 1], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
hold on
scatter3(impedance1_41./10^6, impedance2_41./10^6, impedance3_41./10^6, 50, 'filled', 'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', [1 0.75 0.5], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
grid on
box on
legend('t = 0.05 mm','t = 0.01 mm','Location','Best') % should be adjusted
xlabel('Z_x (MRayl)')
ylabel('Z_y (MRayl)')
zlabel('Z_z (MRayl)')

figure
scatter3(stiffness1_4./10^9, stiffness2_4./10^9, stiffness3_4./10^9, 50, 'filled', 'MarkerFaceColor', [0 0 0.5], 'MarkerEdgeColor', [0.5 0.5 1], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
hold on
scatter3(stiffness1_41./10^9, stiffness2_41./10^9, stiffness3_41./10^9, 50, 'filled', 'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', [1 0.75 0.5], 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5)
grid on
box on
xlim([0,0.5])
ylim([0,0.5])
zlim([0,0.5])
legend('t = 0.05 mm','t = 0.01 mm','Location','Best') % should be adjusted
xlabel('C_{11} (GPa)')
ylabel('C_{22} (GPa)')
zlabel('C_{33} (GPa)')

%%
% close all

% for i = 1: 100
%     data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\Rod_Data\Rod',num2str(i),'.mat']);
%     volumeFraction_1 = [volumeFraction_1; data.VolumeFraction];
%     impedance1_1 = [impedance1_1; data.Impedance(1)];
%     impedance2_1 = [impedance2_1; data.Impedance(2)];
%     impedance3_1 = [impedance3_1; data.Impedance(3)];
%     stiffness1_1 = [stiffness1_1; data.Stiffness(1)];
%     stiffness2_1 = [stiffness2_1; data.Stiffness(2)];
%     stiffness3_1 = [stiffness3_1; data.Stiffness(3)];
% end
% 
% for i = 1: 100
%     data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\Rod_Data\Rod',num2str(i),'.mat']);
%     volumeFraction_12 = [volumeFraction_12; data.VolumeFraction];
%     impedance1_12 = [impedance1_12; data.Impedance(1)];
%     impedance2_12 = [impedance2_12; data.Impedance(2)];
%     impedance3_12 = [impedance3_12; data.Impedance(3)];
%     stiffness1_12 = [stiffness1_12; data.Stiffness(1)];
%     stiffness2_12 = [stiffness2_12; data.Stiffness(2)];
%     stiffness3_12 = [stiffness3_12; data.Stiffness(3)];
% end
% 
% for i = 1: 100
%     data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\Rod_Data\Rod',num2str(i),'.mat']);
%     volumeFraction_13 = [volumeFraction_13; data.VolumeFraction];
%     impedance1_13 = [impedance1_13; data.Impedance(1)];
%     impedance2_13 = [impedance2_13; data.Impedance(2)];
%     impedance3_13 = [impedance3_13; data.Impedance(3)];
%     stiffness1_13 = [stiffness1_13; data.Stiffness(1)];
%     stiffness2_13 = [stiffness2_13; data.Stiffness(2)];
%     stiffness3_13 = [stiffness3_13; data.Stiffness(3)];
% end
% 
% % for i = 1: 100
% %     data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\Rod_Data_Small\Rod',num2str(i),'.mat']);
% %     volumeFraction_S = [volumeFraction_S; data.VolumeFraction];
% %     impedance1_S = [impedance1_S; data.Impedance(1)];
% %     impedance2_S = [impedance2_S; data.Impedance(2)];
% %     impedance3_S = [impedance3_S; data.Impedance(3)];
% %     stiffness1_S = [stiffness1_S; data.Stiffness(1)];
% %     stiffness2_S = [stiffness2_S; data.Stiffness(2)];
% %     stiffness3_S = [stiffness3_S; data.Stiffness(3)];
% % end
% % 
% % for i = 1: 100
% %     data = load(['C:\Users\BrinsonLab\Desktop\Modeling Random Metamaterial\Graph_Data\Rod_Data_New\Rod',num2str(i),'.mat']);
% %     volumeFraction_N = [volumeFraction_N; data.VolumeFraction];
% %     impedance1_N = [impedance1_N; data.Impedance(1)];
% %     impedance2_N = [impedance2_N; data.Impedance(2)];
% %     impedance3_N = [impedance3_N; data.Impedance(3)];
% %     stiffness1_N = [stiffness1_N; data.Stiffness(1)];
% %     stiffness2_N = [stiffness2_N; data.Stiffness(2)];
% %     stiffness3_N = [stiffness3_N; data.Stiffness(3)];
% % end
% 
% figure
% % plot3(impedance1_1./10^6, impedance2_1./10^6, impedance3_1./10^6, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none')
% 
% plot3(impedance1_12./10^6, impedance2_12./10^6, impedance3_12./10^6, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none')
% hold on
% plot3(impedance1_13./10^6, impedance2_13./10^6, impedance3_13./10^6, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none')
% grid on
% legend('Rod, N=5, h=0.05 mm','Rod, N=5, h=0.1 mm','Location','Best')
% xlabel('Z_x (MRayl)')
% ylabel('Z_y (MRayl)')
% zlabel('Z_z (MRayl)')
% 
% figure
% % plot3(stiffness1_1./10^9, stiffness2_1./10^9, stiffness3_1./10^9, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none')
% 
% plot3(stiffness1_12./10^9, stiffness2_12./10^9, stiffness3_12./10^9, 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none')
% hold on
% plot3(impedance1_13./10^6, impedance2_13./10^6, impedance3_13./10^6, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none')
% grid on
% legend('Rod, N=5, h=0.05 mm','Rod, N=5, h=0.1 mm','Location','Best')
% xlabel('C_{11} (GPa)')
% ylabel('C_{22} (GPa)')
% zlabel('C_{33} (GPa)')
