clc;
% Import data files
delimiterIn = ' ';
MembraneMeshZ = importdata('membrane mesh-Z.txt', delimiterIn);
MembraneMeshR = importdata('membrane mesh-R.txt', delimiterIn);
SolutionMeshZ = importdata('solution mesh-Z.txt', delimiterIn);
SolutionMeshR = importdata('solution mesh-R.txt', delimiterIn);

MembraneMeshZ = MembraneMeshZ*1e4; % convert the unit of Z direction of membrane to um 
MembraneMeshR = MembraneMeshR*10; % convert the unit of R direction of membrane to mm
SolutionMeshZ = SolutionMeshZ*10; % convert the unit of Z direction of solution to mm
SolutionMeshR = SolutionMeshR*10; % convert the unit of R direction of solution to mm
% Import concentrations, and converse the unito to Molar
MembraneOx = importdata('membrane_oxidised_speices_concentration@peak.txt', delimiterIn);
MembraneRe = importdata('membrane_reduced_speices_concentration@peak.txt', delimiterIn);
SolutionOx = importdata('solution_oxidised_speices_concentration@peak.txt', delimiterIn);
SolutionRe = importdata('solution_reduced_speices_concentration@peak.txt', delimiterIn);

figure;
subplot(2, 2,1);
surf(MembraneMeshR, MembraneMeshZ, MembraneOx, 'edgecolor', 'none');
hold on;
ax = surf(MembraneMeshR, MembraneMeshZ, MembraneOx, 'edgecolor', 'none');
set(ax, 'Xdata', -1*get(ax, 'Xdata'));
h = colorbar;
xlabel('R / mm');
ylabel('Thickness / \mum');
title('A');
xlabel(h, 'molar');
view(2);
hold off;

subplot(2, 2, 2);
surf(MembraneMeshR, MembraneMeshZ, MembraneRe, 'edgecolor', 'none');
hold on;
ax = surf(MembraneMeshR, MembraneMeshZ, MembraneRe, 'edgecolor', 'none');
set(ax, 'Xdata', -1*get(ax, 'Xdata'));
h = colorbar;
xlabel('R / mm');
ylabel('Thickness / \mum');
title('B');
xlabel(h, 'molar');
view(2);
hold off;

subplot(2, 2, 3);
surf(SolutionMeshR, SolutionMeshZ, SolutionOx, 'edgecolor', 'none');
hold on;
ax = surf(SolutionMeshR, SolutionMeshZ, SolutionOx, 'edgecolor', 'none');
set(ax, 'Xdata', -1*get(ax, 'Xdata'));
h = colorbar;
xlabel('R / mm');
ylabel('Thickness / mm');
title('C');
xlabel(h, 'molar');
view(2);
hold off;

subplot(2, 2, 4);
surf(SolutionMeshR, SolutionMeshZ, SolutionRe, 'edgecolor', 'none');
hold on;
ax = surf(SolutionMeshR, SolutionMeshZ, SolutionRe, 'edgecolor', 'none');
set(ax, 'Xdata', -1*get(ax, 'Xdata'));
h = colorbar;
xlabel('R / mm');
ylabel('Thickness / mm');
title('D');
xlabel(h, 'molar');
view(2);
hold off;

%{
crystalR = crystalmeshZ(1,51);
solutionR = solution1meshX(100,1);
crystalZ = crystalmeshX;
crystalX = crystalmeshZ - crystalR;
crystale = crystalpkce;
crystalh = crystalpkch;

solutionX = solution1meshX(1:100, 1:100);
solutionZ = solution1meshZ(1:100, 1:100);
solutione = solutionpkce(1:100, 1:100);
solutionh = solutionpkch(1:100, 1:100);

figure;
surf(crystalX, crystalZ, crystalh,'edgecolor','none');
hold on;
ax = surf(crystalX, crystalZ, crystalh,'edgecolor','none');
set(ax, 'Xdata',-1*get(ax,'Xdata'));
xlabel('x/cm');
ylabel('z/cm');
axis equal;
%title('Crystal composition at the half of the current peak period')
colorbar;
col1 = caxis;
grid off;
view(2);
saveas(gcf,'reduced species at the half of peak current', 'fig');
saveas(gcf,'reduced species at the half of peak current', 'jpg');

figure;
surf(crystalX, crystalZ, crystale,'edgecolor','none');
hold on;
ax = surf(crystalX, crystalZ, crystale,'edgecolor','none');
set(ax, 'Xdata',-1*get(ax,'Xdata'));
xlabel('x/cm');
ylabel('z/cm');
axis equal;
%title('Crystal composition at the end of the current peak period')
caxis(col1);
colorbar;
grid off;
view(2);
saveas(gcf,'reduced species at the end of peak current', 'fig');
saveas(gcf,'reduced species at the end of peak current', 'jpg');


figure;
surf(solutionX, solutionZ, solutionh,'edgecolor','none');
hold on;
surf(solutionX - 2*solutionR, solutionZ, solutionh, 'edgecolor', 'none');
hold on;
surf(solutionX + 2*solutionR, solutionZ, solutionh, 'edgecolor', 'none');
hold on;
surf(solutionX - 4*solutionR, solutionZ, solutionh, 'edgecolor', 'none');
hold on;
surf(solutionX + 4*solutionR, solutionZ, solutionh, 'edgecolor', 'none');
hold on;
surf(solutionX - 6*solutionR, solutionZ, solutionh, 'edgecolor', 'none');
hold on;
surf(solutionX + 6*solutionR, solutionZ, solutionh, 'edgecolor', 'none');
hold on;
surf(solutionX - 8*solutionR, solutionZ, solutionh, 'edgecolor', 'none');
hold on;
surf(solutionX + 8*solutionR, solutionZ, solutionh, 'edgecolor', 'none');
hold on;
ax = surf(solutionX, solutionZ, solutionh,'edgecolor','none');
set(ax, 'Xdata',-1*get(ax,'Xdata'));
hold on;
ax = surf(solutionX - 2*solutionR, solutionZ, solutionh, 'edgecolor', 'none');
set(ax, 'Xdata', -1*get(ax,'Xdata'));
hold on;
ax = surf(solutionX + 2*solutionR, solutionZ, solutionh, 'edgecolor', 'none');
set(ax, 'Xdata', -1*get(ax,'Xdata'));
hold on;
ax = surf(solutionX - 4*solutionR, solutionZ, solutionh, 'edgecolor', 'none');
set(ax, 'Xdata', -1*get(ax,'Xdata'));
hold on;
ax = surf(solutionX + 4*solutionR, solutionZ, solutionh, 'edgecolor', 'none');
set(ax, 'Xdata', -1*get(ax,'Xdata'));
hold on;
ax = surf(solutionX - 6*solutionR, solutionZ, solutionh, 'edgecolor', 'none');
set(ax, 'Xdata', -1*get(ax,'Xdata'));
hold on;
ax = surf(solutionX + 6*solutionR, solutionZ, solutionh, 'edgecolor', 'none');
set(ax, 'Xdata', -1*get(ax,'Xdata'));
hold on;
ax = surf(solutionX - 8*solutionR, solutionZ, solutionh, 'edgecolor', 'none');
set(ax, 'Xdata', -1*get(ax,'Xdata'));
hold on;
ax = surf(solutionX + 8*solutionR, solutionZ, solutionh, 'edgecolor', 'none');
set(ax, 'Xdata', -1*get(ax,'Xdata'));
xlabel('x/cm');
ylabel('z/cm');
axis equal;
%title('Hydrion distribution at the half of the current peak period')
colorbar;
col2 = caxis;
grid off;
view(2);
saveas(gcf,'hydrion at the half of peak current', 'fig');
saveas(gcf,'hydrion at the half of peak current', 'jpg');

figure;
surf(solutionX, solutionZ, solutione,'edgecolor','none');
hold on;
surf(solutionX - 2*solutionR, solutionZ, solutione, 'edgecolor', 'none');
hold on;
surf(solutionX + 2*solutionR, solutionZ, solutione, 'edgecolor', 'none');
hold on;
surf(solutionX - 4*solutionR, solutionZ, solutione, 'edgecolor', 'none');
hold on;
surf(solutionX + 4*solutionR, solutionZ, solutione, 'edgecolor', 'none');
hold on;
surf(solutionX - 6*solutionR, solutionZ, solutione, 'edgecolor', 'none');
hold on;
surf(solutionX + 6*solutionR, solutionZ, solutione, 'edgecolor', 'none');
hold on;
surf(solutionX - 8*solutionR, solutionZ, solutione, 'edgecolor', 'none');
hold on;
surf(solutionX + 8*solutionR, solutionZ, solutione, 'edgecolor', 'none');
hold on;
ax = surf(solutionX, solutionZ, solutione,'edgecolor','none');
set(ax, 'Xdata',-1*get(ax,'Xdata'));
hold on;
ax = surf(solutionX - 2*solutionR, solutionZ, solutione, 'edgecolor', 'none');
set(ax, 'Xdata', -1*get(ax,'Xdata'));
hold on;
ax = surf(solutionX + 2*solutionR, solutionZ, solutione, 'edgecolor', 'none');
set(ax, 'Xdata', -1*get(ax,'Xdata'));
hold on;
ax = surf(solutionX - 4*solutionR, solutionZ, solutione, 'edgecolor', 'none');
set(ax, 'Xdata', -1*get(ax,'Xdata'));
hold on;
ax = surf(solutionX + 4*solutionR, solutionZ, solutione, 'edgecolor', 'none');
set(ax, 'Xdata', -1*get(ax,'Xdata'));
hold on;
ax = surf(solutionX - 6*solutionR, solutionZ, solutione, 'edgecolor', 'none');
set(ax, 'Xdata', -1*get(ax,'Xdata'));
hold on;
ax = surf(solutionX + 6*solutionR, solutionZ, solutione, 'edgecolor', 'none');
set(ax, 'Xdata', -1*get(ax,'Xdata'));
hold on;
ax = surf(solutionX - 8*solutionR, solutionZ, solutione, 'edgecolor', 'none');
set(ax, 'Xdata', -1*get(ax,'Xdata'));
hold on;
ax = surf(solutionX + 8*solutionR, solutionZ, solutione, 'edgecolor', 'none');
set(ax, 'Xdata', -1*get(ax,'Xdata'));
xlabel('x/cm');
ylabel('z/cm');
axis equal;
%title('Hydrion distribution at the end of the current peak period')
caxis(col2);
colorbar;
grid off;
view(2);
saveas(gcf,'hydrion at the end of peak current', 'fig');
saveas(gcf,'hydrion at the end of peak current', 'jpg');

clear crystalX crystalZ solutione solutionh crystale crystalh ax crystalR
solutionR solutionX solutionZ;
%}
