%% The usual stuff
clear;
clc;
close all;

%% Add modules
mrstModule add ad-blackoil ad-core ad-props 

%% Grid 
nx = 10;
ny = 10;
nz = 20;
x_size = 0.2;
y_size = 0.2;
z_size = 1;
G = cartGrid([nx, ny, nz], [x_size, y_size, z_size]);
G = computeGeometry(G);

%% Set up rock properties
rock = makeRock(G, 100*milli*darcy, .1);

%% Set up fluid
fluid = initSimpleADIFluid('phases', 'WO',...
                                'c', [1e-11,1e-11]/psia,...
                                'n', [2,2],...
                                'mu',[1, 1]*centi*poise,...
                                'rho',[1000, 0.01]);
                            
%% Setting different regions with different rel perms and cap press
region1 = find(G.cells.centroids(:,1)<=x_size/3);
region2 = find(G.cells.centroids(:,1)<=2*x_size/3 & ...
               G.cells.centroids(:,1)>=x_size/3);
region3 = find(G.cells.centroids(:,1)>=2*x_size/3 & ...
               G.cells.centroids(:,1)<=x_size);
rock.regions.saturation = ones(G.cells.num,1);
rock.regions.saturation(region2) = 2;
rock.regions.saturation(region3) = 3;

%% Setting relative permeabilities
krw1 = @(sw)sw.^1;
kro1 = @(so)so.^1;
krw2 = @(sw)sw.^2;
kro2 = @(so)so.^2;
krw3 = @(sw)sw.^3;
kro3 = @(so)so.^3;
Pe1 = 1*kilo;
Pe2 = 10*kilo;
Pe3 = 100*kilo;
pc1 = @(sw)Pe1*sw.^-0.5;
pc2 = @(sw)Pe2*sw.^-0.5;
pc3 = @(sw)Pe3*sw.^-0.5;
fluid.krW = {krw1, krw2, krw3};
fluid.krO = {kro1, kro2, kro3};
fluid.pcOW = {pc1, pc2, pc3};

%% Set up model
gravity reset off;
model = TwoPhaseOilWaterModel(G, rock, fluid);
model = model.validateModel();

%% Boundary conditions
bc = [];
bc = fluxside(bc, G, 'top', 1e-7, 'sat', [0, 1]);
bc = pside(bc, G, 'bottom', 0, 'sat', [0, 1]);

%% Initializing state 
W = [];
state = initResSol(model.G, 1000*psia, [1, 0]);
state.wellSol = initWellSolAD(W, model, state);

%% Solver
solver = NonLinearSolver();

%% Figure
fig1 = figure('Position',[100,100,900,800]);
fig1.Color = 'w';

%% Time loop
dt = 0.01*day;
tmax = 1000*dt;
t = 0;
while t<=tmax
    
    disp(['Time = ',num2str(t/day), ' days'])
    state = solver.solveTimestep(state, dt, model, 'bc', bc);
    
    press = model.getProps(state, 'PhasePressures');
    
    figure(fig1)
    subplot(1,3,1);
    colormap(parula);
    p = plotCellData(G,state.s(:,2));
    p.EdgeAlpha = 0.3;
    colorbar;
    caxis([0,1]);
    set(gca,'FontSize',16);
%     axis equal;
    view(-10, 3);
    xlabel('x');
    ylabel('y');
    title('CO2 Saturation');
%     drawnow;
    
    subplot(1,3,2);
    colormap(parula);
    p = plotCellData(G,press{1}/kilo);
    p.EdgeAlpha = 0.3;
    colorbar;
    set(gca,'FontSize',16);
%     axis equal;
    view(-10, 3);
    xlabel('x');
    ylabel('y');
    title('Water Pressure [kPa]');
    
    subplot(1,3,3);
    colormap(parula);
    p = plotCellData(G,press{2}/kilo);
    p.EdgeAlpha = 0.3;
    colorbar;
    set(gca,'FontSize',16);
%     axis equal;
    view(-10, 3);
    xlabel('x');
    ylabel('y');
    title('CO2 Pressure [kPa]');
    
    figure(fig1)
    drawnow;
    t = t+dt;
    
end