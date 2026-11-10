%% ==============================================================================
% BIOPCL BIOPRINTED MANDIBLE -- MATHEMATICAL MODEL & BIOMECHANICAL SIMULATION
% Anatomical Construct Analysis | Gibson-Ashby Porosity | FEA Nodal Mapping
% ==============================================================================
%
% REQUIRED FINAL BIOPCL MATERIAL SPECIFICATION:
%   Density                       : 1200 kg/m^3 (1.2 g/cm^3)
%   Young's Modulus               : 350 MPa (0.35 GPa)
%   Poisson's Ratio               : 0.30
%   Tensile Yield Strength        : 15 MPa
%   Tensile Ultimate Strength     : 23 MPa
%   Compressive Yield Strength    : 15 MPa
%   Compressive Ultimate Strength : 23 MPa
%
% MATHEMATICAL FRAMEWORK:
%   - Stress field invariance under material change (force-controlled FEA)
%   - Deformation scaling: u_BioPCL = u_ABS * (E_ABS / E_BioPCL)
%   - Gibson-Ashby cellular solids: E(P) = E_0*(1-P)^n
%   - Frost mechanostat using Voigt composite modulus (scaffold + bone)
% ==============================================================================

clear; clc; close all;

fprintf('===========================================================================\n');
fprintf('   BIOPCL BIOPRINTED MANDIBLE -- COMPUTATIONAL BIOMECHANICAL ANALYSIS      \n');
fprintf('===========================================================================\n\n');

%% 1. MATERIAL SPECIFICATIONS (BIOPCL)
BioPCL.Density      = 1200.0;  % kg/m^3
BioPCL.E            = 350.0;   % MPa
BioPCL.nu           = 0.30;
BioPCL.SigmaYield_T = 15.0;    % MPa (tensile yield)
BioPCL.SigmaUlt_T   = 23.0;    % MPa (tensile ultimate)
BioPCL.SigmaYield_C = 15.0;    % MPa (compressive yield)
BioPCL.SigmaUlt_C   = 23.0;    % MPa (compressive ultimate)
BioPCL.G            = BioPCL.E / (2 * (1 + BioPCL.nu)); % Shear modulus (MPa)

E_ABS = 2000.0; % Baseline ABS modulus from original ANSYS simulation (MPa)
E_bone = 15000.0; % Representative cortical bone modulus for Frost mechanostat (MPa)
BVF = 0.30; % Bone volume fraction at scaffold-tissue interface (Cipitria 2012)

%% 2. IMPORT FEA NODAL DATA (ANSYS SIMULATION EXPORT)
% The data file is TAB-separated with columns:
%   Column 1: Node Number (integer)
%   Column 2: Equivalent (von-Mises) Stress in Pascals (Pa)
dataDir = 'data';
stressFile = fullfile(dataDir, 'tabular_stress_data.csv');

stresses_Pa = [];
if exist(stressFile, 'file')
    % Explicitly set delimiter to tab for robustness
    fid = fopen(stressFile, 'r');
    headerLine = fgetl(fid); % skip header
    dataMat = textscan(fid, '%f %f', 'Delimiter', '\t');
    fclose(fid);
    node_ids = dataMat{1};
    stresses_Pa = dataMat{2};
else
    warning('FEA file %s not found. Using synthetic distribution.', stressFile);
    stresses_Pa = lognrnd(1.5, 0.8, 55226, 1) * 1e6;
end

% Convert Pa -> MPa
stresses_MPa = stresses_Pa / 1e6;
numNodes = length(stresses_MPa);

% Mandible Volume (from ANSYS export volume.txt: Total = 2.427e-5 m^3)
volume_m3  = 2.427e-5;
volume_cm3 = volume_m3 * 1e6;
mass_solid_g = (volume_m3 * BioPCL.Density) * 1000.0;

% Deformation values from Total_deformation.xls (ANSYS ABS baseline)
% Max total deformation = 2.4829e-5 m, Average = 5.9981e-6 m
max_def_ABS_mm = 2.4829e-5 * 1000; % Convert m -> mm
avg_def_ABS_mm = 5.9981e-6 * 1000; % Convert m -> mm

% Scale deformations to BioPCL (linear elasticity: u proportional to 1/E)
modulusRatio = E_ABS / BioPCL.E;
max_def_BioPCL_mm = max_def_ABS_mm * modulusRatio;
avg_def_BioPCL_mm = avg_def_ABS_mm * modulusRatio;

%% 3. NODAL SAFETY FACTOR & CLINICAL CLASSIFICATION
SF_Yield = BioPCL.SigmaYield_T ./ max(stresses_MPa, 1e-6);
SF_Ult   = BioPCL.SigmaUlt_T   ./ max(stresses_MPa, 1e-6);

numSafe  = sum(stresses_MPa <= BioPCL.SigmaYield_T);
numYield = sum(stresses_MPa > BioPCL.SigmaYield_T & stresses_MPa <= BioPCL.SigmaUlt_T);
numFail  = sum(stresses_MPa > BioPCL.SigmaUlt_T);

pctSafe  = (numSafe / numNodes) * 100;
pctYield = (numYield / numNodes) * 100;
pctFail  = (numFail / numNodes) * 100;

meanStress   = mean(stresses_MPa);
maxStress    = max(stresses_MPa);
medianStress = median(stresses_MPa);
stdStress    = std(stresses_MPa);
p95Stress    = prctile(stresses_MPa, 95);
p99Stress    = prctile(stresses_MPa, 99);

%% 4. FROST MECHANOSTAT FUNCTION (Voigt Composite Modulus)
% The Frost mechanostat applies to bone tissue growing INTO the scaffold.
% We use a Voigt composite modulus: E_comp = BVF*E_bone + (1-BVF)*E_scaffold
% This gives physiologically meaningful tissue-level microstrain.
computeMechanostat = @(sigma_MPa, P) ...
    (sigma_MPa ./ (BVF * E_bone + (1-BVF) * BioPCL.E * (1-P)^2)) * 1e6;

classifyMechanostat = @(ue) ...
    ternary(ue < 50, 'Disuse', ...
    ternary(ue < 1500, 'Homeostasis', ...
    ternary(ue < 3000, 'Osteoinductive', 'Overload')));

%% 5. DISPLAY CONSOLE REPORT
fprintf('--- 1. BIOPCL MATERIAL & CONSTRUCT SPECIFICATIONS ---------------------\n');
fprintf('  Material        : Bioprinted BioPCL (Polycaprolactone)\n');
fprintf('  Density         : %.0f kg/m^3 (%.2f g/cm^3)\n', BioPCL.Density, BioPCL.Density/1000);
fprintf('  Young''s Modulus : %.0f MPa (0.35 GPa)\n', BioPCL.E);
fprintf('  Poisson''s Ratio : %.2f\n', BioPCL.nu);
fprintf('  Yield Strength  : %.0f MPa (tensile & compressive)\n', BioPCL.SigmaYield_T);
fprintf('  Ultimate Str.   : %.0f MPa (tensile & compressive)\n', BioPCL.SigmaUlt_T);
fprintf('  Shear Modulus   : %.1f MPa\n', BioPCL.G);
fprintf('  Construct Vol.  : %.2f cm^3\n', volume_cm3);
fprintf('  Solid Mass      : %.2f g\n\n', mass_solid_g);

fprintf('--- 2. FEA NODAL STRESS ANALYSIS (%d Nodes) --------\n', numNodes);
fprintf('  Max Von Mises   : %.2f MPa\n', maxStress);
fprintf('  99th Percentile : %.2f MPa\n', p99Stress);
fprintf('  95th Percentile : %.2f MPa\n', p95Stress);
fprintf('  Mean Stress     : %.2f MPa\n', meanStress);
fprintf('  Median Stress   : %.2f MPa\n', medianStress);
fprintf('  Std. Deviation  : %.2f MPa\n', stdStress);
fprintf('  Max Deform. ABS : %.4f mm (%.1f um)\n', max_def_ABS_mm, max_def_ABS_mm*1000);
fprintf('  Max Deform. PCL : %.4f mm (%.1f um)\n', max_def_BioPCL_mm, max_def_BioPCL_mm*1000);
fprintf('  Avg Deform. PCL : %.4f mm (%.1f um)\n', avg_def_BioPCL_mm, avg_def_BioPCL_mm*1000);
fprintf('  Scaling Factor  : E_ABS/E_PCL = %.4f\n\n', modulusRatio);

fprintf('--- 3. CLINICAL SAFETY CLASSIFICATION -----------------------------------\n');
fprintf('  [SAFE]   <= %.0f MPa : %6d nodes (%.2f%%)\n', BioPCL.SigmaYield_T, numSafe, pctSafe);
fprintf('  [YIELD]  %.0f-%.0f MPa  : %6d nodes (%.2f%%)\n', BioPCL.SigmaYield_T, BioPCL.SigmaUlt_T, numYield, pctYield);
fprintf('  [EXCEED] > %.0f MPa   : %6d nodes (%.2f%%)\n\n', BioPCL.SigmaUlt_T, numFail, pctFail);

%% 6. GIBSON-ASHBY SCAFFOLD OPTIMIZATION TABLE
fprintf('--- 4. GIBSON-ASHBY SCAFFOLD OPTIMIZATION TABLE -------------------------\n');
fprintf('  P(%%)   RelRho  E_eff(MPa)   sy(MPa)   Mass(g)   Median ue    Mechanostat\n');
fprintf('  --------------------------------------------------------------------------\n');
porosities_table = [0.0, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50];
for i = 1:length(porosities_table)
    P = porosities_table(i);
    relRho = 1.0 - P;
    E_eff  = BioPCL.E * (relRho^2.0);
    sy_eff = BioPCL.SigmaYield_T * (relRho^1.5);
    m_g    = mass_solid_g * relRho;
    ue     = computeMechanostat(medianStress, P);
    if ue < 50, status = 'Disuse';
    elseif ue < 1500, status = 'Homeostasis';
    elseif ue < 3000, status = 'Osteoinductive';
    else, status = 'Overload';
    end
    fprintf('  %5.0f  %7.2f  %10.1f  %8.2f  %8.2f  %10.0f  %16s\n', ...
            P*100, relRho, E_eff, sy_eff, m_g, ue, status);
end
fprintf('\n');

%% 7. 7-ZONE ANATOMICAL RESPONSE (Bilateral Molar Bite, 600 N)
% Zone definitions: name, area(mm^2), Ixx(mm^4), porosity
zoneNames = {'Condylar Head', 'Condylar Neck', 'Ramus', 'Mandibular Angle', ...
             'Molar Region', 'Premolar Region', 'Symphysis'};
zoneA     = [95, 80, 320, 380, 290, 260, 220];    % mm^2
zoneIxx   = [280, 210, 1850, 2200, 1600, 1400, 1100]; % mm^4
zoneP     = [0.15, 0.10, 0.30, 0.25, 0.20, 0.30, 0.20]; % design porosity

% Physiological loads (bilateral molar bite, 600 N total)
zoneM  = [4800, 6200, 9800, 12500, 15200, 11800, 8500]; % N.mm bending moment
zoneV  = [180, 220, 310, 380, 420, 350, 280];            % N shear
zoneF  = [120, 150, 200, 250, 300, 220, 180];            % N axial

fprintf('--- 5. 7-ZONE ANATOMICAL RESPONSE (Bilateral Molar Bite, 600 N) ---------\n');
fprintf('  %-20s %8s %8s %8s %6s %6s %8s %18s\n', ...
        'Zone', 'sig_VM', 'sig_b', 'tau', 'SF_y', 'SF_u', 'ue', 'Mechanostat');
fprintf('  %s\n', repmat('-', 1, 98));

zone_sigVM = zeros(1, 7);
zone_SF    = zeros(1, 7);
zone_ue    = zeros(1, 7);
zone_sigB  = zeros(1, 7);
zone_tau   = zeros(1, 7);

for k = 1:7
    A = zoneA(k); Ixx = zoneIxx(k); P = zoneP(k);
    E_eff = BioPCL.E * (1-P)^2;
    sy_eff = BioPCL.SigmaYield_T * (1-P)^1.5;
    su_eff = BioPCL.SigmaUlt_T * (1-P)^1.5;

    % c = distance from neutral axis to extreme fiber
    % For elliptical sections: I ~ A*h^2/12, h = sqrt(12*I/A), c = h/2
    c = sqrt(3.0 * Ixx / A);

    sigma_b = (zoneM(k) * c) / Ixx;
    sigma_a = abs(zoneF(k)) / A;
    sigma_total = sigma_b + sigma_a;

    % Shear stress (Jourawski approximation)
    Q_approx = Ixx * 0.40 / c;
    b_width = A / (2.0 * c);
    tau = (zoneV(k) * Q_approx) / (Ixx * b_width);

    sigma_vm = sqrt(sigma_total^2 + 3.0 * tau^2);
    SF_y = sy_eff / max(sigma_vm, 1e-6);
    SF_u = su_eff / max(sigma_vm, 1e-6);
    ue = computeMechanostat(sigma_vm, P);

    if ue < 50, status = 'Disuse';
    elseif ue < 1500, status = 'Homeostasis';
    elseif ue < 3000, status = 'Osteoinductive';
    else, status = 'Overload';
    end

    zone_sigVM(k) = sigma_vm;
    zone_SF(k) = SF_y;
    zone_ue(k) = ue;
    zone_sigB(k) = sigma_b;
    zone_tau(k) = tau;

    fprintf('  %-20s %8.2f %8.2f %8.2f %6.2f %6.2f %8.0f %18s\n', ...
            zoneNames{k}, sigma_vm, sigma_b, tau, SF_y, SF_u, ue, status);
end

[maxVM, maxIdx] = max(zone_sigVM);
fprintf('\n  Critical zone: %s (sigma_VM = %.2f MPa, SF_yield = %.2f)\n\n', ...
        zoneNames{maxIdx}, maxVM, zone_SF(maxIdx));

%% 8. GENERATE VISUALIZATIONS
fprintf('--- 6. GENERATING PUBLICATION FIGURES ------------------------------------\n');

% Figure 1: Nodal Stress Distribution & Safety Classification
fig1 = figure('Name', 'FEA Nodal Stress Analysis', ...
              'Color', [0.02 0.04 0.06], 'Position', [50, 50, 1400, 550]);

subplot(1, 2, 1);
set(gca, 'Color', [0.04 0.08 0.13], 'XColor', [0.8 0.9 1], 'YColor', [0.8 0.9 1]);
hold on; grid on;
histogram(stresses_MPa, 80, 'BinLimits', [0 45], 'FaceColor', [0 0.6 0.8], ...
          'EdgeColor', 'none', 'Normalization', 'count');
xline(BioPCL.SigmaYield_T, '--', 'Color', [0.22 1.00 0.08], 'LineWidth', 2);
xline(BioPCL.SigmaUlt_T, '--', 'Color', [1 0.2 0.4], 'LineWidth', 2);
xline(meanStress, '-', 'Color', [1 0.8 0], 'LineWidth', 1.5);
xline(p95Stress, ':', 'Color', [0.8 0.27 1], 'LineWidth', 1.5);
xlabel('Von Mises Stress (MPa)', 'FontSize', 11, 'Color', [0.8 0.9 1]);
ylabel('Node Count', 'FontSize', 11, 'Color', [0.8 0.9 1]);
title('ANSYS Nodal Stress Distribution (BioPCL Thresholds)', 'FontSize', 12, 'Color', [0 0.83 1]);
legend({sprintf('Yield (%.0f MPa)', BioPCL.SigmaYield_T), ...
        sprintf('Ultimate (%.0f MPa)', BioPCL.SigmaUlt_T), ...
        sprintf('Mean (%.2f MPa)', meanStress), ...
        sprintf('95th pct (%.2f MPa)', p95Stress)}, ...
        'TextColor', [0.8 0.9 1], 'Location', 'northeast', 'FontSize', 8);

subplot(1, 2, 2);
set(gca, 'Color', [0.04 0.08 0.13], 'XColor', [0.8 0.9 1], 'YColor', [0.8 0.9 1]);
hold on; grid on;
cats = categorical({'Safe (<= 15)', 'Yield (15-23)', 'Exceed (> 23)'});
cats = reordercats(cats, {'Safe (<= 15)', 'Yield (15-23)', 'Exceed (> 23)'});
vals = [pctSafe, pctYield, pctFail];
b = bar(cats, vals, 0.5, 'FaceColor', 'flat');
b.CData(1,:) = [0.22 1.00 0.08];
b.CData(2,:) = [1.00 0.80 0.00];
b.CData(3,:) = [1.00 0.20 0.40];
ylim([0 115]);
ylabel('Percentage of Nodes (%)', 'FontSize', 11, 'Color', [0.8 0.9 1]);
title('Structural Safety Classification (BioPCL)', 'FontSize', 12, 'Color', [0 0.83 1]);
for i = 1:3
    text(i, vals(i) + 3, sprintf('%.2f%%\n(%d nodes)', vals(i), ...
         [numSafe, numYield, numFail]'), 'HorizontalAlignment', 'center', ...
         'Color', [0.9 0.95 1], 'FontWeight', 'bold', 'FontSize', 9);
end

fprintf('  [1/3] Nodal Stress Distribution & Safety Classification\n');
exportgraphics(fig1, 'fig1_biopcl_nodal_stress.png', 'Resolution', 180, ...
               'BackgroundColor', [0.02 0.04 0.06]);

% Figure 2: Gibson-Ashby Porosity Analysis
fig2 = figure('Name', 'Gibson-Ashby Scaffold Optimization', ...
              'Color', [0.02 0.04 0.06], 'Position', [80, 80, 1400, 550]);

P_sweep = linspace(0, 0.60, 100);
E_sweep = BioPCL.E * (1 - P_sweep).^2.0;
Sy_sweep = BioPCL.SigmaYield_T * (1 - P_sweep).^1.5;
Su_sweep = BioPCL.SigmaUlt_T * (1 - P_sweep).^1.5;
mass_sweep = mass_solid_g * (1 - P_sweep);

subplot(1, 2, 1);
set(gca, 'Color', [0.04 0.08 0.13], 'XColor', [0.8 0.9 1], 'YColor', [0.8 0.9 1]);
yyaxis left;
plot(P_sweep*100, E_sweep, '-', 'Color', [0 0.83 1], 'LineWidth', 2.5);
ylabel('Effective Modulus E (MPa)', 'Color', [0 0.83 1]);
yyaxis right;
hold on;
plot(P_sweep*100, Sy_sweep, '--', 'Color', [1 0.42 0.21], 'LineWidth', 2.5);
plot(P_sweep*100, Su_sweep, ':', 'Color', [1 0.2 0.4], 'LineWidth', 2.0);
ylabel('Strength (MPa)', 'Color', [1 0.42 0.21]);
xlabel('Bioprinted Scaffold Porosity P (%)');
title('Gibson-Ashby Mechanical Properties vs. Porosity', 'Color', [0 0.83 1]);
legend({'Modulus E(P)', 'Yield Strength', 'Ultimate Strength'}, ...
       'TextColor', [0.8 0.9 1], 'Location', 'northeast', 'FontSize', 8);
grid on;

subplot(1, 2, 2);
set(gca, 'Color', [0.04 0.08 0.13], 'XColor', [0.8 0.9 1], 'YColor', [0.8 0.9 1]);
hold on; grid on;
plot(P_sweep*100, mass_sweep, '-', 'Color', [0.22 1 0.08], 'LineWidth', 2.5);
fill([P_sweep*100, fliplr(P_sweep*100)], [mass_sweep, zeros(size(mass_sweep))], ...
     [0.22 1 0.08], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
xlabel('Bioprinted Scaffold Porosity P (%)');
ylabel('Construct Mass (grams)', 'Color', [0.8 0.9 1]);
title('BioPCL Construct Mass vs. Porosity', 'Color', [0 0.83 1]);

fprintf('  [2/3] Gibson-Ashby Scaffold Optimization\n');
exportgraphics(fig2, 'fig2_biopcl_gibson_ashby.png', 'Resolution', 180, ...
               'BackgroundColor', [0.02 0.04 0.06]);

% Figure 3: 7-Zone Anatomical Response
fig3 = figure('Name', '7-Zone Biomechanical Response', ...
              'Color', [0.02 0.04 0.06], 'Position', [110, 110, 1500, 700]);

zoneColors = [0 0.83 1; 0 0.6 0.8; 1 0.42 0.21; 1 0.8 0; ...
              0.8 0.27 1; 0 1 0.67; 0.22 1 0.08];

subplot(2, 2, 1);
set(gca, 'Color', [0.04 0.08 0.13], 'XColor', [0.8 0.9 1], 'YColor', [0.8 0.9 1]);
hold on; grid on;
for k = 1:7
    bar(k, zone_sigVM(k), 0.55, 'FaceColor', zoneColors(k,:), 'EdgeColor', 'w');
    text(k, zone_sigVM(k)+1, sprintf('%.1f', zone_sigVM(k)), ...
         'HorizontalAlignment', 'center', 'Color', [0.8 0.9 1], 'FontSize', 8);
end
yline(BioPCL.SigmaYield_T, '--g', 'LineWidth', 1.5);
yline(BioPCL.SigmaUlt_T, '--r', 'LineWidth', 1.5);
set(gca, 'XTick', 1:7, 'XTickLabel', cellfun(@(s) s(1:min(8,end)), zoneNames, 'Uni', 0));
ylabel('Von Mises Stress (MPa)', 'Color', [0.8 0.9 1]);
title('Von Mises Stress per Zone', 'Color', [0 0.83 1]);

subplot(2, 2, 2);
set(gca, 'Color', [0.04 0.08 0.13], 'XColor', [0.8 0.9 1], 'YColor', [0.8 0.9 1]);
hold on; grid on;
sfColors = [1 0.2 0.4; 1 0.2 0.4; 1 0.2 0.4; 1 0.2 0.4; ...
            1 0.2 0.4; 1 0.2 0.4; 1 0.2 0.4]; % All SF < 1 -> red
for k = 1:7
    bar(k, zone_SF(k), 0.55, 'FaceColor', sfColors(k,:), 'EdgeColor', 'w');
    text(k, zone_SF(k)+0.02, sprintf('%.2f', zone_SF(k)), ...
         'HorizontalAlignment', 'center', 'Color', [0.8 0.9 1], 'FontSize', 8);
end
yline(1.0, '--r', 'SF = 1.0', 'LineWidth', 1.5);
yline(2.0, '--g', 'SF = 2.0', 'LineWidth', 1.0);
set(gca, 'XTick', 1:7, 'XTickLabel', cellfun(@(s) s(1:min(8,end)), zoneNames, 'Uni', 0));
ylabel('Safety Factor (Yield)', 'Color', [0.8 0.9 1]);
title('Yield Safety Factor per Zone', 'Color', [0 0.83 1]);

subplot(2, 2, 3);
set(gca, 'Color', [0.04 0.08 0.13], 'XColor', [0.8 0.9 1], 'YColor', [0.8 0.9 1]);
hold on; grid on;
bdata = [zone_sigB; abs(zoneF)./zoneA; zone_tau]';
bh = bar(1:7, bdata, 'grouped');
bh(1).FaceColor = [0 0.83 1];
bh(2).FaceColor = [0.22 1.00 0.08];
bh(3).FaceColor = [1 0.42 0.21];
set(gca, 'XTick', 1:7, 'XTickLabel', cellfun(@(s) s(1:min(8,end)), zoneNames, 'Uni', 0));
ylabel('Stress (MPa)', 'Color', [0.8 0.9 1]);
title('Stress Component Breakdown', 'Color', [0 0.83 1]);
legend({'Bending', 'Axial', 'Shear'}, 'TextColor', [0.8 0.9 1], 'FontSize', 8);

subplot(2, 2, 4);
set(gca, 'Color', [0.04 0.08 0.13], 'XColor', [0.8 0.9 1], 'YColor', [0.8 0.9 1]);
hold on; grid on;
barh(1:7, zone_ue, 0.55, 'FaceColor', [1 0.2 0.4], 'EdgeColor', 'w');
xline(1500, '--', 'Color', [1 0.8 0], 'LineWidth', 1.5);
xline(3000, '--', 'Color', [1 0.2 0.4], 'LineWidth', 1.5);
set(gca, 'YTick', 1:7, 'YTickLabel', cellfun(@(s) s(1:min(8,end)), zoneNames, 'Uni', 0));
xlabel('Tissue Microstrain (ue)', 'Color', [0.8 0.9 1]);
title('Frost Mechanostat Classification', 'Color', [0 0.83 1]);

sgtitle(sprintf('BIOPCL MANDIBLE -- 7-ZONE RESPONSE (600 N Bilateral Molar Bite)'), ...
        'Color', [0 0.83 1], 'FontSize', 13, 'FontWeight', 'bold');

fprintf('  [3/3] 7-Zone Biomechanical Response\n');
exportgraphics(fig3, 'fig3_biopcl_zone_analysis.png', 'Resolution', 180, ...
               'BackgroundColor', [0.02 0.04 0.06]);

fprintf('\n===========================================================================\n');
fprintf('  [DONE] All analyses complete. 3 figures saved to working directory.\n');
fprintf('===========================================================================\n');

%% Helper function
function result = ternary(cond, a, b)
    if cond, result = a; else, result = b; end
end
