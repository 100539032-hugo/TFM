% =========================================================================
% EXTRACCIÓN ROBUSTA OPTISTRUCT CON REGISTRO DE PID Y GID
% =========================================================================

clear; clc;

baseFolder   = 'C:\Users\Hugo\OneDrive\Escritorio\sim_part1';
exportFolder = 'C:\Users\Hugo\OneDrive\Escritorio\surrogate_1';
numSims      = 21;

if ~exist(exportFolder, 'dir')
    mkdir(exportFolder);
end

% Tabla con columna adicional 'PID_Disp_Max'
results = table('Size', [numSims, 7], ...
    'VariableTypes', {'string', 'double', 'double', 'uint32', 'uint32', 'double', 'double'}, ...
    'VariableNames', {'Simulacion', 'Masa_Total', 'Desplazamiento_Max', 'Nodo_Disp_Max', 'PID_Disp_Max', 'Energia_Deformacion', 'Factor_Pandeo_BLF'});

for i = 1:numSims
    folderName = sprintf('%03d', i);
    currentFolder = fullfile(baseFolder, folderName);
    results.Simulacion(i) = string(folderName);
    
    if ~exist(currentFolder, 'dir')
        warning('No existe la carpeta: %s', currentFolder);
        continue;
    end
    
    outLines = getFileLines(currentFolder, '*.out');
    pchLines = getFileLines(currentFolder, '*.pch');
    kpiLines = getFileLines(currentFolder, '*.kpi');
    
    % 1. MASA TOTAL
    results.Masa_Total(i) = extractMass(outLines);
    
    % 2. FACTOR DE PANDEO (BLF)
    results.Factor_Pandeo_BLF(i) = extractBLF(outLines);
    
    % 3. ENERGÍA DE DEFORMACIÓN
    results.Energia_Deformacion(i) = extractStrainEnergy(pchLines, outLines);
    
    % 4. DESPLAZAMIENTO MÁXIMO, NODO (GID) Y COMPONENTE (PID)
    [results.Desplazamiento_Max(i), results.Nodo_Disp_Max(i), results.PID_Disp_Max(i)] = extractMaxDispPanels(kpiLines, pchLines);
end

format shortG;
disp(results);

txtPath = fullfile(exportFolder, 'OptiStruct_KPI_Resultados_2.txt');
writetable(results, txtPath, 'Delimiter', '\t');
fprintf('\nExtracción completada. Archivo guardado en: %s\n', txtPath);

% =========================================================================
% FUNCIONES AUXILIARES
% =========================================================================

function lines = getFileLines(folder, pattern)
    files = dir(fullfile(folder, pattern));
    if ~isempty(files)
        content = fileread(fullfile(folder, files(1).name));
        lines = strsplit(content, {'\r\n', '\n'});
    else
        lines = {};
    end
end

function mass = extractMass(lines)
    mass = NaN;
    inTable = false;
    totalMass = 0;
    foundData = false;
    
    for k = 1:length(lines)
        line = strtrim(lines{k});
        if contains(line, 'CENTER OF GRAVITY TABLE', 'IgnoreCase', true)
            inTable = true;
            totalMass = 0;
            foundData = false;
            continue;
        end
        if inTable
            if foundData && (startsWith(line, '---') || isempty(line))
                mass = totalMass;
                return;
            end
            if startsWith(line, '---') || contains(line, 'X-COG') || contains(line, 'Component')
                continue;
            end
            parts = strsplit(line);
            if length(parts) >= 2
                compID   = str2double(parts{1});
                compMass = str2double(parts{2});
                if ~isnan(compID) && ~isnan(compMass)
                    totalMass = totalMass + compMass;
                    foundData = true;
                end
            end
        end
    end
    if foundData, mass = totalMass; end
end

function blf = extractBLF(lines)
    blf = NaN;
    inBlock = false;
    for k = 1:length(lines)
        line = strtrim(lines{k});
        if contains(line, 'Subcase', 'IgnoreCase', true) && contains(line, 'Eigenvalue', 'IgnoreCase', true)
            inBlock = true;
            continue;
        end
        if inBlock
            parts = strsplit(line);
            if length(parts) >= 3
                modeNum = str2double(parts{2});
                eigVal  = str2double(parts{3});
                if modeNum == 1 && ~isnan(eigVal)
                    blf = eigVal;
                    return;
                end
            end
        end
    end
end

function ese = extractStrainEnergy(pchLines, outLines)
    ese = NaN;
    inBlock = false;
    totalVal = 0;
    found = false;
    for k = 1:length(pchLines)
        line = strtrim(pchLines{k});
        if contains(line, '$GROUP STRAIN ENERGIES', 'IgnoreCase', true)
            inBlock = true;
            continue;
        end
        if inBlock
            if found && (startsWith(line, '$') || isempty(line)), break; end
            if startsWith(line, '$') || isempty(line), continue; end
            parts = strsplit(line);
            if length(parts) >= 3
                val = str2double(parts{3});
                if ~isnan(val)
                    totalVal = totalVal + val;
                    found = true;
                end
            end
        end
    end
    if found, ese = totalVal; return; end
    
    for k = 1:length(outLines)
        line = strtrim(outLines{k});
        if contains(line, 'STRAIN ENERGY', 'IgnoreCase', true) || contains(line, 'COMPLIANCE', 'IgnoreCase', true)
            parts = strsplit(line);
            val = str2double(parts{end});
            if ~isnan(val), ese = val; return; end
        end
    end
end

function [maxDisp, maxGID, maxPID] = extractMaxDispPanels(kpiLines, pchLines)
    maxDisp = NaN;
    maxGID  = uint32(0);
    maxPID  = uint32(0);
    
    if ~isempty(kpiLines)
        inBlock = false;
        maxMag = -1;
        found = false;
        panelRowCount = 0;
        
        for k = 1:length(kpiLines)
            line = strtrim(kpiLines{k});
            if contains(line, '$DISPLACEMENTS', 'IgnoreCase', true)
                inBlock = true;
                continue;
            end
            if inBlock
                if startsWith(line, '$') && ~contains(line, '$SUBCASE')
                    break;
                end
                if startsWith(line, '$') || isempty(line) || contains(line, 'PID')
                    continue;
                end
                
                parts = strsplit(line);
                if length(parts) >= 3
                    panelRowCount = panelRowCount + 1;
                    
                    % Detener la lectura al superar los 6 paneles
                    if panelRowCount > 6
                        break;
                    end
                    
                    pidVal = str2double(parts{1}); % Columna 1: PID
                    gidVal = str2double(parts{2}); % Columna 2: GID
                    magVal = str2double(parts{3}); % Columna 3: MAG
                    
                    if ~isnan(magVal) && (magVal > maxMag)
                        maxMag = magVal;
                        if ~isnan(gidVal), maxGID = uint32(gidVal); end
                        if ~isnan(pidVal), maxPID = uint32(pidVal); end
                        found = true;
                    end
                end
            end
        end
        if found
            maxDisp = maxMag;
            return;
        end
    end
    
    % Respaldo .pch
    if ~isempty(pchLines)
        inBlock = false;
        maxMag = -1;
        found = false;
        for k = 1:length(pchLines)
            line = strtrim(pchLines{k});
            if contains(line, '$DISPLACEMENTS', 'IgnoreCase', true)
                inBlock = true;
                continue;
            end
            if inBlock
                if found && startsWith(line, '$') && ~contains(line, '$SUBCASE') && ~contains(line, '$REAL'), break; end
                if startsWith(line, '$') || isempty(line) || startsWith(line, '-CONT-'), continue; end
                parts = strsplit(line);
                if length(parts) >= 5 && (strcmp(parts{2}, 'G') || strcmp(parts{2}, 'C'))
                    nodeID = str2double(parts{1});
                    ux     = str2double(parts{3});
                    uy     = str2double(parts{4});
                    uz     = str2double(parts{5});
                    if ~isnan(ux) && ~isnan(uy) && ~isnan(uz)
                        mag = sqrt(ux^2 + uy^2 + uz^2);
                        if mag > maxMag
                            maxMag = mag;
                            maxGID = uint32(nodeID);
                            maxPID = uint32(0);
                            found = true;
                        end
                    end
                end
            end
        end
        if found, maxDisp = maxMag; end
    end
end