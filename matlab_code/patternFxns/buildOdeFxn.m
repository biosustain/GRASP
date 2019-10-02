function buildOdeFxn(ensemble, kineticFxn, strucIdx)
% Builds the ODE function .m file that is used to simulate the model.
%
% USAGE:
%
%    buildOdeFxn(ensemble, kineticFxn, strucIdx)
%
% INPUTS:
%    ensemble (`struct`):   model ensemble
%    kineticFxn (`char`):	name of the kinetic function
%    strucIdx (`int`):      ID of the model structure
%
% OUTPUT:
%    written .m file with the ODEs 
%
% .. Authors:
%       - Marta Matos	2019 original code 

currentPath = regexp(mfilename('fullpath'), '(.*)/', 'match');
filepath = fullfile(currentPath{1}, '..', '..', 'reactions', [ensemble.description, '_', num2str(strucIdx)], [kineticFxn,'.m']);

fIn = fopen(filepath);
fOut = fopen([filepath(1:(end-2)), '_ode.m'], 'w');

elseCount = 0;
lineIn = fgets(fIn);

while ischar(lineIn)
    
    if startsWith(lineIn, 'function')
        lineIn = strrep(lineIn, '[f,grad]', 'y');
        lineIn = strrep(lineIn, '(x,', '_ode(y,Eref,metsRefConc,');
        fprintf(fOut,lineIn);
    end
    
    if startsWith(lineIn, 'else')
        elseCount = elseCount + 1;
    end
    
    if (startsWith(lineIn, 'v = ') || startsWith(lineIn, 'E = ')) && elseCount == 1
        lineIn = strrep(lineIn, 'x', 'Eref');
        fprintf(fOut,lineIn);
    end

    if startsWith(lineIn, 'm_')
        lineIn = strrep(lineIn, '= x(', '= y(');
        fprintf(fOut,lineIn);
    end

    if ~isempty(regexp(lineIn, 'E(\d+', 'match'))
        
        lineIn = regexprep(lineIn,'E\((\d+),:\) = x\(\d+,:', 'E($1,:) = Eref($1,:');
        
        fprintf(fOut,lineIn);
    end

    if startsWith(lineIn, 'v(') || startsWith(lineIn, 'E(kinInactRxns')
        lineIn = strrep(lineIn,'size(x,', 'size(y,');
        fprintf(fOut,lineIn);
    end
   
    lineIn = fgets(fIn);
end

fprintf(fOut, '\ny = (1./(metsRefConc.*10^6)) .* (Sred*(E.*v));');

fclose(fIn);
fclose(fOut);

end