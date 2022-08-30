
%
% Write the Bezier extraction operators for a NURBS mesh defined
% by the given knots to jem-jive file.
%
% VP Nguyen
% Cardiff University, UK
% Feburary, 2013.

convert2DNurbs
generateIGA2DMesh

C        = bezierExtraction2D(uKnot,vKnot,p,q);

noElems  = size(C,3);

fileName = 'curvedBeam.mesh';

file = fopen(fileName, 'wt');

fprintf(file, '<Nodes>\n');

for i=1:length(controlPts)
   fprintf(file, '  %1d %1d %2.4f', i, controlPts(i,1),controlPts(i,2));
   fprintf(file, ';\n');
end

fprintf(file, '</Nodes>\n\n');

fprintf(file, '<Elements>\n');

for i=1:size(element,1)
   fprintf(file, '  %1d %1d', i-1, element(i,:) );
   fprintf(file, ';\n');
end

fprintf(file, '</Elements>\n\n');

fprintf(file, '<ElementDatabase name="C">\n');

fprintf(file, ' <Column name = "irows" type = "int">\n');

for e=1:noElems
    Ce = C(:,:,e);
    [row,col] = find(Ce);
    fprintf(file, '  %1d ', e-1);
    for i=1:length(row)
        fprintf(file, '%1d ', row(i)-1);
    end
    fprintf(file, ';\n');
end

fprintf(file, ' </Column>\n');

fprintf(file, ' <Column name = "jcols" type = "int">\n');

for e=1:noElems
    Ce = C(:,:,e);
    [row,col] = find(Ce);
    
    fprintf(file, '  %d ', e-1);
    for i=1:length(row)
        fprintf(file, '%1d ', col(i)-1);
    end
    fprintf(file, ';\n');
end

fprintf(file, ' </Column>\n');

fprintf(file, ' <Column name = "values" type = "float">\n');
for e=1:noElems
    Ce = C(:,:,e);
    [row,col,val] = find(Ce);
    
    fprintf(file, '  %d ', e-1);
    for i=1:length(row)
        fprintf(file, '%2.1f ', val(i));
    end
    fprintf(file, ';\n');
end
fprintf(file, ' </Column>\n');

fprintf(file, ' <Column name = "weights" type = "float">\n');

for i=1:size(element,1)
    w = weights(element(i,:));
    fprintf(file, '  %1d ',i-1);
    for i=1:length(w)
        fprintf(file, '%2.4f ', w(i));
    end
    fprintf(file, ';\n');
end

fprintf(file, ' </Column>\n');

fprintf(file, '</ElementDatabase>\n\n');

fprintf(file, '<NodeGroup name="gr1">\n{');

for i=1:length(fixedXNodes1)
   fprintf(file, '  %1d', fixedXNodes1(i));
  
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');
fprintf(file, '<NodeGroup name="gr2">\n{');

for i=1:length(fixedXNodes2)
   fprintf(file, '  %1d', fixedXNodes2(i));
  
end
fprintf(file, '}\n');
fprintf(file, '</NodeGroup>\n');

fclose(file);



