function BH = get_defaultMaterials(matNumber)
%get_defaultMaterials returns the BH curve for some materials
% 
% BH = get_defaultMaterials(matNumber)
% returns the BH curve of the FCSMEK default material matnumber, in the
% format BH = [B H];
%
% supplying a matNumber < 1 returns the number of available materials
% 
% Copyright (c) 2016 Antti Lehikoinen / Aalto University

%NOTE: due to possible confidentiality issues, only two materials are
%included here:
% 1 = construction steel (Ovako 520 L)
% 2 = Electrical steel sheet – Bochum STABOLEC 260-50 A
%
% Feel free to add your own :)


if matNumber < 0
    BH = 2; return;
elseif matNumber == 0
    BH = [0 1; 0 1/(pi*4e-7)]'; return
elseif matNumber <= 2
    HBdata = [           48          0.1           19          0.1
           97          0.2           39          0.2
          145          0.3           58          0.3
          194          0.4           77          0.4
          242          0.5           96          0.5
          291          0.6          116          0.6
          339          0.7          135          0.7
          388          0.8          154          0.8
          436          0.9          173          0.9
          492            1          195            1
          601          1.1          228          1.1
          758          1.2          266          1.2
          980          1.3          328          1.3
         1379          1.4          518          1.4
         2137          1.5         1051          1.5
         3585          1.6         2300          1.6
         6458          1.7         4834          1.7
        10373          1.8         9268          1.8
        15274          1.9        16215          1.9
        23176            2        26966            2];
else
    error('Material not available');
end

matIndex = mod(matNumber - 1, 5) + 1;
BH = HBdata(:, (matIndex-1)*2+[2 1]);

end
