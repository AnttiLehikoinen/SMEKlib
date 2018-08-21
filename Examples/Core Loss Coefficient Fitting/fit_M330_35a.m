%Fitting hysteresis (~fB^2) and eddy-current (~f^2B^2) loss coefficients to
% M330-35a datasheet from
% https://cogent-power.com/cms-data/downloads/m330-35a.pdf
%
% (c) 2018 Antti Lehikoinen \ Aalto University

%Raw data as a string; imported by copy-pasting the table contents from the
%pdf to Excel as text, copy-transposing the resulting datacolumn (still in
%Excel), and finally copy-pasting the data-row to Matlab;

datastr = '0,1	0.02	0.07	33.4	0.05	0.12	0.33	1.43	5.95	0,2	0.08	0.18	43.6	0.20	0.48	1.27	5.40	21.7	0,3	0.17	0.32	50.8	0.41	1.02	2.69	11.0	45.1	0,4	0.28	0.48	57.2	0.67	1.68	4.49	18.3	76.2	0,5	0.40	0.66	63.6	0.97	2.47	6.66	27.2	116	0,6	0.53	0.87	70.4	1.30	3.37	9.19	38.1	167	0,7	0.68	1.11	78.1	1.68	4.39	12.11	51.1	230	0,8	0.84	1.39	87.2	2.10	5.54	15.44	66.4	308	0,9	1.02	1.72	98.7	2.56	6.82	19.22	84.5	403	1,0	1.22	2.12	114	3.07	8.25	23.54	106	517	1,1	1.44	2.63	136	3.64	9.86	28.48	130	654	1,2	1.69	3.35	172	4.29	11.6	34.12	159	803	1,3	2.00	4.56	242	5.07	13.7	40.62	193	1,4	2.40	7.40	428	6.06	16.3	48.24	233	1,5	2.94	17.0	1027	7.40	19.6	57.86	279	1,6	3.67	46.2	2576	8.86	23.2	70.24	335	1,7	4.32	110	5409	1,8	4.73	220	9677';

%parsing data
table = zeros(18, 9); %empty array, the size of the original table
C = regexp(datastr, '[_\t]', 'split'); %splitting string
%going through the data; jumping rows at the T-values that use comma as the
%decimal separator
row = 1; table(row, 1) = str2double( strrep(C{1}, ',', '.') );
col = 2;
for k = 2:numel(C)
    if ~isempty( strfind(C{k}, ',') )
        row = row+1;
        col = 2;
        table(row, 1) = str2double( strrep(C{k}, ',', '.') );
    else
        table(row, col) = str2double( C{k} );
        col = col + 1;
    end
end

%creating arrays for losses (W), flux density (Bs) and frequency (fs)
%values
W = table(:, [2 5 6 7 8 9]);
Bs = repmat( table(:,1), 1, size(W,2));
fs = repmat( [50 100 200 400 1000 2500], size(table,1), 1);

W = W(:,1:4); Bs = Bs(:,1:4); fs = fs(:, 1:4); %limiting frequencies
%between 50 and 400 Hz

%fitting loss coefficients (linear least-squares)
d_hyst = fs.*Bs.^2;
d_eddy = fs.^2 .* Bs.^2;
c = [d_hyst( W>0 ) d_eddy( W>0 )] \ W( W>0 )

%recomputing losses
Wcomp = c(1)*d_hyst + c(2)*d_eddy;

Wcomp( W==0 ) = 0;

figure(1); clf; set(gcf, 'Name', 'Datasheet losses');
imagesc(fs(1,:), Bs(:,1), W); colorbar;

figure(2); clf; set(gcf, 'Name', 'Fitted losses');
imagesc(fs(1,:), Bs(:,1), W); colorbar;

figure(3); clf; set(gcf, 'Name', 'Error');
imagesc(fs(1,:), Bs(:,1), W-Wcomp); colorbar;