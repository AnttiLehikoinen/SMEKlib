function [inds_pp, indsInMatrix] = get_relevantPolePairs(freqs, fm, N_pp)
%get_relevantPolePairs returns relevant pole-pair numbers and their
%indices.
% 
% 
% [inds_pp, indsInMatrix] = get_relevantPolePairs(freqs, fm)
% 
% (c) 2017 Antti Lehikoinen / Aalto University


TOL = 1e-3;

freqs_transf = fm*(0:N_pp);

inds_pp = find( ismembertol(freqs_transf, get_RelevantFrequencies(freqs), TOL) ) - 1;

indsInMatrix = sort([inds_pp(inds_pp==0)+1 inds_pp(inds_pp>0)*2 inds_pp(inds_pp>0)*2+1]);

end