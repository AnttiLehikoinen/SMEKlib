classdef SLContainer < handle
    %SLContainer A SMEKlib container class for key-value pairs.
    % 
    % Has the following methods:
    %
    % value = get(key); returns [] if not found.
    %
    % set(key, value) (overwrites existing).
    %
    % add(key, value) (append to existing; only for arrays).
    %
    % keys() returns the keys.
    %
    % copy() returns a deep copy.
    %
    % (c) 2017 Antti Lehikoinen / Aalto University
    properties
        data
    end
    methods 
        function c = SLContainer()
            c.data = struct();
        end
        
        function val = get(c, name)
            if isfield(c.data, name)
                val = c.data.(name);
            else
                val = [];
            end
        end
        
        function c = set(c, name, d)
            c.data.(name) = d;
        end
        
        function c = add(c, name, d)
            if isfield(c.data, name)
                c.data.(name) = [c.data.(name) toRow(d)];
            else
                c.data.(name) = d;
            end
        end
        
        function names = keys(c)
            names = fieldnames(c.data);
        end
        
        function c2 = copy(c)
            c2 = SLContainer();
            names = c.keys();
            for k = 1:numel(names);
                c2.set(names{k}, c.get(names{k}));
            end
        end
    end
end