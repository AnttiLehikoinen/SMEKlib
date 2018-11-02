%{
[ ] Better handling of materials, e.g. through a class object or similar.
[ ] Better and more general definition of conducting domains.

[ ] Error in replicate_sector_fixed: order of ccl nodes gets mixed
    sometimes.

[ ] Generalize gwrap's loadmesh for several element types, by scanning
entire elements description into an array, and then chopping it into
element types.
[x] Add possibility to use pre-defined triangulation with AirgapTriangulation.  
      [ ] Add higher-order elements here.  
[ ] Add higher-order elements to AGT in general. 
    [ ] 2nd-order works.
[ ] Polynomial interpolation inside an AGT-child class--> better BEMF.
[x] Finish WST2 function.

[x] Fix examples!
    [x] New mesh coarse mesh generator functionality
    [ ] Changed sign of imaginary unit in time-harmonic analysis??

[ ] @folder-style class definitions; otherwise Octave dies
[ ] Add extra periodic nodes to 2nd-order AGT
%}