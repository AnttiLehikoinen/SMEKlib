function M = assemble_BoundaryFluxMatrix(msh, boundaryData, v, M)

M = assemble_IntermeshBoundaryFluxMatrix(msh, msh, [boundaryData; boundaryData], v, M);

end