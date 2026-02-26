function writePOSCAR(filename, comment, scale, lattice_vectors, atom_types, atom_counts, coordinate_type, positions, flags)
    % writePOSCAR creates a POSCAR file with the provided data.
    % 
    % Inputs:
    % - filename: Name of the output POSCAR file.
    % - comment: Comment or title for the POSCAR file (first line).
    % - scale: Scaling factor for the lattice vectors.
    % - lattice_vectors: 3x3 matrix representing lattice vectors.
    % - atom_types: Cell array of atom types.
    % - atom_counts: Array of integers indicating the number of each atom type.
    % - coordinate_type: 'Direct' or 'Cartesian' indicating coordinate system.
    % - positions: Nx3 matrix of atomic positions.
    % - flags (optional): Nx3 cell array of 'T'/'F' strings for selective dynamics.

    % Open the file for writing
    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open file: %s', filename);
    end

    % Write the comment/title line
    fprintf(fid, '%s\n', comment);

    % Write the scaling factor
    fprintf(fid, '%.16f\n', scale);

    % Write the lattice vectors
    for i = 1:3
        fprintf(fid, '%.16f %.16f %.16f\n', lattice_vectors(i, :));
    end

    % Write atom types
    fprintf(fid, '%s ', atom_types{:});
    fprintf(fid, '\n');

    % Write atom counts
    fprintf(fid, '%d ', atom_counts);
    fprintf(fid, '\n');

    % Write selective dynamics if flags are provided
    if nargin > 8 && ~isempty(flags)
        fprintf(fid, 'Selective Dynamics\n');
    end

    % Write coordinate type
    fprintf(fid, '%s\n', coordinate_type);

    % Write atomic positions and optional flags
    num_atoms = size(positions, 1);
    for i = 1:num_atoms
        fprintf(fid, '%.16f %.16f %.16f', positions(i, :));
        if nargin > 8 && ~isempty(flags)
            fprintf(fid, ' %s %s %s', flags{i, :});
        end
        fprintf(fid, '\n');
    end

    % Close the file
    fclose(fid);
end
