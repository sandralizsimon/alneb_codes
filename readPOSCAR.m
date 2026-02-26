function poscar_data = readPOSCAR(filename)
   
    % Inputs: filename - the POSCAR file path
    % Outputs: poscar_data - a structure containing extracted data

    % Open the file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end

    % Initialize the output structure
    poscar_data = struct();

    % Read the header line
    poscar_data.comment = fgetl(fid); % First line is a comment/title

    % Read scaling factor
    scale_line = fgetl(fid);
    poscar_data.scale = str2double(scale_line);

    % Read lattice vectors
    poscar_data.lattice_vectors = zeros(3,3);
    for i = 1:3
        line = fgetl(fid);
        poscar_data.lattice_vectors(i, :) = sscanf(line, '%f %f %f');
    end

    % Read atom types and number of atoms
    atom_types_line = fgetl(fid); % Atom types
    poscar_data.atom_types = strsplit(strtrim(atom_types_line));

    atom_counts_line = fgetl(fid); % Number of atoms
    poscar_data.atom_counts = sscanf(atom_counts_line, '%d');

    % Read selective dynamics or coordinate type
    line = fgetl(fid);
    if startsWith(line, 'S', 'IgnoreCase', true)
        poscar_data.selective_dynamics = true;
        line = fgetl(fid); % Read next line
    else
        poscar_data.selective_dynamics = false;
    end

    % Check for Direct or Cartesian coordinates
    if startsWith(line, 'D', 'IgnoreCase', true)
        poscar_data.coordinate_type = 'Direct';
    elseif startsWith(line, 'C', 'IgnoreCase', true) || startsWith(line, 'K', 'IgnoreCase', true)
        poscar_data.coordinate_type = 'Cartesian';
    else
        error('Unknown coordinate type in POSCAR file');
    end
    
    % Read atomic positions
    num_atoms = sum(poscar_data.atom_counts);

    % Preallocate positions and flags
    poscar_data.positions = zeros(num_atoms, 3); % To store coordinates
    poscar_data.flags = cell(num_atoms, 3);      % To store T/F flags as strings

    for i = 1:num_atoms
        % Read a line from the file
        line = fgetl(fid);

        % Split the line into parts
        data = strsplit(strtrim(line));

        % Extract the first 3 values as atomic positions
        poscar_data.positions(i, :) = str2double(data(1:3));

        % Extract the next 3 values as T/F flags
        poscar_data.flags(i, :) = data(4:6);
    end

    % Close the file
    fclose(fid);
end
