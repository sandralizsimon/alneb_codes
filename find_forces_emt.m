function [energy, forces] = find_forces_emt(positions, box, De, a, re, rc, fixatoms,eval)

folder = append("conf",num2str(eval));
mkdir(folder);

% Copy make_coord.m to the new folder
copyfile('get_emt.py', folder);

copyfile('POSCAR_initial.1', folder);

copyfile('readPOSCAR.m', folder);

copyfile('writePOSCAR.m', folder);


% Change directory to the new folder
cd(folder);

% Making .com file
output_file = 'out.txt';

poscar_data = readPOSCAR('POSCAR_initial.1');  
filename = sprintf('POSCAR.1');
writePOSCAR(filename, poscar_data.comment, poscar_data.scale, poscar_data.lattice_vectors, poscar_data.atom_types, poscar_data.atom_counts, poscar_data.coordinate_type, positions, fixatoms)


system(['python3 get_emt.py POSCAR.1 >> ', output_file]);

forces = load('forces_output.txt');
energy = load('energy_output.txt');

energy = (energy - (36.6849522798661)).* 96.4853; % wrto reactant
forces = forces.* 96.4853;
energy = energy - 800;
cd('..');
end
