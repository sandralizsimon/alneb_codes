function [new_positions, new_flag] = create_new_poscar(len, Force_values, fixatoms, positions)
fixed_atoms = [];
move_atoms = [];
j = 1;
k = 1;
for i = 1:len
    if strcmp(fixatoms{i, 1}, 'F') && strcmp(fixatoms{i, 2}, 'F') && strcmp(fixatoms{i, 3}, 'F')
        fixed_atoms = [fixed_atoms;positions(i,:)];
        f_flag(j,:) = {'F', 'F', 'F'};
        j = j+1;
    else
        move_atoms = [move_atoms;positions(i,:)];
        n_flag(k,:) = {'T', 'T', 'T'};
        k = k+1;
    end

end

new_positions = [move_atoms;fixed_atoms];
new_flag = [n_flag;f_flag];
end