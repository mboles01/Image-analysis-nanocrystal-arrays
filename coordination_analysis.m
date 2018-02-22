function out = coordination_analysis(image5,bonds_filtered,border,scale,lattice);

% Counts particle coordination numbers 

% Written by Michael Boles, University of Chicago March 2014, updated
% February 2018

centroid_numbers = [bonds_filtered(:,1), bonds_filtered(:,4)];
[uniques1,numuniques1] = count_unique(centroid_numbers);
uniques = [uniques1,numuniques1];

% add coordination number to bonds_filtered list:
% (particle no.) (x-coord) (y-coord) (coordination no.) (particle no.) (x-coord) (y-coord) (coordination no.) (bond length)

coord_nos = [];
for i=1:length(bonds_filtered)
    coord_no_1_temp = uniques(find(uniques(:,1) == bonds_filtered(i,1)),2);
    coord_no_2_temp = uniques(find(uniques(:,1) == bonds_filtered(i,4)),2);
    coord_nos_temp = [bonds_filtered(i,1:3), coord_no_1_temp, bonds_filtered(i,4:6), coord_no_2_temp, bonds_filtered(i,7)];
    coord_nos = [coord_nos; coord_nos_temp];
end

% eliminate spurious low-coordination data from particles near the edge of
% the image set by input value "border" (in nanometers)

coord_nos_filter = [];
for i = 1:length(coord_nos)
    if coord_nos(i,2) > border/scale && coord_nos(i,2) < length(image5) - border/scale && coord_nos(i,3) > border/scale && coord_nos(i,3) < length(image5) - border/scale == 1
        coord_no_1_filter_temp = coord_nos(i,4);
    else
        coord_no_1_filter_temp = NaN;
    end
    if coord_nos(i,6) > border/scale && coord_nos(i,6) < length(image5) - border/scale && coord_nos(i,7) > border/scale && coord_nos(i,7) < length(image5)-border/scale == 1
        coord_no_2_filter_temp = coord_nos(i,8);
    else
        coord_no_2_filter_temp = NaN;
    end
    coord_nos_filter_temp = [coord_nos(i,1:3), coord_no_1_filter_temp, coord_nos(i,5:7), coord_no_2_filter_temp, coord_nos(i,9)];
    coord_nos_filter = [coord_nos_filter; coord_nos_filter_temp];
end

% create output matrices describing positions of particles with a given
% coordination number

CN2 = [];
for i = 1:length(coord_nos_filter)
    if coord_nos_filter(i,4) == 2
        CN2_1_temp = coord_nos_filter(i,1:3);
    else 
        CN2_1_temp = [];
    end
    if coord_nos_filter(i,8) == 2
        CN2_2_temp = coord_nos_filter(i,5:7);
    else
        CN2_2_temp = [];
    end
    CN2 = unique([CN2; CN2_1_temp; CN2_2_temp],'rows');
end

CN3 = [];
for i = 1:length(coord_nos_filter)
    if coord_nos_filter(i,4) == 3
        CN3_1_temp = coord_nos_filter(i,1:3);
    else 
        CN3_1_temp = [];
    end
    if coord_nos_filter(i,8) == 3
        CN3_2_temp = coord_nos_filter(i,5:7);
    else
        CN3_2_temp = [];
    end
    CN3 = unique([CN3; CN3_1_temp; CN3_2_temp],'rows');
end

CN4 = [];
for i = 1:length(coord_nos_filter)
    if coord_nos_filter(i,4) == 4
        CN4_1_temp = coord_nos_filter(i,1:3);
    else 
        CN4_1_temp = [];
    end
    if coord_nos_filter(i,8) == 2
        CN4_2_temp = coord_nos_filter(i,5:7);
    else
        CN4_2_temp = [];
    end
    CN4 = unique([CN4; CN4_1_temp; CN4_2_temp],'rows');
end

CN5 = [];
for i = 1:length(coord_nos_filter)
    if coord_nos_filter(i,4) == 5
        CN5_1_temp = coord_nos_filter(i,1:3);
    else 
        CN5_1_temp = [];
    end
    if coord_nos_filter(i,8) == 5
        CN5_2_temp = coord_nos_filter(i,5:7);
    else
        CN5_2_temp = [];
    end
    CN5 = unique([CN5; CN5_1_temp; CN5_2_temp],'rows');
end

CN6 = [];
for i = 1:length(coord_nos_filter)
    if coord_nos_filter(i,4) == 6
        CN6_1_temp = coord_nos_filter(i,1:3);
    else 
        CN6_1_temp = [];
    end
    if coord_nos_filter(i,8) == 6
        CN6_2_temp = coord_nos_filter(i,5:7);
    else
        CN6_2_temp = [];
    end
    CN6 = unique([CN6; CN6_1_temp; CN6_2_temp],'rows');
end

out{1} = CN2; 
out{2} = CN3; 
out{3} = CN4; 
out{4} = CN5; 
out{5} = CN6; 


% create output matrices describing the positions of particles
% participating in bonds of various saturation values: 6-6, 6-5, 6-4, 6-3

bonds_66 = []; bonds_65 = []; bonds_64 = []; bonds_63 = [];
for i = 1:length(coord_nos_filter)
    if coord_nos_filter(i,4) == 6 && coord_nos_filter(i,8) == 6
        bonds_66_temp = [coord_nos_filter(i,:)];
        bonds_66 = [bonds_66; bonds_66_temp];
    elseif coord_nos_filter(i,4) == 6 && coord_nos_filter(i,8) == 5 || coord_nos_filter(i,4) == 5 && coord_nos_filter(i,8) == 6        
        bonds_65_temp = [coord_nos_filter(i,:)];
        bonds_65 = [bonds_65; bonds_65_temp];
    elseif coord_nos_filter(i,4) == 6 && coord_nos_filter(i,8) == 4 || coord_nos_filter(i,4) == 4 && coord_nos_filter(i,8) == 6        
        bonds_64_temp = [coord_nos_filter(i,:)];
        bonds_64 = [bonds_64; bonds_64_temp];
    elseif coord_nos_filter(i,4) == 6 && coord_nos_filter(i,8) == 3 || coord_nos_filter(i,4) == 3 && coord_nos_filter(i,8) == 6        
        bonds_63_temp = [coord_nos_filter(i,:)];
        bonds_63 = [bonds_63; bonds_63_temp];
    else
    end
end

out{6} = bonds_63;
out{7} = bonds_64;
out{8} = bonds_65;
out{9} = bonds_66;
out{10} = coord_nos_filter;

