function out = bond_analysis(lattice,minlength,maxlength,scale)
% Created by Michael Boles, March 2014
% Gratefully acknowledge William Irvine (UChicago Physics) for assistance

% Triangulate the set of points
TRI = delaunay(lattice(:,1),lattice(:,2));

% Reduce the n-by-3 matrix to n-by-2
TRI2 = [TRI(:,[1,2]);TRI(:,[1,3]);TRI(:,[2,3])];

% Sort triangle elements so vertex indices are in ascending order across the column
TRIsort = sort(TRI2,2);

% Eliminate double-counted line segments
TRIunique = unique(TRIsort,'rows');

% Add (x,y)-coordinates to this list of connected points
% Bonds will be a n-by-6 matrix containing:
% [point no.  xval  yval  point no.  xval  yval]

bonds = [];
for i = 1:length(TRIunique)
    xval_point1_temp = lattice(TRIunique(i,1),1);
    yval_point1_temp = lattice(TRIunique(i,1),2);
    xval_point2_temp = lattice(TRIunique(i,2),1);
    yval_point2_temp = lattice(TRIunique(i,2),2);
    bonds_temp = [TRIunique(i,1), xval_point1_temp, yval_point1_temp, TRIunique(i,2), xval_point2_temp, yval_point2_temp];
    bonds = [bonds; bonds_temp];
end

% Measure bond lengths and add to bonds_paired
bond_lengths = scale*sqrt((bonds(:,2)-bonds(:,5)).^2+(bonds(:,3)-bonds(:,6)).^2);
bonds_unfiltered = [bonds, bond_lengths];

bonds_filtered = [];
for i = 1:length(bonds_unfiltered);
    if bonds_unfiltered(i,7) > minlength && bonds_unfiltered(i,7) < maxlength
        bonds_filtered_temp = bonds_unfiltered(i,:);
    else
        bonds_filtered_temp = [];
    end
    bonds_filtered = [bonds_filtered; bonds_filtered_temp];
end

out{1} = bonds_unfiltered;
out{2} = bonds_filtered;
