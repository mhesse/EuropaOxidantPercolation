function [dof_bnd_in,dof_bnd_out] = find_bnd_cells(dof_in,dof_out,dof_f_bnd,D,Grid)
% author: Marc Hesse
% date: 22 Sep 2016
% Description:
% Finds the cells along a boundary defined by a given set of faces and 
% separates them in cell insise and outside the boundary.

% Input:
% dof_in    = dof's of the cells inside the domain
% dof_out   = dof's of the cells outside the domain
% dof_f_bnd = dof's of the faces on the boundary
% D         = N by Nf discrete divergence matrix
% Grid      = structure containing pertinent information about grid

% Output:
% dof_bnd_in  = dof's of the cells along outside of the boundary
% dof_bnd_out = dof's of the cells along inside of the boundary
DD = D(:,dof_f_bnd);
dof_bnd_in = Grid.dof(abs(sum(DD,2))>eps); % these are all cells. Using dof_bnd_in for tmp storage

if length(dof_in) > length(dof_out)
    % Cheaper to compare with exterior cells
    dof_bnd_in = Grid.dof(sum(abs(DD),2)>eps); 
    % Note, these are all cells! Using dof_bnd_in for temporary storage
    [dof_bnd_out,i_bnd_out] = intersect(dof_bnd_in,dof_out); % find exterior cells
    dof_bnd_in(i_bnd_out) = []; % delete exterior cells
else
    % Cheaper to compare with interior cells
    dof_bnd_out = Grid.dof(sum(abs(DD),2)>eps); 
    % Note, these are all cells! Using dof_bnd_out for temporary storage
    [dof_bnd_in,i_bnd_in] = intersect(dof_bnd_out,dof_in); % find interior cells
    dof_bnd_out(i_bnd_in) = []; % delete interior cells
end

% if length(dof_in) > length(dof_out)
%     % Cheaper to compare with exterior cells
%     dof_bnd_in = Grid.dof(abs(sum(DD,2))>eps); 
%     % Note, these are all cells! Using dof_bnd_in for temporary storage
%     [dof_bnd_out,i_bnd_out] = intersect(dof_bnd_in,dof_out); % find exterior cells
%     dof_bnd_in(i_bnd_out) = []; % delete exterior cells
% else
%     % Cheaper to compare with interior cells
%     dof_bnd_out = Grid.dof(abs(sum(DD,2))>eps); 
%     % Note, these are all cells! Using dof_bnd_out for temporary storage
%     [dof_bnd_in,i_bnd_in] = intersect(dof_bnd_out,dof_in); % find interior cells
%     dof_bnd_out(i_bnd_in) = []; % delete interior cells
% end


