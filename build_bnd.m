function [B,N,fn] = build_bnd(Param,Grid,I)
% author: Marc Hesse
% date: 06/09/2015
% Description:
% This function computes the operators and r.h.s vectors for both Dirichlet
% and Neumann boundary conditions. 
% Note:
% N is not created from I the same way B is created from I, because 
% the vector dof_dir contains the columns that must be eliminated rather
% then the colums that are retained in N. If you wanted to do it this way
% you would have to create a new vector
% dof_non_dir = setdiff(dof,dof_dir)
% I suspect that the set-operators are expensive on large vectors, hence
% we simply eliminate the rows.
%
% Input:
% Param = structure containing all information about the physical problem
%         in particular this function needs the fields
%         Param.dof_dir = Nc by 1 column vector containing 
%                         the dof's of the Dirichlet boundary.
%         Param.dof_neu = N by 1 column vector containing 
%                         the dof's of the Neumann boundary.
%         Param.qb      = column vector of prescribed fluxes on Neuman bnd.
% Grid = structure containing all pertinent information about the grid.
% I = identity matrix in the full space
%
% Output:
% B = Nc by N matrix of the Dirichlet constraints
% N = (N-Nc) by N matrix of the nullspace of B
% fn = N by 1 r.h.s. vector of Neuman contributions
%
% Example call:
% >> Grid.xmin = 0; Grid.xmax = 1; Grid.Nx = 10;
% >> Grid = build_grid(Grid);
% >> [D,G,I]=build_ops(Grid);
% >> Param.dof_dir   = Grid.dof_xmin;    % identify cells on Dirichlet bnd
% >> Param.dof_f_dir = Grid.dof_f_xmin;  % identify faces on Dirichlet bnd
% >> Param.dof_neu   = Grid.dof_xmax;    % identify cells on Neumann bnd
% >> Param.dof_f_neu = Grid.dof_f_xmax;  % identify cells on Neumann bnd
% >> Param.qb = 1;                   % set bnd flux
% >> [B,N,fn] = build_bnd(Param,Grid,I);

%% If Neumann condtions are not set explicitly, make them natural
if ~isfield(Param,'dof_neu');   Param.dof_neu   = []; end
if ~isfield(Param,'dof_f_neu'); Param.dof_f_neu = []; end
if ~isfield(Param,'qb');        Param.qb        = []; end

%% Check input format
if isrow(Param.dof_dir)   && length(Param.dof_dir)>1;   error('Param.dof_dir is not a column vector'); end
if isrow(Param.dof_neu)   && length(Param.dof_neu)>1;   error('Param.dof_neu is not a column vector'); end
if isrow(Param.dof_f_dir) && length(Param.dof_f_dir)>1; error('Param.dof_f_dir is a not column vector'); end
if isrow(Param.dof_f_neu) && length(Param.dof_f_neu)>1; error('Param.dof_f_neu is a not column vector'); end
if isfield(Param,'qb') && isrow(Param.qb) && length(Param.qb)>1;        error('Param.qb is not a column vector'); end

%% Check if sufficient BC's are specified
if length(Param.dof_dir) ~= length(Param.g);  error('Insufficient number of constraints specified.'); end
if length(Param.dof_neu) ~= length(Param.qb); error('Insufficient number of boundary fluxes specified.'); end

%% Other checks to implement
% doubling up on constraints
% cell # == face #

%% Dirichlet boundary conditions
B = I(Param.dof_dir,:);
N = I; N(:,Param.dof_dir) = [];

%% Neumann boundary conditions
if isempty(Param.dof_neu)
    fn = spalloc(Grid.N,1,0);                     % allocate sparse zero vector
else
    fn = spalloc(Grid.N,1,length(Param.dof_neu)); % allocate sparse vector
    fn(Param.dof_neu) = Param.qb.*Grid.A(Param.dof_f_neu)./Grid.V(Param.dof_neu);
end
