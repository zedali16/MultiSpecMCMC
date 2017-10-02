function[ai] = matpower(a,alpha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate matrix power
%
%   Input:
%       1) a - a matrix
%       2) alpha - desired power
%   Main Outputs:
%       1) ai - matrix with desired power
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        small = .000001;
        if numel(a) == 1
                ai=a^alpha;
        else
            [p1, ~] = size(a);
            [eve, eva] = eig(a);
            eva = diag(eva);
            eve = eve./(repmat((diag((eve)'*eve).^0.5),1,p1));
            index = 1:p1;
            index = index(eva>small);
            evai = eva;
            evai(index) = (eva(index)).^(alpha);
            ai = eve*diag(evai)*(eve)';
        end
end