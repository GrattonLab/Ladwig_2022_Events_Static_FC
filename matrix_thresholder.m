function [matrix r kden] = matrix_thresholder(matrix,threshold,thresholdtype)
%
% Name:matrix_thresholder.m
% $Revision: 1.2 $
% $Date: 2011/03/08 20:30:16 $
%
% jdp 10/10/10
% 
% This script takes a 2D matrix and thresholds it at a given value or edge
% density, and returns a thresholded matrix. It also returns the precise
% edge density and minimum r value of the thresholded matrix (which will
% differ from the pre-specified kden or r value used to threshold the
% matrix!)
% 
% USAGE: [matrix r kden] = matrix_thresholder(matrix,threshold,thresholdtype)
% USAGE: [matrix r kden] = matrix_thresholder(mat,0.8,'kden')
% USAGE: [matrix r kden] = matrix_thresholder(rmat,0.25,'r')
% 
% Thresholdtype is 'r' or 'kden'
% 
% NOTE:
% This script presumes an undirected (symmetric) matrix !!!!
% 
% It also ignores values on the diagonal (should be zero anyways) !!!!

d=size(matrix);
if size(d,1)>2
    error('Matrix_thresholder only works on 2D matrices (with reason!)');
end

numpossibleedges=d(1)*(d(1)-1)/2;

switch thresholdtype
    case 'r'
        matrix(matrix<threshold)=0;
        kden=nnz(triu(matrix,1))/numpossibleedges;
        notzero=(triu(matrix,1)~=0);
        r=min(matrix(notzero));
        if isempty(r)
            fprintf('No edges passing threshold %f!\n',threshold);
            r=threshold;
        end
    case 'kden'
        if (threshold<0) || (threshold>1)
            error('kden needs to be 0<=kden<=1');
        end
        
        if d(1)<20000 % if it is a "smaller" matrix that won't destroy macthunder's RAM
            edgesleft=ceil(threshold*numpossibleedges);
            matrix=triu(matrix,1);
            [v i]=sort(matrix(:),'descend');
            keptedges=zeros(d(1));
            keptedges(1:edgesleft)=1;
            keptedges=logical(keptedges);
            v(~keptedges)=0;
            matrix(i)=v;
            matrix=reshape(matrix,d);
            matrix=max(matrix,matrix');
            r=v(edgesleft);
            kden=edgesleft/numpossibleedges;
        
        else % if this is a behemoth that will destroy macthunder's RAM with sort operations..
            fprintf('matrix_thresholder: Using alternate kden thresholding for large matrices to spare RAM\n');
            [matrix kden r] = alternatethreshold(matrix,threshold);
        end
    case 'abs'
        matrix(abs(matrix)<threshold)=0;
        kden=nnz(triu(matrix,1))/numpossibleedges;
        r=threshold;
    otherwise
        error('matrix_thresholder: need to use ''r'', ''kden'', or ''abs'' as thresholding switches.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [matrix kden threshold] = alternatethreshold(matrix,kden)

% this function avoids the vectortized operations that kill RAM for big matrices

d=size(matrix);
numedges=d(1)*d(1);
repsneeded=ceil(log2(numedges)); % this is the maximum useful number of cycles when bisecting differences of numedges

% check that matrix can meet the desired kden
beginkden=nnz(matrix)/numedges;
if (beginkden<kden)
    fprintf('matrix_thresholder: matrix starts with kden %f, cannot achieve kden of %f\n',beginkden,kden);
    error('Supply a different kden or a different matrix');
end

% obtain upper and lower bounds of the matrix
amin=min(min(nonzeros(matrix)));
amax=max(max(matrix));

for i=1:repsneeded
    
    % perform the next obvious threshold
    tempthresh=(amax+amin)/2; % set the next threshold
    
    fprintf('\tCycle %d/%d, threshold is %f\t',i,repsneeded,tempthresh);
    
    matrix2=(matrix>=tempthresh); % apply the threshold, create logical array
    kden2=nnz(matrix2)/(numedges); % calculate the edge density
    kdiff=kden2-kden; % how far from our target kden are we?
    
    fprintf('kden is %f, offtarget is %f\n',kden2,kdiff);
    
    if kdiff>0 % if the matrix is too dense still we bump the lower threshold up
        amin=tempthresh;
    elseif kdiff<0 % if the matrix is too sparse we drop the upper threshold down
        amax=tempthresh;
    elseif kdiff==0 % and if we hit the nail on the head.. we're done!
        break
    end
end

threshold=tempthresh;
matrix(matrix<threshold)=0;
kden=nnz(matrix)/numedges;

fprintf('Final result is r=%f\tkden=%f\tedges=%d\n',threshold,kden,nnz(matrix));
