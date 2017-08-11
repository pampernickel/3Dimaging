function mindistances=dist2NearNeighbFastDistMap(CC1, d)
%% This function returns the distance to the nearest neighbour of two 3D array structs created by bwconncomp of equal size.
% Anistropic pixelsizes are possible. This function calls bwdist.m 
 %
 %%

    mindistances=zeros(CC1.NumObjects,1);
    
    %hbar = parfor_progressbar(CC1.NumObjects,'Computing minimal distances');  
    for i=1: CC1.NumObjects
        
        %CC1mask=accumarray(fliplr(CC1.PixelIdxList{i}),true,[nImage,mImage,slPerCH]);
        %CC1mask= full(sparse(pixellist(:, 2), pixellist(:, 1), true,
        %CC1mask = zeros(CC1.ImageSize(1), CC1.ImageSize(2), CC1.ImageSize(3), 'uint8');
        %displayedmask = zeros(CC1.ImageSize(1), CC1.ImageSize(2), CC1.ImageSize(3), 'double');
        pixlist = CC1.PixelIdxList{i};
        dis=min(d(pixlist));
        %d(pixlist) = 1;
        % ismember(labelmatrix(CC1), CC1.PixelIdxList(i));
        
        %d=bwdistsc(CC1mask(:,:,:), pixelsize);

        %% Looking for the smallest distances
        %pixlist = CC2_mask.PixelIdxList{:};
        %%CC1mask(pixlist) = 1;
        %idist=idist+1;
        
        %%
%         workspace
%         whos
%         who
%         size(displayedmask)
%         displayedmask(pixlist)=d(pixlist);
%         size(displayedmask)
%         colormap gray
%         for ii = 1 : CC1.ImageSize(3)
%             imagesc(displayedmask(:,:,ii));
%            %imagesc(d(:,:,ii));
%             pause(.1);
%         end
        %%
        
        
        
        %dis=min(d(pixlist));
        mindistances(i,1)=dis;
        %hbar.iterate(1);
    end
    %close(hbar);   %close progress bar