%% Automated segmentation of DNA from ilastik segmentation
% moritz.kirschmann@zmb.uzh.ch 2015

close all force;
clear all;

verbose=1;
writeSegmentationResults=1;

%% Loading raw data from single channel images
channels = 3;
inputpath = 'L:\img1\';
outputpath = 'L:\img1\3651._-_3651_Merging001_Merged001\output\';
inputBitdepth = 'uint16';
CellsInCH=2;
BlVesInCH=3;
DepthCorrectionBlVes =1; % In channel
DepthCorrectionCells =1; % In channel
sigmaOfGaussianBlur = [2.0, 1.0, 2.0]; %First Cells, second BlVes in pixels
dilRadiusTissue = 77; % in um
threshold= [100, 60, 10000];%First Cells, second BlVes
minimalSegmentSize= [500, 420, 350];  %First Cells, second BlVes

pixelsize= [0.4541, 0.4541, 1.4995]; % x y z in microns

suff= '.tif' % Suffix of image files
folders = dir(inputpath)

myMap = rand(10000, 3); %Generating random color map
myMap(1,:)= [0,0,0] ;

%Iterate over Folders
for go=1 : numel(folders)   
    go
    folders(go).name
    % Check if folder contains images
    if ((folders(go).isdir) & (~strcmp(folders(go).name, '.')) & (~strcmp(folders(go).name, '..')))
        
        FilePath=strcat(inputpath, folders(go).name);
        FilesInFolder = dir(FilePath);
        % FinalImage=TIFFStack(FileTif);
        
        % Counting number of slices in folder
        numFrames=0;
        for file=1 : numel(FilesInFolder)
            %folders(go).name
            %slice
            ImagePath=strcat(FilePath, filesep ,FilesInFolder(file).name);
            [pathstr,name,ext] = fileparts(ImagePath);
            if  strcmpi(ext,suff)
                numFrames=numFrames+1;
                
                info=imfinfo(ImagePath);
                mImage=info(1).Width;
                nImage=info(1).Height;
            end
        end
        
        numFrames
        slPerCH= numFrames/channels;
        ch=zeros(round(nImage), round(mImage), slPerCH, channels, inputBitdepth);
        
        %% starting read in
        'starting read in'
        for slice=1 : numel(FilesInFolder)
            
            slicePath=strcat(FilePath, filesep ,FilesInFolder(slice).name);
            [pathstr,n,ext] = fileparts(slicePath);
            
            if  strcmpi(ext,suff)
                startIndexZ=regexp(n, '_Z')+2;
                endIndexZ=regexp(n, '_C')-1;
                z=cast(str2num(n(startIndexZ:endIndexZ)), 'uint8');
                
                startIndexC=regexp(n, '_C')+2;
                endIndexC=regexp(n, '_T')-1;
                if numel(endIndexC)>1
                    endIndexC=endIndexC(end);
                end
                c=cast(str2num(n(startIndexC:endIndexC)), 'uint8');
                
                ch(:,:,z+1,c+1)=imread(slicePath);
            end
        end
        
        %% show original data
                if verbose == 1
                    for ii = 1 : slPerCH
        
                        imagesc(ch(:,:,ii,1));
                        pause(.01);
                    end
                    for ii = 1 : slPerCH
        
                        imagesc(ch(:,:,ii,2));
                        pause(.01);
                    end
                    for ii = 1 : slPerCH
        
                        imagesc(ch(:,:,ii,3));
                        pause(.01);
                    end
                end
        
        %% Gaussian blurring the first two channels
        h = fspecial('gaussian', 7, sigmaOfGaussianBlur(1));
        ch(:,:,:,CellsInCH) = imfilter(ch(:,:,:,CellsInCH), h);
        
        h = fspecial('gaussian', 7, sigmaOfGaussianBlur(2));
        ch(:,:,:,BlVesInCH) = imfilter(ch(:,:,:,BlVesInCH), h);
        
        if DepthCorrectionBlVes==1
            for ii = 1 : slPerCH
                intensityAverage(ii)=mean(mean(ch(:,:,ii,BlVesInCH)));
                pause(.01);
            end
            DCBlVes=ch(:,:,:,BlVesInCH);
            for ii = 2 : slPerCH
                DCBlVes(:,:,ii)=round(ch(:,:,ii,BlVesInCH)/intensityAverage(ii)*intensityAverage(1));
                pause(.01);
            end
            if verbose == 1
                for ii = 1 : slPerCH
                    colormap default;
                    image(DCBlVes(:,:,ii));
                    pause(.1);
                end
            end
            
            ch(:,:,:,BlVesInCH)=DCBlVes;
        end
        
        
        if DepthCorrectionCells==1
            for ii = 1 : slPerCH
                intensityAverage(ii)=mean(mean(ch(:,:,ii,CellsInCH)));
                pause(.01);
            end
            DCCells=ch(:,:,:,CellsInCH);
            for ii = 2 : slPerCH
                DCCells(:,:,ii)=round(ch(:,:,ii,CellsInCH)/intensityAverage(ii)*intensityAverage(1));
                pause(.01);
            end
            if verbose == 1
                for ii = 1 : slPerCH
                    colormap default;
                    image(DCCells(:,:,ii));
                    pause(.1);
                end
            end
            
            ch(:,:,:,CellsInCH)=DCCells;
        end
        
        %% thresholding
        
        seg(:,:,:,CellsInCH) = ch(:,:,:,CellsInCH) >= threshold(1);
        seg(:,:,:,BlVesInCH) = ch(:,:,:,BlVesInCH) >= threshold(2);
        
        
        %% finding connected components
        CCCells = bwconncomp(seg(:,:,:,CellsInCH), 18);        
        CCBlVes = bwconncomp(seg(:,:,:,BlVesInCH), 18);
        
        %% finding number of voxels per connected component to filter out small objects
        CCCells_props = regionprops(CCCells, 'Area');
        idx = find([CCCells_props.Area] > minimalSegmentSize(1));
        filteredseg1 = ismember(labelmatrix(CCCells), idx);
        
        CCBlVes_props = regionprops(CCBlVes, 'Area');
        idx = find([CCBlVes_props.Area] > minimalSegmentSize(2));
        filteredseg2 = ismember(labelmatrix(CCBlVes), idx);
        
        %% Show sizefiltered CC
        if verbose == 1
            for ii = 1 : slPerCH
                imagesc(filteredseg1(:,:,ii));
                pause(.01);
            end
        end
        if verbose == 1
            for ii = 1 : slPerCH
                imagesc(filteredseg2(:,:,ii));
                pause(.01);
            end
        end
        
        
        %% finding connected components of the remaining segments AND finding Centers
        CC1 = bwconncomp(filteredseg1, 18);
        %%%%%%%%%%%%%%
        CCCells = bwconncomp(filteredseg1(:,:,:), 18);
        
        
        %% finding number of voxels per connected component to filter out small objects
        CCCells_props = regionprops(CCCells, 'Area');
        
        
        
        %%%%%%%%%%%%%%
        
        CellCenters = regionprops(CC1, 'Centroid');
        
        CC2 = bwconncomp(filteredseg2, 18);
        
        CC1mask = zeros([nImage,mImage,slPerCH], 'uint16');
        for i=1: CC1.NumObjects
            pixlist = CC1.PixelIdxList{i};
            CC1mask(pixlist) = i;
        end
        
        
        CC2mask = zeros([nImage,mImage,slPerCH], 'uint16');
        for i=1: CC2.NumObjects
            pixlist = CC2.PixelIdxList{i};
            CC2mask(pixlist) = i;
        end
        CC2maskBin=CC2mask>0;
        %% Save        
        
        if verbose == 1
            colormap(myMap);
            CLim = [0, CC1.NumObjects]
            for ii = 1 : slPerCH
                image(CC1mask(:,:,ii));
                pause(.1);
            end
            for ii = 1 : slPerCH
                image(CC2mask(:,:,ii));
                pause(.1);
            end
        end
        
        
        if writeSegmentationResults == 1
            
            
            for ii = 1 : slPerCH
                imwrite(CC1mask(:,:,ii), strcat(outputpath, folders(go).name, num2str(ii),'_segCells.tif'));
                
            end
            
            
            for ii = 1 : slPerCH
                imwrite(CC2mask(:,:,ii), strcat(outputpath, folders(go).name, num2str(ii),'_segBV.tif'));
                
            end
            
            
            for ii = 1 : slPerCH
                imwrite(ch(:,:,ii,CellsInCH), strcat(outputpath, folders(go).name,num2str(ii), '_Cells_pre_thres.tif'));
                
            end
            
            
            
            for ii = 1 : slPerCH
                imwrite(ch(:,:,ii,BlVesInCH), strcat(outputpath, folders(go).name,num2str(ii), '_BlVes_pre_thres.tif'));
            end
            
        end
                
        'Postfilter'
        CC1.NumObjects
        CC2.NumObjects
        folders(go).nCells=CC1.NumObjects;
        folders(go).BlVes =CC2.NumObjects;
        
        %clear FinalImage seg1 seg2 seg filteredseg1 filteredseg2 filteredseg3 ch1 ch2 ch3 ch
        
        
        %% Finding shortest distance of Channel 1 objects to Channel 2 objects
        tic;
        d=bwdistsc(CC2maskBin, pixelsize);
        mindistances=dist2NearNeighbFastDistMap(CC1, d);
        folders(go).mindists=mindistances
        
        toc
        
        %% Calculating Distances of randomly distributed Balls in the Cells Channels
        
        
        % To only consider areas within tissue a mask is created by finding
        % all foreground in both channels and imclose this mask to fill all
        % the gaps
        allMask=or(CC1mask>0,CC2mask>0);
        
        
        hd = fspecial('disk', round(dilRadiusTissue/pixelsize(1))) >0;
        he = fspecial('disk', round(dilRadiusTissue/pixelsize(1)/2)) >0;
        %ch(:,:,:,CellsInCH) = imfilter(ch(:,:,:,CellsInCH), h);
        tic;
        tissueMask=allMask;
        parfor ii = 1 : slPerCH
            %tissueMask(:,:,ii) = bwmorph(allMask(:,:,ii), h ,'conv');
            tissueMask(:,:,ii) = imdilate(allMask(:,:,ii), hd);
            tissueMask(:,:,ii) = imerode(tissueMask(:,:,ii), he);
        end
        
        
        FilledBleVS =(CC2mask>0);
        h = (fspecial('disk', 30)) >0;
        parfor ii = 1 : slPerCH
            %tissueMask(:,:,ii) = bwmorph(allMask(:,:,ii), h ,'conv');
            FilledBlVesMask(:,:,ii) = imclose(FilledBleVS(:,:,ii), h);
            
        end
        
        toc
        
        tissueWithoutFilledBlVes=tissueMask-FilledBlVesMask;
        
        
        if verbose == 1
            %colormap(myMap);
            %CLim = [0, CC1.NumObjects]
            for ii = 1 : slPerCH
                imagesc(uint8(allMask(:,:,ii)));
                pause(.2);
            end
            for ii = 1 : slPerCH
                imagesc(uint8(tissueMask(:,:,ii)));
                pause(.3);
            end
            for ii = 1 : slPerCH
                ii
                imagesc(uint8(tissueMask(:,:,ii))+uint8(allMask(:,:,ii))+uint8(FilledBlVesMask(:,:,ii))+uint8(ch(:,:,ii,BlVesInCH)));
                pause(.2);
            end
            for ii = 1 : slPerCH
                ii
                imagesc(uint8(tissueWithoutFilledBlVes(:,:,ii)));
                pause(.2);
            end
        end
        if writeSegmentationResults == 1
            
            imwrite(uint16(tissueMask(:,:,1)), strcat(outputpath, folders(go).name, '_tissuemask.tif'));
            for ii = 2 : slPerCH
                imwrite(uint16(tissueMask(:,:,ii)), strcat(outputpath, folders(go).name, '_tissuemask.tif') ,'WriteMode','append');
            end
            
            imwrite(uint16(FilledBlVesMask(:,:,1)), strcat(outputpath, folders(go).name, '_FilledBV.tif'));
            for ii = 2 : slPerCH
                imwrite(uint16(FilledBlVesMask(:,:,ii)), strcat(outputpath, folders(go).name, '_FilledBV.tif') ,'WriteMode','append');
            end
            
        end
        
        %% Randomly distributed cells within tissue
        nRandomCells=CC1.NumObjects;
        RandomCellsInTissue = nan(3,nRandomCells);
        RandomCellsInTissueMask=zeros(nImage,mImage,slPerCH, 'uint16');
        RCITMlin=CCCells;
        for  j=1 : nRandomCells
            while 1
                r=random('unif',0,1,3,1);
                
                RandomCellCoords = (floor((r.').*[nImage, mImage,slPerCH])+[1,1,1]);
                %tissueMask(RandomCellCoords(1),RandomCellCoords(2), RandomCellCoords(3));
                if tissueMask(RandomCellCoords(1),RandomCellCoords(2), RandomCellCoords(3))
                    j
                    RandomCellCoords(:)
                    RandomCellsInTissue(:,j)=[RandomCellCoords(1),RandomCellCoords(2), RandomCellCoords(3)];
                    %RandomCellsInTissueMask=RandomCellsInTissueMask+uint16(make_filled_spheroid(RandomCellCoords(1),RandomCellCoords(2), RandomCellCoords(3),radiusOfCellInMicrons/pixelsize(1),radiusOfCellInMicrons/pixelsize(2),radiusOfCellInMicrons/pixelsize(3), nImage, mImage, slPerCH));
                    Cellmask = zeros([nImage,mImage,slPerCH], 'uint16');
                    
                    pixlist = CC1.PixelIdxList{j};
                    Cellmask(pixlist) = j;
                    Cellmask=circshift(Cellmask, -round(CellCenters(j).Centroid(1))+RandomCellCoords(2), 2 );
                    Cellmask=circshift(Cellmask, -round(CellCenters(j).Centroid(2))+RandomCellCoords(1), 1 );
                    Cellmask=circshift(Cellmask, -round(CellCenters(j).Centroid(3))+RandomCellCoords(3), 3 );
                    
                    
                    
                    % CellmaskCC=bwconncomp((Cellmask(:,:,:)>0), 18);
                    %find
                    RCITMlin.PixelIdxList{j}=find(Cellmask>0);
                    RandomCellsInTissueMask=RandomCellsInTissueMask+Cellmask;
                    break
                    
                end
            end
            
        end
        if verbose == 1
            %colormap(myMap);
        
            for ii = 1 : slPerCH
                %ii
                imagesc(uint8(RandomCellsInTissueMask(:,:,ii))+uint8(tissueMask(:,:,ii))+uint8(CC1mask(:,:,ii)));
                pause(.2);
            end
            for ii = 1 : slPerCH
                %ii
                imagesc(uint8(RandomCellsInTissueMask(:,:,ii))+uint8(tissueMask(:,:,ii)));
                pause(.2);
            end

        end
        if writeSegmentationResults == 1
            imwrite(uint16(RandomCellsInTissueMask(:,:,1)), strcat(outputpath, folders(go).name, '_RandomCellsInTissue.tif'));
            for ii = 2 : slPerCH
                imwrite(uint16(RandomCellsInTissueMask(:,:,ii)), strcat(outputpath, folders(go).name, '_RandomCellsInTissue.tif') ,'WriteMode','append');
            end
        end
                
        tic;
        
        mindistancesRandomInTissue=dist2NearNeighbFastDistMap(RCITMlin, d);
        folders(go).mindistsRCITM=mindistancesRandomInTissue;
        toc
        
        
        %% Randomly distributed cells within tissue not inside BlVes
        nRandomCells=CC1.NumObjects;
        RandomCellsInTissueWithoutBlVes = nan(3,nRandomCells);
        RandomCellsInTissueWithoutBlVesMask=zeros(nImage,mImage,slPerCH, 'uint16');
        RCITMWBlin=CCCells;
        for  j=1 : nRandomCells
            while 1
                r=random('unif',0,1,3,1);
                
                RandomCellCoords = (floor((r.').*[nImage, mImage,slPerCH])+[1,1,1]);
                if tissueWithoutFilledBlVes(RandomCellCoords(1),RandomCellCoords(2), RandomCellCoords(3))
                    j
                    RandomCellCoords(:)
                    RandomCellsInTissue(:,j)=[RandomCellCoords(1),RandomCellCoords(2), RandomCellCoords(3)];
                    %RandomCellsInTissueMask=RandomCellsInTissueMask+uint16(make_filled_spheroid(RandomCellCoords(1),RandomCellCoords(2), RandomCellCoords(3),radiusOfCellInMicrons/pixelsize(1),radiusOfCellInMicrons/pixelsize(2),radiusOfCellInMicrons/pixelsize(3), nImage, mImage, slPerCH));
                    Cellmask = zeros([nImage,mImage,slPerCH], 'uint16');
                    
                    pixlist = CC1.PixelIdxList{j};
                    Cellmask(pixlist) = j;
                    Cellmask=circshift(Cellmask, -round(CellCenters(j).Centroid(1))+RandomCellCoords(2), 2 );
                    Cellmask=circshift(Cellmask, -round(CellCenters(j).Centroid(2))+RandomCellCoords(1), 1 );
                    Cellmask=circshift(Cellmask, -round(CellCenters(j).Centroid(3))+RandomCellCoords(3), 3 );
                    
                    RCITMWBlin.PixelIdxList{j}=find(Cellmask>0);
                    RandomCellsInTissueWithoutBlVesMask=RandomCellsInTissueWithoutBlVesMask+Cellmask;
                    
                    break
                    
                end
            end
            
        end
        if verbose == 1
            for ii = 1 : slPerCH
                imagesc(uint8(RandomCellsInTissueWithoutBlVesMask(:,:,ii))+uint8(tissueWithoutFilledBlVes(:,:,ii))+uint8(CC1mask(:,:,ii)));
                pause(.2);
            end
            for ii = 1 : slPerCH
                imagesc(uint8(RandomCellsInTissueWithoutBlVesMask(:,:,ii))+uint8(tissueWithoutFilledBlVes(:,:,ii))+uint8(tissueMask(:,:,ii)));
                pause(.2);
            end
        end
        
        if writeSegmentationResults == 1
            
            for ii = 1 : slPerCH
                imwrite(uint16(RandomCellsInTissueWithoutBlVesMask(:,:,ii)), strcat(outputpath, folders(go).name, num2str(ii), '_RandomCellsInTissueWithoutBV.tif'));
               
            end
        end
        %% Finding shortest distance of Randomized Cells to
        
        
        tic;
        
        mindistancesRandomInTissueWithoutBlVes=dist2NearNeighbFastDistMap(RCITMWBlin, d);
        folders(go).mindistsRCITMWBl=mindistancesRandomInTissueWithoutBlVes;
        %hist(folders(go).mindistsRCITMWBl);
        toc
                
        %% Writing  results
        ' Writing  results'
        fid = fopen(strcat(outputpath,'Matlab_results_', folders(go).name,'.txt'),'wt');
        
        % if ((folders(go).isdir) & (~strcmp(folders(go).name, '.')) & (~strcmp(folders(go).name, '..')))
        fprintf(fid, '%s\t', folders(go).name);
        fprintf(fid, '\n' );
        fprintf(fid, '%u\t', folders(go).nCells);
        fprintf(fid, '\n' );
        fprintf(fid, '%u\t', folders(go).BlVes);
        fprintf(fid, '\n' );
        
        
        fprintf(fid, 'Mean distances of Cells \n' );
        fprintf(fid, '%u\t', mean(folders(go).mindists(:)));
        fprintf(fid, '\n' );
        
        fprintf(fid, 'All distances of Cells \n' );
        fprintf(fid, '%u\t', folders(go).mindists(:));
        fprintf(fid, '\n' );
        fprintf(fid, '\n' );
        
        fprintf(fid, 'Mean distances of randomized Cells within Tissue \n' );
        fprintf(fid, '%u\t', mean(folders(go).mindistsRCITM(:)));
        fprintf(fid, '\n' );
        
        fprintf(fid, 'All distances of randomized Cells within Tissue  \n' );
        fprintf(fid, '%u\t', folders(go).mindistsRCITM(:));
        fprintf(fid, '\n' );
        fprintf(fid, '\n' );
        
        
        fprintf(fid, 'Mean distances of randomized Cells within Tissue without Bl. Ves \n' );
        fprintf(fid, '%u\t', mean(folders(go).mindistsRCITMWBl(:)));
        fprintf(fid, '\n' );
        
        fprintf(fid, 'All distances of randomized Cells within Tissue without Bl. Ves \n' );
        fprintf(fid, '%u\t', folders(go).mindistsRCITMWBl(:));
        fprintf(fid, '\n' );
        
        fclose(fid);
        
        clear FinalImage seg1 seg2 seg filteredseg1 filteredseg2 filteredseg3 ch1 ch2 ch3 ch slPerCH FilledBleVS FilledBlVesMask Cellmask RandomCellsInTissueMask RandomCellsInTissueWithoutBlVesMask tissueWithoutFilledBlVes CCCells CCCells_props CC1
        clear CC1 CC1mask CC2 CC2mask CC2maskBin CCBlVes CCBlVes_props CCCells CCCells_props CLim CellCenters Cellmask DCBlVes DCCells FilesInFolder FilledBlVesMask FilledBleVS ImagePath
        clear RCITMWBlin RCITMlin RandomCellCoords RandomCellsInTissue RandomCellsInTissueMask RandomCellsInTissueWithoutBlVes allMask
        clear c ch d endIndexC endIndexZ  filteredseg1  filteredseg h hd he intensityAverage mindistances mindistancesRandomInTissue mindistancesRandomInTissueWithoutBlVes
        clear n nImage nRandomCells name numFrames pathstr pixlist r seg slPerCH slice slicePath startIndexC startIndexZ tissueMask tissueWithoutFilledBlVes
        
    end
end


%% Writing textfile
%     fprintf(fid, '%s\t', 'area in pixels which was considered for counting');
%     fprintf(fid, '\n' );
%     fprintf(fid, '%u\t', pararea(i));
%     fprintf(fid, '\n' );
%     fprintf(fid, '%s\t', 'density in objects per squarepixels');
%     fprintf(fid, '\n' );
%     fprintf(fid, '%u\t', parcounts(:,i)/pararea(i));
%     fprintf(fid, '\n' );
%     fprintf(fid, '\n' );



