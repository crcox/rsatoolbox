function [varargout] = fMRISearchlight(fullBrainVols, binaryMasks_nS, models, betaCorrespondence, userOptions)
%
% fMRISearchlight is a function which takes some full brain volumes of data,
% some binary masks and some models and perfoms a searchlight in the data within
% each mask, matching to each of the models.  Saved are native-space r-maps for
% each model.
%
% [rMaps_sS, maskedSmoothedRMaps_sS, searchlightRDMs[, rMaps_nS, nMaps_nS] =]
%                                 fMRISearchlight(fullBrainVols,
%                                                 binaryMasks_nS,
%                                                 models,
%                                                 betaCorrespondence,
%                                                 userOptions)
%
%       fullBrainVols --- The unmasked beta (or t) images.
%               fullBrainVols.(subject) is a [nVoxel nCondition nSession]-sized
%               matrix. The order of the voxels is that given by reshape or (:).
%
%        binaryMasks_nS --- The native- (subject-) space masks.
%               binaryMasks_nS.(subject).(mask) is a [x y z]-sized binary matrix
%               (the same size as the native-space 3D beta images.
%
%        models --- A stack of model RDMs in a structure.
%               models is a [1 nModels] structure with fields:
%                       RDM
%                       name
%                       color
%
%        betaCorrespondence --- The array of beta filenames.
%               betas(condition, session).identifier is a string which referrs
%               to the filename (not including path) of the SPM beta image. (Or,
%               if not using SPM, just something, as it's used to determine the
%               number of conditions and sessions.)
%               Alternatively, this can be the string 'SPM', in which case the
%               SPM metadata will be used to infer this information, provided
%               that userOptions.conditionLabels is set, and the condition
%               labels are the same as those used in SPM.
%
%
%        userOptions --- The options struct.
%                userOptions.analysisName
%                        A string which is prepended to the saved files.
%                userOptions.rootPath
%                        A string describing the root path where files will be
%                        saved (inside created directories).
%                userOptions.subjectNames
%                        A cell array containing strings identifying the subject
%                        names. Defaults to the fieldnames in fullBrainVols.
%                userOptions.maskNames
%                        A cell array containing strings identifying the mask
%                        names. Defaults to the fieldnames of the first subject
%                        of binaryMasks_nS.
%                userOptions.voxelSize
%                        A tripple consisting of the [x y z] dimensions of each
%                        voxel in mm.
%                userOptions.structuralsPath
%                        A string which contains the absolute path to the
%                        location of the structural images and the normalisation
%                        warp definition file. It can contain the following
%                        wildcards which would be replaced as indicated:
%                                [[subjectName]]
%                                        To be replaced with the name of each
%                                        subject where appropriate.
%                userOptions.saveSearchlightRDMs
%                        A flag that indicates whether to store the RDM for
%                        every searchlight. To be clear, this means to store a
%                        nCondtion x nCondition matrix for every voxel, for
%                        every subject, for every model RDM. This is an insane
%                        amount of memory consumption. Use at your own risk.
%                        Default is false.
%
% The following files are saved by this function:
%        userOptions.rootPath/Maps/
%                userOptions.analysisName_fMRISearchlight_Maps.mat
%                        Contains the searchlight statistical maps in struct so
%                        that rMaps_nS.(modelName).(subject).(maskName),
%                        rMaps_sS.(modelName).(subject).(maskName),
%                        maskedSmoothedRMaps_sS.(modelName).(subject).(maskName)
%                        and nMaps_nS.(modelName).(subject).(maskName) contain
%                        the appropriate data.
%        userOptions.rootPath/RDMs/
%                userOptions.analysisName_fMRISearchlight_RDMs.mat
%                        Contains the RDMs for each searchlight so that
%                        searchlightRDMs.(subject)(:, :, x, y, z) is the RDM.
%                        This file is only produced if
%                        userOptions.saveSearchlightRDMs is true.
%        userOptions.rootPath/Details/
%                userOptions.analysisName_fMRISearchlight_Details.mat
%                        Contains the userOptions for this execution of the
%                        function and a timestamp.
%
% Cai Wingfield 2-2010, 3-2010
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

returnHere = pwd; % We'll come back here later

%% Set defaults and check options struct
if ~isfield(userOptions, 'analysisName'), error('fMRISearchlight:NoAnalysisName', 'analysisName must be set. See help'); end%if
if ~isfield(userOptions, 'rootPath'), error('fMRISearchlight:NoRootPath', 'rootPath must be set. See help'); end%if
userOptions = setIfUnset(userOptions, 'subjectNames', fieldnames(fullBrainVols));
userOptions = setIfUnset(userOptions, 'maskNames', fieldnames(binaryMasks_nS.(userOptions.subjectNames{1})));
userOptions = setIfUnset(userOptions, 'saveSearchlightRDMs', false);
if ~isfield(userOptions, 'voxelSize'), error('fMRISearchlight:NoVoxelSize', 'voxelSize must be set. See help'); end%if

% The analysisName will be used to label the files which are eventually saved.
mapsFilename = [userOptions.analysisName, '_fMRISearchlight_Maps.mat'];
RDMsFilename = [userOptions.analysisName, '_fMRISearchlight_RDMs.mat'];
DetailsFilename = [userOptions.analysisName, '_fMRISearchlight_Details.mat'];

promptOptions.functionCaller = 'fMRISearchlight';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'Maps', mapsFilename);
promptOptions.checkFiles(2).address = fullfile(userOptions.rootPath, 'Details', DetailsFilename);

overwriteFlag = overwritePrompt(userOptions, promptOptions);

if overwriteFlag
	
	% Data
	nSubjects = numel(userOptions.subjectNames);
	nMasks = numel(userOptions.maskNames);
	
	searchlightOptions.monitor = false;
	searchlightOptions.fisher = true;
	
	warpFlags.interp = 1;
	warpFlags.wrap = [0 0 0];
	warpFlags.vox = userOptions.voxelSize; % [3 3 3.75]
	warpFlags.bb = [-78 -112 -50; 78 76 85];
	warpFlags.preserve = 0;
	
	fprintf('Shining RSA searchlights...\n');

	for subjectNumber = 1:nSubjects % and for each subject...
	
		tic;%1

		fprintf(['\t...in the brain of subject ' num2str(subjectNumber) ' of ' num2str(nSubjects)]);

		% Figure out which subject this is
		subject = userOptions.subjectNames{subjectNumber};
		
		if ischar(betaCorrespondence) && strcmpi(betaCorrespondence, 'SPM')
			betas = getDataFromSPM(userOptions);
		else
			betas = betaCorrespondence;
		end%if:SPM

		searchlightOptions.nSessions = size(betas, 1);
		searchlightOptions.nConditions = size(betas, 2);

		readFile = replaceWildcards(userOptions.betaPath, '[[subjectName]]', subject, '[[betaIdentifier]]', betas(1,1).identifier);
		subjectMetadataStruct = spm_vol(readFile);
%		subjectMetadataStruct = spawnSPMStruct;
		
		for maskNumber = 1:nMasks % For each mask...
	
			% Get the mask
			maskName = userOptions.maskNames{maskNumber};
			mask = binaryMasks_nS.(subject).(maskName);
			
			% Full brain data volume to perform searchlight on
			singleSubjectVols = fullBrainVols.(subject);

			% Do the searchlight! ZOMG, this takes a while...
      if userOptions.saveSearchlightRDMs
        [rs, ps, ns, searchlightRDMs.(subject)] = rsa.fmri.searchlightMapping_fMRI(singleSubjectVols, models, mask, userOptions, searchlightOptions); % ps are from linear correlation p-values, and so aren't too useful here.
      else
        [rs, ps, ns] = rsa.fmri.searchlightMapping_fMRI(singleSubjectVols, models, mask, userOptions, searchlightOptions); % ps are from linear correlation p-values, and so aren't too useful here.
      end
			
			nMaps_nS.(subject).(maskName) = ns(:,:,:); % How many voxels contributed to the searchlight centred at each point. (Those with n==1 are excluded because the results aren't multivariate.)
			
			for modelNumber = 1:numel(models)
				
				modelName = spacesToUnderscores(models(modelNumber).name);
				
				% Store results in indexed volumes
				rMaps_nS.(modelName).(subject).(maskName) = rs(:,:,:,modelNumber); % r-values for correlation with each model
				
				%% Save native space version
				
				% Write the native-space r-map to a file
				rMapMetadataStruct_nS = subjectMetadataStruct;
				rMapMetadataStruct_nS.fname = fullfile(userOptions.rootPath, 'Maps', [userOptions.analysisName '_rMap_' maskName '_' modelName '_' subject '.img']);
				rMapMetadataStruct_nS.descrip =  'R-map';
				rMapMetadataStruct_nS.dim = size(rMaps_nS.(modelName).(subject).(maskName));
				
				gotoDir(userOptions.rootPath, 'Maps');
				
				rsa.spm.spm_write_vol(rMapMetadataStruct_nS, rMaps_nS.(modelName).(subject).(maskName));
				
				if isfield(userOptions, 'structuralsPath')

					% Write the native-space mask to a file
					maskMetadataStruct_nS = subjectMetadataStruct;
					maskMetadataStruct_nS.fname = fullfile(userOptions.rootPath, 'Maps', [userOptions.analysisName '_nativeSpaceMask_' maskName '_' modelName '_' subject '.img']);
					maskMetadataStruct_nS.descrip =  'Native space mask';
					maskMetadataStruct_nS.dim = size(mask);
					
					rsa.spm.spm_write_vol(maskMetadataStruct_nS, mask);

					% Load in common space warp definition
% 					wildFiles = replaceWildcards(fullfile(userOptions.structuralsPath, ['*' subject '*_seg_sn.mat']), '[[subjectName]]', subject);
                    wildFiles = replaceWildcards(fullfile(userOptions.structuralsPath, ['*_sn.mat']), '[[subjectName]]', subject);
					matchingFiles = dir(wildFiles);
					warpDefFilename = replaceWildcards(fullfile(userOptions.structuralsPath, matchingFiles(1).name), '[[subjectName]]', subject);

					% Warp and write common space r-maps to disk
					spm_write_sn(rMapMetadataStruct_nS,warpDefFilename,warpFlags);

					% Warp and write common space masks to disk
					spm_write_sn(maskMetadataStruct_nS,warpDefFilename,warpFlags);

					% Now read them back in

					% Where are they?
					[warpedPath_rMap, warpedFile_rMap, warpedExt_rMap] = fileparts(rMapMetadataStruct_nS.fname);
					[warpedPath_mask, warpedFile_mask, warpedExt_mask] = fileparts(maskMetadataStruct_nS.fname);

					% Warped versions are prefixed with 'w'
					warpedFile_rMap = ['w' warpedFile_rMap];
					warpedFile_mask = ['w' warpedFile_mask];

					% Read them from the disk
					rMaps_sS.(modelName).(subject).(maskName) = spm_read_vols(spm_vol(fullfile(warpedPath_rMap, [warpedFile_rMap warpedExt_rMap]))); % sS for standard space
					mask_sS = spm_read_vols(spm_vol(fullfile(warpedPath_mask, [warpedFile_mask warpedExt_mask])));

					% Fix the normalisation of the mask
					maskMetadataStruct_sS = spm_vol(fullfile(warpedPath_rMap, [warpedFile_rMap warpedExt_rMap]));
					maskMetadataStruct_sS.fname = fullfile(userOptions.rootPath, 'Maps', [userOptions.analysisName '_commonSpaceMask_' maskName '_' modelName '_' subject '.img']);
					maskMetadataStruct_sS.descrip =  'Common space mask';
					maskMetadataStruct_sS.dim = size(mask_sS);

					maskThreshold = 0.01;
					mask_sS(mask_sS < maskThreshold) = 0;
					mask_sS(isnan(mask_sS)) = 0;
					
					maskMetadataStruct_sS.dim = size(mask_sS);
					
					rsa.spm.spm_write_vol(maskMetadataStruct_sS, mask_sS);

					% Smooth the normalised data

					% Smoothed versions are prefixed with 's'
					smoothedWarpedFile_rMap = ['s' warpedFile_rMap];

					% Smooth it
					smoothingKernel_fwhm = [10 10 10];
					spm_smooth(fullfile(warpedPath_rMap, [warpedFile_rMap warpedExt_rMap]), fullfile(warpedPath_rMap, [smoothedWarpedFile_rMap warpedExt_rMap]), smoothingKernel_fwhm);

					% Read it back in
					smoothedDataMetadataStruct = spm_vol(fullfile(warpedPath_rMap, [smoothedWarpedFile_rMap warpedExt_rMap]));
					smoothedData = spm_read_vols(smoothedDataMetadataStruct);

					% Mask the smoothed data by the sS mask
					maskedData = smoothedData;
					maskedData(mask_sS == 0) = NaN;
					maskedSmoothedRMaps_sS.(modelName).(subject).(maskName) = maskedData;

					% Write it back to disk
					maskedDataMetadataStruct_nS = smoothedDataMetadataStruct;
					maskedDataMetadataStruct_nS.fname = fullfile(userOptions.rootPath, 'Maps', ['msw' userOptions.analysisName '_rMap_' maskName '_' modelName '_' subject '.img']); % 'msw' for 'masked, smoothed, warped'
					maskedDataMetadataStruct_nS.descrip =  'Masked smoothed normalised data';
					maskedDataMetadataStruct_nS.dim = size(maskedData);
					
					rsa.spm.spm_write_vol(maskedDataMetadataStruct_nS, maskedData);
					
				end%if:structuralsPath

			end%for:models
			
			clear fullBrainVolumes rs ps ns;
			
			fprintf(':');

		end%for:maskNumber
		
		t = toc;%1
		fprintf(' [%ds]\n'), ceil(t));
		
	end%for:subjectNumber

	%% Save relevant info

	timeStamp = datestr(now);

	fprintf('Saving searchlight maps to %s' fullfile(userOptions.rootPath, 'Maps', mapsFilename) '\n']);
	gotoDir(userOptions.rootPath, 'Maps');
	if isfield(userOptions, 'structuralsPath')
		save(mapsFilename, 'rMaps_nS', 'rMaps_sS', 'maskedSmoothedRMaps_sS', 'nMaps_nS');
	else
		save(mapsFilename, 'rMaps_nS', 'nMaps_nS');
	end%if
	
	fprintf(['Saving RDMs to ' fullfile(userOptions.rootPath, 'RDMs', RDMsFilename) '\n']);
	gotoDir(userOptions.rootPath, 'RDMs');
  if userOptions.saveSearchlightRDMs
    save(RDMsFilename, 'searchlightRDMs');
  end
	
	fprintf(['Saving Details to ' fullfile(userOptions.rootPath, 'Details', DetailsFilename) '\n']);
	gotoDir(userOptions.rootPath, 'Details');
	save(DetailsFilename, 'timeStamp', 'userOptions');
	
else
	fprintf(['Loading previously saved maps from ' fullfile(userOptions.rootPath, 'Maps', mapsFilename) '...\n']);
	load(fullfile(userOptions.rootPath, 'Maps', mapsFilename));
	fprintf(['Loading previously saved RDMs from ' fullfile(userOptions.rootPath, 'RDMs', RDMsFilename) '...\n']);
	load(fullfile(userOptions.rootPath, 'RDMs', RDMsFilename));
end%if

if nargout == 3
	varargout{1} = rMaps_sS;
	varargout{2} = maskedSmoothedRMaps_sS;
  if userOptions.saveSearchlightRDMs
    varargout{3} = searchlightRDMs;
  end
elseif nargout == 5
	varargout{1} = rMaps_sS;
	varargout{2} = maskedSmoothedRMaps_sS;
  if userOptions.saveSearchlightRDMs
    varargout{3} = searchlightRDMs;
  end
	varargout{4} = rMaps_nS;
	varargout{5} = nMaps_nS;
elseif nargout > 0
	error('0, 3 or 5 arguments out, please.');
end%if:nargout

cd(returnHere); % And go back to where you started

end%function
