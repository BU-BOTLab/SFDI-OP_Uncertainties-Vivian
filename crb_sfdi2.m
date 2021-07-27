function [CRB, sqrtCrbMua, sqrtCrbMusp, rho] = crb_sfdi2(DRC, fxFreqs, op, RdErr, crbFlag, mode)
%
% Code to generate CRBs and/or rCRBs for SFDI. 
% See V.Pera et al., "Optical property uncertainty estimates for spatial frequency domain imaging,"
% to appear in Biomed. Opt. Express. All references to sections
% and equations are from this paper.
%
% For best best results, diffuse reflectance error model (RdErr) should be determined for each 
% instrument and given set of operating and data processing parameters (see paper). We have provided 
% our RdErr for demonstration purposes only.
%
% INPUTS:
% DRC = structure containing the diffuse reflectance cube (Sec. 3.1.2) with the
%   following fields: 
%   [DRC.Mua, DRC.Musp] = meshgrid(1e-4:1e-4:.25, .1:1e-3:3),
%   where Mua = 1e-4:1e-4:.25 (mm^-1) and Musp = .1:1e-3:3 (mm^-1), 
%   mua = absorption coefficient and musp = reduced scattering coefficient.
%   DRC.fx = nx1 vector of spatial frequencies (mm^-1)
%   DRC.M = diffuse reflectance for each combination of mua, musp, and fx specified
%   above.
%
% fxFreqs = nx1 cell array, with each row specifying a frequency
%   combination for which to compute CRB/rCRB.
%
% op = nx2 matrix of mua, musp pairs for which to compute CRB/rCRB.
%
% RdErr = diffuse reflectance error model (Sec. 2); nx1 vector of
%   fractional errors (sigma over Rd) for each spatial frequency in
%   DRC.fx.
%
% crbFlag = 0: compute reduced CRB (rCRB);
%         = 1: compute CRB.
%
% mode = 1: compare several frequency combinations for several pairs of optical properties;
%           generate bar plots (e.g., Figs. 3 and 4).
%      = 2: generate CRB maps for 1 frequency combination(e.g., Fig. 5).
%           Note: could be any number of spatial frequencies, not just a pair.
%
% OUTPUTS:
% CRB = 2x2xn matrix, where n = (number of op pairs x number of frequency combinations).
%   Note: do not output when mode = 2 ("map" mode), as it will be huge and all
%   the info is captured in other 3 variables.
%
% sqrtCrbMua = Absorption coefficient uncertainty: mxn matrix, where m = number of frequency
%   combinations and n = number of op pairs.
%
% sqrtCrbMusp = Reduced scattering coefficient uncertainty: mxn matrix, where 
%   m = number of frequency combinations and n = number of op pairs.
%
% rho = Pearson correlation coefficient: mxn matrix, where m = number of frequency
%   combinations and n = number of op pairs.
%
% ------------------------------------------
% EXAMPLES
%
% Example 1: recreate Fig. 3 in paper
% op = [.005,.73; .005,1.89; .035,.73; .035,1.89];
% [CRB, sqrtCrbMua, sqrtCrbMusp, rho] = crb_sfdi2(DRC, fxFreqs, op, RdErr_wv(:,1), 1, 1);
%
% Example 2: recreate Fig. 4 in paper
% op = [.005,.73; .005,1.89; .035,.73; .035,1.89];
% [CRB, sqrtCrbMua, sqrtCrbMusp, rho] = crb_sfdi2(DRC, fxFreqs2, op, RdErr_wv(:,1), 1, 1);
%
% Example 3: create output similar to Fig. 5 in paper
% opBig = [DRC.Mua(:),DRC.Musp(:)];
% [~, sqrtCrbMua, sqrtCrbMusp, rho] = crb_sfdi2(DRC, {[0,.1]}, opBig, RdErr_wv(:,1), 1, 2);
% Note: this does not reproduce Fig. 5 exactly as this was generated with a
% Monte Carlo-based diffuse reflectance cube; the one provided here is 
% diffusion-based (See Sec. 3.1.2).
%
% V. Pera, 20 Dec 2017
% ------------------------------------------


% COMPUTE CRBs

if mode~=1 && mode~=2
    error('Mode must be equal to 1 or 2.')
end

if crbFlag~=0 && crbFlag~=1
    error('crbFlag must be equal to 0 or 1.')
end

% tic  % uncomment to display execution time

% spacing for gradient calc
dx = diff(DRC.Mua(1,1:2));
dy = diff(DRC.Musp(1:2,1));

FI = nan(2,2); % Fisher information matrix
CRB = nan(2,2,size(op,1)*size(fxFreqs,1));
sqrtCrbMua = nan(size(fxFreqs,1),size(op,1)); % extract CRB info to facilitate plotting
sqrtCrbMusp = nan(size(fxFreqs,1),size(op,1));
rho = nan(size(fxFreqs,1),size(op,1));

count = 1;

for ii=1:size(fxFreqs,1)
    
    disp(['Now working on ',num2str(ii),'/',num2str(size(fxFreqs,1)),' frequency combinations'])
    
    idx = nan(1,size(fxFreqs{ii},2));
    dxM = nan(size(DRC.M,1),size(DRC.M,2),size(fxFreqs{ii},2));
    dyM = nan(size(dxM));
    myCov = zeros(size(fxFreqs{ii},2));
    if crbFlag
        dxCpt = myCov;
        dyCpt = myCov;
    end
    
    for ff = 1:size(fxFreqs{ii},2) % loop over freqs
        [~,idx(ff)] = min(abs(DRC.fx-fxFreqs{ii}(ff))); % find freqs in LUT
        
        % Derivatives for Fisher info matrix (FI)
        [dxM(:,:,ff), dyM(:,:,ff)] = gradient(squeeze(DRC.M(:,:,idx(ff))), dx, dy);        
    end
    
    for jj=1:size(op,1)
                
        dxMpt = [];
        dyMpt = [];
    
        [~,idxMua] = min(abs(DRC.Mua(1,:)-op(jj,1))); % find mua, musp in LUT
        [~,idxMusp] = min(abs(DRC.Musp(:,1)-op(jj,2)));
        
        % Compute components of FI
        for ff = 1:size(fxFreqs{ii},2) % loop over freqs      
            dxMpt = [dxMpt, dxM(idxMusp,idxMua,ff)]; % Eqs. (7a)-(7d)
            dyMpt = [dyMpt, dyM(idxMusp,idxMua,ff)];
            
            % Build covariance matrix
            myCov(ff,ff) = (RdErr(idx(ff)) * DRC.M(idxMusp,idxMua,idx(ff)))^2;
%             myCov(ff,ff) = (RdErr(idx(ff)))^2;
            if crbFlag
                dxCpt(ff,ff) = 2 * RdErr(idx(ff))^2 * DRC.M(idxMusp,idxMua,idx(ff)) * dxM(idxMusp,idxMua,ff);
                dyCpt(ff,ff) = 2 * RdErr(idx(ff))^2 * DRC.M(idxMusp,idxMua,idx(ff)) * dyM(idxMusp,idxMua,ff);
            end
        end        
        
        if ~crbFlag           
            FI(1,1) = dxMpt * (myCov \ dxMpt.');  % FI for rCRB

            FI(1,2) = dxMpt * (myCov \ dyMpt.');

            FI(2,1) = dyMpt * (myCov \ dxMpt.');

            FI(2,2) = dyMpt * (myCov \ dyMpt.');
            
        else
        
            FI(1,1) = dxMpt * (myCov \ dxMpt.') + .5 * trace( myCov\dxCpt*(myCov\dxCpt) );  % FI for CRB

            FI(1,2) = dxMpt * (myCov \ dyMpt.') + .5 * trace( myCov\dxCpt*(myCov\dyCpt) );

            FI(2,1) = dyMpt * (myCov \ dxMpt.') + .5 * trace( myCov\dyCpt*(myCov\dxCpt) );

            FI(2,2) = dyMpt * (myCov \ dyMpt.') + .5 * trace( myCov\dyCpt*(myCov\dyCpt) );           
        end
        
        CRB(:,:,count) = inv(FI);
        
        sqrtCrbMua(ii,jj) = sqrt(CRB(1,1,count));
        sqrtCrbMusp(ii,jj) = sqrt(CRB(2,2,count));
        rho(ii,jj) = CRB(1,2,count) / sqrtCrbMua(ii,jj) / sqrtCrbMusp(ii,jj);
        
        count = count + 1;
        
        r = rem(jj,floor(size(op,1)*.2));
        if r == 0
            disp([num2str(jj/size(op,1)*100),'% complete'])
%             toc  % uncomment to display execution time
        end
    
    end
    
end
% toc  % uncomment to display execution time


% PLOT RESULTS

if mode == 1
     
    % Mua Uncertainty
    figure; bar(100 * sqrtCrbMua./repmat(op(:,1).',size(fxFreqs,1),1));   
    set(gcf,'name','sqrt(CRB) for Mua')
    xlabel('Spatial Frequency Combinations')
    ylabel('sqrt(CRB) for Mua (%)')
    %legend(mylegend)
    
    % Musp Uncertainty
    figure; bar(100 * sqrtCrbMusp./repmat(op(:,2).',size(fxFreqs,1),1));
    set(gcf,'name','sqrt(CRB) for Musp')
    xlabel('Spatial Frequency Combinations')
    ylabel('sqrt(CRB) for Musp (%)')
      
elseif mode == 2
     
    % Mua Uncertainty
    figure; imagesc(DRC.Mua(1,:), DRC.Musp(:,1), reshape(sqrtCrbMua,size(DRC.Mua))./DRC.Mua*100); colorbar; axis('xy')
    xlabel('mua (1/mm)'); ylabel('musp (1/mm)')
    set(gcf,'name','sqrt(CRB) for Mua (%)')
    
    % Musp Uncertainty
    figure; imagesc(DRC.Mua(1,:), DRC.Musp(:,1), reshape(sqrtCrbMusp,size(DRC.Mua))./DRC.Musp*100); colorbar; axis('xy') 
    xlabel('mua (1/mm)'); ylabel('musp (1/mm)')
    set(gcf,'name','sqrt(CRB) for Musp (%)')
    
    % Rho
    figure; imagesc(DRC.Mua(1,:), DRC.Musp(:,1), reshape(rho,size(DRC.Mua))); colorbar; axis('xy') 
    xlabel('mua (1/mm)'); ylabel('musp (1/mm)')
    set(gcf,'name','Rho')
    
end
% toc
