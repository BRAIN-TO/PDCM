function DCM = create_SPM_file_for_DCM(data,nscan,TR,onset,dur,cnam,cutoff,rnam,bins)
nses  = length(data);
ncon  = length(cnam);
%---------------------------------------------------------------------------
        SPM.nscan          = ones(1,nses)*nscan;
        
        % basis functions and timing parameters
        %---------------------------------------------------------------------------
        
        SPM.xBF.name       = 'hrf';
        SPM.xBF.order      = 0;
        SPM.xBF.length     = 32;                % length in seconds
        SPM.xBF.T          = bins;       	% number of time bins per scan
        SPM.xBF.T0         = 1;		% middle slice/timebin
        SPM.xBF.UNITS      = 'secs';           % OPTIONS: 'scans'|'secs' for onsets
        SPM.xBF.Volterra   = 1;                 % OPTIONS: 1|2 = order of convolution
        SPM.xY.RT          = TR;
        SPM.xBF.dt         = SPM.xY.RT/SPM.xBF.T;
        
        
        % Trial specification: Onsets, duration (UNITS) and parameters for modulation
        %---------------------------------------------------------------------------
        for ses=1:nses
           % ci = [1 1 2 2];
            for c=1:ncon
       
                SPM.Sess(ses).U(c).name      = {cnam{c}};
                SPM.Sess(ses).U(c).ons       = onset{ses}{c};
                SPM.Sess(ses).U(c).dur       = dur{ses}{c}; 
                SPM.Sess(ses).U(c).P(1).name = 'none';
            end
            Ui = spm_get_ons(SPM,ses);
            
            
            % Change input format to DCM input format and correct 32 bin offset that is inserted by spm_get_ons
            % (NB: for "normal" DCMs this is corrected for in spm_dcm_ui)
            Us(ses).U.u    = [];   
            Us(ses).U.name = {};
            for  i = 1:length(Ui)
                Us(ses).U.u             = [Us(ses).U.u Ui(i).u(33:end,1)];  % correct 32 bin offset that is inserted (as hard coded value!) by spm_get_ons
                Us(ses).U.name{end + 1} = Ui(i).name{1};
                % any parametric modulators?
                if size(Ui(i).u,2) > 1
                    for j = 2:size(Ui(i).u,2)
                        Us(ses).U.u             = [Us(ses).U.u Ui(i).u(33:end,j)];
                        Us(ses).U.name{end + 1} = Ui(i).P(j-1).name;
                    end
                end
            end
            Us(ses).U.dt = Ui(1).dt;
        end

        
        
        % low frequency confound: high-pass cutoff (secs) [Inf = no filtering]
        %---------------------------------------------------------------------------
        for ses=1:nses
            SPM.xX.K(ses).HParam = cutoff;
            SPM.xX.K(ses).RT = TR;
            SPM.xX.K(ses).row = ones(nscan,1);
            SPM.Sess(ses).C.C = [];
            SPM.Sess(ses).C.name = {};

        end
        SPM.xX.K = spm_filter(SPM.xX.K);
%         

        
                % get cell array of region structures
        %------------------------------------------------------------------
%         n = spm_input('Enter number of regions','+1','r',[],1);
%         for i=1:n,
%             str=sprintf('Region %d',i);
%             xY(i).name = spm_input(['Name for ',str],'+1','s');
%             % Make up spurious VOI info
%             % for compatibility with spm_dcm_display
%             xY(i).xyz = [i i i]'*10;
%             xY(i).XYZmm = [i i i]'*10;
%             xY(i).s=1;
%             xY(i).spec=1;
%         end

Y.y = [];
Y.u = [];
X0 = SPM.xX.K(1).X0;
X0 = [ones(nscan,1),X0];
no_regions = size(data(1).y,2);
for i = 1:no_regions
    %-Compute regional response in terms of first eigenvariate
    %--------------------------------------------------------------------------
    y = data.y(:,i);
    [m,n]   = size(y);
    
    % apply high pass filter
    for j = 1:n
        w = (X0'*y(:,j));
        if size(X0,2)>1
        y(:,j) = y(:,j) - X0(:,2:end)*w(2:end);
        end
    end
   
    [v,s,v] = svd(y'*y);
    s       = diag(s);
    v       = v(:,1);
    u       = y*v/sqrt(s(1));
    
    d       = sign(sum(v));
    u       = u*d;
    v       = v*d;
    u       = u*sqrt(s(1)/n);
    
    %-Set in structure
    %--------------------------------------------------------------------------
    xY(ses).y    = y;
    xY(ses).u    = u;
    xY(ses).v    = v;
    xY(ses).s    = s;
    Y.y          = [Y.y,mean(y,2)];
    Y.u          = [Y.u,u];
end

Y.X0  =  X0;
U.u   = [];
for i = 1:no_regions,
    Y.name{i} = rnam{i};
end
% added by YUexin 20260826
for i = 1:size(cnam,2)
    U.name{i} = cnam{i}; 
end

U.u  = Us(1).U.u;
U.dt = Us(1).U.dt;
Y.Q   = spm_Ce(ones(1,no_regions)*nscan);
Y.dt   = TR;
Y.secs = Y.dt*nscan;

DCM.Y = Y;
DCM.U = U;
