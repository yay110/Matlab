function [U,S,Uss,Sss,PSD] = transformFourier(fps,s,prefiltering,plotting)
% function [U,S,Uss,Sss,PSD] = transformFourier(fps,s,prefiltering,plotting)
%
% PURPOSE:
%    The discrete-space Fourier transform of the signal s in x-space, into 
%    the signal S in the inverse X-space.
%
% INPUT:
%    fps: Frequency of datapoints in x-space, e.g. [Hz] or [1/m]
%
%    s: Vector OR 2D array of the signal.
% 
%    prefiltering: Leave empty (no pre-filtering), use 2 for Gaussian 
%                  windowing, or larger values for hyper-Gaussian window 
%                  filter for larger values. If a vector is specified, the
%                  second element is the standard deviation (default, one 
%                  quarter of the input period).
%
%    plotting: Logical true or false input determines whether the effect of 
%    this function is plotted.
% 
% OUTPUT:
%    U: The inverse space to "u", e.g. frequency [Hz] or k-space [1/meter].
%
%    S: The Fourier transform of "s".
%
%    Uss: The single sided version of "U", e.g. [0:1] rather than [-1:1].
%
%    Sss: The single sided Fourier transform of "s".
%
%    PSD: Power spectral density of "S". Does not yet work for the 2D
%    Fourier transform.
%
% NOTE:
%    Use abs(Sss) to get the single sided amplitude of the Fourier
%    transform of 's', and angle(Sss) for the single sided phase of the
%    Fourier transform of 's'.
%
% CREATED: Martin Verner Kristensen, Aarhus University, October 2011.
% Version 1&2: Tom Vettenburg & Martin Verner Kristensen, St Andrews, June 2013.
% Version 3: Martin Verner Kristensen, St Andrews, June 2014.

    %% Initialisation
    if nargin < 3,	prefiltering = [];	end
    if nargin < 4,	plotting = false;	end
    
    % Ensuring input data consists of column vectors
    if iscell(s)  && max(size(s)) == 2                                      % Checking if input data is a two dimensional cell
        if size(s,1) ~= size(u{1},1), s	= s'; end
        assert(size(s,1) == size(u{1},1) && size(s,2) == size(u{2},1), 'Bad cell [u] or matrix [s] dimensions!!!')
    elseif isnumeric(s) && min(size(s)) == 1                                % Checking if input data is a one dimensional vector
        if size(s,1) == 1, s = s'; end
    end
    
    %% Do the (interpolated) Fourier transform.
    s               = ifftshift(s);
    
    %% Improve output resolution by a factor of two so we get at least as good a resolution compared to zero-padding.
%     s=[s s]; 
%     t=[t t(end)+dt+t];

    %% 1D Fourier Transformation
    if isnumeric(s)
        L               = size(s,1);                                        % Number of measurements
        T               = (L-1)/fps;
        df              = 1/T;
        dt              = 1/fps;                                            % Time between two measurements
        t               = linspace(0,T,L)';
        
        % Multiply input data with a Gaussian to avoid sampling artefacts
        if ~isempty(prefiltering)
            if length(prefiltering)>1
                sigma           = prefiltering(2);
            else
                sigma           = dt*L/4;
            end
            windowWeigths   = L*(1/(sigma*sqrt(2*pi))).*exp(-0.5*((t-t(1+floor(end/2)))/sigma).^abs(prefiltering(1)));
            % Integrates to just below 1 due to the curtailing.
            windowWeigths   = L*windowWeigths./sum(windowWeigths);
            s               = s.*windowWeigths;
        end
        
        % Fourier transform and normalisation
        S               = fftshift(fft(s));                                 % Discrete-space Fourier transform for s
        S               = S/L;                                              % Normalised so that the frequency amplitude corresponds to the input signal amplitude
        
        % Ensuring symmetry and building frequency space
        if ~isodd(L)
            S               = [S; S(1)];                                  	% Makes both sides of the Fourier transform around zero equally long
            U               = fps/2*linspace(-1,1,L+1)';                  	% The inverse space to "u", e.g. frequency [Hz] or k-space [1/meter]
        else
            U               = fps/2*linspace(-1,1,L)';                    	% The inverse space to "u", e.g. frequency [Hz] or k-space [1/meter]
        end

        % Determining the amplitude spectrum
        indss           = 1+floor(size(S,1)/2):size(S,1);
        Sss             = S(indss);
        Sss             = abs([Sss(1); sqrt(2)*Sss(2:end)]);
        Uss             = fps/2*linspace(0,1,length(Sss))';                	% The single sided version of "U", e.g. [0:1] rather than [-1:1].
        
        % Power spectrum
        FT              = fft(s)/fps;                                       % Discrete-time Fourier transform for X.                                          % http://en.wikipedia.org/wiki/Discrete-time_Fourier_transform
        PSD             = fftshift(FT .* conj(FT) / T);                     % Power spectral density
        if ~isodd(L);	PSD = [PSD; PSD(1)]; end
        PSD             = PSD(indss);

    %% 2D Fourier Transformation
    % does not currently work, need to be updated to use the fps input.
    elseif iscell(s)
        x               = u{1};
        y               = u{2};
        Lx              = length(x);                                     	% Number of measurements in x-space
        Ly              = length(y);                                       	% Number of measurements in y-space
        kx              = (Lx-1)/(x(end)-x(1));                             % Frequency of datapoints in x-space [Hz] or [1/m]
        ky              = (Ly-1)/(y(end)-y(1));                           	% Frequency of datapoints in y-space [Hz] or [1/m]
        NFFTx           = 2^nextpow2(Lx);                               	% Next power of 2 from length of x
        NFFTy           = 2^nextpow2(Ly);                               	% Next power of 2 from length of y
        U{1}            = kx/2*linspace(-1,1,NFFTx+1)';                     % The inverse space to "u", e.g. frequency [Hz] or k-space [1/meter]
        U{2}            = ky/2*linspace(-1,1,NFFTy+1)';                     % The inverse space to "u", e.g. frequency [Hz] or k-space [1/meter]
        S               = fftshift(fft2(s,NFFTx,NFFTy)/(Lx*Ly));            % Fourier transform with zero-padding and shift
        S               = [S; S(1,:)];                                    	% Makes both sides of the Fourier transform around zero equally long
        S               = [S, S(:,1)];                                    	% Makes both sides of the Fourier transform around zero equally long
        Uss{1}          = kx/2*linspace(0,1,NFFTx/2+1)';                    % The single sided version of "U", e.g. [-1:1, 0:1] rather than [-1:1, -1:1].
        Uss{2}          = U{2};                                             % The single sided version of "U", e.g. [-1:1, 0:1] rather than [-1:1, -1:1].
        Sss             = 2*S(:,NFFTx/2+1:end);                          	% Single-Sided Spectrum of s
        Sss(NFFTx/2+1,1) = Sss(NFFTx/2+1,1)/2;                              % Correcting for the previous line's mistake of doubling the DC component
%         PSD             = SssA.*conj(SssA)/(nm/f);                       	% Power spectral density, does not yet work
        PSD             = [];%S.*conj(S)/(nmx*nmy);                      	% Power spectral density, does not yet work
    else
        error('[u] must be either a vector or a cell containing two vectors!!!')
    end

    %% Plotting effect of function
    if plotting
        if isnumeric(t)
            % 1D Fourier transform
            if ~isempty(prefiltering)
                figure,plot(t,windowWeigths);
                    set(gca,'FontSize',16)
                    xlabel('time','FontSize',16);
                    ylabel('window weight','FontSize',16);
                    title('Gaussian filter','FontSize',16)
            end
            figure,semilogy(Uss,abs(Sss))
                set(gca,'FontSize',16)
                title('Fourier Spectra','FontSize',16)
            figure,loglog(Uss,PSD)
                set(gca,'FontSize',16)
                title('Power Spectra','FontSize',16)
%                 legend('hide')
        elseif iscell(t)
            % 2D Fourier transform
%             BW = im2bw(abs(S)/max(max(abs(S))), 0.008);
            BW = 1;
            figure,imagesc(U{1},U{2},log10(BW.*abs(S))),colormap(jet)%,colorbar ,%subplot(212), imagesc(log10(abs((fm2)))),colorbar %,'EdgeColor','none'
                set(gca,'FontSize',16)
                title('Fourier Spectra','FontSize',16)
                xlabel('k_x'), ylabel('k_y')
                legend('hide')
        end
    end