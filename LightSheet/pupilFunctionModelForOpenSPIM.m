%
% Utility function for fitLihtSheetModelToRecording and processOpenSPIMVideos
%
function pupil=pupilFunctionModelForOpenSPIM(params,V)
    if numel(params)<9
        params(9)=0;
    end
    pupil=zeros(size(V));
    for termIdx=4:6,
        factor=params(termIdx);
        if termIdx==4,
            factor=factor-params(6);
        end
        pupil=pupil+factor*V.^(termIdx-1);
    end
    pupilAmplitude=1+(params(7)-params(9)).*V + params(8).*V.^2 + params(9).*V.^3; % Can be negative == pi phase change!
    pupil=pupilAmplitude.*exp(2i*pi*pupil);
end