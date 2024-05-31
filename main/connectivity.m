function connMatrix = connectivity(EEG,start,samp,chan,ref,reref,fh,fl,...
    wind,overlap,measure)


    % Fetch the data
    eegData = data(EEG,start,samp,chan,ref);

    % Compute the connectivity matrix
    connMatrix = measure.association(eegData,reref,fh,fl,wind,overlap);

end
