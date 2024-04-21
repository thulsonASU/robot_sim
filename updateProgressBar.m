function updateProgressBar(progress, uiFigure)
    % Create a persistent variable for the progress dialog
    persistent pd

    % If the progress dialog does not exist, create it
    if isempty(pd)
        % Create a new figure using the uifigure function
        % uiFigure = uifigure('Visible', 'on');  % Invisible figure

        pd = uiprogressdlg(uiFigure, 'Title','Please Wait',...
                           'Message','Running Simulation...');
    end

    % Update the progress value
    pd.Value = progress;

    % Update the message
    pd.Message = sprintf('Running Simulation... %.2f%%', progress*100);

    % If the progress is 100%, close the dialog
    if progress >= 1
        close(pd);
        pd = [];
    end
end