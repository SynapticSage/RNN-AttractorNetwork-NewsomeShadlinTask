function saveThis(handle,savedir,filename,type,nestFolder)
    if nargin < 5
        nestFolder='';
    end

    filename=[filename '.' type];
    if isequal(type,'png'); type=''; end;

    savestr = fullfile(savedir,nestFolder,type);
    warning off; mkdir(savestr); warning on;
    savestr=fullfile(savestr,filename);

    warning off; saveFigureFull(handle,savestr); warning on;
end