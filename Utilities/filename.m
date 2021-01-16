function [file,ext] = filename(str)
if iscell(str)
    if nargout > 1
        for i = 1:length(str)
            [file{i},ext{i}] = filename(str{i});
        end
    else
        for i = 1:length(str)
            file{i} = filename(str{i});
        end
    end
    return
end

[~,file,ext] = fileparts(str);

if nargout < 2
    if ~isempty(ext)
        file = [file ext];
    end
end