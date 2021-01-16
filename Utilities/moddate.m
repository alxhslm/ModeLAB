function datenum = moddate(file)
if iscell(file)
    for i = 1:length(file)
        datenum(i) = moddate(file{i});
    end
    return
end
path = fileparts(file);
if isempty(path)
    file = which(file);
end
d = dir(file);
datenum = d.datenum;