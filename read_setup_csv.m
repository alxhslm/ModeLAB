function P = read_setup_csv(file)
% edit(file);
fid  = fopen(file);
C = {};
X = {};
while ~feof(fid)
    line = fgetl(fid);
    row = textscan(line,['%s' repmat('%f',1,50)],'Delimiter',',','CollectOutput',true);
    C(end+1) = row{1};
    X{end+1} = row{2}(~isnan(row{2}));
end
fclose(fid);

for i = 1:length(C)
    if C{i}(1) ~= '#'
        P.(C{i}) = X{i};
    end
end