function P = read_setup_csv(file)

% edit(file);
fid  = fopen(file);
C = {};
X = {};
while ~feof(fid)
    line = fgetl(fid);
    if line(1) ~= '%'
        row = textscan(line,['%s' repmat('%s',1,50)],'Delimiter',',','CollectOutput',true);     
        field = row{1}{1};
        data = row{1}(2:end);
        data = data(~cellfun(@isempty,data));
        if field(1) ~= '#'
            if isempty(data)
                data = [];
            elseif ~isnan(str2double(data{1}))
                data = cellfun(@str2double,data);
            end
            P.(field) = data;
        end
    else
        %skip commented lines
    end
end
fclose(fid);

for i = 1:length(C)
    if C{i}(1) ~= '#'
        P.(C{i}) = X{i};
    end
end