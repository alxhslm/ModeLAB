function P = read_setup_csv(file)

% edit(file);
fid  = fopen(file);
C = {};
X = {};
while ~feof(fid)
    line = fgetl(fid);
    iComment = find(line == '%',1);
    if ~isempty(iComment)
        line = line(1:(iComment-1));
    end
    if ~isempty(line)
        row = textscan(line,['%s' repmat('%s',1,50)],'Delimiter',',','CollectOutput',true);     
        field = row{1}{1};
        data = row{1}(2:end);
        data = data(1:find(~cellfun(@isempty,data),1,'last'));
        data(cellfun(@isempty,data)) = {NaN};
        if field(1) ~= '#'
            if isempty(data)
                data = [];
            elseif ~isnan(str2double(data{1})) || strcmp(data{1},'NaN')
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

if isfield(P,'fL')
    P.fLMode = P.fL;
end
if isfield(P,'fH')
    P.fHMode = P.fH;
end