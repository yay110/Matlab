function[FN DIR pressure] = filenames(FOLDER,TYPE,SIZE,p)
%MVGK mar. 2013
%FN         : List of file names, fn = array.
%DF         : Number of dropped files.
%DIR        : List of directories, DIR.cam1 = string, DIR.cam2 = string.
%FOLDER     : Name of folder containing data.
%TYPE       : Type of file, e.g. 'avi' or 'mat'.
%SIZE       : Size range of file in bytes [min max].
%p          : Extracting pressure variable [true/false].

%% Identifying the correct DATA drives
if nargin < 4,	p = false;	end
if p == false;  pressure  = [];	end

if 0
    drivesList      = getdrives('-nofloppy');                    	%Finds available drives and lists them in a structure

    for i = 1:size(drivesList,2)
        [~,vol]     = system(['vol ', drivesList{i}(1:2)]);
        if strcmp(strcat(vol(end-9:end)),'2470-857A')               %'2470-857A' is the serial number of the external harddisk
            cam1    = [drivesList{i} 'CS_data\Cam1\'];
            cam2    = [drivesList{i} 'CS_data\Cam2\'];
        end
    %     if strcmp(strcat(vol(end-9:end)),'228C-1CDD')               %'228C-1CDD' is the serial number of cam1 disk
    %         cam1    = drivesList{i};
    %     end
    %     if strcmp(strcat(vol(end-9:end)),'129D-5A06')               %'129D-5A06' is the serial number of cam2 disk
    %         cam2    = drivesList{i};
    %     end
    end
    DIR.cam1        = [cam1 FOLDER];
    DIR.cam2        = [cam2 FOLDER];
else
    DIR = [];
end

%% Identifying the correct files
cf              = cd;
cd(FOLDER)
list            = dir;
kk              = 0;

for i = 1:size(list,1)
    dateModified(i)         = 60^2*hour(list(i).date)+60*minute(list(i).date)+second(list(i).date);
    fileName0{i}            = list(i).name;
    fileSize(i)             = list(i).bytes;
end
[dateModified,id]       = sort(dateModified,'ascend');
for ii = 1:size(list,1)
    fileName{ii}            = fileName0{id(ii)};
end
fileSize                = fileSize(id);

for k = 1:size(list,1)
    if size(fileName{k},2) > 3;
        if strcmp(fileName{k}(end-length(TYPE)+1:end), TYPE);
            if fileSize(k) > SIZE(1) && fileSize(k) < SIZE(2)
                kk              = kk+1;
                FN{kk,:}        = fileName{k}(1:end-4);
                if p
                    if k-1<1
                        ind = k+1;
                    elseif k+1>length(dateModified)
                        ind = k-1;
                    else
                        timediffm       = abs(dateModified(k-1)-dateModified(k));
                        timediffp       = abs(dateModified(k+1)-dateModified(k));
                        if timediffm < timediffp
                            ind = k-1;
                        else
                            ind = k+1;
                        end
                    end
                    load([fileName{ind}(1:end-4),'.mat'],'press')
                    pressure(kk)    = press;
                end
            end
        end
    end
end
cd(cf)
end

%% Interesting codes
% FOLDER          = [DATE, ' - Cell stretching, rect. stress\'];

% f           = filesep; %automatic file separator

% [s,r]= dos('wmic logicaldisk get name');
% [s,r] = dos('net use');

% if strfind(vol(2:7), 'Disken')
%     a           = regexp(vol,'(?m)(?<drive>[A-Z]) er (?<label>\w+)$','names');
% else
%     a           = regexp(vol,'(?m)(?<drive>[A-Z]) is (?<label>\w+)$','names');
% end