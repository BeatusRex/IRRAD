clear
source_link='C:\Users\gantu\OneDrive\Belgeler\MATLAB\Master\proje\LED\m10\comma_folder\';
target_link='C:\Users\gantu\OneDrive\Belgeler\MATLAB\Master\proje\LED\m10\';
ngroup=25;LED={'UV';'Blue';'White';'Red';'Green'};LED=sort(LED); nLED=length(LED);
num=1:5; n=length(num);
matfiles=dir(fullfile(source_link, '*.txt'));
nfiles=length(matfiles);
for i=1:nfiles
    names{i,1}=matfiles(i).name;
    fname=names{i};
    if fname
        data=char(textread([source_link, fname], '%s', 'delimiter', '\n', 'whitespace', ''));
        for k=1:size(data, 1)
            f=findstr(data(k, :), ',');
            data(k, f)='.';
        end
        ind=findstr(fname, '.');
        fid=fopen([target_link, fname(1:ind-1), fname(ind:length(fname))], 'w');
        for k=1:size(data, 1)-1
            fprintf(fid, '%s\r\n', data(k, :));
        end
        fprintf(fid, '%s', data(size(data, 1), :));
        fclose(fid);
    end
end
disp('mission completo')
run analiz_m10