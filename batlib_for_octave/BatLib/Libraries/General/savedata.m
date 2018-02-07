function savedata(data)
    if data.path.SaveKey       

%                 system    = data.path.Sytem;
%                 winheader = data.path.WinHeader;
%                 linheader = data.path.LinHeader;
                 filepath  = data.path.Basedata;
%                 filepath = fullfile( changepath(filepath,system,winheader,linheader),data.Title);
                 filepath = fullfile(changepath(filepath,data.path),[data.Title ".mat"]);
                 err = [];
                 e = 0;
                disp(['Saving data:' filepath])
            try   
                save('-v7',filepath,'data');
                e = 0;
            catch err
                disp(err.message)
                e = 1;
            end
            
            if e
                try   
                save(filepath,'data','-v6');
                e = 0;
                catch err
                    disp(err.message)
                end
            end
            
            if e
                try   
                save(filepath,'data','-v7');
                e = 0;
                catch err
                    disp(err.message)
                    e = 1;
                end
            end
            
            if ~e
                disp('Data saved ')
            end
    end