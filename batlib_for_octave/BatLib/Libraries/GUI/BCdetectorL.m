function varargout = BCdetectorL(varargin)
% BCDETECTORL M-file for BCdetectorL.fig
%      BCDETECTORL, by itself, creates a new BCDETECTORL or raises the existing
%      singleton*.
%
%      H = BCDETECTORL returns the handle to a new BCDETECTORL or the handle to
%      the existing singleton*.
%
%      BCDETECTORL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BCDETECTORL.M with the given input arguments.
%
%      BCDETECTORL('Property','Value',...) creates a new BCDETECTORL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BCdetectorL_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BCdetectorL_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BCdetectorL

% Last Modified by GUIDE v2.5 07-Aug-2013 18:34:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BCdetectorL_OpeningFcn, ...
                   'gui_OutputFcn',  @BCdetectorL_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before BCdetectorL is made visible.
function BCdetectorL_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BCdetectorL (see VARARGIN)

% Choose default command line output for BCdetectorL
handles.output = hObject;
clean(handles)

handles.AuxFunPath = [pwd '\Aux Functions'];
addpath(handles.AuxFunPath)

% Update handles structure
load('BCini.mat')
handles.framesize = framesize;
handles.freqaxeslim = freqaxeslim;
handles.windowtype = windowtype;
handles.colormap = colormap;%imadjust(1-gray(128),[0; 1],[],2);
handles.EFT = EFT;
set(handles.labelmenu,'String',labellist);


guidata(hObject, handles);




% UIWAIT makes BCdetectorL wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BCdetectorL_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% 


% --- Executes on slider movement.
function threshslide_Callback(hObject, eventdata, handles)
% hObject    handle to threshslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
alpha = get (handles.threshslide,'Value');
handles.wavinfo(handles.nfile).DetData.Param1 = alpha;
set(handles.edparam1,'String',num2str(alpha,'%3.1f'))

% detect push button script (!!)

handles = getdetparameter(handles);
handles = detect(handles);
%Activation of the detection data;
    n = handles.nfile;
    handles.wavinfo(n).DetDataKey = 1;
%     data = get(handles.audiofiletable,'Data');
%     data{n,3} = true;
%     data{n,4} = char(handles.wavinfo(n).DetData.Label);
%     set(handles.audiofiletable,'Data',data)
%display tasks
handles = displaydetection(handles);
handles = displayenergy(handles);
displaydetectiondata (handles);

guidata(hObject,handles)


% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function threshslide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshslide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function ednamedata_Callback(hObject, eventdata, handles)
% hObject    handle to ednamedata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ednamedata as text
%        str2double(get(hObject,'String')) returns contents of ednamedata as a double


% --- Executes during object creation, after setting all properties.
function ednamedata_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ednamedata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in savebutton.
function savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edparam1_Callback(hObject, eventdata, handles)
% hObject    handle to edparam1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edparam1 as text
%        str2double(get(hObject,'String')) returns contents of edparam1 as a double


% --- Executes during object creation, after setting all properties.
function edparam1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edparam1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Detect.
function Detect_Callback(hObject, eventdata, handles)
% hObject    handle to Detect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = getdetparameter(handles);
handles = detect(handles);
%Activation of the detection data;
    n = handles.nfile;
    handles.wavinfo(n).DetDataKey = 1;
%     data = get(handles.audiofiletable,'Data');
%     data{n,3} = true;
%     data{n,4} = char(handles.wavinfo(n).DetData.Label);
%     set(handles.audiofiletable,'Data',data)
%display tasks
handles = displaydetection(handles);
handles = displayenergy(handles);
displaydetectiondata (handles);

displaypanels (handles)

 handles = displaycurves(handles);

guidata(hObject,handles)


function edfilterlow_Callback(hObject, eventdata, handles)
% hObject    handle to edfilterlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edfilterlow as text
%        str2double(get(hObject,'String')) returns contents of edfilterlow as a double


% --- Executes during object creation, after setting all properties.
function edfilterlow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edfilterlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edwindow_Callback(hObject, eventdata, handles)
% hObject    handle to edwindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edwindow as text
%        str2double(get(hObject,'String')) returns contents of edwindow as a double


% --- Executes during object creation, after setting all properties.
function edwindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edwindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edoverlap_Callback(hObject, eventdata, handles)
% hObject    handle to edoverlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edoverlap as text
%        str2double(get(hObject,'String')) returns contents of edoverlap as a double


% --- Executes during object creation, after setting all properties.
function edoverlap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edoverlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edframe_Callback(hObject, eventdata, handles)
% hObject    handle to edframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edframe as text
%        str2double(get(hObject,'String')) returns contents of edframe as a double


% --- Executes during object creation, after setting all properties.
function edframe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in specbutton.
function specbutton_Callback(hObject, eventdata, handles)
% hObject    handle to specbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clean(handles)
handles = getspecparameter(handles);
displayfftinfo (handles);
[handles,P] = BCspec(handles);
handles = displayspec (handles,P);

handles.wavinfo(handles.nfile).DetDataKey = 0;
n = handles.nfile;
%         data      = get (handles.audiofiletable,'Data');
%         data{n,3} = false;
%         set (handles.audiofiletable,'Data',data);
handles = spec2detdata(handles);
%handles = setdetparameter(handles);
handles = displayenergy (handles);

%Activation of the stored data of the spectrogram;
handles.wavinfo(handles.nfile).SpecDataKey = 1;
displaypanels(handles)

    %eventdata.Indices = [2 n];
    %audiofiletable_CellSelectionCallback(handles.audiofiletable, eventdata, handles)
guidata(hObject,handles)




% --- Executes on button press in selectbutton.
function selectbutton_Callback(hObject, eventdata, handles)
% hObject    handle to selectbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on selection change in audiofilelist.
function audiofilelist_Callback(hObject, eventdata, handles)
% hObject    handle to audiofilelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.nfile = get(handles.audiofilelist,'Value');
clean(handles)
displayaudioinfo(handles);
guidata(hObject,handles);


% Hints: contents = get(hObject,'String') returns audiofilelist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from audiofilelist


% --- Executes during object creation, after setting all properties.
function audiofilelist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to audiofilelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in labelmenu.
function labelmenu_Callback(hObject, eventdata, handles)
% hObject    handle to labelmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns labelmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from labelmenu
 n = handles.nfile;
 if handles.wavinfo(n).DetDataKey
            labelindex = get(handles.labelmenu,'Value');
            labelist = get(handles.labelmenu,'String');
            label = labelist{labelindex};
    handles.wavinfo(n).DetData.LabelList = labelist;
    handles.wavinfo(n).DetData.Label =  label;
    handles.wavinfo(n).DetData.LabelIndex = labelindex;
    
%     data      = get (handles.audiofiletable,'Data');
%     data{n,4} = char(handles.wavinfo(n).DetData.Label);
%     set (handles.audiofiletable,'Data',data);
    %eventdata.Indices = [2 n];
    %audiofiletable_CellSelectionCallback(handles.audiofiletable, eventdata, handles)
    guidata(hObject,handles);
 end



% --- Executes during object creation, after setting all properties.
function labelmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to labelmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edenergyfilter_Callback(hObject, eventdata, handles)
% hObject    handle to edenergyfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edenergyfilter as text
%        str2double(get(hObject,'String')) returns contents of edenergyfilter as a double


% --- Executes during object creation, after setting all properties.
function edenergyfilter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edenergyfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in smscheckbox.
function smscheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to smscheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of smscheckbox
switch get(hObject,'Value')
    case 0
    set(handles.edframe,'Enable','off')
    
    case 1
    set(handles.edframe,'Enable','on')
    
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edfilterlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edfilterlow as text
%        str2double(get(hObject,'String')) returns contents of edfilterlow as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edfilterlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function displayaudioinfo(handles)
        n = handles.nfile;
        samples = num2str(handles.wavinfo(n).samples,'%1.0f samp.');
        time = num2str(handles.wavinfo(n).samples/handles.wavinfo(n).sampling_frequency,...
                       '%3.3f seg.');
        %set(handles.txname,'String',handles.wavinfo(n).name)
        set(handles.txsamples,'String',samples)
        set(handles.txlength,'String',time)
        set(handles.txfs,'String',num2str(handles.wavinfo(n).sampling_frequency/1000,'%3.1f kHz'))
        set(handles.txbits,'String',[num2str(handles.wavinfo(n).bits) ' bits'])
        

    function displayfftinfo (handles)
        n  = handles.nfile;
        fs = handles.wavinfo(n).sampling_frequency;
        dt = handles.wavinfo(n).SpecData.window/fs*1e3;
        df = 1./dt;
        set(handles.txwindow,'String',num2str(dt,'= %2.2f ms'))
        set(handles.txfreqres,'String',num2str(df,'freq. res. = %5.3g kHz'))
        
        
function handles = displayspec (handles,P0)
        n = handles.nfile;
        fs = handles.wavinfo(n).sampling_frequency;
        
        
        flim =[handles.freqaxeslim(1) min(handles.freqaxeslim(2),fs/2)];
        ftick = [0:5:300]*1000;
        fticklabel = cellstr(num2str([0:5:300]'));
        fticklabel(2:2:end)={''};
             
       t = handles.wavinfo(n).SpecData.Time;
       f = handles.wavinfo(n).SpecData.F;
       
        set(handles.specaxes,'Units','pixels');
        pos = get(handles.specaxes,'Position');
                        
        T = t(1):handles.framesize/pos(3):t(end);
        F = flim(1):(flim(2)- flim(1))/(pos(4)-1):flim(2);   
        handles.wavinfo(n).SpecData.imageF = F;
        handles.wavinfo(n).SpecData.imageT = T;
        
        switch nargin 
           case 1
               if isfield(handles.wavinfo(n).SpecData,'image')
                img = handles.wavinfo(n).SpecData.image;
                if size(img,3)==1;
                    ImPdB = ind2rgb(img,handles.colormap);
                else
                   ImPdB=img;                    
                end
               else
                  ImPdB = ones(length(T),length(F),3); 
               end
               kPxx = isfield(handles.wavinfo(n).SpecData,'Pxx');
               if kPxx
                  Pxx = handles.wavinfo(n).SpecData.Pxx;
                  Fxx = f;
               end
           case 2
               Pxx = 10*log10(mean(abs(P0),2));
               Fxx = f;
               a = 10e3;
               K = ceil(size(P0,2)/a);
               if K < 2;
                P = interpn(f,t,P0,F,T','linear');
                P =  10*log10(abs(P));
               else
                   P = zeros(length(F),length(T));
                   iK = ceil((1:size(P0,2))/a);
                   for k = 1:K;
                       it = find(iK==k);
                       iT = find(t(it(1))<=T & T <= t(it(end)));
                       PK = interpn(f,t(it),P0(:,it),F,T(iT)','linear');
                       P(:,iT)=PK;
                   end
                   P =  10*log10(abs(P));
               end
            maxP = max(P(:));
            PdBlim = [mean(P(isfinite(P))) maxP];
            grayPdB = mat2gray(P,PdBlim);
            [iPdB, map] = gray2ind(grayPdB, length(handles.colormap));
            ImPdB = ind2rgb(iPdB,handles.colormap);
            
            handles.wavinfo(n).SpecData.image = iPdB;
            handles.wavinfo(n).SpecData.imagemap = map;
            handles.wavinfo(n).SpecData.PdBlim = PdBlim;
            handles.wavinfo(n).SpecData.Pxx =Pxx;
            kPxx = 1;
            %ImPdB(find(ImPdB<0)) = 0;
        end
        Lcolor = [1 1 1]*.75; 
        alpha  = 5;
        
            
        axes(handles.specaxes)
        handles.specsurf = image(T,F,ImPdB);
        set(handles.specaxes,'YLim',flim,'YDir','normal');
        set(handles.specaxes,'Yticklabel','');
        axes(handles.poweraxes);
        if kPxx
            [Ppeaks, Pvalleys] = peakdet(Pxx,alpha);
            hold on
            plot(handles.poweraxes,Pxx,Fxx,'-b')
            if Ppeaks
            plot(handles.poweraxes,Ppeaks(:,2),Fxx(Ppeaks(:,1)),...
                            '^','Color',Lcolor,'MarkerSize',5)
            end
            hold off
            set(handles.poweraxes,'Xdir','reverse')
        end
        
        %handles.specsurf = surf(handles.specaxes,T,F,P_dB,'EdgeColor','none');
        
        %axis (handles.specaxes,'xy','tight');
        %colormap(cm); 
        %set(handles.specaxes,'View',[0 90]);
        %caxis(handles.specaxes,[mean(P_dB(find(P_dB>-Inf))) max(P_dB(:))])
        
%         ylabel(handles.specaxes,'Freq. (kHz)');
%         set(handles.specaxes,'YLim',flim,'YDir','normal');
%         settimeaxes(handles)
          %set(handles.specaxes,...
%                               'Ytick',ftick,'YTickLabel',ftick/1000 )
        ylabel(handles.poweraxes,'Freq. (kHz)');
        set(handles.poweraxes,'YLim',flim,'YDir','normal');
        settimeaxes(handles)
        
        %ttick = get(handles.specaxes,'Xtick');
        set(handles.poweraxes,...
                               'Ytick',ftick,'YTickLabel',fticklabel)
        % set(handles.specaxes,'Color',cm(1,:));
        %set(handles.specaxes,'Units','pixels')
        
    function settimeaxes(handles)
        n  = handles.nfile;
        Ti = handles.wavinfo(n).SpecData.Time(1);
        Tf = handles.wavinfo(n).SpecData.Time(end);
        if (Tf-Ti)<handles.framesize
        set (handles.timeslider,'Enable','off',...
            'Min',Ti,'Max',Tf,...
            'Value',Ti)
        set(handles.specaxes,'XLim',[Ti Tf])
%         set(handles.threshaxes,'XLim',[Ti Tf])
%         set(handles.threshaxes,'XLim',[Ti Tf])
%         set(handles.threshaxes,'XLim',[Ti Tf])
        else
        set (handles.timeslider,'Enable','on',...
            'Min',Ti,'Max',(Tf-handles.framesize),...
            'Value',Ti)
        set(handles.specaxes,'XLim',[Ti Ti+handles.framesize])
        b = (Tf-Ti)/(Tf-handles.framesize-Ti);
        set(handles.timeslider,'SliderStep',[0.005 b]);
%         set(handles.threshaxes,'XLim',[Ti Ti+handles.framesize])
%         set(handles.threshaxes,'XLim',[Ti Ti+handles.framesize])
%         set(handles.threshaxes,'XLim',[Ti Ti+handles.framesize])
        end
        
        
function handles = displaydetection(handles)
n  = handles.nfile;
D = handles.wavinfo(n).DetData.DetectionSignal;
T = handles.wavinfo(n).DetData.Time;

if handles.wavinfo(n).SpecData.filter
    flim = get(handles.specaxes,'Ylim');
    t = T(1);
    f = handles.wavinfo(n).SpecData.filter(1)/flim(2);
    dt= T(end)-T(1);
    df= diff(handles.wavinfo(n).SpecData.filter)/flim(2);
    
    d = D*df + f;    
    plot(handles.detaxes,T,d,'-w','LineWidth',1)
    
    axes(handles.detaxes)
    rectangle('Position',[t,f,dt,df],'EdgeColor',[0 0.7 0.7])
else
    plot(handles.detaxes,T,D,'-g','LineWidth',1)
end


if isfield(handles.wavinfo(n).DetData,'PhaseTime')
    phasecolor = {'w','g','y','m'};
    axes(handles.detaxes)
    for k = 2:4
        Tph = handles.wavinfo(n).DetData.PhaseTime{k};
        if ~isempty(Tph);
            for i = 1:size(Tph,1)
            t = Tph(i,1);
            f = handles.wavinfo(n).SpecData.filter(1)/flim(2);
            dt= Tph(i,2)-Tph(i,1);
            df= diff(handles.wavinfo(n).SpecData.filter)/flim(2);
            rectangle('Position',[t,f,dt,df],'EdgeColor',phasecolor{k},...
                       'LineWidth',2)        
            end                                  
        end    
    end
end
set(handles.detaxes,'XLim',get(handles.specaxes,'XLim'),'YLim',[0 1])
set(handles.detaxes,'Color','none')
axis (handles.detaxes,'off')



function handles = displayenergy(handles)
       n  = handles.nfile; 
        %elim = [-100 0];
        etick = [-120:5:0];
    
       T = handles.wavinfo(n).SpecData.Time;
       if isempty(handles.wavinfo(n).DetDataKey) || handles.wavinfo(n).DetDataKey == 0 
            E = handles.wavinfo(n).SpecData.E_dB;
            
       elseif handles.wavinfo(n).DetDataKey == 1
           switch handles.wavinfo(n).DetData.EFilterSwitch
               case 0
            E = handles.wavinfo(n).DetData.E_dB;
               case 1
            E = handles.wavinfo(n).DetData.Ef_dB;
           end
       end
        
%        if handles.energyplot
%         set(handles.energyplot,'XData',[],'YData',[])
%        end
       
       handles.energyplot = plot(handles.energyaxes,T,E,'-b');
             
        ylabel(handles.energyaxes,'Energy (dB)');
        xlabel(handles.energyaxes,'Time (sec.)');

        set(handles.energyaxes,'XLim',get(handles.specaxes,'XLim'));
        %ttick = get(handles.specaxes,'Xtick');
        set(handles.energyaxes,'Ytick',etick,'YLim',[max(min(E),-150) max(E)])
        %set(handles.energyaxes,'Ytick',etick,'YLim',[min(E) max(E)])
       %Threshold plot in the energy axes
       handles = displaythreshold(handles);
    function handles = displaycurves(handles)
    if get(handles.eftswitch,'value');
     n = handles.nfile;
     path = handles.wavinfo(n).path;
     fs = handles.wavinfo(n).sampling_frequency; 
     windowtype = handles.windowtype;
     windowlength = handles.wavinfo(n).SpecData.window;
     win = window(windowtype,windowlength);
     ovp =handles.wavinfo(n).SpecData.overlap;
     TCall = handles.wavinfo(n).DetData.TimeCalls;
     
    MethodIndex = handles.wavinfo(n).DetData.MethodIndex;
    switch MethodIndex
        case {2,3};
            
            beta = handles.EFT.Beta;
            Plim = handles.wavinfo(n).SpecData.PdBlim;
            image = double(handles.wavinfo(n).SpecData.image);
            map = handles.wavinfo(n).SpecData.imagemap;
            I = ind2gray(image,map);
            PdB = imlincomb(Plim(2)-Plim(1),I,Plim(1),'double');
            F = handles.wavinfo(n).SpecData.imageF;
            T = handles.wavinfo(n).SpecData.imageT;
            T0 = handles.wavinfo(n).DetData.Tpeak;
            Filt = handles.wavinfo(n).SpecData.filter;
            C = nan(size(T));
     for i = 1:handles.wavinfo(n).DetData.NCalls
         [v iT0] = min(abs(T-T0(i)));
         iT = find(TCall(i,1)<=T & T<= TCall(i,2));
         iF = find(Filt(1)<=F & F<=Filt(2));
         
         maxPdB = max(max(PdB(iF,iT0)));
         handles.EFT.MinPeakHeight  = maxPdB-beta;
         
         P = PdB(iF,iT);
         curve = EFTcurves (P,F(iF),T(iT),handles.EFT,T0(i));
        [c, it, ic] = intersect(T,curve{3});
        C(it)=curve{2};
     end
     set(gcf,'CurrentAxes',handles.detaxes);
     fl = get(handles.specaxes,'Ylim');
     df = fl(2)-fl(1);
     hold on
     plot(T,(C-fl(1))/df,'w','linewidth',2)
     hold off
    end
    end
             
    function handles = displaythreshold (handles)
       n  = handles.nfile;
        %Threshold plot in the energy axes
       T = handles.wavinfo(n).SpecData.Time;
       
       if isempty(get(handles.edparam1,'String')) || isempty(get(handles.edparam1,'String'))
        switch get(handles.detmethodmenu,'Value')
            case 1
                meanE = mean(handles.wavinfo(n).SpecData.E_dB);
                maxE  = max(handles.wavinfo(n).SpecData.E_dB);
                AlphaDefault = meanE + (maxE-meanE)*0.5;
                set(handles.edparam1,'String',num2str(AlphaDefault,'%2.2f'))
            case 2 
                Estd = std(handles.wavinfo(n).SpecData.E_dB)*2;
                set(handles.edparam1,'String',num2str(Estd,'%2.2f'))
                set(handles.edparam2,'String',num2str(50,'%2.2f'))
            case 3
                Estd = std(handles.wavinfo(n).SpecData.E_dB)*2;
                set(handles.edparam1,'String',num2str(Estd,'%2.2f'))
                set(handles.edparam2,'String',num2str(25,'%2.2f'))
        end
        end
       
       
       
       alpha = eval(get(handles.edparam1,'String'));
       elim = get(handles.energyaxes,'YLim');
        
       Lcolor = [1 1 1]*.75; 
       
       if handles.wavinfo(n).DetDataKey
       D2 = handles.wavinfo(n).DetData.DetectionSignal2;
       handles.threshplot =  plot(handles.threshaxes,T,D2,...
           '-','LineWidth', 1.1,'Color',Lcolor);
            switch handles.wavinfo(n).DetData.MethodIndex
                                
                case {1,2,3}
                    hold (handles.threshaxes,'on')
                    plot(handles.threshaxes,T(handles.wavinfo(n).DetData.Peaks(:,1)),handles.wavinfo(n).DetData.Peaks(:,2),...
                        '^','Color',Lcolor,'MarkerSize',5)
                    plot(handles.threshaxes,T(handles.wavinfo(n).DetData.Valleys(:,1)),handles.wavinfo(n).DetData.Valleys(:,2),...
                        'v','Color',Lcolor,'MarkerSize',5)
                    hold (handles.threshaxes,'off')
            end
       else
       handles.threshplot =  plot(handles.threshaxes,[T(1) T(end)],[1 1]*alpha,...
           '-','LineWidth', 1.1,'Color',Lcolor);
       end
       
       set(handles.threshaxes,'XLim',get(handles.specaxes,'XLim'),'YLim',elim,...
           'XTick',[],'YTick',alpha)
       set(handles.threshaxes,'Color','none')
       
        %Set of the slide bar
       if handles.wavinfo(n).DetDataKey
           switch handles.wavinfo(n).DetData.MethodIndex
                case 1
                   if alpha<elim(1) || alpha>elim(2)
                    set (handles.threshslide,'Enable','off')
                   else
                    set (handles.threshslide,'Enable','on')
                    set (handles.threshslide,'Min',elim(1),'Max',elim(2),'Value',alpha)
                   end
                case {2,3}
                    if alpha<=0 || alpha>abs(diff(elim))
                        set (handles.threshslide,'Enable','off')
                    else
                    set (handles.threshslide,'Enable','on')
                    set (handles.threshslide,'Min',0.01,'Max',abs(diff(elim)),'Value',alpha)
                    end
           end
       else
            set (handles.threshslide,'Enable','off')
       end
       
        function displaypanels(handles)
            
            if handles.wavinfo(handles.nfile).SpecDataKey
                    set(handles.detpanel,'Visible','on')
                if handles.wavinfo(handles.nfile).DetDataKey
                    set(handles.histpanel,'Visible','on')
                    set(handles.labelpanel,'Visible','on')
                else
                    set(handles.histpanel,'Visible','off')
                    set(handles.labelpanel,'Visible','off')
                end     
            else
                    set(handles.detpanel,'Visible','off')
                    set(handles.histpanel,'Visible','off')
                    set(handles.labelpanel,'Visible','off')
            end
    function handles = setspecparameter(handles)
        n = handles.nfile;
        %Window
        set(handles.popupwindow,'Value',handles.wavinfo(n).SpecData.guidat{1});
        %Overlap
        set(handles.popupoverlap,'Value',handles.wavinfo(n).SpecData.guidat{2});
        %Filter        
        if handles.wavinfo(n).SpecData.filter == 0;
            set(handles.edfilterlow,'String', '15');
            set(handles.edfilterhigh,'String','150');
            set(handles.edfilterswtich,'Value',0);
            set(handles.edfilterlow,'Enable','off');
            set(handles.edfilterhigh,'Enable','off');  
        else
            set(handles.edfilterlow,'String', num2str(handles.wavinfo(n).SpecData.filter(1)/1e3,'%2.3g'));
            set(handles.edfilterhigh,'String',num2str(handles.wavinfo(n).SpecData.filter(2)/1e3,'%2.3g'));
            set(handles.edfilterswtich,'Value',1);
            set(handles.edfilterlow,'Enable','on');
            set(handles.edfilterhigh,'Enable','on');
        end
        %Spectral Mean Substraction
        switch handles.wavinfo(n).SpecData.SMSswitch
            case 1
                set(handles.smscheckbox,'Value',1);
                set(handles.edframe,'String',num2str(handles.wavinfo(n).SpecData.SMSframe,'%2.3g'));
                set(handles.edframe,'Enable','on');
            case 0
                set(handles.smscheckbox,'Value',0);
                set(handles.edframe,'String','0.1');
                set(handles.edframe,'Enable','off');
        end
    function handles = getspecparameter(handles)
        n = handles.nfile;
    % Sampling frequency
        fs = handles.wavinfo(n).sampling_frequency;
    % Window size
        handles.wavinfo(n).SpecData.guidat{1} = get(handles.popupwindow,'Value');
        a = handles.wavinfo(n).SpecData.guidat{1} + 3;
        handles.wavinfo(n).SpecData.window = 2.^a;
    % Overlap
        handles.wavinfo(n).SpecData.guidat{2} = get(handles.popupoverlap,'Value');
        b = handles.wavinfo(n).SpecData.guidat{2} - 1;
        handles.wavinfo(n).SpecData.overlap = floor( (1-1/(2.^b))* 2.^a );
    % Spectral Mean Substraction
        handles.wavinfo(n).SpecData.SMSframe =  eval(char(get(handles.edframe,'String')));
        handles.wavinfo(n).SpecData.SMSswitch = get(handles.smscheckbox,'Value');
    
    % Filter
    
        switch get(handles.edfilterswtich,'Value')
            case 0
            handles.wavinfo(n).SpecData.filter = 0;
            case 1 
                if isempty(strtrim(get(handles.edfilterlow,'String')))
                    Flow = 0;
                else
                    Flow = eval(char(get(handles.edfilterlow,'String')))*1e3;
                end
                
                if isempty(strtrim(get(handles.edfilterhigh,'String')))
                    Fhigh = 0;
                else
                    Fhigh = eval(char(get(handles.edfilterhigh,'String')))*1e3;
                end
                handles.wavinfo(n).SpecData.filter = [Flow Fhigh] ;
        end
    function handles = getdetparameter(handles)
        n = handles.nfile;
        
        MethodIndex = get(handles.detmethodmenu,'Value');        
        Param1 = eval(get(handles.edparam1,'String'));
        switch MethodIndex
            case {1,2}
            Param2 = eval(get(handles.edparam2,'String'))/100;
            case{3}
            Param2 = eval(get(handles.edparam2,'String'))*1e-3;
        end
        EFilterWindow = eval(get(handles.edenergyfilter,'String'))/1e3;
        EFilterSwitch = get(handles.checkboxefliter,'Value');
        
        handles.wavinfo(n).DetData.MethodIndex = MethodIndex;
        handles.wavinfo(n).DetData.Param1 = Param1;
        handles.wavinfo(n).DetData.Param2 = Param2;
        handles.wavinfo(n).DetData.EFilterSwitch = EFilterSwitch;
        handles.wavinfo(n).DetData.EFilterWindow = EFilterWindow;
        
        
        
    function handles = setdetparameter(handles)
n = handles.nfile;
if isempty(handles.wavinfo(n).DetDataKey) || handles.wavinfo(n).DetDataKey == 0
        switch get(handles.detmethodmenu,'Value')
            case 1
                meanE = mean(handles.wavinfo(n).SpecData.E_dB);
                maxE  = max(handles.wavinfo(n).SpecData.E_dB);
                AlphaDefault = round(meanE + (maxE-meanE)/3);
                set(handles.edparam1,'String',num2str(AlphaDefault,'%2.2f'))
            case 2
                set(handles.edparam1,'String',num2str(10,'%2.2f'))
                set(handles.edparam2,'String',num2str(75,'%2.2f'))
            case 3
                set(handles.edparam1,'String',num2str(10,'%2.2f'))
                set(handles.edparam2,'String',num2str(25,'%2.2f'))                
        end
else
    % Method Type and Parameters
        MethodIndex = handles.wavinfo(n).DetData.MethodIndex;
        alpha       = handles.wavinfo(n).DetData.Param1;
        beta        = handles.wavinfo(n).DetData.Param2;
        
        set(handles.detmethodmenu,'Value',MethodIndex)
                
    switch MethodIndex
        case 1
            set (handles.txalpha,'String','Thresh.(dB)')
            set (handles.txbeta,'Visible','off')
            set (handles.edparam2,'Visible','off')
            set(handles.edparam1,'String',num2str(alpha,'%2.1f'));
            set(handles.edparam2,'String',num2str(beta*100,'%2.1f'));
        case 2
            set (handles.txalpha,'String','Sensibility (dB)')
            set (handles.txbeta,'String','Dur (ms)')
            set (handles.txbeta,'Visible','on')
            set (handles.edparam2,'Visible','on')
            set(handles.edparam1,'String',num2str(alpha,'%2.1f'));
            set(handles.edparam2,'String',num2str(beta*100,'%2.1f'));
        case 3
            set (handles.txalpha,'String','Sensibility (dB)')
            set (handles.txbeta,'String','Dur (ms)')
            set (handles.txbeta,'Visible','on')
            set (handles.edparam2,'Visible','on')
            set(handles.edparam1,'String',num2str(alpha,'%2.1f'));
            set(handles.edparam2,'String',num2str(beta*1000,'%2.1f'));
    end
         
        %Energy Filter
        EFilterWindow = handles.wavinfo(n).DetData.EFilterWindow;
        EFilterSwtich = handles.wavinfo(n).DetData.EFilterSwitch;
        set(handles.edenergyfilter,'String',num2str(1e3*EFilterWindow,'%2.2g'));
        set(handles.checkboxefliter,'Value',EFilterSwtich);        
        switch handles.wavinfo(n).DetData.EFilterSwitch
            case 0
                set(handles.edenergyfilter,'Enable','off')
            case 1
                set(handles.edenergyfilter,'Enable','on')
        end
        %Label
        labelindex    = handles.wavinfo(n).DetData.LabelIndex;
        set(handles.labelmenu,'Value',labelindex);
end
  

    function  [handles , P , E] = BCspec(handles)
     n = handles.nfile;
     fs = handles.wavinfo(n).sampling_frequency;
     maxlength=1e6;
     
        F = [];
        T = [];
        P = [];
        E = [];
    if handles.wavinfo(n).samples>maxlength
            h = waitbar (0,['Reading file: ']);           
            Npieces = ceil(handles.wavinfo(n).samples/maxlength);
           for k = 1:Npieces
                waitbar ((2*k-1)/(2*Npieces),h,['Reading file: ' handles.wavinfo(n).name ' part' num2str(k) '/' num2str(Npieces) ]);            
            y = wavread(handles.wavinfo(n).path, [maxlength*(k-1)+1 min(maxlength*k,handles.wavinfo(n).samples)]);
            
                waitbar((2*k)/(2*Npieces),h,['Calculating spectrogram: ' handles.wavinfo(n).name ' part' num2str(k) '/' num2str(Npieces) ])
            [Tsl,Fsl,Psl,Esl] = BCsms (y,handles);    
              
            F = [Fsl];
            P = [P Psl];
            E = [E Esl];
                if k==1
                    T = [T Tsl];
                else
                T = [T (Tsl+T(end))];
                end
            end
            
    else
                         
              h = waitbar (0,['Reading file: ' handles.wavinfo(n).name]);
              
              y = wavread(handles.wavinfo(n).path);
              
              waitbar(0.5,h,['Calculating spectrogram: ' handles.wavinfo(n).name])
              
              [T,F,P,E] = BCsms (y,handles);
             
            
    end
    waitbar(1,h)
handles.wavinfo(n).SpecData.Time  = T;
handles.wavinfo(n).SpecData.F     = F;
handles.wavinfo(n).SpecData.E_dB  = E;
    close(h)

        function [T,F,P,E_dB] = BCsms (x,handles)
     n = handles.nfile;
     fs = handles.wavinfo(n).sampling_frequency; 
     windowtype = handles.windowtype;
     windowlength = handles.wavinfo(n).SpecData.window;
     win = window(windowtype,windowlength);
     ovp =handles.wavinfo(n).SpecData.overlap;
     blocklength = handles.wavinfo(n).SpecData.SMSframe;
     filter = handles.wavinfo(n).SpecData.filter;
     E_mode = 'broadband';
     SMS_switch = handles.wavinfo(n).SpecData.SMSswitch;
     
[X,F,T,P] = spectrogram(x,win,ovp,[],fs);
dT = T(2)-T(1);
imaxT = length(T);
%blocklength in to samples
    s_block = round (blocklength*fs);

%Numbers of frames (columns) contained in the block

    N = round (s_block/dT);

% Numbers of blocks along the time T
    nmax = ceil(imaxT/N);

    if SMS_switch
    Psms = nan(size(P));
    for    n = 1:nmax
    ind = (n-1)*N+1:min(imaxT,n*N);  
    Pn = P(:,ind);     
    Psms(:,ind) = sms (Pn); 
    end
    % Quitar valores infinitos
    Psms(isnan(Psms))= 0;
    P = Psms;
    end
% Estimación de la energía
            if filter==0
               row = 1:length(F);
            else
               row = find(filter(1)<= F & F <=filter(2)); 
            end
                                         %(revisar  si es max o sum el mejor metodo)
        switch E_mode
            case 'broadband' 
                E = nansum (P(row,:));
            % Normalización de la energía
            %E = E;%/size(P,1)%/((3.^0.5)*pow2(16));  %(revisar si se tiene que dividir por fs también)        
            case 'peak'
                [Ep ind] = max (P(row,:));
                %Fp = F(ind)';   
                E = Ep;
            % Normalización de la energía
                E = E/((3.^0.5)*pow2(16));
        end
    E_dB = 10*log10(E);

            function [X_sms] = sms (X)
   
    X1 = log (abs(X));  L = size(X,2);
    Xp = sum(X1,2)/L*ones(1,L);
    X2 = X1 - Xp;
    
    X_sms = exp(X2);
    
        function handles = spec2detdata(handles)
            n  = handles.nfile;
        handles.wavinfo(n).DetData.E_dB      = handles.wavinfo(n).SpecData.E_dB;
        handles.wavinfo(n).DetData.Time      = handles.wavinfo(n).SpecData.Time;
        size(handles.wavinfo(n).DetData.E_dB)
        size(handles.wavinfo(n).DetData.Time )
                    

    function handles = detect (handles)
    n  = handles.nfile;
    
    E_dB        = handles.wavinfo(n).DetData.E_dB;
    T           = handles.wavinfo(n).DetData.Time;
    Method      = handles.wavinfo(n).DetData.MethodIndex;
    alpha       = handles.wavinfo(n).DetData.Param1;
    beta        = handles.wavinfo(n).DetData.Param2;
    wE          = handles.wavinfo(n).DetData.EFilterWindow;
    fs          = handles.wavinfo(n).sampling_frequency;
    EFswitch    = handles.wavinfo(n).DetData.EFilterSwitch ;
    wt          = handles.wavinfo(n).SpecData.window/fs;
    
    dT = T(2)- T(1);
    
    if EFswitch
    %Filtro de la energía
        w  = round(wE/dT); 
        Ef_dB = smooth(E_dB,w);

%     v = gausswin(wE); a = floor(wE/2);
%     Ef_dB = conv([E_dB 0],v)/sum(v);
%         %Ajustes
%        Ef_dB([1:a (end-a+1):end])=[];
%        Ef_dB([1:a (end-a+1):end])=E_dB([1:a (end-a+1):end]);
%        if length(Ef_dB) ~= length(E_dB)
%            Ef_dB(1)=[];
%        end
        
        E = Ef_dB;
     else
        E = E_dB;
    end

D = zeros (1,length(T));
TCall = [];


    switch Method
    
        case 1
            
            D (find (E>=alpha))= 1;
            Ds = sign(D);   Ds(1) = 0; Ds(end) = 0;
            ind_i = find (diff(Ds)==1)+1;
            ind_f = find (diff(Ds)==-1)+1;
            Ti = T(ind_i);
            Tf = T(ind_f);
            
            TCall = [Ti' Tf'];
            N = length(ind_i);
            
            D2      = ones(1,length(T))*alpha;
            Peaks   = ones(N,2);
            Valleys = ones(N,2);
            
             for j = 1:N
             [emax, im]= max(E(ind_i(j):ind_f(j)));
             imax = ind_i(j)-1 + im; 
             Peaks(j,:)  = [imax(1) emax(1)];
             
             [emin, im]= min(E(ind_i(j):ind_f(j)));
             imin = ind_i(j)-1 + im; 
             Valleys(j,:)   = [imin(1) emin(1)];
             end
            
        case 2
            
            [Peaks,Valleys] = peakdet(E,alpha);
            iP = Peaks(:,1);     
            P = Peaks(:,2);

            if size(Peaks,1)> size(Valleys,1)
               [V2 i2] = min (E(iP(end):end));
               i2 = i2 + iP(end)- 1;
                Valleys = [Valleys; i2 V2];
            end
 
            
            iV = Valleys(:,1); 
            V = Valleys(:,2);
            % Determine High of every peak P
                    h = [];fl = [];
                    %Found floor level for every peak
                    for i = 1:length(P)
                        if  i == 1
                            V0   =  min (E(1:iP(i)));
                            fl(i) = max ([V0 V(i)]); 
                        else
                             fl(i) = max( [V(i-1) V(i)]);  
                         end
                    end    
                        h = abs(P-fl');
                        
            D2 = ones(1,length(T))*min(E);
                        
            for i = 1:length(P)

                        if i==1
                            si = 1; sf = iV(i); 
                        else
                            si = iV(i-1); sf = iV(i);
                        end
                             Eslice = E(si:sf);
                             peak = iP(i) - si + 1;

                        Thresh = P(i)- beta*h(i);   
                       
                        if peak==1
                            ind_i = 1;
                        else
                        ind_i = find( Eslice (1:peak)  < Thresh,1,'last');
                        end
                        ind_f = find( Eslice (peak:end)< Thresh,1,'first') - 1 + peak;
                    
                        ind  = [ind_i ind_f] + si-1;
                        
                        D(ind(1):ind(2))= 1;
                        TCall  = [TCall;  T(ind(1)) T(ind(2))];  
         
                        D2 (si:sf) = fl(i);
                        D2 (ind(1):ind(2))= P(i);
            end
           
        case 3
            [Peaks,Valleys] = peakdet(E,alpha);
            iP = Peaks(:,1);     
            P = Peaks(:,2);
            
            if size(Peaks,1)> size(Valleys,1)
               [V2 i2] = min (E(iP(end):end));
               i2 = i2 + iP(end)- 1;
                Valleys = [Valleys; i2 V2];
            end

            iV = Valleys(:,1); 
            V = Valleys(:,2);
            
            ds = (beta/dT)/2;
            smax = length(T);
            diP = abs(diff([1;iP;smax]));
            
            % Determine High of every peak P
                    h = [];fl = [];
                    %Found floor level for every peak
                    for i = 1:length(P)
                        if  i == 1
                            V0   =  min (E(1:iP(i)));
                            fl(i) = max ([V0 V(i)]); 
                        else
                             fl(i) = max( [V(i-1) V(i)]);  
                         end
                    end
                        
            D2 = ones(1,smax)*min(E);
            
            for i = 1:length(iP)
                
                if diP(i)<2*ds;
                    if i>1
                    si = iV(i-1); 
                    else
                    si = 1; 
                    end
                else
                    si = iP(i)-ceil(ds);
                end
                
                if diP(i+1)<2*ds;
                   sf = iV(i);
                else
                   sf = min(smax,iP(i)+ceil(ds));
                end
                  
                    D(si+1:sf-1)= 1;
                    TCall  = [TCall;  T(si) T(sf)];
                    
                    if i==1
                            s2i = 1; s2f = iV(i); 
                    else
                            s2i = iV(i-1); s2f = iV(i);
                    end
         
                        D2 (s2i:s2f) = V(i);
                        D2 (si+1:sf-1)= P(i);
            end
           
            
    end
    Tpeak = T(Peaks(:,1));
    SCall = round(TCall*fs);
    NCalls = size(TCall,1);
       
       labelindex = get(handles.labelmenu,'Value');
       labelist = get(handles.labelmenu,'String');
%        label = cell(NCalls,1);
%        for L = 1:NCalls
%          label{L} = char(labelist{labelindex});
%        end
       label = labelist{labelindex};
       
     if EFswitch
     handles.wavinfo(n).DetData.Ef_dB = Ef_dB;
     end
     handles.wavinfo(n).DetData.LabelList = labelist;
     handles.wavinfo(n).DetData.Label =  label;
     handles.wavinfo(n).DetData.LabelIndex = labelindex;  
     handles.wavinfo(n).DetData.NCalls = NCalls;
     handles.wavinfo(n).DetData.TimeCalls = TCall;
     handles.wavinfo(n).DetData.SampleCalls = SCall;
     handles.wavinfo(n).DetData.Peaks   = Peaks;
     handles.wavinfo(n).DetData.Valleys = Valleys;
     handles.wavinfo(n).DetData.Tpeak  = Tpeak;
     handles.wavinfo(n).DetData.dTpeak = [Tpeak-wt/2 Tpeak+wt/2];
     handles.wavinfo(n).DetData.DetectionSignal = D;
     handles.wavinfo(n).DetData.DetectionSignal2 = D2;
     
     
function displaydetectiondata (handles)
     n  = handles.nfile;
    
     NCalls = handles.wavinfo(n).DetData.NCalls;
     TCalls = handles.wavinfo(n).DetData.TimeCalls;
     meanTCalls = mean(TCalls(:,2)-TCalls(:,1))*1e3;
     stdTCalls = std(TCalls(:,2)-TCalls(:,1))*1e3;
     
     set(handles.txNdet,'String',num2str(NCalls))
     set(handles.txSdur,'String',num2str(stdTCalls,'%3.1f ms'))
     set(handles.txMdur,'String',num2str(meanTCalls,'%3.1f ms'))
     
     y = (TCalls(:,2)- TCalls(:,1))*1e3;
     x = 0.25:0.5:27.25;
     hist(handles.timehist,y,x)
     set(handles.timehist,'XLim',[0 27.5],'Xtick',0:5:25)
function displaylabel(handles)
        n  = handles.nfile;
        labelindex = handles.wavinfo(n).DetData.LabelIndex;
            set(handles.labelmenu,'Value',labelindex)
        function clean(handles)
            cla(handles.specaxes,'reset');
            cla(handles.poweraxes,'reset');
            cla(handles.energyaxes,'reset');
            cla(handles.detaxes,'reset');
            cla(handles.threshaxes,'reset');
            set(handles.threshaxes,'Color','none')
            axis(handles.threshaxes,'off')
            set(handles.detaxes,'Color','none')
            axis(handles.detaxes,'off')
            set(handles.timeslider,'Enable','off')
        function setall(handles)
            clean(handles)
            displayaudioinfo(handles);
            if handles.wavinfo(handles.nfile).SpecDataKey
                        setspecparameter(handles);
                        displayspec(handles);
                        setdetparameter(handles);
                        displayenergy (handles);
            end

            if handles.wavinfo(handles.nfile).DetDataKey
                        setdetparameter(handles);
                        handles = displaydetection(handles);
                        handles = displayenergy(handles);
                        displaydetectiondata (handles);
                        displaylabel(handles);
            end
            displaypanels (handles);
% --------------------------------------------------------------------
function filemenu_Callback(hObject, eventdata, handles)
% hObject    handle to filemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function axesmenu_Callback(hObject, eventdata, handles)
% hObject    handle to axesmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function labelmanu_Callback(hObject, eventdata, handles)
% hObject    handle to labelmanu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function labellistmenu_Callback(hObject, eventdata, handles)
% hObject    handle to labellistmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function specaxesmenu_Callback(hObject, eventdata, handles)
% hObject    handle to specaxesmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function colormapmenu_Callback(hObject, eventdata, handles)
% hObject    handle to colormapmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function poweraxesmenu_Callback(hObject, eventdata, handles)
% hObject    handle to poweraxesmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function selectmenu_Callback(hObject, eventdata, handles)
% hObject    handle to selectmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function savemenu_Callback(hObject, eventdata, handles)
% hObject    handle to savemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on slider movement.
function timeslider_Callback(hObject, eventdata, handles)
% hObject    handle to timeslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
T1 = get(handles.timeslider,'Value');
    set(handles.specaxes,'XLim',[T1 T1+handles.framesize])
    set(handles.energyaxes,'XLim',[T1 T1+handles.framesize])
    set(handles.threshaxes,'XLim',[T1 T1+handles.framesize])
    set(handles.detaxes,'XLim',[T1 T1+handles.framesize])



% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function timeslider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeslider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function opentool_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to opentool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clean(handles)

[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select a session');

if FilterIndex
file = load([PathName FileName],'BCdetdata');

handles.wavinfo = file.BCdetdata;
handles.nfile = 0;
%set(handles.uipanelpath,'Title',''))
%set(handles.audiofilelist,'String',strrep ({handles.wavinfo.path},[path '\'],''),'Value',1)

%         path_sliced = regexp({handles.wavinfo.path}', '[^\\]*[^\\]', 'match');
%         data = [];
%         for i = 1: length(handles.wavinfo)
%                 data {i,1} = char(path_sliced{i}{end-1});
%                 data {i,2} = char(handles.wavinfo(i).name);
%             if handles.wavinfo(i).DetDataKey
%                 data{i,3} = true;
%                 data{i,4} = handles.wavinfo(i).DetData.Label;
%             else
%                 data {i,3} = false;
%                 data {i,4} = '';
%             end
%         end
        
%         columnformat = {'char','char','logical','char'};
%         columneditable =  [false false true false]; 
%         set(handles.audiofiletable,'Data',data,...
%             'ColumnFormat', columnformat,...
%             'ColumnEditable', columneditable)
        audiolist = regexp({handles.wavinfo.path}','\\[^\\]*\\[^\\]*$', 'match');
        audiolist =  char(cellstr([audiolist{:}]'));
        
        set(handles.audiolistbox,'String',audiolist)
        set(handles.audiolistbox,'Enable','on')
        set(handles.audiolistbox,'Value',1)
        

handles.Session.Path = PathName;
handles.Session.FilenName = FileName;
%displayaudioinfo(handles)
guidata(hObject,handles)
end


% --------------------------------------------------------------------
function FilterIndex = savetool_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to savetool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'wavinfo')
BCdetdata = handles.wavinfo;
% deleteind =[];
% % for i = 1:length(CData)
% %     if isempty(CData(i).TimeCalls)
% %         deleteind = [deleteind i];
% %     end
% % end
% CData(deleteind)=[];
% BCdetIm =[];
% for i = 1:length(BCdetdata)
%     if BCdetdata(i).SpecDataKey
%     BCdetIm(i).image = BCdetdata(i).SpecData.image;
%     BCdetIm(i).imagemap = BCdetdata(i).SpecData.imagemap;
%     BCdetdata(i).SpecData = rmfield(BCdetdata(i).SpecData, 'image');
%     end
%     
% end

if isfield(handles,'Session')
    [FileName,PathName,FilterIndex] = uiputfile('*.mat','Save detection data',...
        [handles.Session.Path '\' handles.Session.FilenName]);
    if FilterIndex
         disp('saving...')
         %uiwait(handles.figure1)
        save([PathName FileName],'BCdetdata','-v6')
        %save([PathName FileName],'BCdetIm','-v6','-append')
         %uiresume(handles.figure1)
         disp('save task finish')
    end

else
    
    DefaultName = 'BCdetdata';
    files = what(pwd);
    if strmatch([DefaultName '.mat'], files.mat)
        i=1;
       while strmatch([DefaultName '-' num2str(i) '.mat'], files.mat,'exact')
           i =i+1;
       end
       DefaultName = [DefaultName '-' num2str(i)];
    end

    [FileName,PathName,FilterIndex] = uiputfile('*.mat','Save detection data',[pwd '\' DefaultName '.mat']);
    if FilterIndex
         disp('saving...')
%         uiwait(handles.figure1)
        save([PathName FileName],'BCdetdata','-v6')
        %save([PathName FileName],'BCdetIm','-v6','-append')
%         uiresume(handles.figure1)
         disp('save task finish')
    end
    handles.Session.Path = PathName;
    handles.Session.FilenName = FileName;
end
guidata(hObject,handles)
end


% --------------------------------------------------------------------
function axestool_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to axestool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

answer = inputdlg({'Time axis length (sec.)','Frequency axis interval (kHz)','Window'},...
                    'Spectrogram axes properties',1,...
                    {num2str(handles.framesize),['[' num2str(handles.freqaxeslim/1e3) ']'],handles.windowtype});
if ~isempty(answer)
handles.framesize = str2num(answer{1});
handles.freqaxeslim = eval(answer{2})*1e3;
handles.windowtype = answer{3};
framesize = handles.framesize;
freqaxeslim = handles.freqaxeslim;
windowtype = handles.windowtype;
save('BCini.mat','framesize','freqaxeslim','windowtype','-append')
guidata(hObject,handles)
end


% --- Executes on selection change in popupwindow.
function popupwindow_Callback(hObject, eventdata, handles)
% hObject    handle to popupwindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = getspecparameter(handles);
displayfftinfo (handles)
guidata(hObject,handles)

% Hints: contents = get(hObject,'String') returns popupwindow contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupwindow


% --- Executes during object creation, after setting all properties.
function popupwindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupwindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupoverlap.
function popupoverlap_Callback(hObject, eventdata, handles)
% hObject    handle to popupoverlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupoverlap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupoverlap
handles = getspecparameter(handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function popupoverlap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupoverlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edfilterhigh_Callback(hObject, eventdata, handles)
% hObject    handle to edfilterhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edfilterhigh as text
%        str2double(get(hObject,'String')) returns contents of edfilterhigh as a double


% --- Executes during object creation, after setting all properties.
function edfilterhigh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edfilterhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in edfilterswtich.
function edfilterswtich_Callback(hObject, eventdata, handles)
% hObject    handle to edfilterswtich (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edfilterswtich

switch get(hObject,'Value')
    case 0
    set(handles.edfilterlow,'Enable','off')
    set(handles.edfilterhigh,'Enable','off')  
    case 1
    set(handles.edfilterlow,'Enable','on')
    set(handles.edfilterhigh,'Enable','on')
end
    


% --- Executes when entered data in editable cell(s) in audiofiletable.
function audiofiletable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to audiofiletable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)



% --- Executes when selected cell(s) is changed in audiofiletable.
function audiofiletable_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to audiofiletable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
if eventdata.Indices
    switch eventdata.Indices(2)
        case {1,2}
            if handles.nfile ~= eventdata.Indices(1)
            handles.nfile = eventdata.Indices(1);
            clean(handles)
            displayaudioinfo(handles);
                if handles.wavinfo(handles.nfile).SpecDataKey
                    setspecparameter(handles);
                    displayspec(handles);
                    setdetparameter(handles);
                    displayenergy (handles);
                end
                
                if handles.wavinfo(handles.nfile).DetDataKey
                    setdetparameter(handles);
                    handles = displaydetection(handles);
                    handles = displayenergy(handles);
                    displaydetectiondata (handles);
                end
     
            guidata(hObject,handles);
            end
            displaypanels (handles);
        case 3
    end
end


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkboxefliter.
function checkboxefliter_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxefliter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxefliter
switch get(hObject,'Value')
    case 0
    set(handles.edenergyfilter,'Enable','off')
    case 1
    set(handles.edenergyfilter,'Enable','on')
end


% --------------------------------------------------------------------
function newsesiontool_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to newsesiontool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clean(handles)
[wavinfo, path] = find_wavfiles;
handles.wavinfo = wavinfo;
handles.nfile = 1;
handles.wavinfo(handles.nfile).SpecDataKey = 0;
handles.wavinfo(handles.nfile).DetDataKey = 0;
set(handles.uipanelpath,'Title',regexp(path,'\\[^\\]*$', 'match'))
%set(handles.audiofilelist,'String',strrep ({handles.wavinfo.path},[path '\'],''),'Value',1)

%         path_sliced = regexp({handles.wavinfo.path}', '[^\\]*[^\\]', 'match');
%         data = [];
%         for i = 1: length(handles.wavinfo)
%         data {i,1} = char(path_sliced{i}{end-1});
%         data {i,2} = char(handles.wavinfo(i).name);
%         data {i,3} = isfield(handles.wavinfo(i),'DetStatus');
%         data {i,4} = '';
%         end
%         columnformat = {'char','char','logical','char'};
%         columneditable =  [false false true false]; 
%         set(handles.audiofiletable,'Data',data,...
%             'ColumnFormat', columnformat,...
%             'ColumnEditable', columneditable)
        audiolist = regexp({handles.wavinfo.path}','\\[^\\]*\\[^\\]*$', 'match');
        audiolist =  char(cellstr([audiolist{:}]'));
        
        set(handles.audiolistbox,'String',audiolist)
        set(handles.audiolistbox,'Enable','on')
        set(handles.audiolistbox,'Value',1)
 %Correct Sampling frequency
 for i = 1:length(handles.wavinfo)
    if handles.wavinfo(i).sampling_frequency<100000;
        handles.wavinfo(i).sampling_frequency = handles.wavinfo(i).sampling_frequency*10;
    end
 end
 
displayaudioinfo(handles)
guidata(hObject,handles)


% --------------------------------------------------------------------
function labellisttool_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to labellisttool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
labellist ={char( get(handles.labelmenu,'String'))};
L = length(get(handles.labelmenu,'String'));

options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='none';

answer = inputdlg('Label list','Set Label List',L,labellist,options);

    if isempty(answer)==0
    newlabellist = cellstr(answer{1});
    set(handles.labelmenu,'String',newlabellist)
    labellist = newlabellist;
    save ('BCini.mat', 'labellist','-append')
    end



function edparam2_Callback(hObject, eventdata, handles)
% hObject    handle to edparam2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edparam2 as text
%        str2double(get(hObject,'String')) returns contents of edparam2 as a double


% --- Executes during object creation, after setting all properties.
function edparam2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edparam2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in detmethodmenu.
function detmethodmenu_Callback(hObject, eventdata, handles)
% hObject    handle to detmethodmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch get(hObject,'Value')
    case 1
        set (handles.txalpha,'String','Thresh.(dB)')
        set (handles.txbeta,'Visible','off')
        set (handles.edparam2,'Visible','off')
    case 2
        set (handles.txalpha,'String','Sensibility (dB)')
        set (handles.txbeta,'String','Depth (%)')
        set (handles.txbeta,'Visible','on')
        set (handles.edparam2,'Visible','on')
    case 3
        set (handles.txalpha,'String','Sensibility (dB)')
        set (handles.txbeta,'String','Dur (ms)')
        set (handles.txbeta,'Visible','on')
        set (handles.edparam2,'Visible','on')
end
 n = handles.nfile;
 %if isempty(get(handles.edparam1,'String')) || isempty(get(handles.edparam1,'String'))
        switch get(handles.detmethodmenu,'Value')
            case 1
                meanE = mean(handles.wavinfo(n).SpecData.E_dB);
                maxE  = max(handles.wavinfo(n).SpecData.E_dB);
                AlphaDefault = meanE + (maxE-meanE)*0.5;
                set(handles.edparam1,'String',num2str(AlphaDefault,'%2.2f'))
            case 2 
                Estd = std(handles.wavinfo(n).SpecData.E_dB)*2;
                set(handles.edparam1,'String',num2str(Estd,'%2.2f'))
                set(handles.edparam2,'String',num2str(50,'%2.2f'))
            case 3
                Estd = std(handles.wavinfo(n).SpecData.E_dB)*2;
                set(handles.edparam1,'String',num2str(Estd,'%2.2f'))
                set(handles.edparam2,'String',num2str(25,'%2.2f'))
        end
 %end


% Hints: contents = cellstr(get(hObject,'String')) returns detmethodmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from detmethodmenu


% --- Executes during object creation, after setting all properties.
function detmethodmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to detmethodmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'Session')
selection = questdlg('Do you want to save changes?',...
      'Close Request Function',...
      'Yes','No','Cancel','Yes'); 
   switch selection, 
      case 'Yes',
        FilterIndex = savetool_ClickedCallback(handles.savetool, [], handles);
        if FilterIndex
        delete(hObject);rmpath(handles.AuxFunPath);
        else
            return
        end
     case 'No'
       delete(hObject);rmpath(handles.AuxFunPath);
       case 'Cancel'
      return 
   end
else
    delete(hObject);rmpath(handles.AuxFunPath);
end
    

% Hint: delete(hObject) closes the figure


% --- Executes on selection change in audiolistbox.
function audiolistbox_Callback(hObject, eventdata, handles)
% hObject    handle to audiolistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.nfile ~= get(hObject,'Value');
    handles.nfile = get(hObject,'Value');
    setall(handles)
guidata(hObject,handles);
end

% Hints: contents = cellstr(get(hObject,'String')) returns audiolistbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from audiolistbox


% --- Executes during object creation, after setting all properties.
function audiolistbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to audiolistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function editdatatool_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to editdatatool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'wavinfo')

        path_sliced = regexp({handles.wavinfo.path}', '[^\\]*[^\\]', 'match');
        tabledata = [];
        inputhandle = handles.audiolistbox;
        labellist = cellstr( get(handles.labelmenu,'String')');
        methodname =  cellstr( get(handles.detmethodmenu,'String')');
        
        hbar = waitbar(0,'Collecting data...');
        N = length(handles.wavinfo);
         for i = 1: N
                tabledata {i,1} = char(path_sliced{i}{end-1});
                tabledata {i,2} = char(handles.wavinfo(i).name);
               
            if handles.wavinfo(i).DetDataKey
                    tabledata{i,3} = true;
                    tabledata{i,4} = handles.wavinfo(i).DetData.Label;
                    tabledata{i,5} = methodname{handles.wavinfo(i).DetData.MethodIndex};
                switch handles.wavinfo(i).DetData.MethodIndex
                    case 1
                    tabledata{i,6} = handles.wavinfo(i).DetData.Param1;
                    case 2
                    tabledata{i,6} = handles.wavinfo(i).DetData.Param1;
                    tabledata{i,7} = handles.wavinfo(i).DetData.Param2*100;
                    case 3
                    tabledata{i,6} = handles.wavinfo(i).DetData.Param1;
                    tabledata{i,7} = handles.wavinfo(i).DetData.Param2*1000;
                end
                    tabledata{i,8} = handles.wavinfo(i).DetData.NCalls;
                if  tabledata{i,8}>1
                    TCalls = handles.wavinfo(i).DetData.TimeCalls ;
                    tabledata{i,9}  = mean(TCalls(:,2)-TCalls(:,1))*1e3;
                    tabledata{i,10} = (std(TCalls(:,2)-TCalls(:,1))*1e3)/tabledata{i,9}*100;
                end
                   
            else
                tabledata {i,3} = false;
                tabledata {i,4} = '';
            end
                    tabledata{i,11} = handles.wavinfo(i).sampling_frequency;
            if handles.wavinfo(i).SpecDataKey
                    tabledata{i,12} = handles.wavinfo(i).SpecData.window;
                    tabledata{i,13} = handles.wavinfo(i).SpecData.overlap/tabledata{i,12}*100;
                    tabledata{i,14} = handles.wavinfo(i).SpecData.window/tabledata{i,11}*1e3;
                    tabledata{i,15} = 1./tabledata{i,14}*1e3;
                    
                if handles.wavinfo(i).SpecData.filter
                    tabledata{i,16} = handles.wavinfo(i).SpecData.filter(1);
                    tabledata{i,17} = handles.wavinfo(i).SpecData.filter(2);
                end
                    
                if handles.wavinfo(i).SpecData.SMSswitch
                    tabledata{i,18} = handles.wavinfo(i).SpecData.SMSframe;
                end
            end
         waitbar(i/N,hbar);   
         end
         close(hbar)
        handles.output = tabledata;
 
        
        
       h = edittable('BCdetector',handles.figure1,tabledata,labellist,...
                    methodname,inputhandle);
       waitfor(h)
       
      input = get(inputhandle,'Userdata');
      
      if input.key
          newdata = input.tabledata;
            for i = 1: length(handles.wavinfo)
                
                if  newdata{i,3} == 1 & handles.wavinfo(i).DetDataKey == 0
                    handles.wavinfo(i).DetDataKey = 0;
                 else
                    handles.wavinfo(i).DetDataKey = newdata{i,3};
                end
                if newdata{i,3}
                    handles.wavinfo(i).DetData.Label =  newdata{i,4};
                end
           end
      end
      
      setall(handles)
      guidata(hObject,handles)
end
% --- Executes on button press in cleardetbutton.
function cleardetbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cleardetbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 n  = handles.nfile;
 if isempty(handles.wavinfo(n).DetDataKey)==0;
      handles.wavinfo(n).DetDataKey = 0;
      setall(handles)
      guidata(hObject,handles)
 end
% --------------------------------------------------------------------
function uipushtool12_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
titles = {'Min Peak Heigth (dB)','Number of Peaks',...
        'Maximum Sweep Rate (kHz/ms)','Maximum Energy Rate (dB/ms)',...
        'Minimum Time Length (ms)','Filter (kHz)'};
values = {num2str(handles.EFT.Beta),num2str(handles.EFT.NumPeaks),...
          num2str(handles.EFT.MaxSweepRate*1e-6),num2str(handles.EFT.MaxEnergyRate*1e-3)...
          num2str(handles.EFT.MinTimeLength*1e3),['[' num2str(handles.EFT.Filter*1e-3) ']']};
    

answer = inputdlg(titles,'EFT cruve',1,values);
if ~isempty(answer)
handles.EFT.Beta = str2num(answer{1});
handles.EFT.NumPeaks = str2num(answer{2});
handles.EFT.MaxSweepRate = str2num(answer{3})*1e6;
handles.EFT.MaxEnergyRate = str2num(answer{4})*1e3;
handles.EFT.MinTimeLength = str2num(answer{5})*1e-3;
handles.EFT.Filter = eval(answer{6})*1e3;
EFT = handles.EFT;
save('BCini.mat','EFT','-append')
guidata(hObject,handles)
end


% --- Executes on button press in eftswitch.
function eftswitch_Callback(hObject, eventdata, handles)
% hObject    handle to eftswitch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of eftswitch

% --------------------------------------------------------------------
function searchtoggletool_OnCallback(hObject, eventdata, handles)
% hObject    handle to searchtoggletool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)set(gcf,'CurrentAxes',handles.specaxes);
n = handles.nfile;
but =1;
t1 = [];
t2 = [];
k = 2;
if isfield(handles.wavinfo(n).DetData,'PhaseTime')
    T = handles.wavinfo(n).DetData.PhaseTime;
else
    T = cell(1,4);
end
while but ==1 || but ==29 || but==28
[t,f,but] = ginput(1);
    flim1 = get(handles.specaxes,'YLim');
    flim2 = handles.wavinfo(n).SpecData.filter;
    tlim = get(handles.specaxes,'XLim');
    if but ==1;
         if tlim(1)<=t & t<=tlim(2)
            if flim1(1)<=f & f<=flim1(2)
                if flim2(1)<=f & f<=flim2(2)
                    if isempty(t1)
                        t1 = t;
                    elseif t>t1
                        t2 = t;
                        T{k} = [T{k}; t1 t2];
                        handles.wavinfo(n).DetData.PhaseTime = T;
                        handles = displaydetection(handles);
                        t1 = [];
                    end
                else
                        but = 0;
                end
            end
         end
    elseif but==29
        eventdata.Key = 'rightarrow';
        figure1_KeyPressFcn(hObject, eventdata, handles)
    elseif but==28
        eventdata.Key = 'leftarrow';
        figure1_KeyPressFcn(hObject, eventdata, handles)
    end
end
guidata(hObject,handles)
set(handles.searchtoggletool,'state','off');
if ~isempty (t2)
    set(handles.approachtoggletool,'state','on');
    %approachtoggletool_OnCallback(hObject, eventdata, handles)
end

% --------------------------------------------------------------------
function searchtoggletool_OffCallback(hObject, eventdata, handles)
% hObject    handle to searchtoggletool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%handles = displaydetection(handles);
guidata(hObject,handles)


% --------------------------------------------------------------------
function approachtoggletool_OnCallback(hObject, eventdata, handles)
% hObject    handle to approachtoggletool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

n = handles.nfile;
but =1;
t1 = [];
t2 = [];
k = 3;
if isfield(handles.wavinfo(n).DetData,'PhaseTime')
    T = handles.wavinfo(n).DetData.PhaseTime;
else
    T = cell(1,4);
end
while  but ==1 || but ==29 || but==28
[t,f,but] = ginput(1);
    flim1 = get(handles.specaxes,'YLim');
    flim2 = handles.wavinfo(n).SpecData.filter;
    tlim = get(handles.specaxes,'XLim');
    if but ==1;
         if tlim(1)<=t & t<=tlim(2)
            if flim1(1)<=f & f<=flim1(2)
                if flim2(1)<=f & f<=flim2(2)
                    if isempty(t1)
                        t1 = t;
                    elseif t>t1
                        t2 = t;
                        T{k} = [T{k}; t1 t2];
                        handles.wavinfo(n).DetData.PhaseTime = T;
                        handles = displaydetection(handles);
                        t1 = [];
                    end
                else
                        but = 0;
                end
            end
         end
    elseif but==29
        eventdata.Key = 'rightarrow';
        figure1_KeyPressFcn(hObject, eventdata, handles)
    elseif but==28
        eventdata.Key = 'leftarrow';
        figure1_KeyPressFcn(hObject, eventdata, handles)
    end
end
guidata(hObject,handles)
set(handles.approachtoggletool,'state','off');
if ~isempty (t2)
    set(handles.terminaltoggletool,'state','on');
    %terminaltoggletool_OnCallback(hObject, eventdata, handles)
end

% --------------------------------------------------------------------
function approachtoggletool_OffCallback(hObject, eventdata, handles)
% hObject    handle to approachtoggletool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%handles = displaydetection(handles);
guidata(hObject,handles)


% --------------------------------------------------------------------
function terminaltoggletool_OnCallback(hObject, eventdata, handles)
% hObject    handle to terminaltoggletool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
n = handles.nfile;
but =1;
t1 = [];
t2 = [];
k = 4;
if isfield(handles.wavinfo(n).DetData,'PhaseTime')
    T = handles.wavinfo(n).DetData.PhaseTime;
else
    T = cell(1,4);
end
while  but ==1 || but ==29 || but==28
[t,f,but] = ginput(1);
    flim1 = get(handles.specaxes,'YLim');
    flim2 = handles.wavinfo(n).SpecData.filter;
    tlim = get(handles.specaxes,'XLim');
    if but ==1;
         if tlim(1)<=t & t<=tlim(2)
            if flim1(1)<=f & f<=flim1(2)
                if flim2(1)<=f & f<=flim2(2)
                    if isempty(t1)
                        t1 = t;
                    elseif t>t1
                        t2 = t;
                        T{k} = [T{k}; t1 t2];
                        handles.wavinfo(n).DetData.PhaseTime = T;
                        handles = displaydetection(handles);
                        t1 = [];
                    end
                else
                        but = 0;
                end
            end
         end
    elseif but==29
        eventdata.Key = 'rightarrow';
        figure1_KeyPressFcn(hObject, eventdata, handles)
    elseif but==28
        eventdata.Key = 'leftarrow';
        figure1_KeyPressFcn(hObject, eventdata, handles)
    end
end
guidata(hObject,handles)
set(handles.terminaltoggletool,'state','off');
if ~isempty (t2)
    set(handles.searchtoggletool,'state','on');
    %searchtoggletool_OnCallback(hObject, eventdata, handles)
end


% --------------------------------------------------------------------
function terminaltoggletool_OffCallback(hObject, eventdata, handles)
% hObject    handle to terminaltoggletool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%handles = displaydetection(handles);
guidata(hObject,handles)


% --------------------------------------------------------------------
function undtoggletool_OnCallback(hObject, eventdata, handles)
% hObject    handle to undtoggletool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)set(gcf,'CurrentAxes',handles.specaxes);
n = handles.nfile;
but =1;

if ~isfield(handles.wavinfo(n).DetData,'PhaseTime')
    handles.wavinfo(n).DetData.PhaseTime = cell(1,4);
else

while  but ==1 || but ==29 || but==28
    [t,f,but] = ginput(1);
    flim = get(handles.specaxes,'YLim');
    flim2 = handles.wavinfo(n).SpecData.filter;
    tlim = get(handles.specaxes,'XLim');
    if but ==1;
         if tlim(1)<=t & t<=tlim(2)
            if flim(1)<=f & f<=flim(2)
                if flim2(1)<=f & f<=flim2(2)
                T0 = handles.wavinfo(n).DetData.PhaseTime;
                    for k = 2:4
                        T = T0{k};
                        if ~isempty(T)
                            for i = 1:size(T,1)
                                if T(i,1)<=t & t<=T(i,2)
                                   T(i,:) = [];
                                   handles.wavinfo(n).DetData.PhaseTime{k} = T;
                                   handles = displaydetection(handles);
                                   guidata(hObject,handles)
                                end                              
                            end    
                        end
                    end
                else
                        but = 0;
                end
            end
         end
    elseif but==29
        eventdata.Key = 'rightarrow';
        figure1_KeyPressFcn(hObject, eventdata, handles)
    elseif but==28
        eventdata.Key = 'leftarrow';
        figure1_KeyPressFcn(hObject, eventdata, handles)
    end
end
end
guidata(hObject,handles)
set(handles.undtoggletool,'state','off');


% --------------------------------------------------------------------
function undtoggletool_OffCallback(hObject, eventdata, handles)
% hObject    handle to undtoggletool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%handles = displaydetection(handles);
guidata(hObject,handles)


% --------------------------------------------------------------------
function dettoggletool_OnCallback(hObject, eventdata, handles)
% hObject    handle to dettoggletool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(gcf,'CurrentAxes',handles.specaxes);
n = handles.nfile;
but =1;
while but ==1
[t,f,but] = ginput(1);
    flim = get(handles.specaxes,'YLim');
    tlim = get(handles.specaxes,'XLim');
    if but ==1;
         if tlim(1)<=t & t<=tlim(2)
            if flim(1)<=f & f<=flim(2)
                if handles.wavinfo(n).DetDataKey
                    T = handles.wavinfo(n).DetData.Time;
                    F = handles.wavinfo(n).SpecData.F;
                    D = handles.wavinfo(n).DetData.DetectionSignal;
                    D2 = handles.wavinfo(n).DetData.DetectionSignal2;
                    TCall = handles.wavinfo(n).DetData.TimeCalls;
                    NCall = handles.wavinfo(n).DetData.NCalls;
                    SCall = handles.wavinfo(n).DetData.SampleCalls;
                    Peaks = handles.wavinfo(n).DetData.Peaks;
                    Valleys = handles.wavinfo(n).DetData.Valleys;
                    fs      = handles.wavinfo(n).sampling_frequency;
                    wt      = handles.wavinfo(n).SpecData.window/fs;
                    [v iT] = min(abs(T-t));
                    iT=iT(1);
                    %[v iF ] = min(abs(F-f));

                    if D(iT)
                        sf = iT+find(~D(iT+1:1:end),1,'first');
                        si = iT-find(~D(iT-1:-1:1),1,'first');
                        D(si:sf)=0;
                        if (si-1)
                        D2(si:sf) = D2(si-1);
                        else
                        D2(si:sf) = D2(sf+1);
                        end
                       for i = 1:size(TCall,1)
                            if TCall(i,1) <= T(iT) & T(iT)<=TCall(i,2)
                                i0 = i;
                            end
                        end
                       TCall(i0,:)=[];
                       SCall(i0,:)=[];
                       Peaks(i0,:)=[];
                       
                       %Valleys(i0,:)=[];
                       
                    else
                        
                            dT = T(2)-T(1);
                            %ds = ceil((handles.wavinfo(n).DetData.Param2/dT)/2);
                            ds = round(8e-3/dT);
                            si = max(iT-ds,1);
                            sf = min(iT+ds,length(D));
%                             si = iT-min(ds,find([D(iT-1:-1:2)   1],1,'first'));
%                             sf = iT+min(ds,find([D(iT+1: 1:end-1) 1],1,'first'));
                        if handles.wavinfo(n).DetData.EFilterSwitch
                            E = handles.wavinfo(n).DetData.Ef_dB(si:sf);
                        else
                            E = handles.wavinfo(n).DetData.E_dB(si:sf);
                        end
                             [pks vly] = peakdet(E,handles.wavinfo(n).DetData.Param1/10);
                            %[pks vly] = peakdet(E,5);
                            if pks
                            iP = pks(:,1)+si-1;     
                            P = pks(:,2);
                            [v imx] = max(P);
                            
                            
                            si2 = max(iP(imx) - ds,1);
                            sf2 = min(iP(imx) + ds,length(D));
                                                
%                             si2 = iP(imx)-min(ds,find([D(iP(imx)-1:-1:2  )   1],1,'first'));
%                             sf2 = iP(imx)+min(ds,find([D(iP(imx)+1: 1:end-1) 1],1,'first'));
                            D(si2+2:sf2-2) = 1;
                            D2(si2+2:sf2-2)= P(imx);
                            
                            Peaks  = [Peaks; iP(imx) P(imx)];
                            TCall  = [TCall;  T(si2) T(sf2)];
                            SCall = round(TCall*fs);
                            end
                            
                            
                    end
                    Tpeak = T(Peaks(:,1));
                    handles.wavinfo(n).DetData.TimeCalls    = TCall;
                    handles.wavinfo(n).DetData.SampleCalls  = SCall;
                    handles.wavinfo(n).DetData.Peaks        = Peaks;
                    handles.wavinfo(n).DetData.Tpeak        = Tpeak;
                    handles.wavinfo(n).DetData.dTpeak       = [Tpeak-wt/2 Tpeak+wt/2];
                    %handles.wavinfo(n).DetData.Valleys      = Valleys;
                    handles.wavinfo(n).DetData.NCalls       = size(TCall,1);
                    handles.wavinfo(n).DetData.DetectionSignal = D;
                    handles.wavinfo(n).DetData.DetectionSignal2 = D2;
                    handles = displaydetection(handles);
                    handles = displayenergy(handles);
                    displaydetectiondata (handles);
                    guidata(hObject,handles)
                end
            end
        end
    end
end
set(handles.dettoggletool,'state','off')


% --------------------------------------------------------------------
function dettoggletool_OffCallback(hObject, eventdata, handles)
% hObject    handle to dettoggletool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = displaydetection(handles);
handles = displayenergy(handles);
handles = displaycurves(handles);
displaydetectiondata (handles);


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
dT = 0.5;
switch eventdata.Key
    case 'rightarrow'
     T0 =  get(handles.timeslider,'Value');
     Tmax = get(handles.timeslider,'Max');
     T = min(T0+dT,Tmax);
     set(handles.timeslider,'Value',T);    
     timeslider_Callback(hObject, eventdata, handles)
     
    case 'leftarrow'
     T0 =  get(handles.timeslider,'Value');
     Tmin = get(handles.timeslider,'Min');
     T = max(T0-dT,Tmin);
     set(handles.timeslider,'Value',T);    
     timeslider_Callback(hObject, eventdata, handles)
end

        


% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key release with focus on figure1 or any of its controls.
function figure1_WindowKeyReleaseFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was released, in lower case
%	Character: character interpretation of the key(s) that was released
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) released
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)
dT = 0.5;
if    gca == handles.detaxes...
   || gca == handles.threshaxes...
   || gca == handles.specaxes...
   || gca == handles.energyaxes
    switch eventdata.VerticalScrollCount
        case 1
         T0 =  get(handles.timeslider,'Value');
         Tmax = get(handles.timeslider,'Max');
         T = min(T0+dT,Tmax);
         set(handles.timeslider,'Value',T);    
         timeslider_Callback(hObject, eventdata, handles)

        case -1
         T0 =  get(handles.timeslider,'Value');
         Tmin = get(handles.timeslider,'Min');
         T = max(T0-dT,Tmin);
         set(handles.timeslider,'Value',T);    
         timeslider_Callback(hObject, eventdata, handles)
    end
end

% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over axes background.
function detaxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to detaxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function uitoggletool11_OnCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%set(gcf,'CurrentAxes',handles.detaxes);
[x, y] = getpts(handles.specaxes)
[x, y] = getline(handles.specaxes)
%get(handles.detaxes,'CurrentPoint')


% --- Executes on mouse press over axes background.
function specaxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to specaxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)