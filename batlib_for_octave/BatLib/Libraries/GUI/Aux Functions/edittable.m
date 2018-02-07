function varargout = edittable(varargin)
% EDITTABLE MATLAB code for edittable.fig
%      EDITTABLE, by itself, creates a new EDITTABLE or raises the existing
%      singleton*.
%
%      H = EDITTABLE returns the handle to a new EDITTABLE or the handle to
%      the existing singleton*.
%
%      EDITTABLE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EDITTABLE.M with the given input arguments.
%
%      EDITTABLE('Property','Value',...) creates a new EDITTABLE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before edittable_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to edittable_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help edittable

% Last Modified by GUIDE v2.5 15-Nov-2011 11:30:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @edittable_OpeningFcn, ...
                   'gui_OutputFcn',  @edittable_OutputFcn, ...
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


% --- Executes just before edittable is made visible.
function edittable_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to edittable (see VARARGIN)

% Choose default command line output for edittable
handles.output = hObject;

dontOpen = false;
mainGuiInput = find(strcmp(varargin, 'BCdetector'));
if (isempty(mainGuiInput)) ...
    || (length(varargin) <= mainGuiInput) ...
    || (~ishandle(varargin{mainGuiInput+1}))
    dontOpen = true;
else
    % Remember the handle, and adjust our position
    handles.BCdet.mainhandle = varargin{mainGuiInput+1};
    handles.BCdet.tabledata = varargin{mainGuiInput+2};
    handles.BCdet.labellist = varargin{mainGuiInput+3};
    handles.BCdet.methodname = varargin{mainGuiInput+4};
    handles.BCdet.inputhandle = varargin{mainGuiInput+5};

        columnformat = {'char','char','logical',handles.BCdet.labellist,handles.BCdet.methodname,...
                        'numeric','numeric','numeric','numeric','numeric',...
                        'numeric','numeric','numeric','numeric','numeric',...
                        'numeric','numeric','numeric'};
        columneditable =  [false false true  true  false,...
                           false false false false false,...
                           false false false false false,...
                           false false false]; 
        columnname =   {'Folder', 'File name', 'Det. Status', 'Label','Method',...
                        'Param1','Param2','Num. Det.','Duration','%Std',...
                        'Fs', 'Window','Overlap','dt','df',...
                        'Filter low', 'Filter High','SMS'};
        
        set(handles.uitable1,...
                        'Data',handles.BCdet.tabledata,...
                        'ColumnName', columnname,...
                        'ColumnFormat', columnformat,...
                        'ColumnEditable', columneditable,...
                        'ColumnWidth','auto',...
                        'RearrangeableColumns','on');
       
    % Position to be relative to parent:
    parentPosition = getpixelposition(handles.BCdet.mainhandle);
    currentPosition = getpixelposition(hObject);  
    % Set x to be directly in the middle, and y so that their tops align.
    newX = parentPosition(1) + (parentPosition(3)/2 - currentPosition(3)/2);
    newY = parentPosition(2) + (parentPosition(4)/2 - currentPosition(4)/2);
    %newY = parentPosition(2) + (parentPosition(4) - currentPosition(4));
    newW = currentPosition(3);
    newH = currentPosition(4);
    
    set(hObject, 'Units','pixels','Position', [newX, newY, newW, newH]);
  %Update handles structure
    guidata(hObject,handles)
end


if dontOpen
   disp('-----------------------------------------------------');
   disp('Improper input arguments. Pass a property value pair') 
   disp('whose name is "changeme_main" and value is the handle')
   disp('to the changeme_main figure, e.g:');
   disp('   x = changeme_main()');
   disp('   changeme_dialog(''changeme_main'', x)');
   disp('-----------------------------------------------------');
else
   uiwait(hObject);
end



 %UIWAIT makes edittable wait for user response (see UIRESUME)
 %uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = edittable_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.figure1;
delete(hObject);



% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
output.tabledata = get(handles.uitable1,'Data');
% Obtain handles using GUIDATA with the caller's handle 
if(ishandle(handles.BCdet.mainhandle))
selection = questdlg('Do you want to save changes?',...
      'Close Request Function',...
      'Yes','No','Cancel','Yes'); 
   switch selection, 
      case 'Yes'
      output.key = 1;
      set(handles.BCdet.inputhandle,'UserData',output)
      uiresume(handles.figure1); 
     case 'No'
       output.key = 0;
       set(handles.BCdet.inputhandle,'UserData',output)
       uiresume(handles.figure1);
       case 'Cancel'
      return 
   end
else
    uiresume(handles.figure1);
end


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected cell(s) is changed in uitable1.
function uitable1_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
