function varargout = RETINA_EXUDATES_DETECTION_CNN_GUI(varargin)
% RETINA_EXUDATES_DETECTION_CNN_GUI M-file for RETINA_EXUDATES_DETECTION_CNN_GUI.fig
%      RETINA_EXUDATES_DETECTION_CNN_GUI, by itself, creates a new RETINA_EXUDATES_DETECTION_CNN_GUI or raises the existing
%      singleton*.
%
%      H = RETINA_EXUDATES_DETECTION_CNN_GUI returns the handle to a new RETINA_EXUDATES_DETECTION_CNN_GUI or the handle to
%      the existing singleton*.
%
%      RETINA_EXUDATES_DETECTION_CNN_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RETINA_EXUDATES_DETECTION_CNN_GUI.M with the given input arguments.
%
%      RETINA_EXUDATES_DETECTION_CNN_GUI('Property','Value',...) creates a new RETINA_EXUDATES_DETECTION_CNN_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RETINA_EXUDATES_DETECTION_CNN_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RETINA_EXUDATES_DETECTION_CNN_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RETINA_EXUDATES_DETECTION_CNN_GUI

% Last Modified by GUIDE v2.5 11-Jan-2017 13:05:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RETINA_EXUDATES_DETECTION_CNN_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @RETINA_EXUDATES_DETECTION_CNN_GUI_OutputFcn, ...
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


% --- Executes just before RETINA_EXUDATES_DETECTION_CNN_GUI is made visible.
function RETINA_EXUDATES_DETECTION_CNN_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RETINA_EXUDATES_DETECTION_CNN_GUI (see VARARGIN)

% Choose default command line output for RETINA_EXUDATES_DETECTION_CNN_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RETINA_EXUDATES_DETECTION_CNN_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
set(handles.hSrcImage, 'XTick',[],'YTick',[],'XTickLabel',[], 'YTickLabel',[],'XGrid','off','YGrid','off','xcolor',get(handles.hSrcImage,'color'),'ycolor',get(handles.hSrcImage,'color'));
set(handles.axes8, 'XTick',[],'YTick',[],'XTickLabel',[], 'YTickLabel',[],'XGrid','off','YGrid','off','xcolor',get(handles.hSrcImage,'color'),'ycolor',get(handles.axes8,'color'));
set(handles.axes9, 'XTick',[],'YTick',[],'XTickLabel',[], 'YTickLabel',[],'XGrid','off','YGrid','off','xcolor',get(handles.hSrcImage,'color'),'ycolor',get(handles.axes9,'color'));
set(handles.axes10, 'XTick',[],'YTick',[],'XTickLabel',[], 'YTickLabel',[],'XGrid','off','YGrid','off','xcolor',get(handles.hSrcImage,'color'),'ycolor',get(handles.axes10,'color'));
% set(handles.axes11, 'XTick',[],'YTick',[],'XTickLabel',[], 'YTickLabel',[],'XGrid','off','YGrid','off','xcolor',get(handles.hSrcImage,'color'),'ycolor',get(handles.axes11,'color'));
set(handles.axes22, 'XTick',[],'YTick',[],'XTickLabel',[], 'YTickLabel',[],'XGrid','off','YGrid','off','xcolor',get(handles.hSrcImage,'color'),'ycolor',get(handles.axes22,'color'));
set(handles.axes23, 'XTick',[],'YTick',[],'XTickLabel',[], 'YTickLabel',[],'XGrid','off','YGrid','off','xcolor',get(handles.hSrcImage,'color'),'ycolor',get(handles.axes23,'color'));
set(handles.axes31, 'XTick',[],'YTick',[],'XTickLabel',[], 'YTickLabel',[],'XGrid','off','YGrid','off','xcolor',get(handles.hSrcImage,'color'),'ycolor',get(handles.axes31,'color'));
set(handles.axes32, 'XTick',[],'YTick',[],'XTickLabel',[], 'YTickLabel',[],'XGrid','off','YGrid','off','xcolor',get(handles.hSrcImage,'color'),'ycolor',get(handles.axes32,'color'));

addpath('./CNN_Toolbox')
% --- Outputs from this function are returned to the command line.
function varargout = RETINA_EXUDATES_DETECTION_CNN_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
global hiddenImage
varargout{1} = hiddenImage;


% --- Executes on button press in hRetinaImage.
function hRetinaImage_Callback(hObject, eventdata, handles)
% hObject    handle to hRetinaImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global testRetinaImg fileTest
[fileTest pathTest] = uigetfile({'*.png';'*.gif';'*.jpg';'*.*'},'Please select an Test Retina Image');
if fileTest==0
    warndlg('Select Input Image');
else   
    testRetinaImg = imread([pathTest fileTest]);
    axes(handles.hSrcImage);
    imshow(testRetinaImg);
    set(handles.editType,'String','');
    temp = uint8(ones(size(testRetinaImg))*255);
    axes(handles.axes23), imshow(temp)
    axes(handles.axes22), imshow(temp)
    axes(handles.axes31), imshow(temp)
    axes(handles.axes32), imshow(temp)
    axes(handles.axes8), imshow(temp)
    axes(handles.axes9), imshow(temp)
    axes(handles.axes10), imshow(temp)   
end



% --- Executes on button press in hEncrypt.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to hEncrypt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in hAddNoise.
function hAddNoise_Callback(hObject, eventdata, handles)
% hObject    handle to hAddNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in hClassification.
function hClassification_Callback(hObject, eventdata, handles)
% hObject    handle to hClassification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global testRetinaImg cnn classRes
[a, C]=max(cnn.layers{cnn.no_of_layers}.outputs, [],1);
if(C ~= 1 & a <0.9)
    disp('Normal Retina');
    set(handles.editType,'String','Normal');
    classRes = 'Normal';
else
    disp('Abnormal Retina');
    set(handles.editType,'String','Abnormal');    
    classRes = 'Abnormal';
end


function editStage_Callback(hObject, eventdata, handles)
% hObject    handle to editStage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editStage as text
%        str2double(get(hObject,'String')) returns contents of editStage as a double


% --- Executes during object creation, after setting all properties.
function editStage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editStage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in hCNNFeatmap.
function hCNNFeatmap_Callback(hObject, eventdata, handles)
% hObject    handle to hCNNFeatmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global testRetinaImg cnn

load Retina_DATA_CNN_model.mat

[cnn, convFeatMaps] = ffcnn(cnn, imresize(rgb2gray(testRetinaImg),[128 128]));
featMap = cell(2,5);
cnt=1;
for f1=1:2
    for f2=1:5
        featMap{f1,f2} = convFeatMaps{cnt};
        cnt=cnt+1;
    end
end
featMap = cell2mat(featMap);
axes(handles.axes23), imshow(featMap,[]);


% --- Executes on button press in hSegment.
function hSegment_Callback(hObject, eventdata, handles)
% hObject    handle to hSegment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global testRetinaImg Ives Iod

R = testRetinaImg(:,:,1);
G = testRetinaImg(:,:,2);
B = testRetinaImg(:,:,3);

% Vessel Segmentation
% Processing Green band
Gc = imclose(G, strel('disk', 10));
Gdiff =( imabsdiff((G), (Gc)));
axes(handles.axes22), imshow(Gdiff,[]);
cumHist = cumulativeHistogram(uint8(Gdiff));
cumHist = bwmorph(cumHist, 'majority', 1);
Ives = bwareaopen(cumHist, 500);
G1 = G;
G1(~Ives) = 0;
B1 = G1 < 120 & G1>0;
B2 = bwareaopen(B1,500);
Ives = bwareaopen(B1,1000);
axes(handles.axes31), imshow(Ives,[]);

% Processing Red band
numIter = 40;
delta_t = 1/7;
kappa = 30;
R1 = imresize(R, 0.25);
Ives1 = imresize(Ives, 0.25);
adImage = anisodiff2D(R1,numIter,delta_t,kappa,Ives1);
adImage = imresize(adImage, size(R));
od = adImage>200;
Iod = logical(od);
%figure(3), imshow((Iod),[])
[labelled nobj] = bwlabel(Iod);
stats = regionprops(labelled, 'area');
if nobj~=0
    for s = 1: nobj
        A(s) = stats(s).Area;
    end    
    maxArea = find(A==max(A));
    labelled(find(labelled~=maxArea)) = 0;
    Iod = logical(labelled);
end
axes(handles.axes32), imshow(Iod,[]);


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in hExudates.
function hExudates_Callback(hObject, eventdata, handles)
% hObject    handle to hExudates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global testRetinaImg classRes Ives Iod  fileTest
[row col nch] = size(testRetinaImg);
Fovea = segmentFovea(testRetinaImg, Ives,Iod);
dim = ndims(testRetinaImg);
if(dim == 3)
    IGray = rgb2gray(testRetinaImg);
else
    IGray = testRetinaImg;
end
se = strel('disk',12);
IGray = imclose(IGray,se);
Threshold = 5;
exudates = exudateslExtract(IGray, Threshold);
axes(handles.axes8), imshow(exudates,[]);
exu1 = logical(exudates);
exu1 = bwareaopen(exu1, 50);
exu1 = bwareaopen(exu1, 100);
[labelled nobj] = bwlabel(exu1);
Iod1 = imdilate(Iod, strel('disk', 30));
labelled1 = labelled;
h = waitbar(0,'Exudates detection in progress...');
for n = 1: nobj
    temp = labelled;
    temp(find(temp~=n)) = 0;
    temp1 = and(logical(temp), logical(Iod1));
    if(length(find(temp1~=0))>0)
        labelled1(find(labelled==n)) = 0;
    end
    waitbar(n/nobj);
end
close(h);
Exudates = imclearborder(labelled1);

exudatesNormal = exudateslExtract(rgb2gray(testRetinaImg), 10);
exudatesNormal = bwareaopen(exudatesNormal, 20);

[labelledNormal nobjNormal] = bwlabel(logical(exudatesNormal));
labelled2 = labelledNormal;
h = waitbar(0,'Exudates detection in progress...');
for n = 1: nobjNormal
    temp = labelledNormal;
    temp(find(temp~=n)) = 0;
    temp = imdilate(temp, strel('disk',3));
    
    temp1 = and(logical(temp), logical(Exudates));
    if(length(find(temp1~=0))==0)
        labelled2(find(labelledNormal==n)) = 0;
    end
    waitbar(n/nobj);
end
close(h);

ExudatesFinal = imfill(logical(labelled2), 'holes');
axes(handles.axes9), imshow(ExudatesFinal,[]);

% Soft Exudates detection
ycbcrI = rgb2ycbcr(testRetinaImg);
BWsoft = histeq(ycbcrI(:,:,1))>100;
BWsoft2 = bwareaopen(BWsoft,50);
BWsoft3 = bwareaopen(BWsoft2,500);
BWsoft4 = BWsoft2-BWsoft3;
SoftExudatesFinal = imfill(logical(BWsoft4), 'holes');

R = testRetinaImg(:,:,1);
G = testRetinaImg(:,:,2);
B = testRetinaImg(:,:,3);
R(find(ExudatesFinal)) = 0;
G(find(ExudatesFinal)) = 255;
B(find(ExudatesFinal)) = 0;

R(find(SoftExudatesFinal)) = 255;
G(find(SoftExudatesFinal)) = 255;
B(find(SoftExudatesFinal)) = 0;

outImg(:,:,1) = R;
outImg(:,:,2) = G;
outImg(:,:,3) = B;

imwrite(outImg,['OutImages\' fileTest(1:length(fileTest)-4) '.png']);
axes(handles.axes10), imshow(uint8(outImg),[]);
PerformenceComparison;

temp1 = and(ExudatesFinal, Fovea);

L1 = length(find(temp1~=0));

L2 = length(find(Fovea~=0));

L3 = length(find(ExudatesFinal~=0));

if(L3~=0)
    if(L1~=0)
        if(L1 > 0.1 * L2)
            disp('Severe DR');
        else if((L1> 0.05*L2) & (L1 < 0.1*L2))
                disp('Moderate DR');
            else
                disp('Mild DR');
            end
        end
    else
        disp('No DR');
    end
else
    disp('No DR');
end
