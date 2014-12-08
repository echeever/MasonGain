function varargout = flowGraphGUI(varargin)
% FLOWGRAPHGUI MATLAB code for flowGraphGUI.fig
%      FLOWGRAPHGUI, by itself, creates a new FLOWGRAPHGUI or raises the existing
%      singleton*.
%
%      H = FLOWGRAPHGUI returns the handle to a new FLOWGRAPHGUI or the handle to
%      the existing singleton*.
%
%      FLOWGRAPHGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FLOWGRAPHGUI.M with the given input arguments.
%
%      FLOWGRAPHGUI('Property','Value',...) creates a new FLOWGRAPHGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before flowGraphGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to flowGraphGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help flowGraphGUI

% Last Modified by GUIDE v2.5 24-May-2013 00:28:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @flowGraphGUI_OpeningFcn, ...
    'gui_OutputFcn',  @flowGraphGUI_OutputFcn, ...
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
end
% End initialization code - DO NOT EDIT


% --- Executes just before flowGraphGUI is made visible.
function flowGraphGUI_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to flowGraphGUI (see VARARGIN)

% Choose default command line output for flowGraphGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes flowGraphGUI wait for user response (see UIRESUME)
% uiwait(handles.mainWindow);

init(handles)
end



% --- Outputs from this function are returned to the command line.
function varargout = flowGraphGUI_OutputFcn(~, ~, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;
end


% --- Executes on mouse press over axes background.
function mainAxes_ButtonDownFcn(~, ~, ~)
end


% --- Executes during object creation, after setting all properties.
function mainTable_CreateFcn(~, ~, ~)
end


% --- Executes on selection change in pathList.
function pathList_Callback(~, ~, handles)
showTraces(handles)
end


% --- Executes during object creation, after setting all properties.
function pathList_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on selection change in loopList.
function loopList_Callback(~, ~, handles)
showTraces(handles)
end


% --- Executes during object creation, after setting all properties.
function loopList_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes when entered data in editable cell(s) in mainTable.
function mainTable_CellEditCallback(~, eventdata, handles)
% hObject    handle to mainTable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

tableData = get(handles.mainTable,'Data');

% only if they just edited the value cell?
if eventdata.Indices(2) == 3 || ...
        size(cell2mat(tableData(eventdata.Indices(1),3)),2) 
    runMason(handles)
end
end


% --------------------------------------------------------------------
function menu_file_Callback(~, ~, ~)
end

% --------------------------------------------------------------------
function menu_new_Callback(~, ~, handles)
clearGraph(handles);
clearTable(handles);
clearLists(handles);
end


% --------------------------------------------------------------------
function menu_open_Callback(~, ~, handles)
[FileName,PathName,FilterIndex] = uigetfile('*.net','Pick a netlist file');

if(FileName)
    f = [PathName FileName];
    loadNet(handles,f)
end
end

% --------------------------------------------------------------------
function menu_save_Callback(~, ~, handles)
[FileName,PathName,FilterIndex] = uiputfile('*.net','Save Netlist','out.net');

if(FileName)
    f = [PathName FileName];
    saveNet(handles,f);
end
end

% --------------------------------------------------------------------
function menu_view_gains_Callback(~, ~, handles)
if strcmp(get(handles.menu_view_gains,'Checked'),'on')
    set(handles.menu_view_gains,'Checked','off')
    set(handles.menu_view_unity,'Enable','off')
else
    set(handles.menu_view_gains,'Checked','on')
    set(handles.menu_view_unity,'Enable','on')
end

runMason(handles)
end


% --------------------------------------------------------------------
function menu_view_unity_Callback(~, ~, handles)

if strcmp(get(handles.menu_view_unity,'Checked'),'on')
    set(handles.menu_view_unity,'Checked','off')
else
    set(handles.menu_view_unity,'Checked','on')
end

runMason(handles)
end


% --------------------------------------------------------------------
function menu_view_labels_Callback(~, ~, handles)
if strcmp(get(handles.menu_view_labels,'Checked'),'on')
    set(handles.menu_view_labels,'Checked','off')
else
    set(handles.menu_view_labels,'Checked','on')
end
end


% --------------------------------------------------------------------
function menu_view_eqn_Callback(~, ~, handles)
if strcmp(get(handles.menu_view_eqn,'Checked'),'on')
    set(handles.menu_view_eqn,'Checked','off')
else
    set(handles.menu_view_eqn,'Checked','on')
end

runMason(handles)
end


% --------------------------------------------------------------------
function menu_saveImage_Callback(~, ~, handles)

tmpFig = figure();
tmpAx = copyobj(handles.mainAxes,tmpFig);
set(tmpAx,'units', 'normalized', 'position', [0.025 0.025 0.95 0.95]);
set(tmpFig,'position',[1,1,1000,450]);
end


% --------------------------------------------------------------------
function menu_view_renderer_Callback(~, ~, ~)
end

% --------------------------------------------------------------------
function menu_view_painters_Callback(~, ~, handles)
setRenderer(handles,'painters')
end

% --------------------------------------------------------------------
function menu_view_opengl_Callback(~, ~, handles)
setRenderer(handles,'opengl')
end

% --------------------------------------------------------------------
function menu_view_zbuffer_Callback(~, ~, handles)
setRenderer(handles,'zbuffer')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
%              START OF FUNCTIONS THAT AREN'T UI CALLBACKS              %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function init(handles)

addpath('MasonsRule');

setRenderer(handles,'painters');

clearGraph(handles);
clearTable(handles);
clearLists(handles);
clearEqn(handles);
end


function setRenderer(handles,renderer)
set(handles.mainWindow,'Renderer',renderer);

set(handles.menu_view_painters,'Checked',on_if_match(renderer,'painters'));
set(handles.menu_view_opengl  ,'Checked',on_if_match(renderer,'opengl'  ));
set(handles.menu_view_zbuffer ,'Checked',on_if_match(renderer,'zbuffer' ));
end

function c = on_if_match(a,b)
if strcmp(a,b)
    c = 'on';
else
    c = 'off';
end
end


function clearGraph(handles)
ax = handles.mainAxes;
cla(ax)
axis(ax,'equal')
set(ax, ...
    'XTick',[], ...
    'YTick',[], ...
    'xcolor','w', ...
    'ycolor','w');

set(ax,'UserData',[])
end


function drawNet(handles,net,names)

clearGraph(handles);
ax = handles.mainAxes;

draw_gains = strcmp(get(handles.menu_view_gains,'Checked'),'on');
draw_unity = strcmp(get(handles.menu_view_unity,'Checked'),'on');

plotHandles=[];
for row=net'
    
    name = names(row(1));
    
    if ~draw_gains || (~draw_unity && str2double(name) == 1)
        name = '';
    end
    
    plotHandles(:,row(1))=drawBranch(row(2),row(3),name,ax);
end

set(handles.mainAxes,'UserData',plotHandles)
end


function handles=drawBranch(x0,x1,label,axes)

% figure out yoff
dst=x1-x0;
if dst==1
    dst=0;
end
yoff=dst*.22;

% the halfway point between ends
mdle=(x1+x0)/2;

% draw the line
[Lx,Ly]=MakeQuadraticBezier(30,[x0,0;mdle,yoff*2;x1,0]);
lh=line( ...
    'Parent',axes,...
    'XData',Lx,...
    'YData',Ly,...
    'LineWidth',2.5,...
    'Color',[0,0,0],...
    'LineSmoothing','on');


% draw the node
Nx=[x0,x1];
Ny=[0,0];

nh=line( ...
    'Parent',axes,...
    'XData',Nx,...
    'YData',Ny,...
    'LineStyle','none',...
    'Color',[0,0,0],...
    'Marker','.',...
    'MarkerSize',24,...
    'LineSmoothing','on');


% draw the arrow
direct=1-2*(x1<x0);

scl=.09;
offx=mdle;
offy=yoff;

Ax=[-.6,-1,1.4,-1 ];
Ay=[  0,.6,  0,-.6];

ah=patch(...
    'Parent',axes,...
    'XData',direct*Ax*scl+offx,...
    'YData',Ay*scl+offy,...
    'LineSmoothing','on');

% draw the label
th=text( ...
    'Position',[mdle,yoff+(.12*direct)],...
    'String',label,...
    'Parent',axes,...
    'FontSize',18,...
    'HorizontalAlignment','center', ...
    'interpreter','tex');

handles=[lh,nh,ah,th];
end


function [X,Y] = MakeQuadraticBezier(np,PP)
% Makes the points needed for a cubic bezier curve
% np   number of points on the curve
% PP   vector of 3 control points
%
%   example PP:  [0,0; 1,1; 2,0]
%

t=linspace(0,1,np);

X=(1-t).^2*PP(1,1)+2*(1-t).*t*PP(2,1)+t.^2*PP(3,1);
Y=(1-t).^2*PP(1,2)+2*(1-t).*t*PP(2,2)+t.^2*PP(3,2);
end


function clearLists(handles)

set(handles.loopList,'ListboxTop',1); % scroll to the top; otherwise poof?
set(handles.loopList,'String','');
set(handles.loopList,'value',1);

% storing corresponding traces in loopList's UserData
set(handles.loopList,'UserData',{[0]});


set(handles.pathList,'ListboxTop',1); % scroll to the top; otherwise poof?
set(handles.pathList,'String','');
set(handles.pathList,'value',1);

% storing corresponding traces in pathList's UserData
set(handles.pathList,'UserData',{[0]});
end


function clearTable(handles)

set(handles.mainTable,'Data',{'','',''});
end


function updateTable(handles,net,names)

tableData = {};

for i = 1:size(net,1)
    tableData(i,:) = {num2str(net(i,2)),num2str(net(i,3)),cell2mat(names(i))};
end
tableData = [tableData;{'','',''}];

set(handles.mainTable,'Data',tableData);
end


function [net,names] = loadFromTable(handles)

tableData = get(handles.mainTable,'Data');

net=[];
names={};
for i = 1:size(tableData,1)
    name = cell2mat(tableData(i,3));
    if size(name,2)
        net = [net ; [str2double(cell2mat(tableData(i,1))),...
            str2double(cell2mat(tableData(i,2)))]];
        names = [names ; name];
    end
end

% add line numbers
net=[(1:size(net,1));net']';
end


function runMason(handles)
% save as temp, load. load runs mason.

tmp_filename = '_tmp_.net';

saveNet(handles,tmp_filename)
loadNet(handles,tmp_filename)

delete(tmp_filename)
end


function loadNet(handles,filename)
% calls mason(), which also loads the net.

[Num,Den,P,L,net,names] = mason(filename,0,0);

drawNet(handles,net,names);

updateTable(handles,net,names);

drawEqn(handles,Num,Den);


% let's talk about these P and L structures.
%
% forward Paths (for path N)
% P{1,N} = {
%   Coeff:[1,2,3],
%   Node:[1,2,3,4]
% }
%
% Loops (for order O)
% L{1,O} = {
%   NumberLoops: 3,
%   Coeff: {[2,8,4,9],[2,8,6,10],[4,9,6,10]}
%   Node:  {[2,3,2,4,5,4],[2,3,2,6,7,6],[4,5,4,6,7,6]}
% }
%
% we're probably mostly interested in the Coeff fields.
%
% the list will enumerate paths and loops -
%  these'll be selectable, so you can highlight.
%

% we'll need to save the paths and loops somewhere
% traces = {
%   [1 2 3] (path)
%   []
%   [4 3 1] (loop)
%   [3 2]
% }

% might as well copy the format that's printed in mason()



% we'll start with Paths

pathListData = {};
pathTraces = {};

pathListData = [pathListData; ' ------- Paths -------'];
pathTraces = [pathTraces; 0];
for i=1:size(P,2)
    pathStruct = cell2mat(P(i));
    pathGainCell = names(P{i}.Coeff);
    pathGainStr=pathGainCell(:);
    pathGain = sprintf('%s*',pathGainStr{:});
    pathGain = pathGain(1:(end-1));  % Remove last "*"
    
    pathListData = [pathListData; ...
        ['P',num2str(i),': ',pathGain]];
    pathTraces = [pathTraces; pathStruct.Coeff];
    
    %     pathStruct = cell2mat(P(i));
    %     pathListData = [pathListData; ...
    %         ['P',num2str(i),': ',num2str(pathStruct.Node)]];
    %     pathTraces = [pathTraces; pathStruct.Coeff];
    
    
    %     pathStruct = cell2mat(P(i));  This is Ames' code
    %     path = pathStruct.Coeff;
    %     pathListData = [pathListData; ['P',num2str(i),': ',num2str(path)]];
    %     pathTraces = [pathTraces; path];
end
pathListData = [pathListData; ' '];
pathTraces = [pathTraces; 0];


set(handles.pathList,'ListboxTop',1); % scroll to the top; otherwise poof?
set(handles.pathList,'String',pathListData);
set(handles.pathList,'value',1); % select the first item

% storing corresponding traces in pathList's UserData
set(handles.pathList,'UserData',pathTraces);


% now loops
loopListData = {};
loopTraces = {};

for order = 1:size(L,2)
    loopsStruct = cell2mat(L(1,order));
    
    if loopsStruct.NumberLoops
        loopListData = [loopListData; ' --- Order ',num2str(order),' Loops ---'];
        loopTraces = [loopTraces; 0];
        
        for i = 1:loopsStruct.NumberLoops
            loopGainCell = names( cell2mat(loopsStruct.Coeff(i)));
            loopGainStr=loopGainCell(:);
            loopGain = sprintf('%s*',loopGainStr{:});
            loopGain = loopGain(1:(end-1));  % Remove last "*"
            
            loopListData = [loopListData; ...
                ['L',num2str(order),num2str(i),': ', loopGain]];
            %             loopListData = [loopListData; ...
            %                 ['L',num2str(order),num2str(i),': ', ...
            %                 num2str(cell2mat(loopsStruct.Node(i)))]];
            loopTraces = [loopTraces; cell2mat(loopsStruct.Coeff(i))];
            %             loop = cell2mat(loopsStruct.Coeff(i));   %Ames' code
            %             loopListData = [loopListData; ['L',num2str(order),num2str(i),': ',num2str(loop)]];
            %             loopTraces = [loopTraces; loop];
        end
        
        loopListData = [loopListData; ' '];
        loopTraces = [loopTraces; 0];
    end
    
    
end

set(handles.loopList,'ListboxTop',1); % scroll to the top; otherwise poof?
set(handles.loopList,'String',loopListData);
set(handles.loopList,'value',1); % select the first item

% storing corresponding traces in loopList's UserData
set(handles.loopList,'UserData',loopTraces);
end

function showTraces(handles)

loopTraces = get(handles.loopList,'UserData');
loopTraceIndex = get(handles.loopList,'Value');
loopTrace = cell2mat(loopTraces(loopTraceIndex));

pathTraces = get(handles.pathList,'UserData');
pathTraceIndex = get(handles.pathList,'Value');
pathTrace = cell2mat(pathTraces(pathTraceIndex));

plotHandles = get(handles.mainAxes,'UserData');


if plotHandles
    
    % intersection of loop and path traces
    iTrace = intersect(loopTrace, pathTrace);
    
    % union of loop and path traces
    uTrace = union(loopTrace, pathTrace);
    uTrace = uTrace(uTrace~=0); %remove zeros
    
    if uTrace
        
        colorBranches(plotHandles,[.8 .8 .8]);        % everything grey
        
        if pathTrace
            colorBranches(plotHandles(:,pathTrace),[0 0 1]); % path = blue
        end
        
        if loopTrace
            colorBranches(plotHandles(:,loopTrace),[1 0 0]); % loop = blue
        end
        
        if iTrace
            colorBranches(plotHandles(:,iTrace),[1 0 1]); % both = magenta
        end
        
        % pull the traced nodes to the top of the uistack
        for h = plotHandles(:,uTrace)
            uistack(h,'top');
        end
        
    else
        colorBranches(plotHandles,'black');  % everything black
    end
end
end

function colorBranches(branches,color)
% a plotHandle list has:
%  1 line
%  2 nodes
%  3 arrow
%  4 text
% e.g. plotHandles(3,1) is the first branch's arrow

set(branches([1,2,4],:),'Color',color) % line, nodes, text
set(branches(3,:),'FaceColor',color)   % arrow face
set(branches(3,:),'EdgeColor',color)   % arrow edge
end

function saveNet(handles,filename)

[net,names] = loadFromTable(handles);

fid = fopen(filename,'w');
for i = 1:size(net,1)
    row = net(i,:);
    name = cell2mat(names(i));
    fprintf(fid,'%d %d %d %s\n',row(1),row(2),row(3),name);
end
fclose(fid);
end

function clearEqn(handles)
ax = handles.eqnAxes;
cla(ax)
%axis(ax,'equal')
set(ax, ...
    'XTick',[], ...
    'YTick',[], ...
    'xcolor','w', ...
    'ycolor','w');

set(ax,'UserData',[])
end

function drawEqn(handles, Num, Den)

clearEqn(handles)

if strcmp(get(handles.menu_view_eqn,'Checked'),'off')
    return
end

s_num = sym(Num);
s_den = sym(Den);

eqn = simple(collect(s_num/s_den,'s'));
eql_latex = ['$' latex(eqn) '$'];

text(...
    'String',eql_latex,...
    'Interpreter','Latex',...
    'Position',[.5,.5],...
    'FontSize',24,...
    'HorizontalAlignment','center',...
    'Parent',handles.eqnAxes)
end


% The following function and subfuctions were originally written by Rob
% Wslton and taken from the MathWorks web site's file exchange.  They were
% adapted by Ames Bielenberg so that it could be used with the GUI wrapper.
function [Num,Den,P,L,Net,Coeff_Names] = mason(NetFile,Start,Stop)
% mason.m
% This function takes a netfile describing a signal flow graph
% with symbolic coefficients and generates an equation representing
% the equivilent term between an independent input node, and dependent
% output node. Please see the *readme* file for a full description, and
% an example netfile.
%
% Author :  Rob Walton
% Organisation : TRLabs and the University of Calgary, Alberta, Canada
% Date  :  January 25th 1999
% Revised : January 20th 2000 (Removed bug that caused the odd loop to be thrown away)
%
% Please email me at <walton@ieee.org> if you find any errors!
%
% USAGE:
%   [Numerator,Denominator] = mason(Netfile,StartNode,StopNode)
%
%   Netfile     - is a string with the name of the netfile to load
%   StartNode   - is the integer number describing the independent input node
%   StopNode    - is the integer number describing the dependent output node
%   Numerator   - is a string containing the equation for the Numerator
%   Denominator - is a string containing the equation for the Denominator


% Modifications by Ames Bielenberg in April 2013:
%  - changed output variables to include P,L,Net,Coeff_Names
%  - default to first and last node if Start and Stop are both 0 (line 65)



%*** Load the the net description file into Net ***
% *** The first column in the file is the coefficient number (These must be in order starting
% at 1). The second and third column are start node and stop node numbers respectively. The last
% column is the name of the coefficient. ** See readme file **

% *** Net will be the first 3 columns of the file. Coeff will be the last column. ***
fid=fopen(NetFile);			% Open the file
if (fid==-1)
    fprintf('\n*** File, %s, not found ***\n\n',NetFile)
    return
end

Net=[];					% Initialize the Numerical Net matrix
line_number=0;				% Line number read from file
Coeff_Names={};				% Initialize cell array of strings

while 1 				% Loop Until end of file, read one line at a time
    line_number=line_number+1;           % Count which line we are on
    x=fscanf(fid,'%d',3);		% Read the three decimal numbers into x
    Coeff=fscanf(fid,'%s\n',1);		% Read the one coefficient name into coeff
    if isempty(x)			% Stop the loop if no data is left in file
        break
    end
    Net(line_number,:)=transpose(x);     % Append one row to bottom of numerical matrix
    Coeff_Names{line_number}= Coeff;     % Append the coefficient name to the coefficient array
end
fclose(fid);				% Remember to close the file!


%*** Determine Number of Coefficients in net file ***
temp=size(Net);				% Determine the number of rows in the net matrix
Number_Coeff=temp(1);


%%%%%%%   ADDED BY AMES   %%%%%%%
%
% default to first and last if Start and Stop are 0.
if Start==0 && Stop==0
    Start = min(min(Net(:,2:3)));
    Stop = max(max(Net(:,2:3)));
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%*** Find all the paths connecting the start and end nodes, along which ***
%*** no node is crossed more than once                                  ***

[PathCoeffList,PathNodeList]=findpaths(Start,Stop,[],[],Net);

% PathCoeffList and PathNodeList are matrixes where each row lists the number of the coefficients or the nodes visited respectively. Each row is padded with zeros to make them the same length.

%fprintf('\n- Path List -\n');
%print_paths(PathCoeffList);

%*** Determine all the first-Order loops ***
LoopCoeffList=[];			% Initialise list of coefficients in for each loop
LoopNodeList=[];			% Initialise list of nodes for each loop found

for index=1:Number_Coeff;      		% Find list of loops originating from each node
    % Get the matrix describing all the loops from at node #index
    [FoundLoops,FoundNodes]=findpaths(index,index,[],[],Net);
    LoopCoeffList=[LoopCoeffList;FoundLoops];	% Append Coefficients of loops
    LoopNodeList=[LoopNodeList;FoundNodes];		% Append nodes visited for each loop
end



% Remove duplicate loops
[LoopCoeffList,LoopNodeList]=RemoveDuplicateLoops(LoopCoeffList,LoopNodeList);

%fprintf('\n\n- 1st Order Loop List -\n');
%print_paths(LoopCoeffList);


%*** Convert to nomenclature used by Pozars RF Book  ***
% P{n} represents the nth path found connecting start to stop
% P{n}.Coeff is an array with a list of the coefficients passed through in order
% P{n}.Nodes is an array listing the nodes passed through in order. Including the end nodes
% NOTE: A cell array is used because the path lengths can be different resulting in different sized arrays



% *** Make a cell array of P the different length paths ***
temp=size(PathCoeffList);			% Determine number of paths
NumberPaths=temp(1);
if (NumberPaths==0);
    fprintf('\n*** There are no paths connecting those nodes ***\n')
    return
end

for index=1:NumberPaths				% Do each path separately
    Coeff=PathCoeffList(index,:);			% Read Coefficients for a path
    P{index}.Coeff=Coeff(1:sum(Coeff>0));		% Strip trailing zeros and put in struct
    Node=PathNodeList(index,:);			% Read node list for a path
    P{index}.Node=[Node(1:sum(Coeff>0)),Stop];    % Append trailing zeros and put in struct.
    % Append the Stop node onto the end of the node list
end


% *** Make a cell array of the first order paths, each row a different order ***
% *** The first column contains the number of paths of that order
% L{LoopOrder}.NumberLoops = number of loops of this order
% L{LoopOrder}.Coeff{n} = Coefficients of nth loop
% L{LoopOrder}.Node{n}  = Nodes of nth loop



temp=size(LoopCoeffList);
NumberLoops=temp(1);				% Determine number of first order paths
L{1}.NumberLoops=NumberLoops;			% Set number of loops in the L{1} struct
for index=1:NumberLoops				% Do each loop seperately
    Coeff=LoopCoeffList(index,:);			% Read coefficients for that loop
    L{1}.Coeff{index}=Coeff(1:sum(Coeff>0));	% Strip Trailing zeros and put in struct
    Node=LoopNodeList(index,:);			% Read Node list for loop
    L{1}.Node{index}=[Node(1:sum(Coeff>0)),Node(1)]; % Strip trailing zeros and put in struct
    % Append Stop node (same as first node in list
end




%*** Determine nth order loops ***
n=1;						% n is the order of loops we are finding

while 1  					% Loop until an order of loop is reached that is empty
    n=n+1;					% Count which order we are on
    L{n}.NumberLoops=0;				% Assume no loops of this order
    
    % compare each first order loop with each n-1th loop. If non touching add to the
    % two combined to the nth loop.
    
    for first=1:L{1}.NumberLoops			% Take each first order loop
        for second=1:L{n-1}.NumberLoops		% Compare with each n-1th loop
            
            if not(AreTouchingLoops(L{1}.Node{first},L{n-1}.Node{second})) % Non Touching loops found
                % Determine if loop is a duplicate
                Duplicate=0;                       % A flag to indicate loop found is duplicate(0=not dup)
                for index=1:L{n}.NumberLoops       % Add this loop if it is not a duplicate entry
                    if IsSameLoop([L{1}.Coeff{first}, L{n-1}.Coeff{second}],L{n}.Coeff{index}) %Duplicate found
                        Duplicate=1;        		% Set the duplicate flag
                    end
                end
                if (Duplicate==0)                       % Add the loop if not a duplicate
                    L{n}.NumberLoops=L{n}.NumberLoops+1;	% Increment the number of loops of that order
                    % For Node and Coeff structs. Append a new array describing the loop of order n found
                    L{n}.Coeff{(L{n}.NumberLoops)}=[L{1}.Coeff{first}, L{n-1}.Coeff{second}];
                    L{n}.Node{(L{n}.NumberLoops)}=[L{1}.Node{first}, L{n-1}.Node{second}];
                end
            end
        end
    end
    
    if (L{n}.NumberLoops==0)			% If no loops of order n where found, then break
        break					% There will be no loops of order n+1
    end
end

% ***** Display File info *****
fprintf('\n-- Network Info --\n')
fprintf('Net File   : ');fprintf(NetFile);fprintf('\n');
fprintf('Start Node : %d\n',Start);
fprintf('Stop Node  : %d\n',Stop);


% ***** Display the paths found ******
fprintf('\n----- Paths -----\n')
for pathn=1:length(P)				% Look at each Path and display it's Coeff numbers
    fprintf('P%d : ',pathn);                       % on a different line
    fprintf('%d ',P{pathn}.Node); %%ECEC was .Coeff
    fprintf('\n');
end

% ***** Display all the loops found *****

for loop_order=1:length(L)-1            	% Look at each loop order (last order has no loops
    fprintf('\n- Order %d Loops -\n',loop_order)  % Print header describing loop order
    for loop_number=1:L{loop_order}.NumberLoops   % Look at each loop of that order
        fprintf('L%d%d : ',loop_order,loop_number)  % Display coefficients on a different line
        fprintf('%d ',L{loop_order}.Coeff{loop_number})
        fprintf('\n')
    end
end


% *******************************************************
% ************ Generate the final equation **************
% *******************************************************


% For now the equations are written in terms of the coefficient number : c#
% the coefficient strings will be substituted later

% Determine Numerator
Num='';				% Assume Numerator is empty to start
for pathn=1:length(P)            % Add Each path and related
    Num=sprintf('%s%s*(1', Num, CoeffToString(P{pathn}.Coeff));            %    Pn*(1 ...
    for order=1:length(L)-1       % Append each order of sums of non-touching loops
        % if order is odd order append a minus, otherwise a plus
        if (rem(order,2)==1)
            Num=sprintf('%s-',Num);
        else
            Num=sprintf('%s+',Num);
        end
        % Append the sum of all the nontouching loops that  don't touch the current path
        Num=[Num,PrintSumsNotTouching(L,order,P{pathn}.Node)];
    end
    Num=sprintf('%s)+',Num);   	% Close the bracket around paths list of sums
end



Num=Num(1:length(Num)-1);       % Remove the extra plus sign on the end.NOTE using /b screws up the symbolic math later

% Determine Denominator
Den='1';			% Denominator always start with a zero
for order=1:length(L)-1		% Add order subtract the sum of each orders loops
    % if order is odd order append a minus, otherwise a plus,
    if (rem(order,2)==1)
        Den=sprintf('%s-',Den);
    else
        Den=sprintf('%s+',Den);
    end
    %Add or subtract all the loops
    % KLUDGE: all the loops not 2 the path with nodes 999999999 are added
    %         That definetly should be all of them!
    Den=[Den,PrintSumsNotTouching(L,order,[9999999 999999])];  %Sums of all the loops of order order
end


fprintf('\nThe variables returned are strings describing\n')
fprintf('the numerator and Denominator of the transfer equation.\n')
fprintf('If you have the symbolic toolbox, use Denominator=sym(Denominator)\n');
fprintf('and Numerator=sym(Numerator) to make these symbolic equations.\n')
fprintf('You can now use simple(Numerator/Denominator) to boil the whole\n')
fprintf('thing down. You could also use simple(Numerator) to simplify the\n')
fprintf('Numerator on it'' own.\n\n')
% ****** Convert to Symbolic and do substitutions *******

for coeff_num=length(Coeff_Names):-1:1;	%Cycle through Coefficient array, substituting each one
    orig=sprintf('c%d',Net(coeff_num,1)); % for each line generate c[Coeff Number] to replace
    Den=strrep(Den,orig,Coeff_Names{coeff_num});	% Replace all the c#s with the strings from net file
    Num=strrep(Num,orig,Coeff_Names{coeff_num});
end % This loop had to count down so there was no risk of C12 being replace by C1
end

%*************************************************************************************************
function Touching=AreTouchingLoops(Nodes1,Nodes2)
%*************************************************************************************************
% This function takes two arrays describing two sets of nodes visited(each padded with zeros).
% Return 1 if they are they are touching loops.

% Determine length of loop arrays with zeros removed
Loop1Length=sum(Nodes1>0);
Loop2Length=sum(Nodes2>0);

for first=1:Loop1Length
    for second=1:Loop2Length
        if (Nodes1(first)==Nodes2(second))
            Touching=1;
            return;
        end
    end
end

Touching=0;
end


%*************************************************************************************************
function StrMult=CoeffToString(Coefficients)
%*************************************************************************************************
% StrMult=CoeffToString(Coefficients)
% Coefficients is an array with coefficients c1,c2..cN

N=length(Coefficients);     			% Get length of string
StrMult=sprintf('c%d',Coefficients(1));		% Start with first coefficient
for n=2:N		% Append each coefficent in list with * before it
    StrMult=[StrMult, sprintf('*c'),sprintf('%d',Coefficients(n))];
end
end


%*************************************************************************************************
function [PathUp,NodesUp]=findpaths(StartNode,StopNode,Path,Nodes,Net)
%*************************************************************************************************
%[PathUp,NodesUp]=findpaths(StartNode,StopNode,Path,Nodes,Net)
%
%Iterative function to find path between StartNode and StopNode. Net is the array with the network
%list in it. Path is the single path to date for a given route through the tree. PathUp is a
%list of all paths terminated below that node that are sucesfull.
%Nodes is the list of nodes tvaersed so far on the way down

% Determine number of coefficients in net
temp=size(Net);
NumberCoeff=temp(1,1);

PathUp=[];
NodesUp=[];

% Terminate branch and return nothing if the Nodes to date contains repetitions.
for index=1:NumberCoeff
    if not(isempty(Nodes))  % Only compare if the the Path has data in it
        if (sum(Nodes==index)>1)
            PathUp=[];
            %	fprintf('Repeated Node : ');
            %	fprintf('%d ',Nodes);
            %	fprintf('\n');
            return
        end
    end
end

% Terminate branch and return path if start and stop nodes are the same
if ((StartNode==StopNode) & (length(Path>1)))
    PathUp=Path;
    NodesUp=Nodes;
    %fprintf('Sucessfull Path : ');
    %fprintf('%d ',Path);
    %fprintf('\n');
    return
end


% Check for all branches leaving StartNode, and call another iteration for them
for index=1:NumberCoeff
    if (StartNode==Net(index,2))
        % Iterate with appended coeff to path and new startnode
        [FoundPath,FoundNodes]=findpaths(Net(index,3),StopNode,[Path,Net(index,1)],[Nodes,StartNode],Net);
        if not(isempty(FoundPath))  %Only append if not empty
            PathUp=[PathUp;[FoundPath,zeros(1,NumberCoeff+1-length(FoundPath))]];
            NodesUp=[NodesUp;[FoundNodes,zeros(1,NumberCoeff+1-length(FoundPath))]];
        end
    end
end
end

%*************************************************************************************************
function Same=IsSameLoop(Loop1,Loop2)
%*************************************************************************************************
% This function takes two arrays describing two loops(Can be padded with zeros if desired).
% Return 1 if they are they describe the same circular loop.

% Determine length of loop arrays with zeros removed
Loop1Length=sum(Loop1>0);		% Count all the non zero terms
Loop2Length=sum(Loop2>0);

%Return 0 if different length since the loops can't be the same!
if (Loop1Length~=Loop2Length)
    Same=0;
    return
end

%They are the same length so see if the contain the same nodes, but in any order.

% sort the nodes and subtract the two vectors. The resulting vector componments will all be zero, only if the lists contains the same values.

if (sum(abs(sort(Loop1)-sort(Loop2)))==0)
    Same=1;				% Loops are the same
else
    Same=0;				% Loops are different
end
end


%*************************************************************************************************
function Str=PrintSumsNotTouching(L,order,Pnodes)
%*************************************************************************************************
% Str=PrintSumsNotTouching(L,path)
% L is the array of structs containing all the first and higher order loops.
% Pnodes is the array of nodes that decibing a path
%
% The function returns a string with the sum off all the loops of order order
% that do not touch the path. The sum is contained within a set of brackets

No_NonTouching=1;  			% Flag set so indacate no nontouching loops found
Str=('(');				% Open the first bracket
for n=1:L{order}.NumberLoops		% Look at each loop of thet order
    if not(AreTouchingLoops(Pnodes,L{order}.Node{n}))   		%The loop doesn't touch the path
        Str=sprintf('%s%s+',Str,CoeffToString(L{order}.Coeff{n}));% So add its coefficients
        No_NonTouching=0;			% Set flag to indacet a nontouching loop was found
    end
end
Str=Str(1:(length(Str)-1));		% Remove the extra plus sign (or open bracket
% if not loops where found
Str=sprintf('%s)',Str);			% Append the closed bracket
%If the sum is zero return zero instead

if No_NonTouching==1			% If no loops foun then return '0' instead
    Str='0';
end
end

%*************************************************************************************************
function [LoopList,NodeList]=RemoveDuplicateLoops(LoopList,NodeList);
%*************************************************************************************************
% [LoopList,NodeList]=RemoveDuplicateLoops(LoopList,NodeList)
% Returns a subset of the LoopList matrix, where each duplicate row has been removed
% This function works on the initial loops description where each loop is a row in the
% the matrix padded with zeros to make them all the same length

temp=size(LoopList);		% Get number of loops
NumberLoops=temp(1);

% Compare each loop with all the loops after it. And remove its duplicates from after it.
first=1;			% The first loop
while (first<=NumberLoops)      % Loop until all the loops have been used as first loops
    second=first+1;		% Start the second loop to compare as being one ahead of first loop
    while (second<=NumberLoops)  % choose the second loop to compare as all the ones after first
        % Remove the extra loop found at the second loop index.
        % NOTE: A for loop was not used since the number of loops is not constant
        if (IsSameLoop(LoopList(first,:),LoopList(second,:))==1) %Remove row at second loop
            LoopList=[LoopList(1:second-1,:);LoopList(second+1:NumberLoops,:)];
            NodeList=[NodeList(1:second-1,:);NodeList(second+1:NumberLoops,:)];
            NumberLoops=NumberLoops-1;       % Decrement the number o loops
        else
            second=second+1;	% Only increment second if no loop was removed
            % Otherwise a new loop is move to second
        end
    end
    first=first+1;				% increment the first loop pointer
end
end
