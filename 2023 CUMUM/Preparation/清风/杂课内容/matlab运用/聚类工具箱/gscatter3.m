function [varargout]=gscatter3(x,y,z,g,clr,mrksle,mrksze,mrkfclr,leg,legloc)
%------------------------------HELP START----------------------------------
% gscatter3 - 3D scatter plot by group

% Copyright:
% V. Salai Selvam
% Associate Professor & Head
% Department of Electronics & Communication Engineering
% Sriram Engineering College, Perumalpattu - 602024

%------------------------------HELP FINIS----------------------------------

if nargin<4,
    disp('Insufficient no of input arguments.');
    return;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isvector(x)||~isvector(y)||~isvector(z),
    disp('x, y and z must be vectors.');
    return;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (length(x)~=length(y))||(length(y)~=length(z))||(length(z)~=length(x)),
    disp('x, y and z must be vectors of same size.');
    return;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=x(:)';
y=y(:)';
z=z(:)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~iscellstr(g),
    if ischar(g),
        g=cellstr(g)';
    else
        g=g(:)';
        g=cellstr(num2str(g'))';
    end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(g)~=length(x),
    disp('If g is a character array, the rows of g must be equal to the length of x, y and z OR');
    disp('if g is not a character array, the length of g must be equal to the length of x, y and z.');
    return;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gu=unique(g);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(gu)>128,
    disp('The function can only handle up to 128 groups.');
    return;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clrstr={'b','g','r','c','m','y','k'};
clrmaps={'jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink','lines'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('clr','var')||isempty(clr),
    if length(gu)<=7,
        clr=clrstr(1,1:length(gu));
    else
        clrmap=colormap([clrmaps{1,1},'(128)']);
        inds=fix(linspace(1,128,length(gu)));
        clr=clrmap(inds,:);
    end;
else
    if iscellstr(clr),
        if ~isvector(clr),
            disp('If clr is a cell array of strings, clr must be a vector.');
            return;
        end;
        clr=clr(:)';
        if (length(clr)~=length(g))&&(length(clr)~=length(gu)),
            disp('If clr is a cell array of strings, the length of clr must be equal to the length of x, y and z.');
            return;
        end;
        if length(clr)==length(g),
            clr=unique(clr);
        end;
        sm=0;
        for i=1:length(clr),
            sm=sm+sum(strcmpi(clr{1,i},clrstr));
        end;
        if sm==0,
            disp('Color string must be one of the following supported types:');
            disp(clrstr);
            return;
        end;
    elseif ischar(clr),
        if ~isvector(clr),
            disp('If clr is a character array, clr must be a vector.');
            return;
        end;
        clr=clr(:)';
        if sum(strcmpi(clr,clrmaps))~=0,
            clrmap=colormap([clr,'(128)']);
            inds=fix(linspace(1,128,length(gu)));
            clr=clrmap(inds,:);
        elseif length(clr)==length(gu),
            clr=cellstr(clr')';
            sm=0;
            for i=1:length(clr),
                sm=sm+sum(strcmpi(clr{1,i},clrstr));
            end;
            if sm==0,
                disp('Color string must be one of the following supported types:');
                disp(clrstr);
                return;
            end;
        elseif length(clr)==length(g),
            clr=unique(clr);
            clr=cellstr(clr')';
            sm=0;
            for i=1:length(clr),
                sm=sm+sum(strcmpi(clr{1,i},clrstr));
            end;
            if sm==0,
                disp('Color string must be one of the following supported types:');
                disp(clrstr);
                return;
            end;
        else
            disp('If clr is a character array, clr must be one of the following supported colormap names:');
            disp(colormaps);
            disp('OR a character array of the following supported color strings of length equal to');
            disp('the length of g or the no. of groups.');
            disp(clrstr);
            return;
        end;
    end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mrkslestr={'+','o','*','.','x','s','d','^','v','<','>','p','h'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('mrksle','var')||isempty(mrksle),
    if length(gu)>length(mrkslestr),
        mrksle=repmat(mrkslestr(1,1),1,length(gu));
    else
        mrksle=mrkslestr(1,1:length(gu));
    end;
else
    if ~iscellstr(mrksle)&&~ischar(mrksle),
        disp('Marker style must be either a cell array of strings or a character array.');
        return;
    end;
    if ~isvector(mrksle),
        disp('If mrksle is a cell array of strings or a character array, mrksle must be a vector.');
        return;
    end;
    if iscellstr(mrksle)||ischar(mrksle),
        if length(mrksle)~=length(g)&&length(mrksle)~=length(gu),
                  mrksle = repmat(mrksle,length(gu),1);  % 清风注释
            % disp('If mrksle is either a cell array of strings or a character array,');
            % disp('the length of mrksle must be either of length of g or of length equal to no of groups.');
            % return;
        end;
    end;
    mrksle=mrksle(:)';
    if length(mrksle)==length(g),
        mrksle=unique(mrksle);
        if ischar(mrksle),
            mrksle=cellstr(mrksle')';
        end;
    else
        if ischar(mrksle),
            mrksle=cellstr(mrksle')';
        end;
    end;
    sm=0;
    for i=1:length(mrksle),
        sm=sm+sum(strcmpi(mrksle{1,i},mrkslestr));
    end;
    if sm==0,
        disp('Marker style must be one of the following supported types:');
        disp(mrkslestr);
        return;
    end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('mrksze','var')||isempty(mrksze),
    mrksze=5;
end;
% 清风注释
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if ~isscalar(mrksze),
%     disp('Marker size must be a scalar.');
%     return;
% end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('mrkfclr','var')||isempty(mrkfclr),
    mrkfclr='none';
else
    if sum(strcmpi(mrkfclr,{'auto','none'}))==0,
        mrkfclr='auto';
    end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(lower(mrkfclr),'auto'),
    mrkfclr=clr;
else
    mrkfclr=repmat({mrkfclr},1,length(gu));
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=zeros(1,length(gu));

mrksze11 = mrksze;  % 清风添加

for i=1:length(gu),
    inds=strcmpi(gu{1,i},g);
    h(1,i)=plot3(x(inds),y(inds),z(inds));
    hold on;
    % 清风添加
     if ~isscalar(mrksze11)
         mrksze = mrksze11(i);
     end
    if ~iscellstr(clr)&&~ischar(clr),
        set(h(1,i),'LineStyle','none','Color',clr(i,:),'Marker',mrksle{1,i},'MarkerSize',mrksze,'MarkerFaceColor',mrkfclr(i,:));
    else
        set(h(1,i),'LineStyle','none','Color',clr{1,i},'Marker',mrksle{1,i},'MarkerSize',mrksze,'MarkerFaceColor',mrkfclr{1,i});
    end;
end;
if ~exist('leg','var')||isempty(leg),
    leg=1;
end;
if ~exist('legloc','var')||isempty(legloc),
    legloc='NorthEastOutside';
end;
if leg==1,
    legend(gu,'Location',legloc);
end;
box on; grid on; hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout~=0,
    varargout={h};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2012, Salai Selvam V
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
