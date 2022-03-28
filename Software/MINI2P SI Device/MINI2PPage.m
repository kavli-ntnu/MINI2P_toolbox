classdef MINI2PPage < dabs.resources.configuration.ResourcePage
    properties
        etLUT
        etScope
        etSystem
        etObjective
    end
    
    methods
        function obj = MINI2PPage(hResource,hParent)
            obj@dabs.resources.configuration.ResourcePage(hResource,hParent);
        end
        
        function makePanel(obj,hParent)
            most.gui.uicontrol('Parent',hParent,'Style','text','RelPosition', [10 25 300 17],'Tag','txLUT','String','Transform Matrix file (leave empty if not exist)','HorizontalAlignment','left');
            obj.etLUT = most.gui.uicontrol('Parent',hParent,'Style','edit','RelPosition', [10 45 290 20],'Tag','etLUT','HorizontalAlignment','left');
            most.gui.uicontrol('Parent',hParent,'RelPosition', [310 45 70 20],'String','Open','Tag','pbLUT','Callback',@(varargin)obj.openFile());
            
            most.gui.uicontrol('Parent',hParent,'Style','text','RelPosition', [10 70 120 17],'Tag','txSystem','String','MINI2P System name?','HorizontalAlignment','left');
            obj.etSystem = most.gui.uicontrol('Parent',hParent,'Style','edit','RelPosition', [10 90 90 17],'Tag','etSystem','String','Which system?','HorizontalAlignment','left');
            
            most.gui.uicontrol('Parent',hParent,'Style','text','RelPosition', [10 115 120 17],'Tag','txScope','String','MINI2P scope name?','HorizontalAlignment','left');
            obj.etScope = most.gui.uicontrol('Parent',hParent,'Style','edit','RelPosition', [10 135 90 17],'Tag','etScope','String','Which MINI2P?','HorizontalAlignment','left');
            
            most.gui.uicontrol('Parent',hParent,'Style','text','RelPosition', [10 160 120 17],'Tag','txObjective','String','Objective type?','HorizontalAlignment','left');
            obj.etObjective = most.gui.uicontrol('Parent',hParent,'Style','edit','RelPosition', [10 180 90 17],'Tag','etObjective','String','Which objective?','HorizontalAlignment','left');
            
        end
        
        function redraw(obj)
            obj.etLUT.String = obj.hResource.transformMatrixDirectory;
            obj.etSystem.String = obj.hResource.system;
            obj.etScope.String = obj.hResource.scope;
            obj.etObjective.String = obj.hResource.objective;
        end
        
        function apply(obj)
            most.idioms.safeSetProp(obj.hResource,'transformMatrixDirectory',obj.etLUT.String);
            most.idioms.safeSetProp(obj.hResource,'system',obj.etSystem.String);
            most.idioms.safeSetProp(obj.hResource,'scope',obj.etScope.String);
            most.idioms.safeSetProp(obj.hResource,'objective',obj.etObjective.String);
            
            obj.hResource.saveMdf();
            obj.hResource.reinit();
        end
        
        function remove(obj)
            obj.hResource.deleteAndRemoveMdfHeading();
        end
        
        function openFile(obj)

            
            [file,path] = uigetfile('*.*','Select Transform Matrix File');
            
            if isnumeric(file)
                return % user cancelled
            end
            
            file = fullfile(path,file);
            obj.etLUT.String = file;
        end
    end
end



% ----------------------------------------------------------------------------
% Version1: (2022.03.22) at Kavli Institute, Trondheim, NTNU
% The MINI2P distortion detection and correction method is from the paper Zong, et.al,Cell,2022
% The MINI2P distortion detection GUI was written by Weijian Zong in the Moser lab
% The SI device of MINI2P was modified from standard device in SI by Weijian Zong in the Moser lab, Mitchell Sandoe and Jacob Franklin in Vidrio Technologies, LLC
% Updated version of this device and the user manual can be found in https://github.com/kavli-ntnu/MINI2P_toolbox
% ----------------------------------------------------------------------------
% Copyright (C) 2021 Vidrio Technologies, LLC
% 
% ScanImage (R) 2021 is software to be used under the purchased terms
% Code may be modified, but not redistributed without the permission
% of Vidrio Technologies, LLC
% 
% VIDRIO TECHNOLOGIES, LLC MAKES NO WARRANTIES, EXPRESS OR IMPLIED, WITH
% RESPECT TO THIS PRODUCT, AND EXPRESSLY DISCLAIMS ANY WARRANTY OF
% MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.
% IN NO CASE SHALL VIDRIO TECHNOLOGIES, LLC BE LIABLE TO ANYONE FOR ANY
% CONSEQUENTIAL OR INCIDENTAL DAMAGES, EXPRESS OR IMPLIED, OR UPON ANY OTHER
% BASIS OF LIABILITY WHATSOEVER, EVEN IF THE LOSS OR DAMAGE IS CAUSED BY
% VIDRIO TECHNOLOGIES, LLC'S OWN NEGLIGENCE OR FAULT.
% CONSEQUENTLY, VIDRIO TECHNOLOGIES, LLC SHALL HAVE NO LIABILITY FOR ANY
% PERSONAL INJURY, PROPERTY DAMAGE OR OTHER LOSS BASED ON THE USE OF THE
% PRODUCT IN COMBINATION WITH OR INTEGRATED INTO ANY OTHER INSTRUMENT OR
% DEVICE.  HOWEVER, IF VIDRIO TECHNOLOGIES, LLC IS HELD LIABLE, WHETHER
% DIRECTLY OR INDIRECTLY, FOR ANY LOSS OR DAMAGE ARISING, REGARDLESS OF CAUSE
% OR ORIGIN, VIDRIO TECHNOLOGIES, LLC's MAXIMUM LIABILITY SHALL NOT IN ANY
% CASE EXCEED THE PURCHASE PRICE OF THE PRODUCT WHICH SHALL BE THE COMPLETE
% AND EXCLUSIVE REMEDY AGAINST VIDRIO TECHNOLOGIES, LLC.
% ----------------------------------------------------------------------------
