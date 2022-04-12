classdef MINI2PWidget < dabs.resources.widget.Widget
    properties
        pbDistortionDetection
    end
    
    methods
        function obj = MINI2PWidget(hResource,hParent)
            obj@dabs.resources.widget.Widget(hResource,hParent);
            
            try
                obj.redraw();
            catch ME
                most.ErrorHandler.logAndReportError(ME);
            end
        end
        
        function delete(obj)
        end
    end
    
    methods
        function makePanel(obj,hParent)
            hFlow = most.gui.uiflowcontainer('Parent',hParent,'FlowDirection','TopDown','margin',0.001);
            most.gui.uicontrol('Parent',hFlow,'String','Distortion Detection','Callback',@(varargin)obj.distDetect);
            most.gui.uicontrol('Parent',hFlow,'String','Distortion Correction','Callback',@(varargin)obj.hResource.distortionCorrection);
        end
        
        function redraw(obj)
           
        end
        
        function distDetect(varargin)
            matlab.apputil.run('MINI2PDistortionDetection');
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
