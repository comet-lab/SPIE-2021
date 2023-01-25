classdef EndoWrist < handle
    %ENDOWRIST Container class to hold a endoscope/ wrist combination
    %   Detailed explanation goes here
    % 
    %   Author: Jesse F d'Almeida <jfdalmeida@wpi.edu>
    %
    %   Last revision: 5/29/2020
    
    properties
        endo        % Endoscope object
        wrist       % Wrist object
        
        % mesh models
        robotModel  % complete model
        endoModel   % all endoscope meshes
        wristModel  % all wrist meshes
        
        % Transformations
        transformations
        pose
     
    end
    
    methods
        function self = EndoWrist(IID, IOD, OID, OOD, g_inner, g_outer, n, cutouts)
            %ENDOWRIST Construct an instance of this class
            %   Creates instances of endoscope and wrist given config
            self.endo = Endoscope();
            self.wrist = Wrist(IID, IOD, OID, OOD, g_inner, g_outer, n, cutouts);
        end
        
        function fwkine(self, q, baseTransform)
            %FWKINE creates the homogeneous transformation matrix for base
            %to tip
            %   q [kappa a theta dz dl phi tau] = configuration
            %       k [1/m] = curvature
            %       theta [rad] = base rotation
            %       dz [m] = base translation
            %       dl [m] = tendon displacement
            %       phi [rad] = wrist base rotation
            %       tau [m] = wrist advancement
            %
            %  baseTransform [4x4 matrix] *optional* homogenous transformation matrix
            %   where the endoscope begins
            
            % default baseTransform of none
            if ~exist('baseTransform', 'var')
                baseTransform = eye(4);
            end
            
            % parse the configuration
            q_endo = q(1:3);
            q_wrist = q(4:6);
            
            self.endo.fwkine(q_endo, baseTransform);     % fwkin endo
            
            wristbase = self.endo.wristT;              % wrist base
            self.wrist.fwkine(q_wrist, wristbase);      % fwkin wrist
            
            % add all transformations
            all_t = cat(3, self.endo.transformations, self.wrist.transformations);
            self.transformations = all_t;
            
            all_p = [self.endo.pose self.wrist.pose];
            self.pose = all_p(:,end);
        end
        
        function robotModel = makePhysicalModel(self)
            % ROBOTMODEL generate mesh models for the endoscope and wrist
            
            % create meshes
            endoModel = self.endo.makePhysicalModel();
            wristModel = self.wrist.makePhysicalModel();
            
            self.endoModel = endoModel;
            self.wristModel = wristModel;
            
            % add meshes to model
            robotModel.surface.Xe = endoModel.surface.X;
            robotModel.surface.Ye = endoModel.surface.Y;
            robotModel.surface.Ze = endoModel.surface.Z;
            
            robotModel.surface.Xw = wristModel.surface.X;
            robotModel.surface.Yw = wristModel.surface.Y;
            robotModel.surface.Zw = wristModel.surface.Z;
            
            robotModel.surface.X = [endoModel.surface.X(:); wristModel.surface.X(:)];
            robotModel.surface.Y = [endoModel.surface.Y(:); wristModel.surface.Y(:)];
            robotModel.surface.Z = [endoModel.surface.Z(:); wristModel.surface.Z(:)];
            
           
            %return
            self.robotModel = robotModel;
        end
    end
end

