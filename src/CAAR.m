classdef CAAR
    %CAAR A class representing a concentric agonist-antagonist robot
    %   This class stores information and performs important calculations
    %   for a CAAR. This class can perform kinematics and statics
    %   calculations.

    properties
        IID;                % [mm]  inner tube inner diameter
        IOD;                % [mm]  inner tube outer diameter
        OID;                % [mm]  outer tube inner diameter
        OOD;                % [mm]  outer tube outer diameter
        n;                  % []    number of notches
        h;                  % [mm]  height of notches
        g_inner             % [mm]  depth of inner tube notches
        g_outer             % [mm]  depth of outer tube notches
        orientations;       % [rad] orientation of notches
        c;                  % [mm]  distance of uncut sections between notches
    end

    methods

        function obj = CAAR(IID, IOD, OID, OOD, n, h, g_inner, g_outer, orientations, c)
            obj.IID = IID;
            obj.IOD = IOD;
            obj.OID = OID;
            obj.OOD = OOD;
            obj.n = n;
            obj.h = h;
            obj.g_inner = g_inner;
            obj.g_outer = g_outer;
            obj.orientations = orientations;
            obj.c = c;
        end

        function deltaL = displacementNeeded(obj, angle)
            % DISPLACEMENTNEEDED Calculates the inner tube displacement
            % needed for this CAAR to achieve the given bending angle
            %       deltaL = displacementNeeded(angle) calculates the
            %       distance that the inner tube needs to be pulled back
            %       wrt the outer tube to achieve the given bending angle
            %       in radians.
            
            % Initialize deltaL
            deltaL = 0;

            % Find the neutral bending plane of the inner tube
            % Formulas taken from Oliver-Butler et al, Swaney et al
            ro_inner = obj.IOD;
            ri_inner = obj.IID;
            phio_inner = 2 * acos((obj.g_inner - ro_inner) / ro_inner);
            phii_inner = 2 * acos((obj.g_inner - ro_inner) / ri_inner);
            yo_inner = (4 * ro_inner * (sin(phio_inner/2))^3) / (3 * (phio_inner - sin(phio_inner)));
            yi_inner = (4 * ri_inner * (sin(phii_inner/2))^3) / (3 * (phii_inner - sin(phii_inner)));
            Ao_inner = (ro_inner^2 * (phio_inner - sin(phio_inner))) / (2);
            Ai_inner = (ri_inner^2 * (phii_inner - sin(phii_inner))) / (2);
            y_inner = (yo_inner * Ao_inner - yi_inner * Ai_inner) / (Ao_inner - Ai_inner);

            % Find the neutral bending plane of the outer tube
            % Formulas taken from Oliver-Butler et al, Swaney et al
            ro_outer = obj.OOD;
            ri_outer = obj.OID;
            phio_outer = 2 * acos((obj.g_outer - ro_outer) / ro_outer);
            phii_outer = 2 * acos((obj.g_outer - ro_outer) / ri_outer);
            yo_outer = (4 * ro_outer * (sin(phio_outer/2))^3) / (3 * (phio_outer - sin(phio_outer)));
            yi_outer = (4 * ri_outer * (sin(phii_outer/2))^3) / (3 * (phii_outer - sin(phii_outer)));
            Ao_outer = (ro_outer^2 * (phio_outer - sin(phio_outer))) / (2);
            Ai_outer = (ri_outer^2 * (phii_outer - sin(phii_outer))) / (2);
            y_outer = (yo_outer * Ao_outer - yi_outer * Ai_outer) / (Ao_outer - Ai_outer);

            % Relationship between total bending angle, inner tube
            % displacement, and neutral bending plane.
            % Formula taken from Oliver-Butler et al
            deltaL = angle * (y_inner + y_outer);

        end

        function HTMs = fkine(obj, q)
            % FKINE Calculate the forward kinematics of a wrist given arc parameters q
            %       HTMs = fkine(q) calculates the set of all intermediate
            %       homogeneous transformation matrices for the wrist given arc
            %       parameters q. These matrices are returned in the form
            %       [T0uncut, T1cut, T1uncut, ... , Tncut, Tnuncut].
            %       Note that fkine(q) assumes that all notches have an
            %       orientation of 0 radians on the tube.


            % Extract the arc parameters from the input vector q
            % These represent the 3 DoF of the tube
            deltaL  = q(1) / obj.n; % inner tube displacement [mm]
            % Positive displacement means inner
            % tube is *pulled*
            alpha   = q(2); % rotation at base of tube about z axis [deg]
            tau     = q(3); % linear translation of entire body along z axis [mm]


            % Initialize a matrix to store all the transformation matrices.
            % This matrix is initialized with the given rotation alpha
            % about the z axis and the given translation tau along the z
            % axis. This first 4x4 matrix is T0uncut, or the matrix
            % representing the transformation of the 0th uncut section of
            % the tube.
            T0uncut = [cosd(alpha) -sind(alpha) 0 0;
                sind(alpha) cosd(alpha) 0 0;
                0 0 1 tau;
                0 0 0 1];
            allHTMs = T0uncut;


            % These 2 formulas taken from original CAAR paper, which
            % doesn't help much for figuring out the transform per notch
            %             % Find the total angle that the CAAR bends
            %             theta = deltaL / (y_outer + y_inner);
            %
            %             % Find the total arc length of the CAAR
            %             s = (obj.n * obj.h(1) - )

            % For each notch in the tube, calculate the trnasformation
            % matrix of the cut notch and the uncut section immediately
            % following the notch
            for j=1:obj.n

                % Uncut height of the notch
                h_j = obj.h;

                %                 ***From last year's HW 2***
                %                 % Formulas taken from Swaney et al to find y, the distance
                %                 % from the center of the tube to its neutral bending plane
                %                 phio = 2 * acos((obj.widths(j) - ro) / ro);
                %                 phii = 2 * acos((obj.widths(j) - ro) / ri);
                %                 yo = (4 * ro * (sin(phio/2))^3) / (3 * (phio - sin(phio)));
                %                 yi = (4 * ri * (sin(phii/2))^3) / (3 * (phii - sin(phii)));
                %                 Ao = (ro^2 * (phio - sin(phio))) / (2);
                %                 Ai = (ri^2 * (phii - sin(phii))) / (2);
                %                 y = (yo * Ao - yi * Ai) / (Ao - Ai);
                %
                %                 % Calculate the curvature and arc length using the formulas
                %                 % given in lecture
                %                 K = (deltaL) / (h_j * (ri + y) - deltaL * y);
                %                 l_j = h_j / (1 + K * y);

                % Same idea as a single notched-tube, but now there's 2
                % tubes

                % Find the neutral bending plane of the inner tube
                % Formulas taken from Oliver-Butler et al, Swaney et al
                ro_inner = obj.IOD;
                ri_inner = obj.IID;
                phio_inner = 2 * acos((obj.g_inner - ro_inner) / ro_inner);
                phii_inner = 2 * acos((obj.g_inner - ro_inner) / ri_inner);
                yo_inner = (4 * ro_inner * (sin(phio_inner/2))^3) / (3 * (phio_inner - sin(phio_inner)));
                yi_inner = (4 * ri_inner * (sin(phii_inner/2))^3) / (3 * (phii_inner - sin(phii_inner)));
                Ao_inner = (ro_inner^2 * (phio_inner - sin(phio_inner))) / (2);
                Ai_inner = (ri_inner^2 * (phii_inner - sin(phii_inner))) / (2);
                y_inner = (yo_inner * Ao_inner - yi_inner * Ai_inner) / (Ao_inner - Ai_inner);

                % Find the neutral bending plane of the outer tube
                % Formulas taken from Oliver-Butler et al, Swaney et al
                ro_outer = obj.OOD;
                ri_outer = obj.OID;
                phio_outer = 2 * acos((obj.g_outer - ro_outer) / ro_outer);
                phii_outer = 2 * acos((obj.g_outer - ro_outer) / ri_outer);
                yo_outer = (4 * ro_outer * (sin(phio_outer/2))^3) / (3 * (phio_outer - sin(phio_outer)));
                yi_outer = (4 * ri_outer * (sin(phii_outer/2))^3) / (3 * (phii_outer - sin(phii_outer)));
                Ao_outer = (ro_outer^2 * (phio_outer - sin(phio_outer))) / (2);
                Ai_outer = (ri_outer^2 * (phii_outer - sin(phii_outer))) / (2);
                y_outer = (yo_outer * Ao_outer - yi_outer * Ai_outer) / (Ao_outer - Ai_outer);


                theta_j = deltaL / (y_inner + y_outer);
                l_j = obj.h - y_outer * theta_j;
                K_j = theta_j / l_j;


                %                 % Here is the disconnect. How can we get from y_inner and
                %                 % y_outer to K (curvature of the notch) and l (arc length
                %                 % of the central axis of the notch)
                %
                %
                %                 % Formulas taken from second CAAR paper.
                %                 K_j = K_oj / (1 + y_outer * K_oj);       % Curvature of this notch
                %                 l_j = h_j / (1 - y_outer * K_j);       % Arc length of this notch


                % Robot independent kinematics below
                % Now that we have found the arc parameters, we can create
                % the HTMs. These arc -> task transforms are independent of
                % the implementation of the robot.

                % Use the template given in lecture to create the
                % transformation matrix of the cut section. Handle the edge
                % case where K=0 by treating the section as uncut.
                if abs(K_j) > 1e-6
                    Tcut = [cos(K_j*l_j) 0 sin(K_j*l_j) (1-cos(K_j*l_j))/K_j;
                        0 1 0 0;
                        -sin(K_j*l_j) 0 cos(K_j*l_j) sin(K_j*l_j)/K_j;
                        0 0 0 1];
                else
                    Tcut = [1 0 0 allHTMs(1, 8*j-5);
                        0 1 0 allHTMs(2, 8*j-5);
                        0 0 1 allHTMs(3, 8*j-5);
                        0 0 0 1];
                end

                % Calculate the formula of the uncut section. This section
                % is not curved, so its transformation is just a
                % translation along the previous z axis
                Tuncut = [1 0 0 0;
                    0 1 0 0;
                    0 0 1 obj.c;
                    0 0 0 1];

                % Append the transformation matrices calculated in this
                % iteration of the loop onto the matrix of all HTMs
                allHTMs = [allHTMs, Tcut, Tuncut];
            end

            HTMs = allHTMs;

        end


        function maxStrain = calcMaxStrain(obj, q)
            % CALCMAXSTRAIN Calculates the maximum strain of the wrist with given arc parameters q
            %
            %   maxStrain = calcMaxStrain(q) calculates the maximum strain
            %   (given as a decimal between 0 and 1) that the wrist will
            %   experience when subjected to the arc parameters q


            % Extract the arc parameters from the input vector q0
            deltaL  = q(1) / obj.n; % tendon displacement [mm]
            alpha   = q(2); % rotation at base of tube about z axis [deg]
            tau     = q(3); % linear translation of entire body along z axis [mm]

            % Calculate the inner and outer radius of the tube using the
            % inner and outer diameter
            ro_outer = obj.OOD / 2;
            ri_outer = obj.OID / 2;
            ro_inner = obj.IOD / 2;
            ri_inner = obj.IID / 2;

            % Initialize the maximum strain to 0
            maxStrain = [0, 0];

            % Calculate the strain at each notch in the tube
            for i=1:obj.n

                % Uncut height of the notch
                h_j = obj.h;

                %                 % Formulas taken from Swaney et al to find y, the distance
                %                 % from the center of the tube to its neutral bending plane
                %                 phio = 2 * acos((obj.widths(i) - ro) / ro);
                %                 phii = 2 * acos((obj.widths(i) - ro) / ri);
                %                 yo = (4 * ro * sin(phio/2)^3) / (3 * (phio - sin(phio)));
                %                 yi = (4 * ri * sin(phii/2)^3) / (3 * (phii - sin(phii)));
                %                 Ao = (ro^2 * (phio - sin(phio))) / (2);
                %                 Ai = (ri^2 * (phii - sin(phii))) / (2);
                %                 y = (yo * Ao - yi * Ai) / (Ao - Ai);

                % Find the neutral bending plane of the inner tube
                % Formulas taken from Oliver-Butler et al, Swaney et al
                phio_inner = 2 * acos((obj.g_inner - ro_inner) / ro_inner);
                phii_inner = 2 * acos((obj.g_inner - ro_inner) / ri_inner);
                yo_inner = (4 * ro_inner * (sin(phio_inner/2))^3) / (3 * (phio_inner - sin(phio_inner)));
                yi_inner = (4 * ri_inner * (sin(phii_inner/2))^3) / (3 * (phii_inner - sin(phii_inner)));
                Ao_inner = (ro_inner^2 * (phio_inner - sin(phio_inner))) / (2);
                Ai_inner = (ri_inner^2 * (phii_inner - sin(phii_inner))) / (2);
                y_inner = (yo_inner * Ao_inner - yi_inner * Ai_inner) / (Ao_inner - Ai_inner);

                % Find the neutral bending plane of the outer tube
                % Formulas taken from Oliver-Butler et al, Swaney et al
                phio_outer = 2 * acos((obj.g_outer - ro_outer) / ro_outer);
                phii_outer = 2 * acos((obj.g_outer - ro_outer) / ri_outer);
                yo_outer = (4 * ro_outer * (sin(phio_outer/2))^3) / (3 * (phio_outer - sin(phio_outer)));
                yi_outer = (4 * ri_outer * (sin(phii_outer/2))^3) / (3 * (phii_outer - sin(phii_outer)));
                Ao_outer = (ro_outer^2 * (phio_outer - sin(phio_outer))) / (2);
                Ai_outer = (ri_outer^2 * (phii_outer - sin(phii_outer))) / (2);
                y_outer = (yo_outer * Ao_outer - yi_outer * Ai_outer) / (Ao_outer - Ai_outer);

                %                 % Calculate the curvature using the formula given in
                %                 % lecture
                %                 K = deltaL / (h_j * (ro + y) - deltaL * y);

                theta_j = deltaL / (y_inner + y_outer);
                l_j = obj.h - y_outer * theta_j;
                K_j = theta_j / l_j;

                % Calculate the strain using the formula given in lecture
                strain_outer = K_j * (ro_outer - y_outer) / (1 + K_j * y_outer);
                strain_inner = K_j * (ro_inner - y_inner) / (1 + K_j * y_inner);

                % Update the max strain if applicable
                if abs(strain_outer) > abs(maxStrain(1))
                    maxStrain(1) = strain_outer;
                end

                % Update the max strain if applicable
                if abs(strain_inner) > abs(maxStrain(2))
                    maxStrain(2) = strain_inner;
                end

            end

            % If the max strain is greater than 0.08, issue a warning. 0.08
            % is the max strain that nitinol can experience before it
            % deforms permanently
            if maxStrain(1) > 0.08 || maxStrain(2) > 0.08
                disp("WARNING: The maximum strain for this arm is greater than the recoverable strain of 8%. Reconsider your design.");
            end

        end

        function q = IK(obj, X, accuracy, posbend)
            % IK Calculates the inverse kinematics of the CAAR given a
            % desired task space position X (x, y, z) to a given accuracy.
            % The positive z axis is considered to point downward into the
            % throat.
            %       q = IK(X, accuracy) calculates the joint space
            %       configuration (deltaL, alpha, tau) required to achieve
            %       a given task space position (within the given accuracy).
            %       Assumes all notches have the same geometry.
            %       Calculates joint configuration by starting from home
            %       configuration and incrementing the joint variables
            %       using proportional control until the robot reaches the
            %       desired position.
            %       Bends in the positive direction if posbend is true, and
            %       in the negative direction if posbend is false.

            % Extract the coordinates from the input vector
            x_desired = X(1);
            y_desired = X(2);
            z_desired = X(3);

            % Convert to polar (actually cylindrical) coords
            r_desired = sqrt(x_desired .^ 2 + y_desired .^ 2);
            theta_desired = rad2deg(atan2(y_desired, x_desired));

            % Initialize the output vector to the default configuration
            deltaL = 0;
            alpha = 0;
            tau = 0;
            q = [deltaL; alpha; tau];

            % Of the 3 task-space DoF (r, theta, z),
            %       r depends on deltaL
            %       theta is equivalent to alpha
            %       z depends on deltaL and tau

            % We focus on achieving desired r first so we can solve for
            % deltaL.

            % Get the task space position based on the current joint space
            % config. The current r will be the x position that FK outputs
            HTMs_curr = fkine(obj, q);
            X_curr = getXfromHTMs(obj, HTMs_curr);
            r_curr = X_curr(1);

            % Define K_p of sorts
            Gain_deltaL = 0.2;

            % Adjust for the negative bending case if needed
            if posbend == false
                alpha = rad2deg(pi);
                Gain_deltaL = -1 * Gain_deltaL;
            end

            % Calculate error in r
            r_error = r_desired - r_curr;

            % Iterate until r_error has a small enough magnitude
            while abs(r_error) > accuracy

                % Update value of tau based on P control
                deltaL = deltaL + Gain_deltaL * r_error;

                % Check FK again to get the updated r value
                q = [deltaL; alpha; tau];
                HTMs_curr = fkine(obj, q);
                X_curr = getXfromHTMs(obj, HTMs_curr);
                r_curr = X_curr(1);

                % Calculate the new r error
                r_error = r_desired - r_curr;

            end

            % Now that we know how to achieve the desired r through
            % bending (deltaL), we can turn our attention to z. Now we
            % solve for tau.

            % Get the task space position based on the current joint space
            % config. The current z will be the z position that FK outputs
            HTMs_curr = fkine(obj, q);
            X_curr = getXfromHTMs(obj, HTMs_curr);
            z_curr = X_curr(3);

            % Calculate error in z
            z_error = z_desired - z_curr;

            % Define K_p of sorts
            Gain_tau = 0.35;

            % Iterate until z_error has a small enough magnitude
            while abs(z_error) > accuracy

                % Update value of tau based on P control
                tau = tau + Gain_tau * z_error;

                % Check FK again to get the updated r value
                q = [deltaL; alpha; tau];
                HTMs_curr = fkine(obj, q);
                X_curr = getXfromHTMs(obj, HTMs_curr);
                z_curr = X_curr(3);

                % Calculate the new r error
                z_error = z_desired - z_curr;

            end

            % We have solved for deltaL and tau. Now we solve for alpha.
            % If we are bending in the positive direction, alpha is the
            % same as theta_desired. If we are bending in the negative
            % direction, alpha is pi plus theta_desired.

            alpha = alpha + theta_desired;

            % We don't want unnecessary bending, so make sure alpha is
            % between -pi and pi.

            if alpha > 180
                alpha = alpha - 360;
            elseif alpha < -180
                alpha = alpha + 360;
            end

            % Assemble the final values of deltaL, alpha, and tau into an
            % output vector q

            q = [deltaL, alpha, tau];

        end

        function X = getXfromHTMs(obj, HTMs)
            % GETXFROMHTMS Extracts the task space position in Cartesian
            %              coordinates of the distal end of the final
            %              "link" as given by a series of intermediate HTMs.

            % Get the product of all the HTMs
            total_HTM = getTotalHTM(obj, HTMs);

            % Extract the Cartesian coords of the distal end
            X = total_HTM(1:3, 4);

        end

        function total_HTM = getTotalHTM(obj, HTMs)
            % GETTOTALHTM Takes in a series of HTMs and multiplies them
            %               together to get one total HTM

            % Figure out how many HTMs we were given
            size_HTMs = size(HTMs);
            num_HTMs = size_HTMs(2) / 4;

            % Initialize a matrix to hold the result
            total_HTM = eye(4);

            % Multiply the HTMs together
            for i=1:num_HTMs
                total_HTM = total_HTM * HTMs(:, 4*i-3:4*i);
            end

        end


%         function HTMs = fkine2(obj, q)
%             % FKINE Calculate the forward kinematics of a wrist given arc parameters q
%             %       HTMs = fkine(q) calculates the set of all intermediate
%             %       homogeneous transformation matrices for the wrist given arc
%             %       parameters q. These matrices are returned in the form
%             %       [T0uncut, T1cut, T1uncut, ... , Tncut, Tnuncut].
%             %       Note that fkine2(q) does NOT assume that all notches have an
%             %       orientation of 0 radians on the tube.
% 
%             % Extract the arc parameters from the input vector q
%             deltaL  = q(:, 1); % tendon displacement [mm]
%             alpha   = q(1, 2); % rotation at base of tube about z axis [deg]
%             tau     = q(1, 3); % linear translation of entire body along z axis [mm]
% 
%             % Consolidate the array of notch orientations into an array of
%             % the unique orientations.
%             uniqueOrientations = unique(obj.orientations);
% 
%             % Count how many times each of the orientations appears
%             for i = 1:length(uniqueOrientations)
%                 countOfOrientations(i,1) = sum(obj.orientations==uniqueOrientations(i)); % number of times each unique value is repeated
%             end
% 
%             % Create expandedDeltaL: An array of length obj.n holding the
%             % exact tendon displacement for each notch. This accounts for
%             % the multiple wires and the distribution of deltaL among
%             % notches of the same orientation.
%             for i=1:obj.n
%                 index = find(uniqueOrientations==obj.orientations(i));
%                 expandedDeltaL(i) = deltaL(index) / countOfOrientations(index);
%             end
% 
%             % Use the tube's inner and outer diameter to calculate its
%             % inner and outer radius
%             ro = obj.OOD / 2;
%             ri = obj.IID / 2;
% 
%             % Initialize a variable 'allHTMs' to store a matrix of
%             % homogeneous transformation matrices. Start it off by
%             % translating in the z axis based on tau and rotating about the
%             % z axis based on alpha
%             allHTMs = [cosd(alpha) -sind(alpha) 0 0;
%                 sind(alpha) cosd(alpha) 0 0;
%                 0 0 1 tau;
%                 0 0 0 1];
% 
%             % Keep track of the previous notch's orientation. This way, we
%             % know how much we have to rotate to get to the current notch.
%             prevOrientation = 0;
% 
%             % Calculate the transformation matrix (both cut and uncut) for
%             % each notch in the tube
%             for i=1:obj.n
% 
%                 % Uncut height of the notch
%                 h = obj.heights(i);
% 
%                 % Formulas taken from Swaney et al to find y, the distance
%                 % from the center of the tube to its neutral bending plane
%                 phio = 2 * acos((obj.widths(i) - ro) / ro);
%                 phii = 2 * acos((obj.widths(i) - ro) / ri);
%                 yo = (4 * ro * sin(phio/2)^3) / (3 * (phio - sin(phio)));
%                 yi = (4 * ri * sin(phii/2)^3) / (3 * (phii - sin(phii)));
%                 Ao = (ro^2 * (phio - sin(phio))) / (2);
%                 Ai = (ri^2 * (phii - sin(phii))) / (2);
%                 y = (yo * Ao - yi * Ai) / (Ao - Ai);
% 
%                 % Calculate the curvature using the formula given in
%                 % lecture
%                 K = expandedDeltaL(i) / (h * (ri + y) - expandedDeltaL(i) * y);
% 
%                 % Calculate the arc length using the formula given in
%                 % lecture
%                 s = h / (1 + K * y);
% 
%                 % What is the current orientation of the notch? To be
%                 % compared with the orientation of the previous notch
%                 currOrientation = obj.orientations(i);
% 
%                 % Create a transformation matrix (to be incorporated in
%                 % Tcut) that handles any changes in orientation since the
%                 % last notch
%                 Torient = [cos(currOrientation - prevOrientation) -sin(currOrientation - prevOrientation) 0 0;
%                     sin(currOrientation - prevOrientation) cos(currOrientation - prevOrientation) 0 0;
%                     0 0 1 0
%                     0 0 0 1];
% 
%                 prevOrientation = currOrientation;
% 
%                 % Calculate the transformation matrix of the cut section.
%                 % Be sure to account for the orientation of the notch; do
%                 % this by multiplying by Torient.
% 
%                 % If K is not zero (+/- a little tolerance), then calculate
%                 % the transformation matrix using the formula given in
%                 % lecture and in provided papers
%                 if abs(K) > 1e-4
%                     Tcut = Torient * [cos(K*s) 0 sin(K*s) (1-cos(K*s))/K;
%                         0 1 0 0;
%                         -sin(K*s) 0 cos(K*s) sin(K*s)/K;
%                         0 0 0 1];
% 
%                     % If K is zero, then the notch is not bent, so we can just
%                     % translate along the previous frame's z axis. No rotation
%                     % required
%                 else
%                     Tcut = Torient * [1 0 0 allHTMs(1, 8*i-5);
%                         0 1 0 allHTMs(2, 8*i-5);
%                         0 0 1 allHTMs(3, 8*i-5);
%                         0 0 0 1];
%                 end
% 
%                 % Also calculate the transformation matrix of the uncut
%                 % section. All that's required here is a translation along
%                 % the previous frame's z axis.
%                 Tuncut = [1 0 0 0;%Tcut(1,3) * obj.distToNotches(i);
%                     0 1 0 0;%Tcut(2,3) * obj.distToNotches(i);
%                     0 0 1 obj.distToNotches(i);
%                     0 0 0 1];
% 
%                 % Append the matrices for the cut and uncut sections to the
%                 % matrix of all matrices.
%                 allHTMs = [allHTMs, Tcut, Tuncut];
%             end
% 
%             % Set the return value to the matrix of all transformation
%             % matrices. Of the form
%             % [T0uncut, T1cut, T1uncut, ..., TNcut, TNuncut]
%             HTMs = allHTMs;
% 
%         end

    end
end

