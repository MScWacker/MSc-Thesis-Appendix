function [Q_loss_tot, T_x_tot] = frame_dragV5(R, L_s, beta, n_f, n_parts, d_f, theta_0, omega, v_ref, h_ref, psi)
%V5 with frame sides partitioned into multiple parts 

% R is vector of segment radi, top to bottom
% L_s is vector of segment length, top to bottom
% n_f is the number of corners
% d_f is the frame diameter
% theta is position on ring
% theta_0 is start orientation of ring
% delta_theta is relative position of tethers
% x cylindrical coordinates along tether axis
% omega is rotational speed (assumed to be equal throughout the TRPT)
% n_parts number of pieces each frame segment is cut into

%Calculate the TRPT frame drag for a given timestep
Cdf = 1.2;      % frame drag coefficient
Cff = 0.04;     % frame skin friction drag coefficient
Clf = 0;        % frame lift coefficient
rho = 1.225;    % air density

%% Position of each node in global coordinates
delta_theta = (0:(2*pi)/n_f:2*pi-((2*pi)/n_f))';
radius = zeros(n_f, length(R));
theta = zeros(n_f, length(R));
x = zeros(n_f, length(R));
cb = cos(beta);
sb = sin(beta);

% position of each corner node in cylindrical coordinates
for i = 1:length(R)
    radius(:, i) = R(i);
    theta(:, i) = delta_theta + theta_0(i);
    if i == length(R) % if the bottom ring
        x(:, i) = 0;
    else
        x(:, i) = sum(L_s(i:end));
    end
end
       
% transform corner points to cartisian coordinates in global ref-frame
Rx = x.*cb - radius.*cos(theta).*sb;
Ry = radius.*sin(theta);
Rz = x.*sb + radius.*cos(theta).*cb;     % +R(end) is added at wind speed calculations for raised ground station

%% Frame vectors from segment to segment
x_f = zeros(n_f, length(R));
y_f = zeros(n_f, length(R));
z_f = zeros(n_f, length(R));

% Find all frame vectors within one segment
for i = 1:n_f-1
    x_f(i, :) = Rx(i,:) - Rx(i+1,:);
    y_f(i, :) = Ry(i,:) - Ry(i+1,:);
    z_f(i, :) = Rz(i,:) - Rz(i+1,:);
end
x_f(end, :) = Rx(end,:) - Rx(1,:);
y_f(end, :) = Ry(end,:) - Ry(1,:);
z_f(end, :) = Rz(end,:) - Rz(1,:);

% frame lengths
l_f = sqrt(x_f.^2 + y_f.^2 + z_f.^2);

% unit vector along each frame length
ex_f = x_f./l_f;
ey_f = y_f./l_f;
ez_f = z_f./l_f;

%% Drag on each segment  
Q_loss_tot = zeros(1, length(R));
T_x_tot = zeros(1, length(R));
  
for i = 1:length(R) % for each segment
    % partition frame in n_parts parts
    l_fp =  l_f(1, i)/n_parts;      
    n_p = (0.5:n_parts-0.5)';  
    
    for j = 1:n_f   % for each frame side
        % find central point of each frame part
        Fx = Rx(j, i) - n_p.*l_fp.*ex_f(j, i);
        Fy = Ry(j, i) - n_p.*l_fp.*ey_f(j, i);
        Fz = Rz(j, i) - n_p.*l_fp.*ez_f(j, i);
        % find the alignment of the fram side
        ef = ones(n_parts,3).*[ex_f(j, i) ey_f(j, i) ez_f(j, i)];

        %Transform from cart. wind to cart. ring reference frame
        %Rfx = Fx.*cb + Fz.*sb;
        Rfy = Fy;
        Rfz = -Fx.*sb + Fz.*cb;
        
%         if (~isreal(Rfy) || ~isreal(Rfy))
%             disp('Frame drag has complex number')
%             Rfy = imag(Rfy);
%             Rfz = real(Rfz);
%         end

        %find azimuth position (zeta) and radius of the selected points on the
        %frame
        zeta = atan2(Rfy, Rfz);
        Rf = sqrt(Rfy.^2 + Rfz.^2);
        cz = cos(zeta);
        sz = sin(zeta);
    
        % find the wind speed at height of each frame part point
        Vw_f = v_ref.*((Fz + R(end))./h_ref).^psi; % Ground station height is added to get correct wind speed 

        % apparent wind in wind reference frame
        Vax = -omega.*Rf.*sb.*sz + Vw_f;
        Vay = -omega.*Rf.*cz;
        Vaz = omega.*Rf.*cb.*sz;

        VaN = [Vax Vay Vaz];
        eVa = VaN./vecnorm(VaN, 2, 2);
        alpha = acos(dot(ef, eVa, 2));  % angle of attack
        Va_mag = vecnorm(VaN, 2, 2);    % magnitude of appearant wind speed

        % Find the magnitude of the aerodynamic force components
        Fdt_mag = 0.5*rho*Cdf*d_f.*l_fp.*Va_mag.^2.*sin(alpha).^2;
        Fdr_mag = 0.5*rho*Cff*pi*d_f.*l_fp.*Va_mag.^2.*cos(alpha).^2;
        Flt_mag = 0.5*rho*Clf*d_f.*l_fp.*Va_mag.^2.*sin(alpha).^2;

        % Find the direction of the aerodynamic force components
        etest = cross(eVa, ef, 2);
        eFlt = etest./vecnorm(etest, 2, 2);
        eFlt(isnan(eFlt)) = 0;
        
        etest = cross(ef, eFlt, 2);
        eFdt = etest./vecnorm(etest, 2, 2);
        eFdt(isnan(eFdt)) = 0;

        % Find the force vectors and total aerodynamic drag force
        Fdt = Fdt_mag.*eFdt;
        Fdr = Fdr_mag.*ef;
        Flt = Flt_mag.*eFlt;
        Fdrag = Fdt + Fdr + Flt;
  
        % Transform force into the frame reference frame
        FdragX = Fdrag(:, 1).*cb + Fdrag(:, 3).*sb;    % in Axis of rotation x
        FdragY = Fdrag(:, 1).*sb.*sz + Fdrag(:, 2).*cz - Fdrag(:, 3).*cb.*sz;   % in Tangential axis
        % FdragZ = -Fdrag(:, 1).*sb.*cz + Fdrag(:, 2).*sz + Fdrag(:, 3).*cb.*cz;
    
        %Find the tension contribution
        T_x_tot(i) = T_x_tot(i)+sum(FdragX); % total axial tension contribution on rotational axis
        % T_z_tot(i) = T_z_tot(i)+sum(FdragZ); % total radial tension within disc    

        %Find the torque loss
        Q_loss_tot(i) = Q_loss_tot(i)+sum(FdragY.*Rf); % drag loss from each frame part
    end
end
if any(isnan(Q_loss_tot))
    disp('Some calculated forces are labeled as nan');
end
if n_f<3  % When frame only has two corners, then there are less frames
    Q_loss_tot = Q_loss_tot/2;
end
Q_loss_tot = real(Q_loss_tot);
T_x_tot = real(T_x_tot);

end