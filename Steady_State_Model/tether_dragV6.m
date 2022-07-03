function [Q_loss_tot, T_x_tot, Node_pos, theta, d_t] = tether_dragV6(R, L_s, beta, n_t, n_parts, theta_0, omega, v_ref, h_ref, psi, T)
% V4 where tethers are partitioned into n_parts pieces
% V5 is returning the node positions
% V6 calculated tether diameter from tensile forces

% R is vector of segment radi, top to bottom
% L_s is vector of segment length, top to bottom
% n_t is number of tethers per segment
% theta is position on ring
% theta_0 is start orientation of ring
% delta_theta is relative position of tethers
% x cylindrical coordinates along tether axis
% omega is rotational speed (assumed to be equal throughout the TRPT)


%Calculate the tether drag for a given timestep
Cdt = 1.2;      % tether drag coefficient
Cft = 0.04;     % tether skin friction drag coefficient
Clt = 0;        % tether lift coefficient
rho = 1.225;    % air density
sigma = 249.8e6;  % Tensile strength if kite surf line https://www.liros.com/userdata/files/LIROS_rope_catalogue_2022_2023_eng_210_x_297_mm_low_res.pdf
sf = 5;         % safety factor;
%% Position of each node in global coordinates
delta_theta = (0:(2*pi)/n_t:2*pi-((2*pi)/n_t))';
radius = zeros(n_t, length(R));
theta = zeros(n_t, length(R));
x = zeros(n_t, length(R));
cb = cos(beta);
sb = sin(beta);

% position of each node in cylindrical coordinates
for i = 1:length(R)
    if i == length(R) % if the bottom ring
        radius(:, i) = R(i);
        theta(:, i) = delta_theta + theta_0(i);
        x(:, i) = 0;
    else
        radius(:, i) = R(i);
        theta(:, i) = delta_theta + theta_0(i);
        x(:, i) = sum(L_s(i:end));
    end
end
        
% transform to cartisian coordinates
Rx = x.*cb - radius.*cos(theta).*sb;
Ry = radius.*sin(theta);
Rz = x.*sb + radius.*cos(theta).*cb; % +R(end) for raised ground station is added later

Node_pos = zeros(length(R)*n_t,3);
% Node_con = zeros(length(R)*n_t);
% count = 0;
for i = 1:length(R)
    for j = 1:n_t
        Node_pos((i-1)*n_t+j,:) = [Rx(j,i) Ry(j,i) Rz(j,i)+R(end)];
    end
end

%% Tether vectors from segmen to segment
x_t = zeros(n_t, length(L_s));
y_t = zeros(n_t, length(L_s));
z_t = zeros(n_t, length(L_s));
for i = 1:length(L_s)
    x_t(:, i) = Rx(:, i) - Rx(:, i+1);
    y_t(:, i) = Ry(:, i) - Ry(:, i+1);
    z_t(:, i) = Rz(:, i) - Rz(:, i+1);
end

% tether lengths
l_t = sqrt(x_t.^2 + y_t.^2 + z_t.^2);

% find tether thickness
    % angle between tether and rotation axis is dot product of their
    % normalised vectors which can be simplified to cos(gamma) = l_s/l_t
Ft = T(1:end-1)./n_t.*l_t(1,:)./L_s';
D_t = 2* sqrt(Ft*sf/(sigma*pi));
d_t = max(D_t);

% unit vector along each tether
ex_t = x_t./l_t;
ey_t = y_t./l_t;
ez_t = z_t./l_t;

%% Drag on each tether  
Q_loss_tot = zeros(1, length(R));
T_x_tot = zeros(1, length(R));

n_p = (0.5:n_parts-0.5)';   % Each tether half is split into n_parts pieces

for i = 1:length(R) % for each section
    for j = 1:n_t   % for each tether
        if i == 1 %if the top ring
            l_tp =  l_t(j,i)/2/n_parts;
            % half the section below the ring
            Tx = Rx(j, i) - n_p.*l_tp.*ex_t(j, i);
            Ty = Ry(j, i) - n_p.*l_tp.*ey_t(j, i);
            Tz = Rz(j, i) - n_p.*l_tp.*ez_t(j, i);
            et = ones(n_parts,3).*[ex_t(j, i) ey_t(j, i) ez_t(j, i)];

        elseif i == length(R) % if the last ring
            l_tp =  l_t(j,i-1)/2/n_parts;
            % half the section above the ring
            Tx = Rx(j, i) + n_p.*l_tp.*ex_t(j, i-1);
            Ty = Ry(j, i) + n_p.*l_tp.*ey_t(j, i-1);
            Tz = Rz(j, i) + n_p.*l_tp.*ez_t(j, i-1);
            et = ones(n_parts,3).*[ex_t(j, i-1) ey_t(j, i-1) ez_t(j, i-1)];

        else
            %half the section above ring
            l_tp =  l_t(j,i-1)/2/n_parts;
            Tx = Rx(j, i) + n_p.*l_tp.*ex_t(j, i-1);
            Ty = Ry(j, i) + n_p.*l_tp.*ey_t(j, i-1);
            Tz = Rz(j, i) + n_p.*l_tp.*ez_t(j, i-1);
            et = ones(n_parts,3).*[ex_t(j, i-1) ey_t(j, i-1) ez_t(j, i-1)];

            % half the section below the ring
            l_tp =  l_t(j,i)/2/n_parts;
            Tx = [Tx; Rx(j, i) - n_p.*l_tp.*ex_t(j, i)];
            Ty = [Ty; Ry(j, i) - n_p.*l_tp.*ey_t(j, i)];
            Tz = [Tz; Rz(j, i) - n_p.*l_tp.*ez_t(j, i)];
            et = [et; ones(n_parts,3).*[ex_t(j, i) ey_t(j, i) ez_t(j, i)]];
        end
    
 
        % Transform from wind to ring reference frame
%         Rtx = Tx.*cb + Tz.*sb;
        Rty = Ty;
        Rtz = -Tx.*sb + Tz.*cb;

        if (~isreal(Rty) || ~isreal(Rty))
            disp('Tether drag has complex number')
            Rty = imag(Rty);
            Rtz = real(Rtz);
        end
            
        % find azimuth position (zeta) and radius of the selected points
        zeta = atan2(Rty, Rtz);
        Rt = sqrt(Rty.^2 + Rtz.^2);
        cz = cos(zeta);
        sz = sin(zeta);

        %% Improved tether drag
            
        % find the wind speed at height of each tether point
        Vw_t = v_ref.*((Tz + R(end))./h_ref).^psi;

        % apparent wind in wind reference frame
        Vax = -omega.*Rt.*sb.*sz + Vw_t;
        Vay = -omega.*Rt.*cz;
        Vaz = omega.*Rt.*cb.*sz;

        VaN = [Vax Vay Vaz];
        eVa = VaN./vecnorm(VaN, 2, 2);
        alpha = acos(dot(et, eVa, 2));  % angle of attack
        Va_mag = vecnorm(VaN, 2, 2);    % magnitude of appearant wind speed

        % Find the magnitude of the aerodynamic force components
        Fdt_mag = 0.5*rho*Cdt*d_t.*l_tp.*Va_mag.^2.*sin(alpha).^2;
        Fdr_mag = 0.5*rho*Cft*pi*d_t.*l_tp.*Va_mag.^2.*cos(alpha).^2;
        Flt_mag = 0.5*rho*Clt*d_t.*l_tp.*Va_mag.^2.*sin(alpha).^2;

        % Find the direction of the aerodynamic force components
        etest = cross(eVa, et, 2);
        eFlt = etest./vecnorm(etest, 2, 2);
        eFlt(isnan(eFlt)) = 0;

        etest = cross(et, eFlt, 2);
        eFdt = etest./vecnorm(etest, 2, 2);
        eFdt(isnan(eFdt)) = 0;

        % Find the force vectors and total aerodynamic drag force
        Fdt = diag(Fdt_mag)*eFdt;
        Fdr = diag(Fdr_mag)*et;
        Flt = diag(Flt_mag)*eFlt;
        Fdrag = Fdt + Fdr + Flt;

        % Transform force into the tether reference frame
        FdragX = Fdrag(:, 1).*cb + Fdrag(:, 3).*sb;
        FdragY = Fdrag(:, 1).*sb.*sz + Fdrag(:, 2).*cz - Fdrag(:, 3).*cb.*sz;   % Tangential axis
        % FdragZ = -Fdrag(:, 1).*sb.*cz + Fdrag(:, 2).*sz + Fdrag(:, 3).*cb.*cz;

        %Find the tension contribution
        T_x_tot(i) = T_x_tot(i)+sum(FdragX); % total axial tension contribution on rotational axis
        % T_z_tot(i) = T_z_tot(i)+sum(FdragZ); % total radial tension within disc    

        %Find the torque loss
        Q_loss = FdragY.*Rt; % loss from each tether
        Q_loss_tot(i) = Q_loss_tot(i)+sum(Q_loss); % total loss from all tethers around the segment
    end
end
Q_loss_tot = real(Q_loss_tot);
T_x_tot = real(T_x_tot);
d_t = real(d_t);
end