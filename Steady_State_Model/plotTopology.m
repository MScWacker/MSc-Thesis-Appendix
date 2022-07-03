function plotTopology(Node_pos, R, L_s, l_b, n_r, beta, theta, txt)
global n_R n_L n_t n_f

f_con = [(1:n_f) 1];
t_con = (0:n_L)*n_t;

figure('Name', 'Daisy sketch', 'Color', 'w')  
    plot3(Node_pos(:,1),Node_pos(:,2),Node_pos(:,3),'o','color','k','MarkerSize',3)
    hold on;
    for i = 0:n_R-1
        plot3(Node_pos(i*n_f+f_con,1),Node_pos(i*n_f+f_con,2),Node_pos(i*n_f+f_con,3),'color','r') 
    end   
    for i = 1:n_t
        plot3(Node_pos(i+t_con,1),Node_pos(i+t_con,2),Node_pos(i+t_con,3),'color','b')
    end   
    
    n_b = 3;
    l_c = l_b/5;  % Chord length of blade [m]
    dtheta = 0:2*pi/n_b:2*pi-2*pi/n_b; 
    for j = 1:n_r
        % in radia coordinates
        x_r = sum(L_s(j:end))*ones(4,1);
        y_r = [l_c; -l_c; -l_c; l_c]; 
        z_r = [R(j)+l_b/2; R(j)+l_b/2; R(j)-l_b/2; R(j)-l_b/2];

        for i = 1:n_b
            xyz_w = [x_r y_r z_r]*transMat(beta, dtheta(i)) + [zeros(4,1), zeros(4,1), R(end)*ones(4,1)];
            fill3([xyz_w(:,1);xyz_w(1,1)]', [xyz_w(:,2);xyz_w(1,2)]', [xyz_w(:,3);xyz_w(1,3)]','b')  
        end
    end
    % plot lines from rotor towards lifter
%     conPoint = [sum(L_s)+L_s(1)/2 0 0; sum(L_s)+L_s(1) 0 0]*transMat(beta, 0) + [0, 0, R(end)];
%     for i = 1:n_t
%         plot3([Node_pos(i,1);conPoint(:,1)],[Node_pos(i,2);conPoint(:,2)],[Node_pos(i,3);conPoint(:,3)],'color','b')
%     end 
    title(txt, 'FontWeight','normal');
    axis equal; 
%     annotation('textbox', [0.25 0.4 0.3 0.3], 'String',txt,'FitBoxToText','on');
    hold off;
%     axis off; 
    set(gca,'fontsize',14)
    
end
% 
% %Tension losses
% figure
%     scatter(1:n_r,T_rot,'filled')
%     hold on;
%     plot(1:length(R),T_fDrag,'-o')
%     plot(1:length(R),T_tDrag,'-o')
%     plot(1:length(R),T,'-o')
%     legend('Rotor Induced Tension', 'Frame Induced Tension', 'Tether Induced Tension','Net Axial Tension')
%     xlabel('TRPT Polygon index []')
%     ylabel('Axial Tension [N]')
%     title('Axial Tension through TRPT')
%     grid on;
%     hold off;

% %Torque losses
% figure
%     scatter(1:n_r,Q_rot,'filled')
%     hold on;
%     plot(1:length(R),Q_fDrag,'-o')
%     plot(1:length(R),Q_tDrag,'-o')
%     plot(1:length(R),Q,'-o')
%     scatter(length(R),Q(end)*eta_gen,'filled')
%     lgd = legend('Rotor Induced Torque', 'Frame Drag Loss', 'Tether Drag Loss','Net Transmitted Torque','Generator output');
%     lgd.Location = 'east';
%     xlabel('TRPT Polygon index []')
%     ylabel('Torque [Nm]')
%     title('Torque transmission through TRPT')
%     grid on;
%     hold off;
%     
% % System Topology
% figure
%     plot(cos(beta)*[0 cumsum(flip(L_s))],sin(beta)*[0 cumsum(flip(L_s))]+R(end),'--x','color','k')
%     hold on;
%     plot(cos(beta)*[0 cumsum(flip(L_s))]-sin(beta)*flip(R), sin(beta)*[0 cumsum(flip(L_s))]+cos(beta)*flip(R)+R(end),'-x','color','k')
%     plot(cos(beta)*[0 cumsum(flip(L_s))]+sin(beta)*flip(R), sin(beta)*[0 cumsum(flip(L_s))]-cos(beta)*flip(R)+R(end),'-x','color','k')
%     for i = 1:n_r
%         plot([cos(beta)*sum(L_s(i:end))+sin(beta)*(R(i)+l_b/2), cos(beta)*sum(L_s(i:end))+sin(beta)*(R(i)-l_b/2)], ...
%             [sin(beta)*sum(L_s(i:end))-cos(beta)*(R(i)+l_b/2)+R(end), sin(beta)*sum(L_s(i:end))-cos(beta)*(R(i)-l_b/2)+R(end)],'color','b') 
%         plot([cos(beta)*sum(L_s(i:end))-sin(beta)*(R(i)+l_b/2), cos(beta)*sum(L_s(i:end))-sin(beta)*(R(i)-l_b/2)], ...
%             [sin(beta)*sum(L_s(i:end))+cos(beta)*(R(i)+l_b/2)+R(end), sin(beta)*sum(L_s(i:end))+cos(beta)*(R(i)-l_b/2)+R(end)],'color','b')     
%     end
%     ylabel('Height [m]')
%     ylabel('Length [m]')
%     title(['System Topology with ',num2str(n_t),' Tethers and ',num2str(n_f),'-sided Polygons'])
%     axis equal