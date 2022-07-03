function plotOutput(WEI, WEI_g, WEI_f, WEI_t, WEI_b, WEI_l, POW, P_rot, P_pot, eta_trpt, eta_gen)

C_p = P_rot / P_pot; % Coefficient of performance of rotor with tilted swept area
power = 100*[1-C_p, C_p*(1-eta_trpt), C_p*eta_trpt*(1-eta_gen), C_p*eta_trpt*eta_gen];
weight = 100*[WEI_l, WEI_b, WEI_t, WEI_f, WEI_g]/WEI;


disp(['Power output is ', num2str(POW), ' W']) ;
% disp(['Power output should be ', num2str(P_pot * C_p * eta_trpt * eta_gen), ' W']);
disp(['Material weight is ', num2str(WEI), ' kg']);
disp(['Power density is ', num2str(POW/WEI), ' W/kg']);

% power = P_pot*[C_p, C_p*eta_trpt, C_p*eta_trpt*eta_gen*eta_trpt, 1-C_p*eta_trpt*eta_gen*eta_trpt];


figure('Name', 'Daisy sketch', 'Color', 'w')
tiledlayout(2,1)
    nexttile;
    %bar( categorical({['Total of ', num2str(P_pot/1000) ,' kW']}),power,'stacked')
    barh(categorical({'Power Breakdown'}), power,'stacked')
    legend("P_{loss, rotor}","P_{loss, TRPT}","P_{loss, gen}","P_{out}","Location","eastoutside")
    xlabel('P [%]');
    xlim([0,100]);
    grid minor;
    set(gca,'fontsize',12)
    nexttile;
    %bar( categorical({['Total of ', num2str(WEI) ,' kg']}), weight,'stacked',"Location","south")
    barh(categorical({'Mass Breakdown'}), weight,'stacked')
    legend("m_{lifter}","m_{rotor}","m_{tether}","m_{frame}","m_{ground}","Location","eastoutside")
    xlabel('m [%]');
    xlim([0,100]);
    grid minor;
    set(gca,'fontsize',12)

end

