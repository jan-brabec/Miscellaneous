r = 1e-5; %radius of circular restriction
th = 0:pi/50:2*pi;
xunit = r * cos(th);
yunit = r * sin(th);

D0 = 2e-9; %bulk diffusion coefficient

p = 1e5; %number of particles
x = [zeros(p,1) zeros(p,1)]; %restricted diffusion
y = [zeros(p,1)-10e-5 zeros(p,1)]; %free diffusion

T = 20*1e-3; % max time
t = linspace(0,T,101); %time
dt = t(2) - t(1);
n = round(T / dt);

dx = sqrt(2*D0*dt); %stepping

for c = 1:n
    
    %Simulation step
    x = x + randn(size(x,1),2) * dx;
    y = y + randn(size(x,1),2) * dx;
    
    if c == 1
        x_init = x(:,1); %record along x
        y_init = y(:,1); %record along x
    end
    
    %Apply restriction
    if (1)
        ind = find(sqrt(x(:,1).^2 + x(:,2).^2) > r);
        
        for i = 1:numel(ind)
            theta = atan2(x(ind(i),2),x(ind(i),1));
            ro = sqrt(x(ind(i),1).^2+x(ind(i),2).^2) - r; %overshoot
            tmpx = (r-ro) * cos(theta); %place on the opposite side
            tmpy = (r-ro) * sin(theta);
            
            x(ind(i),1) = tmpx;
            x(ind(i),2) = tmpy;
        end
        
        ind = find(sqrt(x(:,1).^2+x(:,2).^2) > r); %those that
        x(ind,1) = 0;
        x(ind,2) = 0;
    end
    
    %Compute
    x_act = x(:,1); %record along x, restricted diffusion
    y_act = y(:,1); %record along x, free diffusion
    
    msd_x(c) = mean((x_act-x_init).^2); %mean square displacement
    ADC_x(c) = 1/2 * msd_x(c)/t(c); %apparent diffusion coeffiecient
    
    msd_y(c) = mean((y_act-y_init).^2);
    ADC_y(c) = 1/2 * msd_y(c)/t(c); %apparent diffusion coeffiecient
    
    %Plot
    
    % h = figure;
    % set(h,'WindowStyle','docked')
    
    subplot(1,3,1)
    plot(x(:,1),x(:,2),'.','Color',[156 97 20]./255,'Markersize', 0.01)
    hold on
    plot(xunit, yunit,'Color',[179 200 186]./255,'Linewidth',4);
    
    hold on
    plot(y(:,1),y(:,2),'.','Color',[4 0 108]./255,'Markersize', 0.01)
    xlim([-15e-5 2e-5])
    ylim([-10e-5 7e-5])
    axis equal
    axis off
    drawnow;
    hold off
    
    subplot(1,3,2)
    plot(t(1:c)*1e3,msd_x(1:c)*1e12,'Color',[156 97 20]./255,'Linewidth',5)
    hold on
    %     plot(t(1:c)*1e3,2*D0.*t(1:c)*1e12,'Color',[156 97 20]./255) %unit test
    plot(t(1:c)*1e3,msd_y(1:c)*1e12,'Color',[4 0 108]./255,'Linewidth',5)
    hold on
    ylim([0, r^2/2]*1e12)
    xlim([0 max(t)*1e3])
    %     xlabel('time [ms]')
    %     ylabel(['\langle\Delta x^2 \rangle [\mum^2]'])
    set(gca,'YTicklabels',[])
    set(gca,'XTicklabels',[])
    drawnow;
    hold off
    
    subplot(1,3,3)
    plot(t(1:c)*1e3,ADC_x(1:c)*1e9,'Color',[156 97 20]./255,'Linewidth',5)
    hold on
    plot(t(1:c)*1e3,smooth(ADC_y(1:c),3)*1e9,'Color',[4 0 108]./255,'Linewidth',5)
    ylim([0 1.1*D0*1e9])
    xlim([0 max(t)*1e3])
    %     ylabel('ADC [\mum^2/ms]')
    %     xlabel('time [ms]')
    set(gca,'YTicklabels',[])
    set(gca,'XTicklabels',[])
    drawnow;
    hold off;
    
    
    %   saveas(h,sprintf('FIG%d.png',c));
    
    % GIF files created by https://gifmaker.me/
    
end



