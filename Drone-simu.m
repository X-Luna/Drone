clear;
tsm = 0.0001;
ts = 0.01;
t = 0;
T = 800;


vp = 10;
vp_theta = 0;
x0 = 0;
x = x0;
y0 = 0;
y = y0;


d0 = 50;
vt0 = 1;
vt = vt0;
vt_theta = -45;
xt0 = -300;
xt = xt0;
yt0 = 600;
yt = yt0;

w = 60;
z_past = d0;



w_max = 30;
wd_max = 90;



f = 3;
w0 = 0;
w1 = 0;


w_input = 0;

 
i = 0;
ii = 0;
t_record = [];
x_record = [];
xt_record = [];
y_record = [];
yt_record = [];

z_record = [];

while(t < T)
    
    t = t + tsm; 
    i = i + 1;
    
    vt = vt0 + 0.1*sin(t);
    vt_theta = -45 + 0.2*t;
    

    xt = xt + vt*tsm*cosd(vt_theta);     
    yt = yt + vt*tsm*sind(vt_theta); 
    x = x + vp*tsm*cosd(vp_theta);
    y = y + vp*tsm*sind(vp_theta);
    vp_theta = vp_theta + w*tsm;
    vp_theta = mod(vp_theta+180,360)-180;
    
    z = sqrt((x - xt)^2+(y - yt)^2);
    
    w0 = (w1+2*f*pi*tsm*w_input)/(1+2*f*pi*tsm); 
    if(-wd_max*tsm <= w0-w && w0-w<= wd_max*tsm)
        w = w0;
    else
        if(w0-w > wd_max*tsm)
            w = w + wd_max*tsm;
        else
            w = w - wd_max*tsm;
        end
    end
    w1 = w;
    
    if(mod(t,ts)<tsm)
        
        ii = ii + 1;
        
            k = 0.01;
            z_dot = (z-z_past)/ts;
            vi = sqrt(vp^2-vt^2*sin(vt_theta-vp_theta)^2);
            if (vp_theta>180)
                if (z>=d0)
                %w_input = rad2deg(k*z_dot + vp/d0);
                    w_input = rad2deg(vp*(k*z_dot + vp/d0)/vi);
                else
                %w_input = rad2deg(k*d0*z_dot/z + vp/d0);
                    w_input = rad2deg(vp*(k*d0*z_dot/z + vp/d0)/vi);
                end
            else    
                if (z>=d0)
                %w_input = rad2deg(k*z_dot + vp/d0);
                    w_input = rad2deg(vp*(-k*z_dot - vp/d0)/vi);
                else
                %w_input = rad2deg(k*d0*z_dot/z + vp/d0);
                    w_input = rad2deg(vp*(-k*d0*z_dot/z - vp/d0)/vi);
                end
            end
            
            z_past = z;
        
        t_record(ii) = t;
        x_record(ii) = x;
        y_record(ii) = y;
        xt_record(ii) = xt;
        yt_record(ii) = yt;
        
        z_record(ii) = z;
    end
    
end

t0 = 0.01;



plot(xt_record(t0*100:end),yt_record(t0*100:end),'bo','linewidth',2);
hold on
plot(x_record(t0*100:end),y_record(t0*100:end),'r','linewidth',2);



%plot(t_record(t0*100:end),z_record(t0*100:end),'b','linewidth',2);

%xlabel('x/m','fontsize',22,'fontname','宋体');
%ylabel('y/m','fontsize',22,'fontname','宋体');
xlabel('t/s','fontsize',22,'fontname','宋体');
ylabel('ρ/m','fontsize',22,'fontname','宋体');
%set(gca,'FontSize',20);
%legend('\fontsize{20}目标轨迹','\fontsize{20}无人机轨迹');


%axis equal
grid on