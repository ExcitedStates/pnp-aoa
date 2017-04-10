function pnp_num(sig, kk, cc, pore_size)

%% quick setting

switch nargin
    case 0,
        kk = 0;
        sig = 1;
        cc = 10;
	pore_size = 1.3;
    case 1,
        kk = 0;
        cc = 10;
	pore_size = 1.3;
    case 2,
        cc = 10;
	pore_size = 1.3;
    case 3,
	pore_size = 1.3;
end

disp(['(' num2str(sig) 'x, ' num2str(kk) ' rate, ' num2str(cc) 'nM, '...
    num2str(pore_size) ' nm )']);
        

%% sim setting

close all;

real_time_movie = 1;
reaction = 0;

n_display = 1000;

%% constants

k_b = 1.38064852e-23; % Boltzmann constant in  m2 kg s-2 K-1
e_c = 1.602e-19; % elementary charge in C
avo = 6.022e23;
r = 8.31446;
faraday = avo * e_c;


%% parameters

% simulation
dt = 1e-12; % s
dx = 0.1e-9; % m

l_t = 5e6 * dt;
l_x = 22.1e-9;
l_y = 60.1e-9;

% experiment

kr = kk; % 12.8 day^-1
d_m = 1.64e-9; % m^2/s
eps =  80 * 8.854e-12; % solvent permittivity in C/V/m
eps_o =  80 * 8.854e-12;
tmp = 300;
z = [1 -1]; % valence charge
k_m = 130e-9 * 1000; 

%% membrane


x_s = 1;

y_pc = (l_y - 20e-9/x_s); % center of the protein is 20 nm above the bottom

% l_px = 1.3e-9;
l_px = pore_size *1e-9;
l_py = 4.5e-9;

l_px = l_px / x_s;
l_py = l_py / x_s;

y_pt = round( (y_pc-l_py/2)/dx );
y_pb = round( (y_pc+l_py/2)/dx );

x_pl = round( (l_x-l_px)/2/dx );
x_pr = round( ((l_x-l_px)/2+l_px)/dx );

% domain 1
% C( y_pt:y_pb, 1:end, : ) = 0;
% domain 2
% C( y_pt:y_pb, 1:end, : ) = 0;
%% init

n_i = length(z);

n_t = ceil(l_t/dt);
n_x = ceil(l_x/dx);
n_y = ceil(l_y/dx);

x_c = round( (n_x+1)/2 );
y_c = round( (n_y+1)/2 );

% analytical

[C_pls, C_min, Ey0, V] = pnp_ana(cc*1e-9, sig, dx, n_y);

% conecentration
C = zeros( n_y, n_x, n_i);
if sig > 0
    C(1:(y_pt-1),:,1) = 1000*repmat( fliplr(C_pls(2:(y_pt-0))).', [1 n_x]);
    C(1:(y_pt-1),:,2) = 1000*repmat( fliplr(C_min(2:(y_pt-0))).', [1 n_x]);
    C(y_pt:y_pb,(x_pl+1):(x_pr-1),1) = 1000*repmat( C_pls(1), [length(y_pt:y_pb) length((x_pl+1):(x_pr-1))]);
    C(y_pt:y_pb,(x_pl+1):(x_pr-1),2) = 1000*repmat( C_min(1), [length(y_pt:y_pb) length((x_pl+1):(x_pr-1))]);
    C((y_pb+1):end,:,1) = 1000*repmat( (C_pls(2:(size(C((y_pb+1):end,:,1),1)+1))).', [1 n_x]);
    C((y_pb+1):end,:,2) = 1000*repmat( (C_min(2:(size(C((y_pb+1):end,:,2),1)+1))).', [1 n_x]);
else
    C(:) = 1e-6 * cc;
    % domain 1
    C( y_pt:y_pb, 1:x_pl, : ) = 0;
    % domain 2
    C( y_pt:y_pb, x_pr:end, : ) = 0;
    
end


% C(:,:,2) = 1000*repmat( fliplr(C_min(2:(end-1))).', [1 n_x]);

% C(C~=0) = 10e-9 * 1000;

c_bulk_pls = C(1,:,1);
c_bulk_min = C(1,:,2);

% C(2:end,:,1) = c_bulk_pls;
% C(2:end,:,2) = c_bulk_min;

% fluxes
Jx = zeros( n_y, n_x+1, n_i);
Jy = zeros( n_y+1, n_x, n_i);
% E-fields
Ex = zeros( n_y, n_x+1);
Ey = zeros( n_y+1, n_x);
% Ey = -1000 * repmat(fliplr(Ey0).', [1 n_x]);
% Ey(1:(end-1),:) = 0;

% ICs

% C(1, :) = 100;
% C(131,107) = 1e0;




%% generate E-fields profile

rho = zeros(size(C(:,:,1)));

sigma = -e_c /25; % custom

rho( y_pt, 1:x_pl ) = 1;
rho( y_pb, 1:x_pl ) = 1;
rho( y_pt, x_pr:end ) = 1;
rho( y_pb, x_pr:end ) = 1;
rho( y_pt:y_pb, x_pl ) = 1;
rho( y_pt:y_pb, x_pr ) = 1;


% [X1, Y1] = meshgrid((0.5*dx):dx:(l_x+0.5*dx), dx:dx:l_y);
% [X2, Y2] = meshgrid(dx:dx:l_x, (0.5*dx):dx:(l_y+0.5*dx));


% for u = y_pt:y_pb
%     disp(['u = ' int2str(u)]);
%     for v = 1:(3*n_x)
%         
%         uu = u;
%         vv = v - n_x;
%         
%         R3 = ((vv-X1).^2 + (uu-Y1).^2 ).^(3/2);
%         Ex = Ex + 1/(4*pi*eps) * rho3(u,v) * (X1-vv) ./ R3;
%         
%         R3 = ((vv-X2).^2 + (uu-Y2).^2 ).^(3/2);
%         Ey = Ey + 1/(4*pi*eps) * rho3(u,v) * (Y2-uu) ./ R3;
%         
%         
%     end
% end



% Ex = sigma*(1/dx)^2 * Ex;
% Ey = sigma*(1/dx)^2 * Ey;

% Ex = max(abs(1000*Ey0(:)))/max(abs(Ex(:))) *  Ex;
% Ey = max(abs(1000*Ey0(:)))/max(abs(Ey(:))) *  Ey;

if sig > 0
    Ey(1:y_pt,:) = -repmat( fliplr(Ey0(1:y_pt)).', [1 n_x]);
    Ey((y_pb+1):end,:) = repmat( (Ey0(1:size(Ey((y_pb+1):end,:),1))).', [1 n_x]);
    Ex( y_pt:y_pb, (x_pl+1):x_pr )...
        = repmat(Ey0(1:length((x_pl+1):x_pr))-fliplr(Ey0(1:length((x_pl+1):x_pr))),...
        [length(y_pt:y_pb) 1]);
end




%% recation

R = zeros(size(C(:,:,1)));

%% detectors

convg_sim = zeros(1, n_t);
mass_sim = zeros(1, n_t);
chrg_sim = zeros(1, n_t);

c_dtr = zeros(1, n_t);

%% iterations

% Ex = E_fx;
% Ey = E_fy;


y2 = linspace(dx,l_y, size(C,1));
y1 = linspace(0.5*dx,l_y+0.5*dx, size(C,1)+1);
tic;

C0 = C;

if real_time_movie>0
    figure(999);
end
for t = 1:n_t
    
    C_pv = C;
    
    % update P
%     P2 = zeros(n_y, n_x);
%     for k = 1:n_i
%         P2 = P2 - faraday/eps * z(k) * C(:,:,k);
%     end
%     
%     for u = 1:max_itr
%         P_pv = P;
%         
%         P( 2:(end-1), 2:(end-1) )...
%             = 0.25 * ( P(2:(end-1), 1:(end-2)) + P(2:(end-1),3:end)...
%             + P(1:(end-2),2:(end-1)) + P(3:end,2:(end-1)) -dx^2 * P2 );
% %         P(1,:) = P(2,:);
%         P(1,:) = 0;% P(2,:);
%         P( (end-1):(end), : ) = v_z;
% 
%         P(:,1) = P(:,(end-1));
%         P(:,end) = P(:,2);
%         
%         conv_er = sum( abs( abs(P(:)) - abs(P_pv(:)) ) ) / sum( abs(P(:)) );
%         if conv_er < max_err || sum(P(:)) == 0
%             break;
%         end
%     end
%     disp(['ER = ' num2str(conv_er) ' @ step ' int2str(t)]);
    
    % update E 
    for k = 1:n_i
        Ex = Ex - (dt*faraday/eps) * z(k) * Jx(:,:,k);
        Ey = Ey - (dt*faraday/eps) * z(k) * Jy(:,:,k);
    end
    

   
    
    % PNP
   
    % update J
    
    for k = 1:n_i
        
        Jx(:,:,k) = - d_m/dx * diff([C(:,end,k) C(:,:,k) C(:,1,k)], 1, 2) ...
            + 1*(faraday/r/tmp) * d_m * z(k) * 0.5 * ( [C(:,end,k) C(:,:,k)] + [C(:,:,k) C(:,1,k)] ) .* Ex(:,:);       
        Jy(2:(end-1),:,k) = - d_m/dx * diff(C(:,:,k), 1, 1)  ...
            + 1*(faraday/r/tmp) * d_m * z(k) * 0.5 * ( C(1:(end-1),:,k) + C(2:end,:,k) ) .* Ey(2:(end-1),:);
    end
    

    % compensation by Jy
    Jy(1,:,1) = - d_m/dx * ( C(1,:,1) - c_bulk_pls)  ...
            + 1*(faraday/r/tmp) * d_m * z(k) * 0.5 * ( c_bulk_pls + C(1,:,1) ) .* Ey(1,:);
    Jy(1,:,2) = - d_m/dx * ( C(1,:,2) - c_bulk_min)  ...
            + 1*(faraday/r/tmp) * d_m * z(k) * 0.5 * ( c_bulk_min + C(1,:,2) ) .* Ey(1,:);    
    
    Jx(y_pt:y_pb, x_pl+1,:) = 0;
    Jx(y_pt:y_pb, x_pr,:) = 0;
    
    Jy(y_pt, [1:x_pl x_pr:end],:) = 0;
    Jy(y_pb+1, [1:x_pl x_pr:end],:) = 0;
    
    
%     Jy(1,:,:) = 0;
    Jy(end,:,:) = 0;
    

%     
    
    
    
    % update C
    
    C = C - dt/dx * ( diff(Jx, 1, 2) + diff(Jy, 1, 1) );
    R(434,111) =  kr * dt * (C(434,111,1)/k_m) / (6e23) / dx^3;
    C(:,:,1) = C(:,:,1) - R;
%     C(1,:,:) = c_bulk;
    
%     C(1,:,1) = C(1,:,1) - sum(R(:))/n_x;
    
    
    % monitors

    c_conv_er = sum(abs( abs(C(:)) - abs(C_pv(:)) )) / sum(abs(C(:)));
    convg_sim(t) = c_conv_er;
    
    mass_sim(t) = sum( C(:) );
    chrg_sim(t) = sum( sum( C(:,:,1) - C(:,:,2) ) );
    c_dtr(t) = C(434,111,1);
    
    if mod(t,n_display) == 0
        disp(['CONVG = ' num2str(c_conv_er) ' @ step ' int2str(t)]);
        t_now = toc;
        eta = t_now/t * (n_t-t);
        
        if c_conv_er == 0
            disp(['early convergence occured at ' int2str(t) 'th step,'...
                ' or simulation time ' num2str(1e9*t*dt) 'ns.']);
            mass_sim((t+1):end) = mass_sim(t);
            chrg_sim((t+1):end) = chrg_sim(t);
            c_dtr((t+1):end) = c_dtr(t);
	    break;
        else
            disp([num2str(round(1000*t/n_t)/10) '% done, ' num2str(round(eta)) 's to go...']);
        end
        
        if real_time_movie>0
            plot(y2, 0.001*C0(:,111,1), 'k', y2, 0.001*C(:,111,1));
%             ylim([1e-9 19e-9]);


    %         imagesc( C(:,:,1) );
    %         caxis([0 2]);
        %     colorbar;
            title(['t = ' int2str(t) '/' int2str(n_t)]);
%             ylim([v_z 0])
    %         colormap(jet);
    %         hold on;
            pause(0.01);
        %     print('-dpng', ['pot_sc050_' int0str(t,4) '.png']);
        end
    end
        
    
    
    
end

hold off;
toc






P = -cumsum( dx * Ey(1:end,:) );


save(['ana_num_' int0str(cc,3) 'nm_k' int0str(kk, 2) '_sc_sig' int0str(1000*sig, 6) '_pore' int0str(10*pore_size, 2)  '.mat'],...
    'C0','C','P','Ex','Ey','Jx','Jy',...
    'convg_sim','mass_sim','chrg_sim','c_dtr','-v7.3');
% 

end
