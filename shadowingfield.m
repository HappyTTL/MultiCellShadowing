%function [M,D_R,D_theta] = shadowingfield(n, r_min, r_max, theta,R0)
for count = 3:9
    n = 50;
    r_min = 20;
    r_max = 1500;
    r = 3000;
    theta = pi/6;
    R0 = 6;
    F_theta = n;
    D_theta = ceil(F_theta*2*pi/theta);
    F_R = 3*n;
    D_R = ceil(F_R*10*log10(r_max/r_min)/R0);
    
    %%%%%%%%%%%%%%%%%%%%%%%Fast Shadowing Fields Generation
    
    Z = randn(D_theta, D_R+F_R-1);
    W = zeros(D_theta, D_R+F_R-1);
    M = zeros(D_theta, D_R);
    for j =1 : D_R+F_R-1
        W(1,j) = sum(Z(1:F_theta,j));
    end
    for i = 1 : D_theta-1
        i_star = mod(i+F_theta-1, D_theta)+1;
        for j =1 : D_R+F_R-1
            W(i+1,j) = W(i,j)-Z(i,j)+Z(i_star,j);
        end
    end
    for i = 1 : D_theta
        M(i,1) = sum(W(i,1:F_R));
    end
    for j = 1 : D_R-1
        for i = 1 : D_theta
            M(i,j+1) = M(i,j)-W(i,j)+W(i,j+F_R);
        end
    end
    M = M./sqrt(F_theta*F_R);
    
    % A = zeros(r,r);
    r_db_min = 0;
    r_db_max = 10*log10(r_max/r_min);
    % for theta = 0:2*pi/D_theta:(2*pi-0.1)
    %     for r = r_db_min:(r_db_max-r_db_min)/D_R:(r_db_max-0.1)
    %         rho = theta;
    %         R = 10.^(0.1.*r)*r_min;
    %         [X,Y] = pol2cart(rho,R);
    %         A(round(X)+r/2,round(Y)+r/2) = M(round(theta*D_theta/2/pi)+1,round((r-r_db_min)*D_R/(r_db_max-r_db_min))+1);
    %     end
    % end
    B = zeros(r,r);
    a = 0;
    b = 0;
    for i = 1:r
        for j = 1:r
            [rho,R] = cart2pol(i-r/2,j-r/2);
            for m = 1 : D_theta
                for n = 1 : D_R
                    if rho > -pi+(m-1)*2*pi/D_theta && rho <= -pi+m*2*pi/D_theta && 10*log10(R/r_min)>r_db_min+(n-1)*(r_db_max-r_db_min)/D_R && 10*log10(R/r_min)<r_db_min+n*(r_db_max-r_db_min)/D_R
                        B(i,j)=M(m,n);
                    end
                end
            end
        end
    end
    filename = ['ShadowField_30_6_' num2str(count) '.mat'];
    save(filename, 'B');
end
% figure(1);
% imagesc(B);
% colormap(jet);
%end