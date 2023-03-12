clc;clear;
format long
%% ==================== Lambert's Orbit Estimation ===================== %%
% Variables:
%
% a             = Semi-Major Axis (km)
% A             = Assumption parameter
% c             = C(z) class of Stumpff Functions
% c_z           = C(z) class of Stumpff Functions array
% del_theta     = Change in true anamoly (degree)
% delta_t       = Elapsed Time between successive vectors (sec)
% e             = Eccentricity Vector of the orbit
% e_m           = Magnitude of the eccentricity vector
% f,g,dg        = Lagrange Coefficients
% f_z           = Shape Parameter function for newton iteration
% h             = Specific Angular Momentum vector (km^2/s)
% h_mag         = Magnitude of specific angular momentum (km^2/s)
% incl          = Inclination of the orbit (degree)
% K             = Vector Node along z-axis
% mu            = Gravitational Parameter (km^3/sec^2)
% N             = Node Line Vector
% N_M           = Magnitide of Node vector
% ohm           = Right Ascension of Ascending node (degrees)
% omega         = Argument of Perigee (degrees)
% r_a,r_p       = Apogee and Perigee Radii (km)
% R1 (I,J,K)    = Components of position vector 1 (km)
% R1_M          = Magnitude of vector R1 (km)
% R1_V          = Position Vector 1 (km)
% R2 (I,J,K)    = Components of position vector 2 (km)
% R2_M          = Magnitude of vector R2 (km)
% R2_V          = Position Vector 2 (km)
% R_cross       = Cross Product of R1 and R2 vectors
% s             = S(z) class of Stumpff Functions
% s_z           = S(z) class of Stumpff Functions array
% theta         = True Anamoly (degrees)
% V1_M          = Velocity magnitude of velocity vector V1 (km/s)
% V1_V          = Velocity vector corresponding to R1 (km/s)
% V2_V          = Velocity vector corresponding to R2 (km/s)
% VR            = Radial velocity vector (km/s)
% y_z           = Stumpff function subsitution assumption
% z             = Trajectory Shape Parameter
%%
% LP165P gravity model 
mu = 4902.80105600000; %km3/s2
Rm = 1738; %km
%R_EARTH'e gerek yok, R_MOON'u bir kere tanýmlamak yeterli 
%R_EARTH = 6371; % Radius of Earth
%R_MOON = 1738; % Radius of Moon

%Initial position and velocity (km, km/s)
R1_V       = [1838.2;0;0]; 
V1_initial = [0;0;1.633148794269431]; %100 km dairesel yörünge için gerekli hýz
%V1_initial = [1;1;1.303148794269431];

%Final position(on the surface of Moon) (km)
R2_V = [Rm;Rm;Rm];

%Time of flight(ToF) 
tof_minute = 30; %minute
tof_second = tof_minute*60; %sec 

%Aþaðýda iþaretlediðim bloðu ayrý bir fonksiyon olarak yazabilirsin. 
%Hem kod çok uzun olmaz, takip edilebilir olur hem de daha anlaþýlýr olur
%Örn:
%[V1_V, V2_V] = lambert(R1_V,R2_V,tof_second)  
%--------------------------------------------------------------------------
R1_M = norm(R1_V);

R2_M = norm(R2_V);

% Estimating nature of the trajectory and true anamoly:
R_cross     = cross(R1_V,R2_V);
del_theta   = acosd(dot(R1_V,R2_V)/(R1_M*R2_M));
incl        = acosd((R_cross(3))/(R1_M*R2_M*sind(del_theta)));

if cosd(incl) > 0
    if R_cross(3) >= 0
        del_theta = del_theta;
    else %R_cross(3) < 0
        del_theta = 360 - del_theta;
    end
else %cosd(incl) < 0
    if R_cross(3) < 0
        del_theta = del_theta;
    else %R_cross(3) >= 0
        del_theta = 360 - del_theta;
    end
end

A = sind(del_theta)*sqrt((R1_M*R2_M)/(1-cosd(del_theta)));

%z kullanýcý tarafýndan girilmemeli
%z'nin baþlangýç deðeri aþaðýdaki þekilde hesaplanabilir(Curtis)
%Aþaðýdaki hesap F(z)'nin iþaret deðiþtirdiði z'yi verecek
z     = -100;
value = -100;
while value < 0
    
    if z > 0
        s = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
    elseif z < 0
        s = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
    else
        s = 1/6;
    end
    
    %Defining C(z) class of stumpff functions:
    if z > 0
        c = (1 - cos(sqrt(z)))/z;
    elseif z < 0
        c = (cosh(sqrt(-z)) - 1)/(-z);
    else
        c = 1/2;
    end
    
    y_z   = R1_M + R2_M + A.*((z*s-1)/(sqrt(c)));
    f_z   = (((y_z)/(c))^1.5)*s + A*sqrt(y_z) - sqrt(mu)*tof_second;
    value = f_z; 
    z = z + 0.1;
end

%{
%  F(z) and F'(z) to predict the trajectory:
z       = (-10:0.1:10);
s_z     = [];
c_z     = [];

%source code
for i=1:201
    % Defininig S(z)of stumpff functions:
    if z(i) > 0
        s = (sqrt(z(i)) - sin(sqrt(z(i))))/(sqrt(z(i)))^3;
    elseif z(i) < 0
        s = (sinh(sqrt(-z(i))) - sqrt(-z(i)))/(sqrt(-z(i)))^3;
    else
        s = 1/6;
    end
    s_z     = [s_z,s];
    
    %Defining C(z) class of stumpff functions:
    if z(i) > 0
        c = (1 - cos(sqrt(z(i))))/z(i);
    elseif z(i) < 0
        c = (cosh(sqrt(-z(i))) - 1)/(-z(i));
    else
        c = 1/2;
    end
    c_z     = [c_z,c];
end
%here 

y_z     = R1_M + R2_M + A.*((z.*s_z-1)./(sqrt(c_z)));
f_z     = (((y_z)./(c_z)).^1.5).*s_z + A.*sqrt(y_z) - sqrt(mu)*tof_second;


figure(1)
hold on
grid on
plot(z,f_z,':','Linewidth',2)
xlabel('z','Fontweight','bold')
ylabel('F(z)','Fontweight','bold')
hold off

z       = input('\nEnter the value of z at inflection point: ');
%}

ratio   = 1; 
% source code
while abs(ratio) > 1.e-10

    % Defininig S(z) class of stumpff functions:
    if z > 0
        s = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
    elseif z < 0
        s = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
    else
        s = 1/6;
    end
    
    %Defining C(z) class of stumpff functions:
    if z > 0
        c = (1 - cos(sqrt(z)))/z;
    elseif z < 0
        c = (cosh(sqrt(-z)) - 1)/(-z);
    else
        c = 1/2;
    end
    %here
    
    %related formulas from the book
    y_z     = R1_M + R2_M + A.*((z.*s-1)./(sqrt(c)));
    f_z     = (((y_z)./(c)).^1.5).*s + A.*sqrt(y_z) - sqrt(mu)*tof_second;

    if z == 0
        df_z    = (sqrt(2)/40)*((y_z)^1.5) + (A/8)*(sqrt(y_z) + A*sqrt(1/2*y_z));
    else
        a1      = (1/(2*z))*(c-(1.5*(s/c))) + 0.75*(s^2/c);
        a2      = (3*(s/c)*sqrt(y_z)) + A*sqrt(c/y_z);
        df_z    = (((y_z/c)^1.5)*a1) + ((A/8)*a2);
    end

    ratio       = f_z/df_z;
    z           = z - ratio;
end

%related formulas from the book
% Calculating lagrange functions:
y_z     = R1_M + R2_M + A*((z*s-1)/(sqrt(c)));
f       = 1 - (y_z/R1_M);
g       = A*sqrt(y_z/mu);
dg      = 1 - (y_z/R2_M);

% Calculating velocity vectors with reference to position vectors:
V1_V    = (1/g)*(R2_V - f*R1_V);
V2_V    = (1/g)*(dg*R2_V - R1_V);
%--------------------------------------------------------------------------

% orbital elements :
V1_M    = norm(V1_V);
VR      = dot(R1_V,V1_V)/R1_M;
h       = cross(R1_V,V1_V);
h_mag   = norm(h);
incl    = acosd(h(3)/h_mag); 
K       = [0 0 1];
N       = cross(K,h);
N_M     = norm(N);

% RA of ascending node:
if N(2) >= 0
    ohm = acosd(N(1)/N_M);
else %N(2) < 0
    ohm = 360 - acosd(N(1)/N_M);
end

% Eccentricity:
e       = (1/mu)*(((V1_M^2)-(mu/R1_M))*R1_V - (R1_M*VR*(V1_V)));
e_m     = norm(e);

% Argument of Perigee
if e(3) >= 0
    omega   = acosd(dot(N,e)/(N_M*e_m));
else %e(3) < 0
    omega   = 360 - acosd(dot(N,e)/(N_M*e_m));
end

% True Anamoly
if VR >= 0
    theta = acosd(dot(e,R1_V)/(R1_M*e_m));
else %VR < 0
    theta = 360 - acosd(dot(e,R1_V)/(R1_M*e_m));
end

% Semi-major axis:
r_p     = (h_mag^2/mu)*(1/(1 + e_m));
r_a     = (h_mag^2/mu)*(1/(1 - e_m));
a       = 0.5*(r_p + r_a);

%Aþaðýdaki kýsma gerek yok(Minumum deltav'yi veren TOF'u bulmak için gerekli Loop
%lambertin baþýna yazman gerek)
%{
% Define the number of revolutions
nrev = 0;
% Calculate the distance between the starting and ending positions
r1r2 = norm(R1_V - R2_V);

% Calculate the transfer angle
c = r1r2 / 2;
s = sqrt((R_EARTH + R_MOON)^2 - r1r2^2) / 2;
transfer_angle = 2 * atan(sqrt((R_EARTH + R_MOON)^2 - r1r2^2) / r1r2);

% Calculate the mean anomaly at departure
M0 = transfer_angle - e * s;

% Calculate the time of flight
tf = M0 * sqrt(a^3 / mu) * (1 + nrev);

% Calculate the true anomaly at departure
v0 = 2 * atan(sqrt((1-e)/(1+e)) * tan(M0/2));

% Calculate the departure velocity
v1 = sqrt(mu/ a) .* (-sin(v0) + e .* sin(M0));

% Calculate the arrival velocity
v2 = sqrt(mu / a) .* (cos(v0) - e .* cos(M0));

% Define the time of flight range
tof_range = 1:100; % In seconds

% Iterate over different time of flights and find the minimum deltav
for i = 1:length(tof_range)
% Calculate the orbital elements and trajectory using the Lambert algorithm
% Calculate the deltav
deltav = norm(v1 - v2);
% Initialize the minimum deltav and corresponding time of flight
min_deltav = inf;
min_tof = 0;
% Update the minimum deltav and corresponding time of flight
if deltav < min_deltav
    min_deltav = deltav;
    min_tof = tof_range(i);
end
end
%}
% Print the results

%Maneuver calculation (km/s)
%delv = norm(V1_V-V1_initial) ;

% Echoing results from the lambert's problem:
clc

fprintf('--------------------------- Initial Input ---------------------------')
fprintf('\n Position Vector (R1)        : [%g, %g, %g] (km)', R1_V(1), R1_V(2), R1_V(3))
fprintf('\n Position Vector (R2)        : [%g, %g, %g] (km)', R2_V(1), R2_V(2), R2_V(3))
fprintf('\n Elapsed time                : %g seconds',tof_second)
fprintf('\n\n-------------------- Orbital Elements and Nature --------------------')
fprintf('\n Velocity Vector (V1)        : [%g, %g, %g] (km/s)', V1_V(1), V1_V(2), V1_V(3))
fprintf('\n Semi-Major Axis             : %g km', a)
fprintf('\n Specific Angular Momentum   : %g km^2/s', h_mag)
fprintf('\n Inclination of Orbit        : %g degrees', incl)
fprintf('\n RA of Ascending Node        : %g degrees', ohm)
fprintf('\n Argument of Perigee         : %g degrees', omega)
fprintf('\n True Anomaly                : %g degrees', theta)
fprintf('\n Eccentricity of Orbit       : %g', e_m)
% Print the results
% fprintf('Minimum deltav: %f m/s\n', min_deltav);
% fprintf('Corresponding time of flight: %d seconds\n', min_tof);


%fprintf('\n DeltaV       : %g', delv)


if e_m == 0
    fprintf('\n\n The orbit is circular.')
elseif e_m == 1
    fprintf('\n\n The orbit is parabolic.')
elseif e_m > 1
    fprintf('\n\n The orbit is hyperbolic.')
else 
    fprintf('\n\n The orbit is elliptic.')
end

if incl > 0
    if incl < 90
        fprintf('\n The given orbit is a prograde trajectory.')
    else
        fprintf('\n The given orbit is a retrograde trajectory.')
    end
else
end
fprintf('\n---------------------------------------------------------------------')