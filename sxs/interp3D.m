% this function interpolate the scattering LUT for a specific input
% combination, linear interpolation
% this is to suit the scattering LUT whose dimension is 5
% X, Y, Z: range/coordinate of each inputs
% Xq, Yq, Zq: specific value/ query points

function Vq = interp3D(X,Y,Z,V,Xq,Yq,Zq)

X_start = X(1);
X_res = X(2)-X(1);

Y_start = Y(1);
Y_res = Y(2)-Y(1);

Z_start = Z(1);
Z_res = Z(2)-Z(1);

% find bin index of query points
X_idx = ceil((Xq-X_start)/X_res);
Y_idx = ceil((Yq-Y_start)/Y_res);
Z_idx = ceil((Zq-Z_start)/Z_res);

X1 = X(X_idx);
Y1 = Y(Y_idx);
Z1 = Z(Z_idx);

% locate the data
V_1_1_1 = V(Y_idx,X_idx,Z_idx,:,:);
V_1_1_1 = squeeze(V_1_1_1);

V_1_2_1 = V(Y_idx,X_idx+1,Z_idx,:,:);
V_1_2_1 = squeeze(V_1_2_1);

V_2_1_1 = V(Y_idx+1,X_idx,Z_idx,:,:);
V_2_1_1 = squeeze(V_2_1_1);

V_2_2_1 = V(Y_idx+1,X_idx+1,Z_idx,:,:);
V_2_2_1 = squeeze(V_2_2_1);

V_1_1_2 = V(Y_idx,X_idx,Z_idx+1,:,:);
V_1_1_2 = squeeze(V_1_1_2);

V_1_2_2 = V(Y_idx,X_idx+1,Z_idx+1,:,:);
V_1_2_2 = squeeze(V_1_2_2);

V_2_1_2 = V(Y_idx+1,X_idx,Z_idx+1,:,:);
V_2_1_2 = squeeze(V_2_1_2);

V_2_2_2 = V(Y_idx+1,X_idx+1,Z_idx+1,:,:);
V_2_2_2 = squeeze(V_2_2_2);

% 2D interploation
d_X = Xq-X1;
d_Y = Yq-Y1;

V_1_1 = (V_1_1_1*(X_res-d_X)+d_X*V_1_2_1)/X_res;
V_1_2 = (V_2_1_1*(X_res-d_X)+d_X*V_2_2_1)/X_res;
V_1 = (V_1_1*(Y_res-d_Y)+d_Y*V_1_2)/Y_res;

V_2_1 = (V_1_1_2*(X_res-d_X)+d_X*V_1_2_2)/X_res;
V_2_2 = (V_2_1_2*(X_res-d_X)+d_X*V_2_2_2)/X_res;
V_2 = (V_2_1*(Y_res-d_Y)+d_Y*V_2_2)/Y_res;

% 3D interpolation
d_Z = Zq-Z1;
Vq = (V_1*(Z_res-d_Z)+d_Z*V_2)/Z_res;






















