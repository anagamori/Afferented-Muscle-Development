function [] = gridSpeedUpFcn(x)

% Initialize the grid and initial and final points
nx1 = x; nx2 =  x;
x1l =  0; x1u = 100;
x2l =  0; x2u = 100;

x1 = linspace(x1l,x1u,nx1+1);
x2 = linspace(x2l,x2u,nx2+1);

limsf1 = 1:nx1+1; limsf2 = 1:nx2+1;

% Initalize other variables
t = 1;
sigmax1 = 0.5; sigmax2 = 1;
sigma = t * [sigmax1^2 0; 0 sigmax2^2];
invSig = inv(sigma);
detSig = det(sigma);

expF = [1 0; 0 1];
n = size (expF, 1);
gausThresh = 10;

small = 0; subs = []; vals = [];

% Iterate through all possible initial
% and final positions and calculate
% the values of exponent and out
% if exponent > gausThresh.

for i1 = 1:nx1+1
    for i2 = 1:nx2+1
        for f1 = limsf1
            for f2 = limsf2

                % Initial and final position
                xi = [x1(i1) x2(i2)]';
                xf = [x1(f1) x2(f2)]';

                exponent = 0.5 * (xf - expF * xi)'...
                    * invSig * (xf - expF * xi);

                if exponent > gausThresh
                    small = small + 1;
                else
                    out = 1 / (sqrt((2 * pi)^n * detSig))...
                        * exp(-exponent);
                    subs = [subs; i1 i2 f1 f2];
                    vals = [vals; out];
                end

            end
        end
    end
end